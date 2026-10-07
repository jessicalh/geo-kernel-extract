package downloadtests

import (
	"archive/zip"
	"bytes"
	"context"
	"crypto/tls"
	"encoding/json"
	"encoding/pem"
	"errors"
	"fmt"
	"io"
	"net"
	"net/http"
	"net/http/httptest"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"
	"sync"
	"testing"
	"time"
)

const (
	payloadBytes = 2 << 20
	prefixBytes  = 512 << 10
	etagV1       = `"fixture-v1"`
	etagV2       = `"fixture-v2"`
	lastModified = "Mon, 02 Feb 2026 12:00:00 GMT"
)

type requestRecord struct {
	Range    string `json:"range"`
	IfRange  string `json:"if_range"`
	Protocol string `json:"protocol"`
	Status   int    `json:"status"`
	ETag     string `json:"etag"`
	Modified string `json:"last_modified"`
}

type transferCase struct {
	records  []requestRecord
	release  chan struct{}
	released bool
}

type peer struct {
	t       *testing.T
	payload []byte
	changed []byte
	bundle  []byte
	// Transfer handlers and the Qt control requests access these cases concurrently.
	mu    sync.Mutex
	cases map[string]*transferCase
}

func payload(version int) []byte {
	data := make([]byte, payloadBytes)
	for i := range data {
		data[i] = byte((i*17 + i/251 + (version-1)*73) % 251)
	}
	return data
}

func newPeer(t *testing.T) *peer {
	t.Helper()
	p := &peer{t: t, payload: payload(1), changed: payload(2), cases: make(map[string]*transferCase)}
	var buffer bytes.Buffer
	archive := zip.NewWriter(&buffer)
	for _, file := range []struct {
		name string
		data []byte
	}{{"run.LGS", []byte("{}")}, {"data.bin", p.payload}} {
		entry, err := archive.CreateHeader(&zip.FileHeader{Name: file.name, Method: zip.Store})
		if err != nil {
			t.Fatal(err)
		}
		if _, err := entry.Write(file.data); err != nil {
			t.Fatal(err)
		}
	}
	if err := archive.Close(); err != nil {
		t.Fatal(err)
	}
	p.bundle = buffer.Bytes()
	return p
}

func (p *peer) routes() http.Handler {
	mux := http.NewServeMux()
	mux.HandleFunc("GET /files/{id}/{mode}", p.transfer)
	mux.HandleFunc("GET /records/{id}", func(w http.ResponseWriter, r *http.Request) {
		p.mu.Lock()
		records := []requestRecord{}
		if c := p.cases[r.PathValue("id")]; c != nil {
			records = append(records, c.records...)
		}
		p.mu.Unlock()
		w.Header().Set("Content-Type", "application/json")
		if err := json.NewEncoder(w).Encode(records); err != nil {
			p.t.Errorf("write request records: %v", err)
		}
	})
	mux.HandleFunc("POST /release/{id}", func(w http.ResponseWriter, r *http.Request) {
		p.mu.Lock()
		c := p.cases[r.PathValue("id")]
		if c != nil && !c.released {
			close(c.release)
			c.released = true
		}
		p.mu.Unlock()
		if c == nil {
			http.Error(w, "unknown transfer", http.StatusNotFound)
			return
		}
		w.WriteHeader(http.StatusNoContent)
	})
	return mux
}

func (p *peer) transfer(w http.ResponseWriter, r *http.Request) {
	id, mode := r.PathValue("id"), r.PathValue("mode")
	p.mu.Lock()
	c := p.cases[id]
	if c == nil {
		c = &transferCase{release: make(chan struct{})}
		p.cases[id] = c
	}
	attempt := len(c.records) + 1
	p.mu.Unlock()

	body, etag := p.payload, etagV1
	validator, modified := etagV1, ""
	if mode == "date-validator" {
		etag, validator, modified = "", lastModified, lastModified
	}
	if strings.HasPrefix(mode, "bundle-") {
		body = p.bundle
	}
	if mode == "changed" && attempt > 1 {
		body, etag = p.changed, etagV2
	}
	offset, status := 0, http.StatusOK
	if value := r.Header.Get("Range"); value != "" {
		var err error
		if !strings.HasPrefix(value, "bytes=") || !strings.HasSuffix(value, "-") {
			p.t.Errorf("%s: unexpected Range %q", id, value)
			http.Error(w, "invalid test range", http.StatusBadRequest)
			return
		}
		offset, err = strconv.Atoi(strings.TrimSuffix(strings.TrimPrefix(value, "bytes="), "-"))
		if err != nil || offset <= 0 || offset >= len(body) {
			p.t.Errorf("%s: invalid resume offset %q", id, value)
			http.Error(w, "invalid resume offset", http.StatusBadRequest)
			return
		}
		if r.Header.Get("If-Range") != validator {
			p.t.Errorf("%s: resume lost its validator: If-Range=%q", id, r.Header.Get("If-Range"))
		}
		status = http.StatusPartialContent
	}
	if mode == "changed" || mode == "ignored-range" {
		offset, status = 0, http.StatusOK
	}
	if attempt == 2 {
		switch mode {
		case "unavailable":
			status = http.StatusServiceUnavailable
		case "unsatisfiable":
			status = http.StatusRequestedRangeNotSatisfiable
		case "wrong-etag":
			etag = etagV2
		}
	}
	record := requestRecord{r.Header.Get("Range"), r.Header.Get("If-Range"), r.Proto, status, etag, modified}
	p.mu.Lock()
	c.records = append(c.records, record)
	p.mu.Unlock()
	p.t.Logf("%s attempt=%d protocol=%s Range=%q If-Range=%q status=%d ETag=%s",
		id, attempt, record.Protocol, record.Range, record.IfRange, status, etag)

	if etag != "" {
		w.Header().Set("ETag", etag)
	}
	if modified != "" {
		w.Header().Set("Last-Modified", modified)
		w.Header().Set("Date", "Mon, 02 Mar 2026 12:00:00 GMT")
	}
	w.Header().Set("Accept-Ranges", "bytes")
	w.Header().Set("Content-Type", "application/octet-stream")
	if status == http.StatusServiceUnavailable || status == http.StatusRequestedRangeNotSatisfiable {
		if status == http.StatusRequestedRangeNotSatisfiable {
			w.Header().Set("Content-Range", fmt.Sprintf("bytes */%d", len(body)))
		}
		w.Header().Set("Content-Length", "0")
		w.WriteHeader(status)
		return
	}
	if status == http.StatusPartialContent {
		start, total := offset, len(body)
		if attempt == 2 && mode == "wrong-offset" {
			start++
		}
		if attempt == 2 && mode == "wrong-total" {
			total++
		}
		w.Header().Set("Content-Range", fmt.Sprintf("bytes %d-%d/%d", start, len(body)-1, total))
	}
	w.Header().Set("Content-Length", strconv.Itoa(len(body)-offset))
	w.WriteHeader(status)

	if attempt == 1 && mode != "slow" && mode != "bundle-complete" {
		if err := writeFlush(w, body[:prefixBytes]); err != nil {
			p.t.Logf("%s prefix interrupted: %v", id, err)
			return
		}
		select {
		case <-c.release:
			p.t.Logf("%s: forcing disconnect after acknowledged prefix", id)
			if r.ProtoMajor == 1 {
				connection, _, err := http.NewResponseController(w).Hijack()
				if err != nil {
					p.t.Errorf("%s: hijack interrupted transfer: %v", id, err)
					return
				}
				// Reset TCP directly: tls.Conn.Close would send a graceful close_notify first.
				tcp := connection.(*tls.Conn).NetConn().(*net.TCPConn)
				lingerErr := tcp.SetLinger(0)
				closeErr := tcp.Close()
				if err := errors.Join(lingerErr, closeErr); err != nil {
					p.t.Errorf("%s: reset interrupted transfer: %v", id, err)
				}
				return
			}
			// net/http's documented sentinel resets the HTTP/2 stream.
			panic(http.ErrAbortHandler)
		case <-r.Context().Done():
			p.t.Logf("%s: held response cancelled: %v", id, r.Context().Err())
			return
		}
	}
	if mode == "slow" {
		// Pacing is intentional: the transfer outlasts the injected inactivity timeout.
		ticker := time.NewTicker(25 * time.Millisecond)
		defer ticker.Stop()
		for offset < len(body) {
			end := min(offset+(32<<10), len(body))
			if err := writeFlush(w, body[offset:end]); err != nil {
				p.t.Logf("%s slow response interrupted: %v", id, err)
				return
			}
			offset = end
			select {
			case <-ticker.C:
			case <-r.Context().Done():
				return
			}
		}
		return
	}
	if err := writeFlush(w, body[offset:]); err != nil {
		// Malformed 206 tests deliberately reject headers and close the body.
		p.t.Logf("%s response interrupted: %v", id, err)
	}
}

func writeFlush(w http.ResponseWriter, data []byte) error {
	if _, err := w.Write(data); err != nil {
		return err
	}
	return http.NewResponseController(w).Flush()
}

func TestQtResumableDownloads(t *testing.T) {
	executable := os.Getenv("H5READER_DOWNLOAD_TEST")
	if executable == "" {
		t.Fatal("Set H5READER_DOWNLOAD_TEST to h5reader_resumable_download_tests, or run its CTest target")
	}
	for _, http2 := range []bool{false, true} {
		t.Run("HTTP2="+strconv.FormatBool(http2), func(t *testing.T) {
			p := newPeer(t)
			server := httptest.NewUnstartedServer(p.routes())
			server.EnableHTTP2 = http2
			server.StartTLS()
			defer server.Close()
			certificate := filepath.Join(t.TempDir(), "ca.pem")
			if err := os.WriteFile(certificate,
				pem.EncodeToMemory(&pem.Block{Type: "CERTIFICATE", Bytes: server.Certificate().Raw}), 0600); err != nil {
				t.Fatal(err)
			}
			ctx, cancel := context.WithTimeout(t.Context(), 90*time.Second)
			defer cancel()
			command := exec.CommandContext(ctx, executable, "-o", "-,txt")
			command.Env = append(os.Environ(), "QT_FORCE_STDERR_LOGGING=1",
				"H5READER_DOWNLOAD_URL="+server.URL+"/", "H5READER_DOWNLOAD_CA="+certificate,
				"H5READER_DOWNLOAD_HTTP2="+strconv.FormatBool(http2))
			command.Stdout = os.Stdout
			command.Stderr = os.Stderr
			if err := command.Run(); err != nil {
				t.Fatalf("Reader Qt tests: %v (context: %v)", err, ctx.Err())
			}
			want := "HTTP/1.1"
			if http2 {
				want = "HTTP/2.0"
			}
			p.mu.Lock()
			defer p.mu.Unlock()
			if len(p.cases) == 0 {
				t.Fatal("Qt tests did not make any fixture requests")
			}
			for id, c := range p.cases {
				for _, record := range c.records {
					if record.Protocol != want {
						t.Errorf("%s: negotiated %s, want %s", id, record.Protocol, want)
					}
				}
			}
		})
	}
}

func TestBundleFixture(t *testing.T) {
	p := newPeer(t)
	archive, err := zip.NewReader(bytes.NewReader(p.bundle), int64(len(p.bundle)))
	if err != nil {
		t.Fatal(err)
	}
	if len(archive.File) != 2 {
		t.Fatalf("bundle contains %d files, want 2", len(archive.File))
	}
	for i, want := range [][]byte{[]byte("{}"), p.payload} {
		file, err := archive.File[i].Open()
		if err != nil {
			t.Fatal(err)
		}
		data, readErr := io.ReadAll(file)
		closeErr := file.Close()
		if readErr != nil || closeErr != nil {
			t.Fatalf("read fixture: %v; close: %v", readErr, closeErr)
		}
		if !bytes.Equal(data, want) {
			t.Fatalf("wrong bytes in %s", archive.File[i].Name)
		}
	}
}

func TestPeerDisconnect(t *testing.T) {
	for _, http2 := range []bool{false, true} {
		t.Run("HTTP2="+strconv.FormatBool(http2), func(t *testing.T) {
			p := newPeer(t)
			server := httptest.NewUnstartedServer(p.routes())
			server.EnableHTTP2 = http2
			server.StartTLS()
			defer server.Close()
			client := server.Client()
			client.Timeout = 5 * time.Second
			response, err := client.Get(server.URL + "/files/probe/disconnect")
			if err != nil {
				t.Fatal(err)
			}
			defer response.Body.Close()
			wantProtocol := 1
			if http2 {
				wantProtocol = 2
			}
			if response.ProtoMajor != wantProtocol {
				t.Fatalf("negotiated %s, want HTTP/%d", response.Proto, wantProtocol)
			}
			prefix := make([]byte, prefixBytes)
			if _, err := io.ReadFull(response.Body, prefix); err != nil {
				t.Fatal(err)
			}
			if !bytes.Equal(prefix, p.payload[:prefixBytes]) {
				t.Fatal("incorrect prefix before disconnect")
			}
			release, err := client.Post(server.URL+"/release/probe", "application/octet-stream", nil)
			if err != nil {
				t.Fatal(err)
			}
			if err := release.Body.Close(); err != nil {
				t.Fatal(err)
			}
			if release.StatusCode != http.StatusNoContent {
				t.Fatalf("release returned HTTP %d", release.StatusCode)
			}
			remainder, err := io.ReadAll(response.Body)
			if err == nil || errors.Is(err, context.DeadlineExceeded) {
				t.Fatalf("expected immediate transport failure, got %v", err)
			}
			if len(remainder) != 0 {
				t.Fatalf("received %d bytes after acknowledged prefix", len(remainder))
			}
			t.Logf("peer disconnect observed: %v", err)
		})
	}
}
