#pragma once
#include <cmath>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

extern "C" {
#include "xdrfile.h"
#include "xdrfile_xtc.h"
}

namespace fs = std::filesystem;

// xdrfile API uses non-const char* — helper to cast safely.
inline char* cpath(const fs::path& p) {
    thread_local std::string buf;
    buf = p.string();
    return buf.data();
}

// Coordinates for a single frame (in nm, as stored in .xtc).
struct XtcFrame {
    float              time;
    int                step;
    int                natoms;
    std::vector<float> x; // natoms * 3, interleaved [x0,y0,z0, x1,y1,z1, ...]
    float              box[3][3];
};

// Read a specific frame from an .xtc file by matching target time (ps).
// Reads sequentially from the start — call once per file with sorted times
// for efficiency. For our use case (~1000 frames, <15 MB), this is fast.
inline XtcFrame read_xtc_frame(const fs::path& xtc_path, double target_time,
                               double tol = 0.5)
{
    int natoms = 0;
    if (read_xtc_natoms(cpath(xtc_path), &natoms) != exdrOK)
        throw std::runtime_error("Cannot read natoms from " + xtc_path.string());

    XDRFILE* xd = xdrfile_open(cpath(xtc_path), "r");
    if (!xd)
        throw std::runtime_error("Cannot open " + xtc_path.string());

    XtcFrame result{};
    result.natoms = natoms;
    result.x.resize(natoms * 3);
    auto* coords = reinterpret_cast<rvec*>(result.x.data());

    int   step;
    float time, prec;
    matrix box;

    while (read_xtc(xd, natoms, &step, &time, box, coords, &prec) == exdrOK) {
        if (std::fabs(time - target_time) <= tol) {
            result.time = time;
            result.step = step;
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    result.box[i][j] = box[i][j];
            xdrfile_close(xd);
            return result;
        }
        if (time > target_time + tol) break;
    }
    xdrfile_close(xd);
    throw std::runtime_error("Frame at t=" + std::to_string(target_time) +
                             " not found in " + xtc_path.string());
}

// Read multiple frames from one .xtc (times must be sorted ascending).
// Much more efficient than calling read_xtc_frame() repeatedly.
inline std::vector<XtcFrame> read_xtc_frames(
    const fs::path& xtc_path,
    const std::vector<double>& sorted_times,
    double tol = 0.5)
{
    if (sorted_times.empty()) return {};

    int natoms = 0;
    if (read_xtc_natoms(cpath(xtc_path), &natoms) != exdrOK)
        throw std::runtime_error("Cannot read natoms from " + xtc_path.string());

    XDRFILE* xd = xdrfile_open(cpath(xtc_path), "r");
    if (!xd)
        throw std::runtime_error("Cannot open " + xtc_path.string());

    std::vector<XtcFrame> results;
    results.reserve(sorted_times.size());

    std::vector<float> buf(natoms * 3);
    auto* coords = reinterpret_cast<rvec*>(buf.data());
    int   step;
    float time, prec;
    matrix box;
    size_t tidx = 0;

    while (tidx < sorted_times.size() &&
           read_xtc(xd, natoms, &step, &time, box, coords, &prec) == exdrOK)
    {
        if (std::fabs(time - sorted_times[tidx]) <= tol) {
            XtcFrame f;
            f.time   = time;
            f.step   = step;
            f.natoms = natoms;
            f.x.assign(buf.begin(), buf.end());
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    f.box[i][j] = box[i][j];
            results.push_back(std::move(f));
            ++tidx;
        }
        if (time > sorted_times.back() + tol) break;
    }
    xdrfile_close(xd);

    if (results.size() != sorted_times.size())
        throw std::runtime_error("Only found " + std::to_string(results.size()) +
                                 "/" + std::to_string(sorted_times.size()) +
                                 " frames in " + xtc_path.string());
    return results;
}

// Read ALL frames from an .xtc file sequentially.
// Returns every frame from start to EOF.  For 50ns of metadynamics
// with 10ps output interval this is ~5000 frames — a few seconds of I/O.
inline std::vector<XtcFrame> read_all_xtc_frames(const fs::path& xtc_path)
{
    int natoms = 0;
    if (read_xtc_natoms(cpath(xtc_path), &natoms) != exdrOK)
        throw std::runtime_error("Cannot read natoms from " + xtc_path.string());

    XDRFILE* xd = xdrfile_open(cpath(xtc_path), "r");
    if (!xd)
        throw std::runtime_error("Cannot open " + xtc_path.string());

    std::vector<XtcFrame> results;
    std::vector<float> buf(natoms * 3);
    auto* coords = reinterpret_cast<rvec*>(buf.data());
    int   step;
    float time, prec;
    matrix box;

    while (read_xtc(xd, natoms, &step, &time, box, coords, &prec) == exdrOK) {
        XtcFrame f;
        f.time   = time;
        f.step   = step;
        f.natoms = natoms;
        f.x.assign(buf.begin(), buf.end());
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                f.box[i][j] = box[i][j];
        results.push_back(std::move(f));
    }
    xdrfile_close(xd);

    if (results.empty())
        throw std::runtime_error("No frames in " + xtc_path.string());

    return results;
}
