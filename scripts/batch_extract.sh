#!/usr/bin/env bash
#
# batch_extract.sh -- Run nmr_extract for every protein in consolidated/.
#
# Usage:
#   batch_extract.sh OUTPUT_DIR [PROTEIN_ID]
#
# First argument:   output base directory (required).
# Second argument:  a single protein ID to process (optional; all if omitted).
#
# Each protein produces two calls:
#   1. --mutant  (WT vs ALA comparison, delta arrays)
#   2. --orca    (ALA standalone features)
#
# Proteins are skipped if OUTPUT_DIR/PROT/done.marker already exists.
# On success, done.marker is written; on failure, the error is logged and
# the script continues to the next protein.
#
# All file paths passed to nmr_extract are fully qualified -- no internal
# directory scanning.  When multiple _nmr.out files match the glob (reruns),
# the lexicographically last one is chosen (latest timestamp in the filename).

set -euo pipefail

readonly CONSOLIDATED="/shared/2026Thesis/consolidated"
readonly NMR_EXTRACT="/shared/2026Thesis/nmr-shielding/build/nmr_extract"

# ── argument handling ────────────────────────────────────────────────
if [[ $# -lt 1 ]]; then
    echo >&2 "Usage: $0 OUTPUT_DIR [PROTEIN_ID]"
    exit 1
fi

readonly OUTPUT_BASE="$(realpath "$1")"
readonly SINGLE_PROTEIN="${2:-}"

if [[ ! -x "${NMR_EXTRACT}" ]]; then
    echo >&2 "FATAL: nmr_extract not found or not executable at ${NMR_EXTRACT}"
    exit 1
fi

if [[ ! -d "${CONSOLIDATED}" ]]; then
    echo >&2 "FATAL: consolidated directory not found at ${CONSOLIDATED}"
    exit 1
fi

mkdir -p "${OUTPUT_BASE}"

# ── helpers ──────────────────────────────────────────────────────────

# Pick the lexicographically last matching _nmr.out file (latest timestamp).
# Returns empty string if no match.
pick_nmr() {
    local pattern="$1"
    local files
    files=( ${pattern} )   # intentional word-splitting on glob
    if [[ "${files[0]}" == "${pattern}" ]]; then
        # glob matched nothing -- the literal pattern string was returned
        echo ""
        return
    fi
    # Sort and take last (timestamps in filenames are YYYYMMDD_HHMMSS, so
    # lexicographic order == chronological order).
    printf '%s\n' "${files[@]}" | sort | tail -n1
}

# ── build protein list ───────────────────────────────────────────────
if [[ -n "${SINGLE_PROTEIN}" ]]; then
    if [[ ! -d "${CONSOLIDATED}/${SINGLE_PROTEIN}" ]]; then
        echo >&2 "FATAL: protein directory not found: ${CONSOLIDATED}/${SINGLE_PROTEIN}"
        exit 1
    fi
    proteins=( "${SINGLE_PROTEIN}" )
else
    proteins=()
    for d in "${CONSOLIDATED}"/*/; do
        proteins+=( "$(basename "$d")" )
    done
fi

total=${#proteins[@]}
succeeded=0
failed=0
skipped=0
fail_list=()

echo >&2 "batch_extract: ${total} protein(s) to consider, output -> ${OUTPUT_BASE}"

# ── main loop (sequential) ───────────────────────────────────────────
for prot in "${proteins[@]}"; do
    dir="${CONSOLIDATED}/${prot}"
    out="${OUTPUT_BASE}/${prot}"

    # Skip if already done.
    if [[ -f "${out}/done.marker" ]]; then
        echo >&2 "[SKIP] ${prot} -- done.marker exists"
        skipped=$((skipped + 1))
        continue
    fi

    echo >&2 "[START] ${prot}"
    t_start=${SECONDS}

    # ── Resolve all input files ──────────────────────────────────────
    wt_xyz="${dir}/${prot}_WT.xyz"
    wt_prmtop="${dir}/${prot}_WT.prmtop"
    wt_nmr="$(pick_nmr "${dir}/${prot}_WT_*_nmr.out")"

    ala_xyz="${dir}/${prot}_ALA.xyz"
    ala_prmtop="${dir}/${prot}_ALA.prmtop"
    ala_nmr="$(pick_nmr "${dir}/${prot}_ALA_*_nmr.out")"

    # Check minimum requirements: both xyz, both prmtop, at least ALA nmr.
    missing=""
    [[ ! -f "${wt_xyz}" ]]    && missing="${missing} WT.xyz"
    [[ ! -f "${wt_prmtop}" ]] && missing="${missing} WT.prmtop"
    [[ ! -f "${ala_xyz}" ]]   && missing="${missing} ALA.xyz"
    [[ ! -f "${ala_prmtop}" ]] && missing="${missing} ALA.prmtop"
    [[ -z "${ala_nmr}" ]]     && missing="${missing} ALA_nmr.out"

    if [[ -n "${missing}" ]]; then
        echo >&2 "[FAIL] ${prot} -- missing files:${missing}"
        failed=$((failed + 1))
        fail_list+=( "${prot}" )
        continue
    fi

    mkdir -p "${out}/wt" "${out}/ala"

    # ── Call 1: --mutant (WT vs ALA) ─────────────────────────────────
    # Requires WT nmr; skip mutant call if absent, but still run ALA.
    mutant_ok=true
    if [[ -n "${wt_nmr}" ]]; then
        if ! "${NMR_EXTRACT}" --mutant \
                --wt-xyz "${wt_xyz}" \
                --wt-prmtop "${wt_prmtop}" \
                --wt-nmr "${wt_nmr}" \
                --ala-xyz "${ala_xyz}" \
                --ala-prmtop "${ala_prmtop}" \
                --ala-nmr "${ala_nmr}" \
                --output "${out}/wt" \
                2>&1; then
            echo >&2 "[FAIL] ${prot} -- nmr_extract --mutant returned non-zero"
            mutant_ok=false
        fi
    else
        echo >&2 "[WARN] ${prot} -- no WT nmr.out, skipping --mutant call"
        mutant_ok=false
    fi

    # ── Call 2: --orca (ALA standalone) ──────────────────────────────
    ala_ok=true
    if ! "${NMR_EXTRACT}" --orca \
            --xyz "${ala_xyz}" \
            --prmtop "${ala_prmtop}" \
            --nmr "${ala_nmr}" \
            --output "${out}/ala" \
            2>&1; then
        echo >&2 "[FAIL] ${prot} -- nmr_extract --orca returned non-zero"
        ala_ok=false
    fi

    elapsed=$(( SECONDS - t_start ))

    if ${mutant_ok} && ${ala_ok}; then
        date -Iseconds > "${out}/done.marker"
        echo >&2 "[OK]   ${prot} -- ${elapsed}s"
        succeeded=$((succeeded + 1))
    else
        echo >&2 "[FAIL] ${prot} -- ${elapsed}s (mutant_ok=${mutant_ok}, ala_ok=${ala_ok})"
        failed=$((failed + 1))
        fail_list+=( "${prot}" )
    fi
done

# ── summary ──────────────────────────────────────────────────────────
echo >&2 ""
echo >&2 "===== batch_extract complete ====="
echo >&2 "  succeeded: ${succeeded}"
echo >&2 "  failed:    ${failed}"
echo >&2 "  skipped:   ${skipped}"
echo >&2 "  total:     ${total}"

if [[ ${#fail_list[@]} -gt 0 ]]; then
    echo >&2 ""
    echo >&2 "Failed proteins:"
    printf >&2 '  %s\n' "${fail_list[@]}"
fi

# Exit non-zero only if everything failed.
if [[ ${succeeded} -eq 0 && ${failed} -gt 0 ]]; then
    exit 1
fi
