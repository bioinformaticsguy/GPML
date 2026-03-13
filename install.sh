#!/usr/bin/env bash
# install.sh
#
# Downloads and sets up all external databases required by the GPML pipeline.
# Safe to re-run: already-downloaded files are resumed, already-verified files
# are skipped, already-extracted files are not overwritten.
#
# Usage:
#   bash install.sh [--data-dir <path>] [--skip-extraction]
#
# Options:
#   --data-dir <path>      Base directory for all database downloads.
#                          Default: Data/  (relative to this script's location)
#   --skip-extraction      Download and verify checksums only; do not unzip.
#                          Useful when downloading on a login node and
#                          extracting later on a compute node.

set -euo pipefail

# ---------------------------------------------------------------------------
# Resolve script directory so paths are correct regardless of where the
# script is invoked from.
# ---------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

DATA_DIR="${SCRIPT_DIR}/Data"
SKIP_EXTRACTION=false
LOG_FILE=""          # set after DATA_DIR is finalised

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

while [[ $# -gt 0 ]]; do
    case "$1" in
        --data-dir)
            [[ -n "${2:-}" ]] || { echo "ERROR: --data-dir requires a value." >&2; exit 1; }
            DATA_DIR="$(realpath -m "$2")"   # resolve to absolute path
            shift 2
            ;;
        --skip-extraction)
            SKIP_EXTRACTION=true
            shift
            ;;
        -h|--help)
            sed -n '2,14p' "${BASH_SOURCE[0]}"   # print the header comment
            exit 0
            ;;
        *)
            echo "ERROR: Unknown argument: $1" >&2
            echo "Run 'bash install.sh --help' for usage." >&2
            exit 1
            ;;
    esac
done

# ---------------------------------------------------------------------------
# Logging — writes to both stdout and a log file
# ---------------------------------------------------------------------------

mkdir -p "${DATA_DIR}"
LOG_FILE="${DATA_DIR}/install.log"

log() {
    local msg="[$(date '+%Y-%m-%d %H:%M:%S')] $*"
    echo "${msg}"
    echo "${msg}" >> "${LOG_FILE}"
}

die() {
    log "ERROR: $*"
    exit 1
}

# ---------------------------------------------------------------------------
# Pre-flight: required tools
# ---------------------------------------------------------------------------

check_tool() {
    local tool="$1"
    command -v "${tool}" >/dev/null 2>&1 || \
        die "'${tool}' is not available. On HPC, try: module load ${tool}"
}

log "=== GPML install.sh starting ==="
log "Script dir : ${SCRIPT_DIR}"
log "Data dir   : ${DATA_DIR}"
log "Log file   : ${LOG_FILE}"
log ""
log "--- Checking required tools ---"

check_tool curl
check_tool java
check_tool unzip
check_tool md5sum

log "All required tools found."
log ""

# ---------------------------------------------------------------------------
# Helper: download a single file
#
# Usage: download_file <dest_dir> <url> [large]
#   large=true  uses '--http1.1 -C -' for resume support on large files
# ---------------------------------------------------------------------------

download_file() {
    local dest_dir="$1"
    local url="$2"
    local large="${3:-false}"
    local filename
    filename="$(basename "${url}")"
    local dest="${dest_dir}/${filename}"

    if [[ -f "${dest}" ]]; then
        log "  Found locally, will resume if incomplete: ${filename}"
    else
        log "  Downloading: ${filename}"
    fi

    if [[ "${large}" == "true" ]]; then
        curl --http1.1 -C - -O --output-dir "${dest_dir}" "${url}" \
            || die "Download failed: ${url}"
    else
        curl --http1.1 -O --output-dir "${dest_dir}" "${url}" \
            || die "Download failed: ${url}"
    fi
}

# ---------------------------------------------------------------------------
# Helper: verify an MD5 checksum
#
# Handles two .md5 file formats:
#   1. Plain hash only:          "d41d8cd98f00b204e9800998ecf8427e"
#   2. Hash + filename (md5sum): "d41d8cd98f00b204e9800998ecf8427e  filename"
# ---------------------------------------------------------------------------

verify_md5() {
    local file="$1"       # absolute path to file being verified
    local md5_file="${file}.md5"

    [[ -f "${md5_file}" ]] || die "MD5 file not found: ${md5_file}"
    [[ -f "${file}" ]]     || die "File to verify not found: ${file}"

    local expected actual
    expected=$(awk '{print $1}' "${md5_file}" | tr '[:upper:]' '[:lower:]')
    actual=$(md5sum "${file}" | awk '{print $1}' | tr '[:upper:]' '[:lower:]')

    if [[ "${expected}" == "${actual}" ]]; then
        log "  Checksum OK: $(basename "${file}")"
    else
        die "Checksum mismatch for $(basename "${file}")
         Expected : ${expected}
         Got      : ${actual}
         The file may be corrupt or incomplete. Delete it and re-run."
    fi
}

# ---------------------------------------------------------------------------
# Helper: check available disk space
#
# Usage: check_disk_space <directory> <required_gb>
# ---------------------------------------------------------------------------

check_disk_space() {
    local dir="$1"
    local required_gb="$2"
    local available_kb available_gb

    available_kb=$(df -P "${dir}" | awk 'NR==2 {print $4}')
    available_gb=$(( available_kb / 1024 / 1024 ))

    if (( available_gb < required_gb )); then
        die "Insufficient disk space in ${dir}.
         Required  : ~${required_gb} GB
         Available : ~${available_gb} GB
         Free up space or specify a different --data-dir."
    else
        log "  Disk space OK: ~${available_gb} GB available (need ~${required_gb} GB)."
    fi
}

# ===========================================================================
# DATABASE: dbNSFP 5.3.1a
# ===========================================================================

log "=== Installing dbNSFP 5.3.1a ==="

DBNSFP_DIR="${DATA_DIR}/dbNSFP"
DBNSFP_BASE_URL="https://dist.genos.us/academic/e55b09"

mkdir -p "${DBNSFP_DIR}"

# ---- Disk space check -------------------------------------------------------
# zip (~50 GB) + grch38 (~50 GB) + gene.gz (a few GB) + extracted zip (~50 GB)
log "--- Checking disk space (need ~200 GB) ---"
check_disk_space "${DBNSFP_DIR}" 200

# ---- Step 1: Download small files first (checksums, index) ------------------
log "--- Step 1/4: Downloading checksum and index files ---"

download_file "${DBNSFP_DIR}" "${DBNSFP_BASE_URL}/dbNSFP5.3.1a.zip.md5"
download_file "${DBNSFP_DIR}" "${DBNSFP_BASE_URL}/dbNSFP5.3.1a_grch38.gz.md5"
download_file "${DBNSFP_DIR}" "${DBNSFP_BASE_URL}/dbNSFP5.3.1a_grch38.gz.tbi"

# ---- Step 2: Download large files (resume-capable) --------------------------
log "--- Step 2/4: Downloading large database files (this will take a while) ---"

download_file "${DBNSFP_DIR}" "${DBNSFP_BASE_URL}/dbNSFP5.3.1a.zip"         true  # ~50 GB
download_file "${DBNSFP_DIR}" "${DBNSFP_BASE_URL}/dbNSFP5.3_gene.gz"        true  # gene annotations
download_file "${DBNSFP_DIR}" "${DBNSFP_BASE_URL}/dbNSFP5.3.1a_grch38.gz"   true  # ~50 GB

# ---- Step 3: Verify checksums -----------------------------------------------
log "--- Step 3/4: Verifying checksums ---"

verify_md5 "${DBNSFP_DIR}/dbNSFP5.3.1a.zip"
verify_md5 "${DBNSFP_DIR}/dbNSFP5.3.1a_grch38.gz"

# ---- Step 4: Extract --------------------------------------------------------
if [[ "${SKIP_EXTRACTION}" == "true" ]]; then
    log "--- Step 4/4: Skipping extraction (--skip-extraction flag set) ---"
else
    log "--- Step 4/4: Extracting dbNSFP5.3.1a.zip ---"
    # -n: never overwrite existing files (safe to re-run)
    unzip -n "${DBNSFP_DIR}/dbNSFP5.3.1a.zip" -d "${DBNSFP_DIR}" \
        || die "Extraction failed."
    log "  Extraction complete."

    # Confirm the Java search tool is present
    if [[ -f "${DBNSFP_DIR}/search_dbNSFP53a.class" ]]; then
        log "  search_dbNSFP53a.class found — ready to use."
    else
        log "  WARNING: search_dbNSFP53a.class not found after extraction."
        log "           Contents of ${DBNSFP_DIR}:"
        ls "${DBNSFP_DIR}" | sed 's/^/             /' | tee -a "${LOG_FILE}"
    fi
fi

log ""
log "=== dbNSFP 5.3.1a installation complete ==="
log ""
log "  Database location : ${DBNSFP_DIR}"
log "  To query dbNSFP   : cd ${DBNSFP_DIR} && java search_dbNSFP53a -i <input> -o <output> -v hg38"
log ""

# ===========================================================================
# Add future database installs below this line, following the same pattern.
# ===========================================================================

log "=== All installations complete. Full log: ${LOG_FILE} ==="
