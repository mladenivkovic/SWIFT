#! /bin/bash

set -u 

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SWIFT_ROOT="$(dirname "$SCRIPT_DIR")"

AUTOMATION_DIR="${SWIFT_ROOT}/HydroAutomation"
BUILD_LOG_DIR="${AUTOMATION_DIR}/build_logs"
BINARY_DIR="${AUTOMATION_DIR}/binaries"
RESULTS_DIR="${AUTOMATION_DIR}/results"
MAKE_JOBS=32

mkdir -p "$BUILD_LOG_DIR" "$BINARY_DIR" "$RESULTS_DIR" 


log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*"
}

log "Test message, automation directories are ready."