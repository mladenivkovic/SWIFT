#! /bin/bash

set -u #Treat unset variables as errors


#Get the swift directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SWIFT_ROOT="$(dirname "$SCRIPT_DIR")"

AUTOMATION_DIR="${SWIFT_ROOT}/HydroAutomation"
BUILD_LOG_DIR="${AUTOMATION_DIR}/build_logs"
BINARY_DIR="${AUTOMATION_DIR}/binaries"
RESULTS_DIR="${AUTOMATION_DIR}/results"
MAKE_JOBS=32

#Optional variables:
# FORCE_REBUILD=1 rebuild binaries if already stashed
# FRESH_IC=1 delete existing inital conditions before each run, forcing regeneration

# dimension -> example directory relative to SWIFT_ROOT
declare -A EXAMPLES=(
    [2]="examples/HydroTests/SodShock_2D"
    [3]="examples/HydroTests/SodShock_3D"
)

mkdir -p "$BUILD_LOG_DIR" "$BINARY_DIR" "$RESULTS_DIR" 

#logs a time based message
log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*"
}

#Loads the exact modules SWIFT needs to configure/build
#Prefer the saved module with my_swift
#Gives the modules to those who don't have my_swift saved
load_modules() {
    module purge

    if module restore my_swift > /dev/null 2>&1; then
        log "Modules restored from my_swift"
    else
        log "Failed to restore modules from my_swift, loading default modules"
        module load intel_comp/2024.2.0 compiler-rt tbb compiler mpi
        module load ucx/1.17.0
        module load parallel_hdf5/1.14.4
        module load fftw/3.3.10
        module load parmetis/4.0.3-64bit
        module load gsl
        module load sundials/5.8.0_c8_single
    fi
}

#Configures and compiles SWIFT for each dimension
#stasges the binary as binaries/swift_<dim>d.
#skips rebuild if the binary already exists unless FORCE_REBUILD=1
build_binary() {
    local dim="$1"
    local binary_path="${BINARY_DIR}/swift_${dim}d"

    if [ -e "$binary_path" ] && [ "${FORCE_REBUILD:-0}" != "1" ]; then
        log "Binary for ${dim}D already exists at ${binary_path}, skipping build."
        return 0
    fi

    cd "$SWIFT_ROOT" || return 1
    load_modules

    log "Configuring SWIFT for hydro-dimension=${dim}..."
    if ! ./configure --with-hydro-dimension="${dim}" > "${BUILD_LOG_DIR}/configure_${dim}d.log" 2>&1; then
        log "Configuration failed for ${dim}D. Check ${BUILD_LOG_DIR}/configure_${dim}d.log for details."
        return 1
    fi

    log "Building SWIFT (${dim}D) with make -j${MAKE_JOBS}... (this takes a while)"
    if ! make -j"${MAKE_JOBS}" > "${BUILD_LOG_DIR}/make_${dim}d.log" 2>&1; then
        log "ERRPR: make failed for ${dim}D. See ${BUILD_LOG_DIR}/make_${dim}d.log"
        return 1
    fi

    cp "${SWIFT_ROOT}/swift" "$binary_path"
    log "Stashed ${dim}D binary at ${binary_path}"
    return 0
    
}

#Copies the stashed binary for a given dimension into place as ./swift
#since each run.sh referes to ./swift as as a relative path
activate_binary() {
    local dim="$1"
    local binary_path="${BINARY_DIR}/swift_${dim}d"

    if [ ! -e "$binary_path" ]; then
        log "ERROR: no stashed binary for ${dim}D at ${binary_path}"
        return 1
    fi

    cp "$binary_path" "${SWIFT_ROOT}/swift"
    log "Actiavted ${dim}D binary as ${SWIFT_ROOT}/swift"
    return 0
}


#Runs a single examples run.sh for the given dimension
#By default reuses any initial conditions
#FRESH_IC=1 to force regeneration
run_example() {
    local dim="$1"
    local example_rel="${EXAMPLES[$dim]}"
    local example_dir="${SWIFT_ROOT}/${example_rel}"
    local example_name
    example_name="$(basename "$example_rel")"

    if [ ! -d "$example_dir" ]; then
        log "ERROR: example directory not found: ${example_dir}"
        return 1
    fi

    cd "$example_dir" || return 1

    if [ "${FRESH_IC:-0}" = "1" ]; then
        log "FRESH_IC=1 set, removing exisiting initial conditions for ${example_name}..."
        rm -f sodShock.hdf5
    fi

    log "Running ${example_name} (${dim}D).."
    if ! ./run.sh > run_automation.log 2>&1; then
        log "ERROR: run.sh faield for ${example_name}. Last 20 lines of log:"
        tail -n 20 run_automation.log
        return 1
    fi

    log "Run succeeded for ${example_name}."
    cd "$SWIFT_ROOT" || return 1

    return 0
}

activate_binary 2
run_example 2