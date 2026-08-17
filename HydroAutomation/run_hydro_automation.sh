#! /bin/bash

set -u #Treat unset variables as errors
set -o pipefail # tee has different exit codes to bash


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

KEEP_FILES=(SodShock.png output.log used_parameters.yml unused_parameters.yml)

SCHEMES=(sphenix minimal)

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
    local scheme="$2"
    local binary_path="${BINARY_DIR}/swift_${scheme}_${dim}d"

    if [ -e "$binary_path" ] && [ "${FORCE_REBUILD:-0}" != "1" ]; then
        log "Binary for ${scheme}/${dim}D already exists at ${binary_path}, skipping build."
        return 0
    fi

    cd "$SWIFT_ROOT" || return 1
    load_modules

    log "Configuring SWIFT for hydro-dimension=${dim}, hydro=${scheme}..."
    if ! ./configure --with-hydro-dimension="${dim}" --with-hydro="${scheme}"> "${BUILD_LOG_DIR}/configure_${scheme}_${dim}d.log" 2>&1; then
        log "Configuration failed for ${scheme}/${dim}D. Check ${BUILD_LOG_DIR}/configure_${scheme}_${dim}d.log for details."
        return 1
    fi

    log "Building SWIFT (${scheme}/${dim}D) with make -j ${MAKE_JOBS}... (this takes a while)"
    if ! make -j"${MAKE_JOBS}" > "${BUILD_LOG_DIR}/make_${scheme}_${dim}d.log" 2>&1; then
        log "ERROR: make failed for ${scheme}/${dim}D. See ${BUILD_LOG_DIR}/make_${scheme}_${dim}d.log"
        return 1
    fi

    cp "${SWIFT_ROOT}/swift" "$binary_path"
    log "Stashed ${scheme}/${dim}D binary at ${binary_path}"
    return 0
    
}

#Copies the stashed binary for a given dimension into place as ./swift
#since each run.sh referes to ./swift as as a relative path
activate_binary() {
    local dim="$1"
    local scheme="$2"
    local binary_path="${BINARY_DIR}/swift_${scheme}_${dim}d"

    if [ ! -e "$binary_path" ]; then
        log "ERROR: no stashed binary for ${scheme}/${dim}D at ${binary_path}"
        return 1
    fi

    cp "$binary_path" "${SWIFT_ROOT}/swift"
    log "Actiavted ${scheme}/${dim}D binary as ${SWIFT_ROOT}/swift"
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
    if ! ./run.sh 2>&1 | tee run_automation.log; then
        log "ERROR: run.sh faield for ${example_name}. Last 20 lines of log:"
        tail -n 20 run_automation.log
        return 1
    fi

    log "Run succeeded for ${example_name}."
    cd "$SWIFT_ROOT" || return 1

    return 0
}

#Copies the keep files for one example into a folder under results
collect_outputs() {
    local example_dir="$1"
    local example_name="$2"
    local ts
    ts="$(date +%Y%m%d_%H%M%S)"
    local out_dir="${RESULTS_DIR}/${example_name}_${ts}"

    cd "$example_dir" || return 1
    mkdir -p "$out_dir"
    log "Collecting outputs for ${example_name} into ${out_dir}"

    for f in "${KEEP_FILES[@]}"; do
        if [ -e "$f" ]; then
            cp -f "$f" "$out_dir"/
        else
            log "Expected file '${f}' Not found, skipping"
        fi
    done

    cp -f run_automation.log "$out_dir"/
    cd "$SWIFT_ROOT" || return 1
    return 0

}

# Extract a value from a top-level YAML section (e.g. yaml_get file.yml Snapshots basename)
yaml_get() {
    local file="$1" section="$2" key="$3"
    awk -v section="$section" -v key="$key" '
        $0 ~ "^"section":[ \t]*$" { in_section=1; next }
        in_section && /^[^ \t#]/ { in_section=0 }
        in_section && $0 ~ "^[ \t]+"key":" {
            line=$0
            sub("^[ \t]+"key":[ \t]*", "", line)
            sub(/[ \t]*#.*/, "", line)
            gsub(/^[ \t"'"'"']+|[ \t"'"'"']+$/, "", line)
            print line; exit
        }
    ' "$file"
}

cleanup_example() {
    local example_dir="$1" example_name="$2"
    cd "$example_dir" || return 1
    log "Cleaning up run generated files for ${example_name}"

    # Names that never vary between examples
    local -a static_junk=("dependency_graph_*.csv" "task_level_*.txt" "timesteps.txt" "*.xmf")
    for pattern in "${static_junk[@]}"; do
        find . -maxdepth 2 -name "$pattern" -exec rm -f {} \;
    done

    # Snapshot/statistics/restart names come from each yml's own basename —
    # loop over all .yml files, since a dir can hold more than one (res variants etc.)
    local yml
    for yml in *.yml; do
        [ -f "$yml" ] || continue
        local snap_base snap_subdir stats_base restart_subdir
        snap_base=$(yaml_get "$yml" "Snapshots" "basename")
        snap_subdir=$(yaml_get "$yml" "Snapshots" "subdir")
        stats_base=$(yaml_get "$yml" "Statistics" "basename")
        restart_subdir=$(yaml_get "$yml" "Restarts" "subdir")
        stats_base="${stats_base:-statistics.txt}"
        restart_subdir="${restart_subdir:-restart}"

        local snap_dir="."; [ -n "$snap_subdir" ] && snap_dir="$snap_subdir"
        [ -n "$snap_base" ] && find "$snap_dir" -maxdepth 1 -name "${snap_base}_[0-9]*.hdf5" -exec rm -f {} \;
        [ -f "$stats_base" ] && rm -f "$stats_base"
        [ -d "$restart_subdir" ] && rm -rf "$restart_subdir"
    done

    rm -f run_automation.log
    log "Directory listing for ${example_name} after cleanup:"
    ls -la
    cd "$SWIFT_ROOT" || return 1
}

# =======MAIN========

SUMMARY_FILE="${AUTOMATION_DIR}/run_summary_$(date +%Y%m%d_%H%M%S).csv"
echo "dimension,scheme,build,activate,run,status" > "$SUMMARY_FILE"

FAILED_COMBOS=()

for dim in 2 3; do
    for scheme in "${SCHEMES[@]}"; do
        log "=== Processing ${scheme} / ${dim}D ==="
        example_dir="${SWIFT_ROOT}/${EXAMPLES[$dim]}"
        example_name="$(basename "${EXAMPLES[$dim]}")"

        build_status="ok"
        activate_status="skipped"
        run_status="skipped"

        if ! build_binary "$dim" "$scheme"; then
            build_status="FAILED"
            FAILED_COMBOS+=("${scheme}/${dim}D (build)")
            echo "${dim},${scheme},${build_status},${activate_status},${run_status},FAILED" >> "$SUMMARY_FILE"
            continue
        fi

        if ! activate_binary "$dim" "$scheme"; then
            activate_status="FAILED"
            FAILED_COMBOS+=("${scheme}/${dim}D (activate)")
            echo "${dim},${scheme},${build_status},${activate_status},${run_status},FAILED" >> "$SUMMARY_FILE"
            continue
        fi
        activate_status="ok"

        if ! run_example "$dim"; then
            run_status="FAILED"
            FAILED_COMBOS+=("${scheme}/${dim}D (run)")
            echo "${dim},${scheme},${build_status},${activate_status},${run_status},FAILED" >> "$SUMMARY_FILE"
            continue
        fi
        run_status="ok"
        
        collect_outputs "$example_dir" "${example_name}_${scheme}"
        cleanup_example "$example_dir" "${example_name}_${scheme}"

        echo "${dim},${scheme},${build_status},${activate_status},${run_status},OK" >> "$SUMMARY_FILE"

    done
done

log "=== Done ==="
log "Summary written to ${SUMMARY_FILE}"
if [ "${#FAILED_COMBOS[@]}" -eq 0 ]; then
    log "All scheme/dimension combinations completed successfully"
else
    log "The following combinations failed: ${FAILED_COMBOS[*]}"
    exit 1
fi
