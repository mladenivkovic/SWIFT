#! /bin/bash

# Optional variables:
# FORCE_REBUILD=1 rebuild binaries if already stashed
# FRESH_IC=1 delete existing inital conditions before each run, forcing regeneration

set -u # Treat unset variables as errors
set -o pipefail # tee has different exit codes to bash

# Get the swift directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SWIFT_ROOT="$(dirname "$SCRIPT_DIR")"

AUTOMATION_DIR="${SWIFT_ROOT}/HydroAutomation"
BUILD_LOG_DIR="${AUTOMATION_DIR}/build_logs"
BINARY_DIR="${AUTOMATION_DIR}/binaries"
RESULTS_DIR="${AUTOMATION_DIR}/results"
MAKE_JOBS=32

# Add the schemes for compilation here
SCHEMES=(sphenix)

# example name -> relative path. Dimension is inferred from the folder's
# own _1D/_2D/_3D suffix rather than tracked separately, so adding a new
# example later just means adding one line here.
declare -A EXAMPLES=(
    # [GreshoVortex_2D]="examples/HydroTests/GreshoVortex_2D"
    # [GreshoVortex_3D]="examples/HydroTests/GreshoVortex_3D"
    # [InteractingBlastWaves_1D]="examples/HydroTests/InteractingBlastWaves_1D"
    # [KelvinHelmholtzGrowthRate_2D]="examples/HydroTests/KelvinHelmholtzGrowthRate_2D"
    # [KelvinHelmholtzGrowthRate_3D]="examples/HydroTests/KelvinHelmholtzGrowthRate_3D"
    # [KelvinHelmholtz_2D]="examples/HydroTests/KelvinHelmholtz_2D"
    # [Noh_1D]="examples/HydroTests/Noh_1D"
    # [Noh_2D]="examples/HydroTests/Noh_2D"
    # [Noh_3D]="examples/HydroTests/Noh_3D"
    # [SedovBlast_1D]="examples/HydroTests/SedovBlast_1D"
    # [SedovBlast_2D]="examples/HydroTests/SedovBlast_2D"
    # [SedovBlast_3D]="examples/HydroTests/SedovBlast_3D"
    # [SodShockSpherical_2D]="examples/HydroTests/SodShockSpherical_2D"
    # [SodShockSpherical_3D]="examples/HydroTests/SodShockSpherical_3D"
    [SodShock_1D]="examples/HydroTests/SodShock_1D"
    [SodShock_2D]="examples/HydroTests/SodShock_2D"
    [SodShock_3D]="examples/HydroTests/SodShock_3D"
    # [SodShock_BCC_3D]="examples/HydroTests/SodShock_BCC_3D"
    # [SquareTest_2D]="examples/HydroTests/SquareTest_2D"
    # [ToroTest2_1D]="examples/HydroTests/ToroTest2_1D"
    # [ToroTest2_2D]="examples/HydroTests/ToroTest2_2D"
    # [ToroTest2_3D]="examples/HydroTests/ToroTest2_3D"
    # [VacuumSpherical_2D]="examples/HydroTests/VacuumSpherical_2D"
    # [VacuumSpherical_3D]="examples/HydroTests/VacuumSpherical_3D"
    # [Vacuum_1D]="examples/HydroTests/Vacuum_1D"
    # [Vacuum_2D]="examples/HydroTests/Vacuum_2D"
    # [Vacuum_3D]="examples/HydroTests/Vacuum_3D"
)

# Infers hydro dimension (1/2/3) from an example's own folder-name suffix.
get_example_dimension() {
    local name="$1"
    case "$name" in
        *_1D) echo 1 ;;
        *_2D) echo 2 ;;
        *_3D) echo 3 ;;
        *) echo "" ;;
    esac
}

# Needed directories
mkdir -p "$BUILD_LOG_DIR" "$BINARY_DIR" "$RESULTS_DIR" 

#logs a time based message
log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*"
}

# Configures and compiles SWIFT for each dimension
# Stashes the resulting binaries as binaries/swift_<scheme>_<dim>d
# Skips rebuild if the binary already exists unless FORCE_REBUILD=1
build_binary() {
    local dim="$1"
    local scheme="$2"
    local binary_path="${BINARY_DIR}/swift_${scheme}_${dim}d"

    if [ -e "$binary_path" ] && [ "${FORCE_REBUILD:-0}" != "1" ]; then
        log "Binary for ${scheme}/${dim}D already exists at ${binary_path}, skipping build."
        return 0
    fi

    cd "$SWIFT_ROOT" || return 1


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

# Copies the stashed binary for a given dimension and scheme into place as ./swift
# since each run.sh referes to ./swift as as a relative path
activate_binary() {
    local dim="$1"
    local scheme="$2"
    local binary_path="${BINARY_DIR}/swift_${scheme}_${dim}d"

    if [ ! -e "$binary_path" ]; then
        log "ERROR: no stashed binary for ${scheme}/${dim}D at ${binary_path}"
        return 1
    fi

    cp "$binary_path" "${SWIFT_ROOT}/swift"
    sync
    local src_hash dst_hash
    src_hash=$(md5sum "$binary_path" | awk '{print $1}')
    dst_hash=$(md5sum "${SWIFT_ROOT}/swift" | awk '{print $1}')
    if [ "$src_hash" != "$dst_hash" ]; then
        log "ERROR: copied binary hash mismatch for ${scheme}/${dim}D (expected ${src_hash}, got ${dst_hash})"
        return 1
    fi

    log "Activated ${scheme}/${dim}D binary as ${SWIFT_ROOT}/swift (verified ${dst_hash})"
    return 0
}


# Runs a single example's run.sh
# By default reuses any initial conditions FRESH_IC=1 to force regeneration
# run marker is used to track the files outputted by the example
run_example() {
    local example_dir="$1"
    local example_name="$2"

    if [ ! -d "$example_dir" ]; then
        log "ERROR: example directory not found: ${example_dir}"
        return 1
    fi

    cd "$example_dir" || return 1

    if [ "${FRESH_IC:-0}" = "1" ]; then
        log "FRESH_IC=1 set, removing existing initial conditions for ${example_name}..."
        find . -maxdepth 1 -name "*.hdf5" ! -name "glass*" -exec rm -f {} \;
    fi

    touch .run_marker

    log "Running ${example_name}..."
    if ! ./run.sh 2>&1 | tee run_automation.log; then
        log "ERROR: run.sh failed for ${example_name}. Last 20 lines of log:"
        tail -n 20 run_automation.log
        cd "$SWIFT_ROOT" || return 1
        return 1
    fi

    log "Run succeeded for ${example_name}."
    cd "$SWIFT_ROOT" || return 1
    return 0
}

# Moves whatever run-output files remain after cleanup (i.e. anything newer
# than the .run_marker dropped at the start of run_example, and not deleted
# as junk) into a timestamped results folder. Works for any example,
# regardless of what it names its plot/log/parameter files.
collect_outputs() {
    local example_dir="$1"
    local dim="$2"
    local scheme="$3"
    local example_name="$4"
    local ts
    ts="$(date +%Y%m%d_%H%M%S)"
    local out_dir="${RESULTS_DIR}/${dim}D/${scheme}/${example_name}_${ts}"

    cd "$example_dir" || return 1

    if [ ! -e .run_marker ]; then
        log "WARNING: no .run_marker found for ${example_name}, skipping collection"
        return 1
    fi

    mkdir -p "$out_dir"
    log "Collecting outputs for ${example_name} into ${out_dir}"

    local moved_any=0
    while IFS= read -r -d '' f; do
        mv "$f" "$out_dir"/
        moved_any=1
    done < <(find . -maxdepth 1 -type f -newer .run_marker ! -name ".run_marker" -print0)

    if [ "$moved_any" -eq 0 ]; then
        log "  (note: nothing newer than the run marker was found to collect)"
    fi

    rm -f .run_marker
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
    # Prepare the working directory
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
# For every scheme and for each dimension
# 1. Build/reuse swift binary
# 2. Activate it as ./swift
# 3. Run every example whose suffix matches the dimension
# 4. On success clean and collect the example outputs
# A per combination row is added to SUMMARY_FILE
# Any failures are collected into FAILED_COMBOS for final report
# ====================
SUMMARY_DIR="${RESULTS_DIR}/summaries"
mkdir -p "$SUMMARY_DIR"
SUMMARY_FILE="${SUMMARY_DIR}/run_summary_$(date +%Y%m%d_%H%M%S).csv"
echo "scheme,dimension,example,build,activate,run,status" > "$SUMMARY_FILE"

FAILED_COMBOS=()

for scheme in "${SCHEMES[@]}"; do
    for dim in 1 2 3; do
        build_status="ok"
        if ! build_binary "$dim" "$scheme"; then
            build_status="FAILED"
        else
            if ! activate_binary "$dim" "$scheme"; then
                build_status="FAILED (activate)"
            fi
        fi

        for example_name in "${!EXAMPLES[@]}"; do
            example_dim="$(get_example_dimension "$example_name")"
            [ "$example_dim" = "$dim" ] || continue

            example_dir="${SWIFT_ROOT}/${EXAMPLES[$example_name]}"
            log "=== Processing ${scheme} / ${example_name} ==="

            if [ "$build_status" != "ok" ]; then
                FAILED_COMBOS+=("${scheme}/${example_name} (build/activate)")
                echo "${scheme},${dim},${example_name},${build_status},skipped,skipped,FAILED" >> "$SUMMARY_FILE"
                continue
            fi

            if ! run_example "$example_dir" "${example_name}_${scheme}"; then
                FAILED_COMBOS+=("${scheme}/${example_name} (run)")
                echo "${scheme},${dim},${example_name},ok,ok,FAILED,FAILED" >> "$SUMMARY_FILE"
                continue
            fi

            cleanup_example "$example_dir" "${example_name}_${scheme}"
            collect_outputs "$example_dir" "$dim" "$scheme" "$example_name"
            echo "${scheme},${dim},${example_name},ok,ok,ok,OK" >> "$SUMMARY_FILE"
        done
    done
done

log "=== Done ==="
log "Summary written to ${SUMMARY_FILE}"
if [ "${#FAILED_COMBOS[@]}" -eq 0 ]; then
    log "All scheme/example combinations completed successfully"
else
    log "The following combinations failed: ${FAILED_COMBOS[*]}"
    exit 1
fi