#!/usr/bin/env bash
set -euo pipefail

SCRATCH_DIR=${SCRATCH_DIR:-/data/numerik/people/fchegini/project/fatemeh/MicroCard/current_papers/BDDC_petsc_SDC_2023/gatevariables/20Nov2023/BDDC/poisson_new_scratch}
REPORT_REPRO_DIR=${REPORT_REPRO_DIR:-$(cd "$(dirname "$0")/.." && pwd)}
MANIFEST="$REPORT_REPRO_DIR/figure3_quantization/run_manifest.csv"
RUN_NODELIST=${RUN_NODELIST:-htc-cmp007,htc-cmp008}
NODE_ARGS=()
if [[ -n "$RUN_NODELIST" ]]; then
  NODE_ARGS=(--nodelist="$RUN_NODELIST")
fi
DEPENDENCY_ARGS=()
if [[ -n "${BUILD_JID:-}" ]]; then
  DEPENDENCY_ARGS=(--dependency="afterok:${BUILD_JID}")
fi
mkdir -p "$REPORT_REPRO_DIR/figure3_quantization" "$REPORT_REPRO_DIR/logs" "$REPORT_REPRO_DIR/results"

if [[ ! -f "$SCRATCH_DIR/run_poisson_mpi_2nodes_8ranks.sbatch" ]]; then
  echo "Missing launcher: $SCRATCH_DIR/run_poisson_mpi_2nodes_8ranks.sbatch" >&2
  exit 2
fi

if [[ ! -f "$MANIFEST" ]]; then
  printf 'submitted_at,job_id,bits,subdomains,global_elements,mpi_ranks,threads,lifting_mode,zigzag,huffman,bitlength,min_interface_size,threshold_interface_types,exclude_corners,command\n' > "$MANIFEST"
fi

for bits in 4 5 6 7 8 16 24 27 32 64; do
  export_string="ALL,BITS=$bits,SUBRES=16,ELMRES=32,MPI_RANKS=8,THREADS=1,CPUS_PER_RANK=1,ITER=20,TOL=1e-10,COMPRESSION=1,GRAPH_LIFTING=0,GRAPH_LIFTING_AFTER_QUANTIZATION=0,ZIGZAG=0,HUFFMAN=0,BITLENGTH=0,HISTOGRAM=0,MIN_INTERFACE_SIZE=0,THRESHOLD_INTERFACE_TYPES=4,EXCLUDE_CORNERS=1,CASE_TAG=figure3-bits-${bits}"
  job_id=$(sbatch --parsable "${DEPENDENCY_ARGS[@]}" --partition=big "${NODE_ARGS[@]}" --nodes=2 --ntasks=8 --cpus-per-task=1 \
    --export="$export_string" \
    --job-name="poisson-f3-b${bits}" \
    --output="$REPORT_REPRO_DIR/logs/figure3-b${bits}-%j.out" \
    --error="$REPORT_REPRO_DIR/logs/figure3-b${bits}-%j.err" \
    "$SCRATCH_DIR/run_poisson_mpi_2nodes_8ranks.sbatch" \
    )
  printf '%s,%s,%s,256,262144,8,1,0,0,0,0,0,4,1,%s\n' \
    "$(date -Iseconds)" "$job_id" "$bits" \
    "BITS=$bits SUBRES=16 ELMRES=32 MPI_RANKS=8 THREADS=1 CPUS_PER_RANK=1 ITER=20 TOL=1e-10 COMPRESSION=1 LIFTING_MODE=0 ZIGZAG=0 HUFFMAN=0 BITLENGTH=0 HISTOGRAM=0 MIN_INTERFACE_SIZE=0 THRESHOLD_INTERFACE_TYPES=4 EXCLUDE_CORNERS=1" >> "$MANIFEST"
  echo "Figure 3 bits=$bits job=$job_id"
done

