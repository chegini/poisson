#!/usr/bin/env bash
set -euo pipefail

SCRATCH_DIR=${SCRATCH_DIR:-/data/numerik/people/fchegini/project/fatemeh/MicroCard/current_papers/BDDC_petsc_SDC_2023/gatevariables/20Nov2023/BDDC/poisson_new_scratch}
REPORT_REPRO_DIR=${REPORT_REPRO_DIR:-$SCRATCH_DIR/reproducibility}
OUT=${1:-metrics.csv}
shift || true

printf 'job_id,state,elapsed,minimum_interface_size,threshold_interface_types,raw_bytes,transmitted_bytes,savings_percent,ratio\n' > "$OUT"

for jid in "$@"; do
  mapfile -t log_files < <(
    find "$SCRATCH_DIR/logs" "$REPORT_REPRO_DIR/logs" \
      -maxdepth 1 -type f -name "*${jid}*" 2>/dev/null
  )
  if [[ ${#log_files[@]} -eq 0 ]]; then
    printf '%s,,,,,,,,\n' "$jid" >> "$OUT"
    continue
  fi

  log_text=$(cat "${log_files[@]}")
  state=$(sacct -n -X -j "$jid" --format=State | awk 'NR==1{print $1}')
  elapsed=$(sacct -n -X -j "$jid" --format=Elapsed | awk 'NR==1{print $1}')
  min_size=$(printf '%s\n' "$log_text" | grep -m1 'Minimum interface size:' | sed -E 's/.*Minimum interface size: *//')
  mask=$(printf '%s\n' "$log_text" | grep -m1 'Threshold interface types:' | sed -E 's/.*Threshold interface types: *//')
  raw=$(printf '%s\n' "$log_text" | grep -m1 'raw sent bytes:' | sed -E 's/.*raw sent bytes: *//')
  transmitted=$(printf '%s\n' "$log_text" | grep -m1 'transmitted sent bytes:' | sed -E 's/.*transmitted sent bytes: *//')
  savings=$(printf '%s\n' "$log_text" | grep -m1 'savings:' | sed -E 's/.*savings: *//; s/%.*//')
  ratio=$(printf '%s\n' "$log_text" | grep -m1 'ratio:' | sed -E 's/.*ratio: *//')
  printf '%s,%s,%s,%s,%s,%s,%s,%s,%s\n' "$jid" "$state" "$elapsed" "$min_size" "$mask" "$raw" "$transmitted" "$savings" "$ratio" >> "$OUT"
done

