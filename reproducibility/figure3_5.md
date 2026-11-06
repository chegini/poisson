cd /data/numerik/people/fchegini/project/fatemeh/MicroCard/current_papers/BDDC_petsc_SDC_2023/gatevariables/20Nov2023/BDDC/poisson_new_scratch

archive="archive/diagnostics_2026-10-24"
mkdir -p "$archive"

mv logs/* "$archive"/ 2>/dev/null || true
mv output/* "$archive"/ 2>/dev/null || true

Then recreate clean directories:

mkdir -p logs output reproducibility
------------------------------------------------------------------------------------------
Your next step is compilation.

cd /data/numerik/people/fchegini/project/fatemeh/MicroCard/current_papers/BDDC_petsc_SDC_2023/gatevariables/20Nov2023/BDDC/poisson_new_scratch

bash -n run_poisson_mpi_2nodes_8ranks.sbatch
bash -n reproducibility/figure3_quantization/submit_sweep.sh
bash -n reproducibility/figure5_subdomains/submit_sweep.sh
bash -n reproducibility/common/collect_metrics.sh

Then submit the build:

build_jid=$(sbatch --parsable \
  --partition=small \
  --nodelist=htc-cmp023 \
  build_poisson.sbatch)

echo "$build_jid"

Monitor it:

squeue -j "$build_jid"

Check the final status:

sacct -j "$build_jid" \
  --format=JobID,State,ExitCode,DerivedExitCode,Elapsed,NodeList

Continue only if the main job is:

COMPLETED  0:0

Then run the container executable check:

srun \
  --partition=big \
  --nodes=1 \
  --nodelist=htc-cmp011 \
  --ntasks=1 \
  --time=00:02:00 \
  bash -lc '
    podman run --rm --pull=never \
      -v "$PWD:/BDDC" \
      -v /data/numerik/people/fchegini/kaskade7_new:/kaskade7 \
      -v /home/htc/fchegini/opencarp_dependencies/petsc_install:/petsc:ro \
      -e LD_LIBRARY_PATH=/petsc/lib \
      -w /BDDC \
      localhost/kaskade7-mpi-build-environment:latest \
      ./bddc --help
  '

If the help text appears, run the small 8-rank validation:

validation_jid=$(sbatch --parsable \
  --partition=big \
  --nodes=1 \
  --nodelist=htc-cmp011 \
  --ntasks=8 \
  --ntasks-per-node=8 \
  --cpus-per-task=1 \
  --export=ALL,MPI_RANKS=8,THREADS=1,CPUS_PER_RANK=1,BITS=16,SUBRES=4,ELMRES=64,ITER=4,TOL=1e-10,COMPRESSION=0,GRAPH_LIFTING=0,GRAPH_LIFTING_AFTER_QUANTIZATION=0,ZIGZAG=0,HUFFMAN=0,BITLENGTH=0,HISTOGRAM=0,MIN_INTERFACE_SIZE=0,THRESHOLD_INTERFACE_TYPES=4,EXCLUDE_CORNERS=1,CASE_TAG=validation-small \
  run_poisson_mpi_2nodes_8ranks.sbatch)

echo "$validation_jid"

sacct -j "$validation_jid" \
  --format=JobID,State,ExitCode,DerivedExitCode,Elapsed,NodeList

Only after this validation succeeds, create a clean campaign directory and submit Figure 3 and Figure 5:

export SCRATCH_DIR="$PWD"
export REPORT_REPRO_DIR="$SCRATCH_DIR/reproducibility/run_$(date +%Y%m%d_%H%M%S)"
export RUN_NODELIST=htc-cmp011,htc-cmp013
export BUILD_JID="$build_jid"

mkdir -p "$REPORT_REPRO_DIR"

bash reproducibility/figure3_quantization/submit_sweep.sh
bash reproducibility/figure5_subdomains/submit_sweep.sh

Use:

squeue -u "$USER"

------------------------------------------------------------------------------------------
Below is the complete reproducibility workflow. It uses:

8 MPI ranks
1 CPU per MPI rank
1 OpenMP thread per rank
2 nodes
4 ranks per node

This means one CPU per running thread/MPI rank. It does not allocate one CPU per subdomain; each rank may process multiple subdomains.

Before running, copy these updated local files to the cluster:

run_poisson_mpi_2nodes_8ranks.sbatch
reproducibility/figure3_quantization/submit_sweep.sh
reproducibility/figure5_subdomains/submit_sweep.sh
reproducibility/common/collect_metrics.sh

chmod +x run_poisson_mpi_2nodes_8ranks.sbatch
chmod +x reproducibility/figure3_quantization/submit_sweep.sh
chmod +x reproducibility/figure5_subdomains/submit_sweep.sh
chmod +x reproducibility/common/collect_metrics.sh

The launcher now:

--cpus-per-task=1
THREADS=1
OMP_NUM_THREADS=1

and mounts the correct MPI installation:

-v /home/htc/fchegini/opencarp_dependencies/petsc_install:/petsc:ro
-e LD_LIBRARY_PATH=/petsc/lib

## A. Enter the cluster project

cd /data/numerik/people/fchegini/project/fatemeh/MicroCard/current_papers/BDDC_petsc_SDC_2023/gatevariables/20Nov2023/BDDC/poisson_new_scratch

Check the launcher:

grep -nE 'cpus-per-task|THREADS|petsc|LD_LIBRARY_PATH|PMI_RANK|PMI_SIZE|OMP_NUM_THREADS' \
  run_poisson_mpi_2nodes_8ranks.sbatch

## B. Compile the executable

build_jid=$(sbatch --parsable \
  --partition=small \
  --nodelist=htc-cmp023 \
  build_poisson.sbatch)

echo "$build_jid"

Monitor:

squeue -j "$build_jid"

Check the result:

sacct -j "$build_jid" \
  --format=JobID,State,ExitCode,DerivedExitCode,Elapsed,NodeList

Continue only if the main build job is COMPLETED with exit code 0:0.

## C. Check the executable inside the container

srun \
  --partition=big \
  --nodes=1 \
  --nodelist=htc-cmp011 \
  --ntasks=1 \
  --time=00:02:00 \
  bash -lc '
    podman run --rm --pull=never \
      -v "$PWD:/BDDC" \
      -v /data/numerik/people/fchegini/kaskade7_new:/kaskade7 \
      -v /home/htc/fchegini/opencarp_dependencies/petsc_install:/petsc:ro \
      -e LD_LIBRARY_PATH=/petsc/lib \
      -w /BDDC \
      localhost/kaskade7-mpi-build-environment:latest \
      ./bddc --help
  '

The command must print the bddc help text.

## D. Small MPI validation

validation_jid=$(sbatch --parsable \
  --partition=big \
  --nodes=1 \
  --nodelist=htc-cmp011 \
  --ntasks=8 \
  --ntasks-per-node=8 \
  --cpus-per-task=1 \
  --export=ALL,MPI_RANKS=8,THREADS=1,CPUS_PER_RANK=1,BITS=16,SUBRES=4,ELMRES=64,ITER=4,TOL=1e-10,COMPRESSION=0,GRAPH_LIFTING=0,GRAPH_LIFTING_AFTER_QUANTIZATION=0,ZIGZAG=0,HUFFMAN=0,BITLENGTH=0,HISTOGRAM=0,MIN_INTERFACE_SIZE=0,THRESHOLD_INTERFACE_TYPES=4,EXCLUDE_CORNERS=1,CASE_TAG=validation-small \
  run_poisson_mpi_2nodes_8ranks.sbatch)

echo "$validation_jid"

Check:

sacct -j "$validation_jid" \
  --format=JobID,State,ExitCode,DerivedExitCode,Elapsed,NodeList

This should complete successfully before submitting the full sweeps.

## E. Prepare a clean reproducibility directory

Use a new directory for each complete campaign:

export SCRATCH_DIR=/data/numerik/people/fchegini/project/fatemeh/MicroCard/current_papers/BDDC_petsc_SDC_2023/gatevariables/20Nov2023/BDDC/poisson_new_scratch

export REPORT_REPRO_DIR="$SCRATCH_DIR/reproducibility/run_$(date +%Y%m%d_%H%M%S)"

mkdir -p "$REPORT_REPRO_DIR"

Set the known working nodes:

export RUN_NODELIST=htc-cmp011,htc-cmp013
export BUILD_JID="$build_jid"

## F. Submit Figure 3

Figure 3 uses:

SUBRES=16
ELMRES=32
256 subdomains
512 × 512 global resolution
8 MPI ranks
1 CPU per rank
1 thread per rank
quantization bits: 4, 5, 6, 7, 8, 16, 24, 27, 32, 64

Submit:

bash reproducibility/figure3_quantization/submit_sweep.sh

The manifest will be written to:

$REPORT_REPRO_DIR/figure3_quantization/run_manifest.csv

## G. Submit Figure 5

Figure 5 keeps the global resolution fixed at 512 × 512:

| Subdomains | SUBRES | ELMRES |
|---:|---:|---:|
| 16 | 4 | 128 |
| 64 | 8 | 64 |
| 256 | 16 | 32 |
| 1024 | 32 | 16 |
| 4096 | 64 | 8 |

Submit:

bash reproducibility/figure5_subdomains/submit_sweep.sh

This submits 40 jobs:

5 subdomain configurations × 8 bit configurations

## H. Monitor all jobs

squeue -u "$USER"

Do not use:

squeue -u "$JOB_ID"

For one job:

squeue -j "$JOB_ID"

For all Figure 3 jobs:

awk -F, 'NR > 1 {print $2}' \
  "$REPORT_REPRO_DIR/figure3_quantization/run_manifest.csv" |
while read -r jid; do
  sacct -n -X -j "$jid" \
    --format=JobID,State,ExitCode,Elapsed,NodeList
done

For all Figure 5 jobs:

awk -F, 'NR > 1 {print $2}' \
  "$REPORT_REPRO_DIR/figure5_subdomains/run_manifest.csv" |
while read -r jid; do
  sacct -n -X -j "$jid" \
    --format=JobID,State,ExitCode,Elapsed,NodeList
done

## I. Collect compression metrics

Figure 3:

mapfile -t jobs < <(
  awk -F, 'NR > 1 && $2 ~ /^[0-9]+$/ {print $2}' \
    "$REPORT_REPRO_DIR/figure3_quantization/run_manifest.csv"
)

bash reproducibility/common/collect_metrics.sh \
  "$REPORT_REPRO_DIR/results/figure3_metrics.csv" \
  "${jobs[@]}"

Figure 5:

mapfile -t jobs < <(
  awk -F, 'NR > 1 && $2 ~ /^[0-9]+$/ {print $2}' \
    "$REPORT_REPRO_DIR/figure5_subdomains/run_manifest.csv"
)

bash reproducibility/common/collect_metrics.sh \
  "$REPORT_REPRO_DIR/results/figure5_metrics.csv" \
  "${jobs[@]}"

Inspect the results:

column -s, -t "$REPORT_REPRO_DIR/results/figure3_metrics.csv"
column -s, -t "$REPORT_REPRO_DIR/results/figure5_metrics.csv"

A valid result row should contain:

COMPLETED
elapsed time
raw_bytes
transmitted_bytes
savings_percent
ratio

Rows with empty metrics indicate failed or incomplete jobs and must not be used in the figures.

## J. Preserve the complete record

Keep these files together:

figure3_quantization/run_manifest.csv
figure5_subdomains/run_manifest.csv
results/figure3_metrics.csv
results/figure5_metrics.csv
logs/
run_poisson_mpi_2nodes_8ranks.sbatch
figure3_quantization/submit_sweep.sh
figure5_subdomains/submit_sweep.sh
common/collect_metrics.sh

Also record the build job ID:

echo "$build_jid" > "$REPORT_REPRO_DIR/build_job_id.txt"

The final report should use only completed jobs from the manifests and corresponding metrics CSV files.