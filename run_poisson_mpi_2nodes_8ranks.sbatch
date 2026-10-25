#!/bin/bash
#SBATCH --job-name=poisson-2n8r
#SBATCH --partition=big
#SBATCH --nodes=2
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --output=logs/poisson-2nodes-8ranks-%j.out
#SBATCH --error=logs/poisson-2nodes-8ranks-%j.err

# Multi-node launch: Slurm starts one container per MPI task. Do not call
# run_case_mpi.sh here because it would start a nested mpirun on one node.
# The application is built against the cluster MPICH installation. Use PMI2,
# which is the Slurm integration supported by this MPICH Hydra runtime.
set -euo pipefail

BASE=/data/numerik/people/fchegini/project/fatemeh/MicroCard/current_papers/BDDC_petsc_SDC_2023/gatevariables/20Nov2023/BDDC
SCRATCH_DIR="$BASE/poisson_new_scratch"
KASKADE7=/data/numerik/people/fchegini/kaskade7_new
IMAGE=localhost/kaskade7-mpi-build-environment:latest
ARCHIVE=/data/numerik/people/fchegini/kaskade7-mpi-build-environment.tar
MPI_IMPL=${MPI_IMPL:-pmi2}
HOST_MPI=/home/htc/fchegini/opencarp_dependencies/petsc_install
MPI_RANKS=${MPI_RANKS:-${SLURM_NTASKS}}
THREADS=${THREADS:-1}
CPUS_PER_RANK=${CPUS_PER_RANK:-1}
BITS=${BITS:-22}
MIN_INTERFACE_SIZE=${MIN_INTERFACE_SIZE:-0}
THRESHOLD_INTERFACE_TYPES=${THRESHOLD_INTERFACE_TYPES:-4}
EXCLUDE_CORNERS=${EXCLUDE_CORNERS:-1}
SUBRES=${SUBRES:-4}
ELMRES=${ELMRES:-64}
ITER=${ITER:-20}
TOL=${TOL:-1e-10}
COMPRESSION=${COMPRESSION:-1}
GRAPH_LIFTING=${GRAPH_LIFTING:-0}
GRAPH_LIFTING_AFTER_QUANTIZATION=${GRAPH_LIFTING_AFTER_QUANTIZATION:-1}
HUFFMAN=${HUFFMAN:-1}
BITLENGTH=${BITLENGTH:-1}
ZIGZAG=${ZIGZAG:-1}
HISTOGRAM=${HISTOGRAM:-0}
INTERFACE_TYPES=${INTERFACE_TYPES:-7}
CRISSCROSS=${CRISSCROSS:-1}
CASE_TAG=${CASE_TAG:-multi-node}
RUN_NAME=${RUN_NAME:-multi-node-ranks-${MPI_RANKS}-bits-${BITS}-job-${SLURM_JOB_ID}}
RUN="output/multi-node/${RUN_NAME}"
export BASE SCRATCH_DIR KASKADE7 IMAGE ARCHIVE HOST_MPI RUN BITS MIN_INTERFACE_SIZE \
  THRESHOLD_INTERFACE_TYPES EXCLUDE_CORNERS THREADS CPUS_PER_RANK \
  SUBRES ELMRES ITER TOL COMPRESSION GRAPH_LIFTING GRAPH_LIFTING_AFTER_QUANTIZATION \
  HUFFMAN BITLENGTH ZIGZAG HISTOGRAM INTERFACE_TYPES CRISSCROSS CASE_TAG

cd "$SCRATCH_DIR"
HOST_WORKDIR="$PWD"
export HOST_WORKDIR
mkdir -p logs "$RUN"

echo "Job ID: $SLURM_JOB_ID"
echo "Nodes: $(scontrol show hostnames "$SLURM_JOB_NODELIST" | tr "\n" " ")"
echo "Tasks: $SLURM_NTASKS"
echo "MPI mode: $MPI_IMPL"
echo "Ranks per node: ${SLURM_NTASKS_PER_NODE:-unknown}"
echo "Bits: $BITS"
echo "Minimum interface size: $MIN_INTERFACE_SIZE"
echo "Threshold interface types: $THRESHOLD_INTERFACE_TYPES"
echo "CPU binding: cores"
echo "CPUs per rank: $CPUS_PER_RANK"

# Keep one stdout/stderr file per MPI task. Otherwise Slurm merges output from
# all ranks, making concurrent solver and compression messages unreadable.
srun --mpi="$MPI_IMPL" --ntasks="$SLURM_NTASKS" \
  --ntasks-per-node="$SLURM_NTASKS_PER_NODE" \
  --cpus-per-task="$CPUS_PER_RANK" \
  --time=04:00:00 \
  --kill-on-bad-exit=0 \
  --output="logs/poisson-multi-%j-rank%t.out" \
  --error="logs/poisson-multi-%j-rank%t.err" \
  bash -lc '
    : "${PMI_FD:?PMI_FD was not provided by Slurm}"
    eval "exec 3<&${PMI_FD}"
    podman run --rm --pull=never \
      --network=host \
      --ipc=host \
      --pid=host \
      --preserve-fds=1 \
      --env-host \
      -e PMI_FD=3 \
      -e PMI_RANK="$PMI_RANK" \
      -e PMI_SIZE="$PMI_SIZE" \
      -v "$HOST_MPI:/petsc:ro" \
      -e PATH=/petsc/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
      -e LD_LIBRARY_PATH=/kaskade7/libs:/petsc/lib:/petsc/lib64:/usr/local/lib64:/usr/local/lib \
      -e OMP_NUM_THREADS="$THREADS" \
      -v "$KASKADE7:/kaskade7" \
      -v "$HOST_WORKDIR:/BDDC" \
      -w /BDDC \
      "$IMAGE" \
      ./bddc \
        --mpi 1 \
        --compression "$COMPRESSION" \
        --compressionBits "$BITS" \
        --compressionMinInterfaceSize "$MIN_INTERFACE_SIZE" \
        --compressionThresholdInterfaceTypes "$THRESHOLD_INTERFACE_TYPES" \
        --graphLifting "$GRAPH_LIFTING" \
        --graphLiftingAfterQuantization "$GRAPH_LIFTING_AFTER_QUANTIZATION" \
        --huffman "$HUFFMAN" \
        --bitlength "$BITLENGTH" \
        --zigzag "$ZIGZAG" \
        --threads "$THREADS" \
        --subres "$SUBRES" \
        --elmres "$ELMRES" \
        --iter "$ITER" \
        --tol "$TOL" \
        --timing 0 \
        --vtk 0 \
        --histogram "$HISTOGRAM" \
        --sigma 1.0 \
        --interfacetypes "$INTERFACE_TYPES" \
        --crisscross "$CRISSCROSS" \
      --bddcExcludeCornersFromCompression "$EXCLUDE_CORNERS"
  '

if command -v sacct >/dev/null 2>&1; then
  sacct -j "$SLURM_JOB_ID" --format=JobID,State,Elapsed,TotalCPU,AllocCPUS,MaxRSS
fi
