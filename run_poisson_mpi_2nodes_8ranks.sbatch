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
KASKADE7=/data/numerik/people/fchegini/kaskade7_new
IMAGE=localhost/kaskade7-mpi-build-environment:latest
ARCHIVE=/data/numerik/people/fchegini/kaskade7-mpi-build-environment.tar
MPI_IMPL=${MPI_IMPL:-pmi2}
HOST_MPI=/home/htc/fchegini/opencarp_dependencies/petsc_install
RUN="output/multi-node/sbatch-${SLURM_JOB_ID}"
export BASE KASKADE7 IMAGE ARCHIVE HOST_MPI RUN

cd "$BASE/poisson_new_scratch"
mkdir -p logs "$RUN"

echo "Job ID: $SLURM_JOB_ID"
echo "Nodes: $(scontrol show hostnames "$SLURM_JOB_NODELIST" | tr "\n" " ")"
echo "Tasks: $SLURM_NTASKS"
echo "MPI mode: $MPI_IMPL"

srun --ntasks="$SLURM_NNODES" --ntasks-per-node=1 bash -lc '
  if ! podman image exists '"$IMAGE"'; then
    podman load -i '"$ARCHIVE"'
  fi
  podman image exists '"$IMAGE"'
'

# Keep one stdout/stderr file per MPI task. Otherwise Slurm merges output from
# all ranks, making concurrent solver and compression messages unreadable.
srun --mpi="$MPI_IMPL" --ntasks="$SLURM_NTASKS" --ntasks-per-node=4 \
  --kill-on-bad-exit=1 \
  --output="logs/poisson-2nodes-8ranks-%j-rank%t.out" \
  --error="logs/poisson-2nodes-8ranks-%j-rank%t.err" \
  bash -lc '
    set -euo pipefail
    : "${PMI_FD:?PMI_FD was not provided by Slurm}"
    eval "exec 3<&${PMI_FD}"
    exec podman run --rm --pull=never \
      --network=host \
      --ipc=host \
      --pid=host \
      --preserve-fds=1 \
      --env-host \
      -e PMI_FD=3 \
      -v "$HOST_MPI:$HOST_MPI:ro" \
      -e PATH=$HOST_MPI/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
      -e LD_LIBRARY_PATH=/kaskade7/libs:$HOST_MPI/lib:$HOST_MPI/lib64:/usr/local/lib64:/usr/local/lib \
      -e OMP_NUM_THREADS=1 \
      -v "$KASKADE7:/kaskade7" \
      -v "$BASE:/BDDC" \
      -w /BDDC/poisson_new_scratch \
      "$IMAGE" \
      ./bddc \
        --mpi 1 \
        --compression 1 \
        --compressionBits 22 \
        --graphLifting 0 \
        --graphLiftingAfterQuantization 1 \
        --huffman 1 \
        --bitlength 1 \
        --zigzag 1 \
        --threads 1 \
        --subres 4 \
        --elmres 64 \
        --iter 20 \
        --tol 1e-10 \
        --vtk 0 \
        --histogram 0 \
        --sigma 1.0 \
        --interfacetypes 7 \
        --crisscross 1 \
        --bddcExcludeCornersFromCompression 1 \
        --prefix "/BDDC/poisson_new_scratch/$RUN"
  '
