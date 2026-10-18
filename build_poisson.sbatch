#!/usr/bin/env bash
#SBATCH --job-name=build-poisson
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --output=logs/build-poisson-%j.out
#SBATCH --error=logs/build-poisson-%j.err

set -euo pipefail

BASE=/data/numerik/people/fchegini/project/fatemeh/MicroCard/current_papers/BDDC_petsc_SDC_2023/gatevariables/20Nov2023/BDDC
IMAGE=localhost/kaskade7-mpi-build-environment:latest
IMAGE_ARCHIVE=/data/numerik/people/fchegini/kaskade7-mpi-build-environment.tar
HOST_MPI=/home/htc/fchegini/opencarp_dependencies/petsc_install
BUILD_TMP="$BASE/poisson_new_scratch/build-tmp-${SLURM_JOB_ID:-manual}"
MPI_CXX="$HOST_MPI/bin/mpicxx"
MPI_INCLUDE_FLAGS="-I$HOST_MPI/include -DMUMPS_MPI_H"
MPI_FORCE_INCLUDE="-include $HOST_MPI/include/mpi.h"
MPI_LIBRARY_PATH="$HOST_MPI/lib:$HOST_MPI/lib64"

cd "$BASE/poisson_new_scratch"
mkdir -p "$BUILD_TMP"

echo "Job ID: ${SLURM_JOB_ID:-unknown}"
echo "Node: $(hostname)"
echo "Start: $(date)"

if ! podman image exists "$IMAGE"; then
  podman load -i "$IMAGE_ARCHIVE"
fi

echo "MPI compiler: $HOST_MPI/bin/mpicxx"

podman run --rm --pull=never \
  -e PATH=$HOST_MPI/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
  -e LD_LIBRARY_PATH=/kaskade7/libs:$MPI_LIBRARY_PATH:/usr/local/lib64:/usr/local/lib \
  -e CPATH=$HOST_MPI/include \
  -e CPLUS_INCLUDE_PATH=$HOST_MPI/include \
  -e LIBRARY_PATH=$MPI_LIBRARY_PATH \
  -e HOST_MPI=$HOST_MPI \
  -e MPI_CXX=$MPI_CXX \
  -e MPI_INCLUDE_FLAGS="-I$HOST_MPI/include -DMUMPS_MPI_H" \
  -e MPI_FORCE_INCLUDE="$MPI_FORCE_INCLUDE" \
  -e TMPDIR=/BDDC/poisson_new_scratch/$(basename "$BUILD_TMP") \
  -v "$HOST_MPI:$HOST_MPI:ro" \
  -v /data/numerik/people/fchegini/kaskade7_new:/kaskade7 \
  -v "$BASE:/BDDC" \
  -w /BDDC/poisson_new_scratch \
  "$IMAGE" \
  bash -lc '
    set -euo pipefail
    test -x "$MPI_CXX"
    test -f "$HOST_MPI/include/mpi.h"
    make clean KASKADE7=/kaskade7 MPI=1 MPICXX="$MPI_CXX" MPI_INCLUDE_FLAGS="-I$HOST_MPI/include -DMUMPS_MPI_H" MPI_FORCE_INCLUDE="$MPI_FORCE_INCLUDE"
    make -j1 KASKADE7=/kaskade7 MPI=1 MPICXX="$MPI_CXX" MPI_INCLUDE_FLAGS="-I$HOST_MPI/include -DMUMPS_MPI_H" MPI_FORCE_INCLUDE="$MPI_FORCE_INCLUDE"
    test -x bddc
    ls -lh bddc bddc-add
  '

echo "Finished: $(date)"
