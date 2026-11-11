#!/usr/bin/env bash
#SBATCH --partition=big
#SBATCH --constraint=Gold6342
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --job-name=p2d-interface
#SBATCH --output=output/poisson-interface-%j.log
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from poisson_new_scratch}"
export LANG=C LC_ALL=C LC_CTYPE=C
poisson_sif=${poisson_sif:-/data/numerik/people/fchegini/containers/kaskade-robin-20261105/kaskade7-robin.sif}
poisson_kaskade=${poisson_kaskade:-/data/numerik/people/fchegini/kaskade7_new}
# subres is per axis: 16 gives 256 subdomains in 2D.
subres=${POISSON_SUBRES:-16}
elmres=${POISSON_ELMRES:-32}
bits=${POISSON_BITS:-16}
# One launcher task reserves the CPUs used by its local MPI ranks.
ranks=${SLURM_CPUS_PER_TASK:-4}
result="/BDDC/poisson_new_scratch/output/poisson-interface-${SLURM_JOB_ID}"
printf 'Fresh Poisson 2D capture: subres=%s elmres=%s bits=%s ranks=%s\n' \
  "$subres" "$elmres" "$bits" "$ranks"
# The existing wrapper rebuilds bddc in this checkout and compiles the replay
# helper. mpirun runs inside one container on this single-node allocation.
srun --mpi=none --nodes=1 --ntasks=1 --ntasks-per-node=1 \
  --cpus-per-task="$ranks" \
  apptainer exec --cleanenv \
  --bind "$poisson_kaskade:/kaskade7" \
  --bind "$(dirname "$PWD"):/BDDC" \
  --pwd /BDDC/poisson_new_scratch "$poisson_sif" \
  /bin/bash tests/run_poisson_frozen_order_check.sh \
  --subres "$subres" --elmres "$elmres" --bits "$bits" \
  --ranks "$ranks" --iterations 1000 --tolerance 1e-7 \
  --output "$result"
printf 'Results: %s/output/poisson-interface-%s\n' "$PWD" "$SLURM_JOB_ID"
