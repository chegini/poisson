#!/bin/bash
#SBATCH --partition=big
#SBATCH --nodes=2
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --job-name=poisson-check
#SBATCH --output=output/poisson-check-%j.log
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from poisson_new_scratch}"
export LANG=C LC_ALL=C LC_CTYPE=C
poisson_sif=${poisson_sif:-/data/numerik/people/fchegini/containers/kaskade-robin-20261105/kaskade7-robin.sif}
poisson_kaskade=${poisson_kaskade:-/data/numerik/people/fchegini/kaskade7_new}
exec python3 -u tests/run_poisson_solver_check.py \
  --ranks "$SLURM_NTASKS" --expected-nodes "$SLURM_JOB_NUM_NODES" \
  --slurm-container "$poisson_sif" --kaskade "$poisson_kaskade" \
  --check-routing --check-local-bypass "$@"
