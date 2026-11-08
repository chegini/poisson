#!/usr/bin/env bash
# Build inside the Kaskade Open MPI image; does not launch MPI or submit jobs.
set -euo pipefail
cd "$(dirname "$0")/.."
kaskade_dir=${KASKADE7:-/kaskade7}
mpi_dir=/usr/lib64/mpi/gcc/openmpi4
export PATH="$mpi_dir/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin"
export LD_LIBRARY_PATH="$mpi_dir/lib64:$mpi_dir/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export LANG=C LC_ALL=C OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
check_dir=$(mktemp -d /tmp/poisson-check.XXXXXX)
mpicxx -std=c++20 -O2 -Wall -Wextra -Werror tests/test_poisson_validation.cpp -o "$check_dir/validation"
"$check_dir/validation"
mpicxx -std=c++20 -O2 -Wall -Wextra -Werror -I . tests/test_poisson_fgmres.cpp -o "$check_dir/fgmres"
"$check_dir/fgmres"
mpicxx -std=c++20 -O2 -Wall -Wextra -Werror -DKASKADE_HAVE_MPI -I . \
  tests/test_poisson_fgmres_mpi.cpp -o "$check_dir/fgmres-mpi"
# Match the existing Kaskade/EMI transfer-fixture compile flags.
mpicxx -std=c++20 -O2 -I "$kaskade_dir" \
  "$kaskade_dir/tests/bddcTransferOrder/bddcTransferOrder.cpp" -o "$check_dir/transfer"
"$check_dir/transfer"
python3 -m unittest discover -s tests -p 'test_*.py'
for dep in bddc.d bddc-add.d; do
  if [[ ! -f "$dep" ]]; then touch "$dep"; fi
done
if ! make -W bddc.cpp -o "$kaskade_dir/libs/libkaskade.a" KASKADE7="$kaskade_dir" \
  MPI=1 MPICXX=mpicxx MPI_FORCE_INCLUDE="-include $mpi_dir/include/mpi.h" \
  MPI_INCLUDE_FLAGS=-DMUMPS_MPI_H bddc; then
  tail -80 gccerr.txt
  exit 1
fi
