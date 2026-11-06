#!/usr/bin/env bash
# Run inside the Open MPI Kaskade container, with /kaskade7 and /BDDC mounted.
set -euo pipefail
cd "$(dirname "$0")/.."
kaskade_dir=${KASKADE7:-/kaskade7}
mpi_dir=/usr/lib64/mpi/gcc/openmpi4
if [[ ! -x "$mpi_dir/bin/mpicxx" ]]; then
  echo 'Run this check inside the Kaskade Open MPI container.' >&2
  exit 2
fi
export PATH="$mpi_dir/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin"
export LD_LIBRARY_PATH="$mpi_dir/lib64:$mpi_dir/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 LANG=C LC_ALL=C
# This validation is single-node. Avoid probing unavailable fabric devices.
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp
check_dir=$(mktemp -d /tmp/poisson-transform-order.XXXXXX)
mpicxx -std=c++20 -O2 -I "$kaskade_dir" \
  "$kaskade_dir/tests/bddcTransferOrder/bddcTransferOrder.cpp" \
  -o "$check_dir/transfer-order"
"$check_dir/transfer-order"
PYTHONDONTWRITEBYTECODE=1 python3 tests/test_poisson_diagnostic.py
# Preserve dependency files; force the application TU to instantiate updated
# transfer headers. Existing Kaskade library objects are not rebuilt by this check.
for dep in bddc.d bddc-add.d; do
  if [[ -f "$dep" ]]; then cp -- "$dep" "$check_dir/$dep"; else touch "$dep"; fi
done
if ! make -W bddc.cpp -o "$kaskade_dir/libs/libkaskade.a" \
    KASKADE7="$kaskade_dir" MPI=1 MPICXX=mpicxx \
    MPI_FORCE_INCLUDE="-include $mpi_dir/include/mpi.h" \
    MPI_INCLUDE_FLAGS='-DMUMPS_MPI_H' bddc; then
  if [[ -f gccerr.txt ]]; then tail -80 gccerr.txt; fi
  exit 1
fi
# Defaults are a small two-rank, 16-bit check. Additional CLI arguments can
# override bits, ranks, resolution and tolerance. No jobs are submitted here.
python3 tests/run_poisson_quantization_diagnostic.py \
  --bits 16 --ranks 2 --subres 4 --elmres 8 \
  --transform-orders before after --mpi-args=--allow-run-as-root "$@"
