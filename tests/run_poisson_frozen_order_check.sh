#!/usr/bin/env bash
# Run inside the existing Open MPI container. Single node, two ranks by default.
set -euo pipefail
cd "$(dirname "$0")/.."
kaskade_dir=${KASKADE7:-/kaskade7}
mpi_dir=/usr/lib64/mpi/gcc/openmpi4
export PATH="$mpi_dir/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin"
export LD_LIBRARY_PATH="$mpi_dir/lib64:$mpi_dir/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 LANG=C LC_ALL=C
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp
check_dir=$(mktemp -d /tmp/poisson-frozen-order.XXXXXX)
mpicxx -std=c++20 -O2 -I "$kaskade_dir" tests/compare_frozen_transfers.cpp \
  -o "$check_dir/replay"
"$check_dir/replay" "$check_dir/fixture" "$check_dir/fixture-replay" 16 --fixture
PYTHONDONTWRITEBYTECODE=1 python3 tests/summarize_frozen_transfers.py "$check_dir/fixture-replay"
PYTHONDONTWRITEBYTECODE=1 python3 tests/test_frozen_transfer_summary.py
bash tests/run_poisson_transform_order_check.sh \
  --subres 4 --elmres 64 --bits 16 \
  --frozen-replay-executable "$check_dir/replay" "$@"
