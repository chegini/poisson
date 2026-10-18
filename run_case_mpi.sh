#!/usr/bin/env bash
set -euo pipefail

BITS=${BITS:-22}
LIFTING_MODE=${LIFTING_MODE:-2}
ZIGZAG=${ZIGZAG:-1}
HUFFMAN=${HUFFMAN:-1}
BITLENGTH=${BITLENGTH:-1}
MPI_RANKS=${MPI_RANKS:-1}
MPI_LAUNCHER=${MPI_LAUNCHER:-/usr/lib64/mpi/gcc/openmpi4/bin/mpirun}
MPI_ALLOW_ROOT=${MPI_ALLOW_ROOT:-0}

MPI_ARGS=(-np "$MPI_RANKS")
if [[ "$MPI_ALLOW_ROOT" == 1 ]]; then
  MPI_ARGS+=(--allow-run-as-root)
fi

THREADS=${THREADS:-2}
SUBRES=${SUBRES:-4}
ELMRES=${ELMRES:-64}
ITER=${ITER:-20}
TOL=${TOL:-1e-10}
VTK=${VTK:-0}
HISTOGRAM=${HISTOGRAM:-0}
SIGMA=${SIGMA:-1.0}
INTERFACE_TYPES=${INTERFACE_TYPES:-7}
CRISSCROSS=${CRISSCROSS:-1}
EXECUTABLE=${EXECUTABLE:-./bddc}
EXCLUDE_CORNERS=${EXCLUDE_CORNERS:-0}

case "$EXCLUDE_CORNERS" in
  0|1) ;;
  *) printf 'Invalid EXCLUDE_CORNERS=%s; use 0 or 1\n' "$EXCLUDE_CORNERS" >&2; exit 2 ;;
esac

case "$LIFTING_MODE" in
  0)
    LIFTING=0
    AFTER_QUANTIZATION=0
    LIFTING_NAME="no-graph-lifting"
    ;;
  1)
    LIFTING=0
    AFTER_QUANTIZATION=1
    LIFTING_NAME="graph-lifting-after-quantization"
    ;;
  2)
    LIFTING=1
    AFTER_QUANTIZATION=0
    LIFTING_NAME="graph-lifting-before-quantization"
    ;;
  *)
    printf 'Invalid LIFTING_MODE=%s\n' "$LIFTING_MODE" >&2
    printf 'Use: 0 (none), 1 (after quantization), or 2 (before quantization)\n' >&2
    exit 2
    ;;
esac

CASE="mpi-ranks-${MPI_RANKS}_${LIFTING_NAME}_bits-${BITS}_corners-${EXCLUDE_CORNERS}_zigzag-${ZIGZAG}_huffman-${HUFFMAN}_bitlength-${BITLENGTH}"
if [[ "$HISTOGRAM" == 1 ]]; then
  CASE="${CASE}_histogram"
fi
RUN="output/poisson-study/${CASE}"
LOG="logs/poisson-${CASE}.log"

mkdir -p "$RUN" "$(dirname "$LOG")"

"$MPI_LAUNCHER" "${MPI_ARGS[@]}" "$EXECUTABLE" \
  --mpi 1 \
  --compression 1 \
  --compressionBits "$BITS" \
  --graphLifting "$LIFTING" \
  --graphLiftingAfterQuantization "$AFTER_QUANTIZATION" \
  --huffman "$HUFFMAN" \
  --bitlength "$BITLENGTH" \
  --zigzag "$ZIGZAG" \
  --threads "$THREADS" \
  --subres "$SUBRES" \
  --elmres "$ELMRES" \
  --iter "$ITER" \
  --tol "$TOL" \
  --vtk 0 \
  --histogram "$HISTOGRAM" \
  --sigma "$SIGMA" \
  --interfacetypes "$INTERFACE_TYPES" \
  --crisscross "$CRISSCROSS" \
  --bddcExcludeCornersFromCompression "$EXCLUDE_CORNERS" \
  --prefix "$RUN" \
  2>&1 | tee "$LOG"
