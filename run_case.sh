#!/usr/bin/env bash
set -euo pipefail

BITS=${BITS:-22}
LIFTING_MODE=${LIFTING_MODE:-before_quantization}
ZIGZAG=${ZIGZAG:-1}
HUFFMAN=${HUFFMAN:-1}
BITLENGTH=${BITLENGTH:-1}

EXECUTABLE=${EXECUTABLE:-./bddc}
MPI=${MPI:-0}
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

case "$LIFTING_MODE" in
  none)
    LIFTING=0
    AFTER_QUANTIZATION=0
    LIFTING_NAME="no-graph-lifting"
    ;;
  before_quantization)
    LIFTING=1
    AFTER_QUANTIZATION=0
    LIFTING_NAME="graph-lifting-before-quantization"
    ;;
  after_quantization)
    LIFTING=0
    AFTER_QUANTIZATION=1
    LIFTING_NAME="graph-lifting-after-quantization"
    ;;
  *)
    printf 'Invalid LIFTING_MODE=%s\n' "$LIFTING_MODE" >&2
    printf 'Use: none, before_quantization, or after_quantization\n' >&2
    exit 2
    ;;
esac

CASE="${LIFTING_NAME}_bits-${BITS}_zigzag-${ZIGZAG}_huffman-${HUFFMAN}_bitlength-${BITLENGTH}"
RUN="output/poisson-study/${CASE}"
LOG="logs/poisson-${CASE}.log"

mkdir -p "$RUN" "$(dirname "$LOG")"

"$EXECUTABLE" \
  --mpi "$MPI" \
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
  --vtk "$VTK" \
  --histogram "$HISTOGRAM" \
  --sigma "$SIGMA" \
  --interfacetypes "$INTERFACE_TYPES" \
  --crisscross "$CRISSCROSS" \
  --prefix "$RUN" \
  2>&1 | tee "$LOG"
