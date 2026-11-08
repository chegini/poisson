# Poisson BDDC scratch example

For the current geometric `bddc` solver matrix and cluster commands, see
[Current solver, transform-order and MPI validation](#current-solver-transform-order-and-mpi-validation).
The dated audit below records the earlier implementation.

## Compression audit and validation (31 October 2026)

The old `bddc` and `bddc-add` used `uint16_t` quantization symbols: requests above
16 bits, including the scripts' default 22, were silently clamped to 16 by
Kaskade. The geometric driver also used float32 transfer slots. Their default
solver was **steepest descent**, despite some historical documentation calling
it CG. `--tol` compared an internal correction/energy quantity, and exhaustion
of the iteration limit could still exit successfully. Historical bit-depth
labels and convergence claims must therefore be checked against their source
and executable versions; old 24/32/48/64 labels do not establish those precisions.

Both drivers now use **uint32 symbols and double transfer slots**, validate
requested precisions in **1..32**, and reject 48/64 rather than silently clamp.
This changes numerical/transport behavior and invalidates direct comparisons
of old byte ratios or results without rerunning matched configurations.
Double slots keep every uint32 quantization symbol exactly representable;
scale/min metadata are still float32, so 32-bit symbols do not guarantee 32-bit
reconstruction accuracy. Quantization-only transports each symbol in an
eight-byte slot; the full codec packs variable-length streams.

The options `--useCg 1`, `--restrictionPasses N`,
`--restrictCompressionBits N` and `--prolongateCompressionBits N` are available
in both executables. Directional precisions default to `-1` (inherit
`--compressionBits`). Extra passes require compression; each pass exchanges
the remaining restriction reconstruction error while the solution stays fixed.
Only one coarse solve and prolongation follow those exchanges.

Legacy solver/stopping defaults remain **SD**, one restriction pass, and
`--relativeTolerance 0` (use `--tol`). New diagnostics explicitly select CG and
`--relativeTolerance 1e-7`. In true-residual mode, the read-only monitor forms
one continuous global solution by averaging shared copies, then evaluates
`b-A*x` with the original assembled local matrices and RHS. It counts each
physical global DOF once and sums each owner contribution once. Zero RHS uses
the absolute residual. A separate consistency check requires
`||local copies - continuous copies||_2 / max(1,||x||_2) <= relativeTolerance`;
this prevents projection from concealing a discontinuous iterate. No diagnostic
projection is written back into the solver. Failure or iteration exhaustion
now returns a nonzero exit code. `convergence.csv` distinguishes solver calls
from accepted updates, and `solution.csv` stores full-precision global values.

The validated full pipeline uses reversible integer graph lifting **after**
quantization (`--graphLifting 0 --graphLiftingAfterQuantization 1`), now the
default transform configuration. At the time of this audit, pre-quantization
lifting was rejected because its neighbor-wise reconstruction did not match nodal
published-residual feedback. That restriction has since been removed after the
producer-owned reconstruction fix; the current campaign tests `before`, `after`
and `both` explicitly. To disable transforms, set both lifting options to zero.
Use the current runner below for matched transform-order comparisons.

Inside the same OpenMPI container as EMI, rebuild without deleting result files:

```bash
cd /BDDC/poisson_new_scratch
make -j1 -W bddc.cpp -W bddc-add.cpp KASKADE7=/kaskade7 MPI=1 MPICXX=mpicxx \
  MPI_FORCE_INCLUDE="-include /usr/lib64/mpi/gcc/openmpi4/include/mpi.h" \
  MPI_INCLUDE_FLAGS="-I/usr/lib64/mpi/gcc/openmpi4/include -DMUMPS_MPI_H"
```

Start with a small matched 2-D Poisson problem (16 subdomains, sigma=1):

```bash
python3 tests/run_poisson_quantization_diagnostic.py \
  --bits 8 12 16 24 32 --restriction-passes 1 \
  --full-compression --ranks 8 --subres 4 --elmres 8 \
  --tolerance 1e-7 --mpi-args="--allow-run-as-root"
```

This runs an assembled UMFPACK direct reference on one rank, uncompressed CG on
the requested ranks, then quantization-only and integer-lifting/Huffman cases.
The mesh, coefficient, constraints, thread counts and tolerance match in all
cases. Failed compressed cases remain recorded while the sweep continues;
reference/baseline failure stops the sweep. The diagnostic does **not** assert
that any particular bit depth converges. Compare two-pass 12-bit and three-pass
8-bit only afterward, retaining exactly the same `subres`, `elmres` and `sigma`.

Each new output directory retains exact commands, binary hash, environment,
logs, convergence histories, global solutions, metrics and a summary. Raw global
solution errors are compared against direct and uncompressed references; no
mean offset is removed and no hidden solution-error threshold is applied.
Net byte savings compare **total** logical payload bytes against uncompressed
CG, including all restriction passes and changed iteration counts. Same-rank
delivery is included; diagnostic reductions, MPI framing, allgather replication
and codebook training are excluded. These are validation timings, not pure MPI
communication timings. Codec-only timers are not instrumented and are no longer
printed as if zero represented a measured cost.

The shell runners additionally expose `COMPRESSION`, `USE_CG`, `REL_TOL`,
`RESTRICT_BITS`, `PROLONGATE_BITS` and `RESTRICTION_PASSES`. `run_case_mpi.sh`
uses `mpirun` from PATH by default; use `MPI_LAUNCHER` to override it and
`MPI_EXTRA_ARGS="-bootstrap fork"` only for an MPI implementation requiring it.
The new diagnostic leaves historical Slurm campaigns and archived data alone;
their old 48/64-bit cases and pre-quantization settings must be revised before
submitting new jobs.

Local checks (no real Poisson or MPI runs):

```bash
c++ -std=c++20 -Wall -Wextra -pedantic tests/test_poisson_validation.cpp -o /tmp/test_poisson_validation
/tmp/test_poisson_validation
python3 -m unittest discover -s tests -p 'test_poisson_diagnostic.py' -v
bash -n run_case.sh run_case_mpi.sh
```

The remaining sections give individual executable options and cluster context.

This directory contains two Poisson examples:

- `bddc`: geometric Poisson problem with BDDC;
- `bddc-add`: algebraic/additive Poisson variant with optional refinement and
  decomposition controls.

The executables accept command-line options in the form:

```bash
./bddc --option value
./bddc-add --option value
```

Use the built-in help for the exact options compiled into the executable:

```bash
./bddc --help
./bddc-add --help
```

## Build

Inside the Kaskade MPI container, build the examples with:

```bash
make -j1 -W bddc.cpp -W bddc-add.cpp KASKADE7=/kaskade7 MPI=1 MPICXX=mpicxx
```

### Cluster workflow

The Slurm build script compiles the MPI-enabled Poisson executables inside
the container:

```bash
sbatch build_poisson.sbatch
```

The build is needed after changing Poisson source files, changing Kaskade
headers, or when `bddc` is missing. Reconnecting to the cluster or running on
a different login or compute node does not by itself require a rebuild because
the project is stored on the shared filesystem.

After a successful build, submit the two benchmark runs independently:

```bash
sbatch run_poisson_serial_8threads.sbatch
sbatch run_poisson_mpi_8ranks.sbatch
```

The first job runs one MPI process with eight Kaskade threads. The second runs
eight MPI ranks with one thread per rank. Logs are written to
`logs/poisson-serial8-<jobid>.out` and `logs/poisson-mpi8-<jobid>.out`.

## Common options

These options are available in both executables.

| Option | Default | Meaning |
| --- | ---: | --- |
| `--prefix` | `.` | Output directory for VTK files |
| `--timing` | `1` | Print timing information |
| `--vtk` | `1` | Write VTK output |
| `--threads` | `-1` | Thread count; `0` sequential, negative uses the default |
| `--subres` | `2` | Subdomains per coordinate direction |
| `--elmres` | `2` | Elements per subdomain and coordinate direction |
| `--iter` | `16` | Maximum BDDC iterations |
| `--tol` | `1e-10` | Legacy internal norm tolerance; active when `relativeTolerance=0` |
| `--relativeTolerance` | `0` | True assembled relative residual target; `0` selects legacy stopping |
| `--useCg` | `0` | `1` selects CG, `0` selects steepest descent |
| `--crisscross` | `1` | Use criss-cross triangulation |

Example:

```bash
./bddc \
  --threads 2 \
  --subres 2 \
  --elmres 2 \
  --iter 100 \
  --tol 1e-10 \
  --timing 1 \
  --vtk 0 \
  --prefix output/poisson_bddc
```

`--threads 0` runs sequentially. A negative value uses the Kaskade default
thread configuration. With `--subres 4`, the two-dimensional problem contains
`4^2 = 16` subdomains. Increasing `--elmres` increases the local problem size.

## Coarse interface constraints

The `--interfacetypes` option selects coarse constraints using bit flags:

| Value | Constraint |
| ---: | --- |
| `1` | Corners |
| `2` | Edges |
| `4` | Faces |
| `7` | Corners, edges, and faces |

Examples:

```bash
# Corners only
./bddc --interfacetypes 1

# Corners and edges
./bddc --interfacetypes 3

# All available interface types
./bddc --interfacetypes 7
```

## Diffusion coefficient

The `bddc` executable supports:

```bash
--sigma 1.0
```

This sets the diffusion coefficient in the top-right part of the unit square.
It can be used to test coefficient jumps and heterogeneous diffusion.

```bash
./bddc --sigma 100.0 --iter 100 --tol 1e-10
```

## MPI options

| Option | Default | Meaning |
| --- | ---: | --- |
| `--mpi` | `0` | Enable owner-only MPI BDDC |

Serial execution does not require `mpirun`:

```bash
./bddc --mpi 0 --vtk 0
```

MPI execution requires the MPI-enabled executable and launcher:

```bash
mpirun --allow-run-as-root -np 2 ./bddc \
  --mpi 1 \
  --vtk 0 \
  --threads 2 \
  --subres 2 \
  --elmres 2 \
  --iter 100 \
  --tol 1e-10 \
  --timing 1 \
  --prefix output/poisson_bddc_mpi
```

When MPI is enabled, subdomains are assigned to MPI ranks using owner-only
storage. Each rank owns only its assigned subdomains and exchanges interface
information with the other ranks. VTK output is disabled for MPI runs because
`--vtk 1` is not supported together with `--mpi 1`.

## Compression options

These options are available in both Poisson executables.

| Option | Default | Meaning |
| --- | ---: | --- |
| `--compression` | `0` | Enable compressed BDDC transfer |
| `--compressionBits` | `16` | Quantization alphabet precision, 1..32 |
| `--restrictCompressionBits` | `-1` | Restriction precision override; `-1` inherits common precision |
| `--prolongateCompressionBits` | `-1` | Prolongation precision override; `-1` inherits common precision |
| `--restrictionPasses` | `1` | Restriction exchanges per compressed solver update, 1..16 |
| `--compressionMinInterfaceSize` | `0` | Minimum scalar interface size for Huffman coding; `0` encodes all interfaces |
| `--compressionThresholdInterfaceTypes` | `4` | Interface type receiving the size threshold: `2=edge`, `4=face`, `0=all` |
| `--graphLifting` | `0` | Floating graph lifting before quantization; set integer-after mode to `0` |
| `--graphLiftingAfterQuantization` | `1` | Apply reversible integer graph lifting after quantization |
| `--huffman` | `1` | Enable Huffman coding |
| `--bitlength` | `1` | Enable bit-length/tail encoding |
| `--histogram` | `0` | Write compression symbol histograms below `--prefix` |

Plain BDDC uses:

```bash
--compression 0
```

The full compression pipeline uses:

```bash
--compression 1 \
--compressionBits 16 \
--graphLifting 0 --graphLiftingAfterQuantization 1 \
--huffman 1 \
--bitlength 1
```

The transfer pipeline is:

```text
interface vector (or remaining restriction error)
-> quantization
-> optional reversible integer graph lifting after quantization
-> bit-length encoding
-> Huffman coding
-> packed payload
```

The receiver applies the inverse operations. Compression variants can be
tested as follows:

`run_case_mpi.sh` includes `SUBRES`, `ELMRES`, and the optional `CASE_TAG` in
the output directory name. This prevents parameter sweeps from overwriting
one another.

```bash
# Quantization only
--compression 1 --compressionBits 16 \
--graphLifting 0 --graphLiftingAfterQuantization 0 --huffman 0 --bitlength 0

# Quantization plus Huffman coding
--compression 1 --compressionBits 16 \
--graphLifting 0 --graphLiftingAfterQuantization 0 --huffman 1 --bitlength 0

# Full compression pipeline
--compression 1 --compressionBits 16 \
--graphLifting 0 --graphLiftingAfterQuantization 1 --huffman 1 --bitlength 1

# Reversible graph lifting after quantization
--compression 1 --compressionBits 16 \
--graphLifting 0 --graphLiftingAfterQuantization 1 \
--huffman 1 --bitlength 1
```

To export the empirical symbol distributions used by the compression setup,
add `--histogram 1`. The output directory is the value of `--prefix`; it is
created automatically for histogram output. The run writes CSV files for
restriction and prolongation:

```text
restrict_raw.csv          quantized symbols before bit-length coding
restrict_zigzag.csv       zigzag/unsigned symbol representation
restrict_bitlength.csv    symbols supplied to the Huffman codebook
prolongate_raw.csv
prolongate_zigzag.csv
prolongate_bitlength.csv
metadata.txt              active transform and bit-length settings
```

Example for the full pipeline:

```bash
./bddc \
  --mpi 0 \
  --compression 1 \
  --compressionBits 16 \
  --graphLifting 0 --graphLiftingAfterQuantization 1 \
  --huffman 1 \
  --bitlength 1 \
  --histogram 1 \
  --threads 2 \
  --subres 8 \
  --elmres 64 \
  --iter 20 \
  --tol 1e-10 \
  --vtk 0 \
  --timing 1 \
  --prefix output/poisson-study/histograms/full16
```

For a direct comparison, use separate output directories:

```bash
# No graph lifting
./bddc --mpi 0 --compression 1 --compressionBits 22 \
  --graphLifting 0 --graphLiftingAfterQuantization 0 \
  --huffman 1 --bitlength 1 --zigzag 1 --vtk 0 \
  --iter 20 --prefix output/poisson-study/no-lifting

# Graph lifting after quantization
./bddc --mpi 0 --compression 1 --compressionBits 22 \
  --graphLifting 0 --graphLiftingAfterQuantization 1 \
  --huffman 1 --bitlength 1 --zigzag 1 --vtk 0 \
  --iter 20 --prefix output/poisson-study/after-quantization
```

The histogram is collected from the first shared-codebook training exchange.
With the current Poisson executable, histogram collection is meaningful when
`--huffman 1` is enabled. The current quantized type is unsigned, so the
zigzag file is an explicit identity-stage record; signed quantized types would
show the signed-to-unsigned zigzag mapping there.

Test the precision/communication trade-off with:

```bash
--compressionBits 8
--compressionBits 12
--compressionBits 16
--compressionBits 24
```

## `bddc`-specific options

| Option | Default | Meaning |
| --- | ---: | --- |
| `--sigma` | `1.0` | Diffusion coefficient in the top-right region |

## `bddc-add`-specific options

| Option | Default | Meaning |
| --- | ---: | --- |
| `--refine` | `6` | Number of uniform mesh refinements |
| `--decompTol` | `1000.0` | Relative weak diagonal-dominance threshold |

Example:

```bash
./bddc-add \
  --refine 6 \
  --decompTol 1000 \
  --subres 2 \
  --elmres 2 \
  --iter 100 \
  --tol 1e-10
```

The `--refine` option changes the algebraic problem size before decomposition
and is important for performance and scalability studies.

## BDDC iterations versus CG

The Poisson `bddc` executable does not expose a standalone CG mode such as
`--solver cg`. The option:

```bash
--iter 100
```

controls BDDC iterations, not standalone CG iterations. A pure Poisson+CG
baseline requires a separate CG executable or solver implementation.

## Reproducible experiment commands

### Poisson + BDDC

```bash
./bddc \
  --mpi 0 \
  --compression 0 \
  --threads 2 \
  --subres 2 \
  --elmres 2 \
  --iter 100 \
  --tol 1e-10 \
  --vtk 0 \
  --timing 1 \
  --prefix output/poisson_bddc
```

### Poisson + BDDC + compression

```bash
./bddc \
  --mpi 0 \
  --compression 1 \
  --compressionBits 16 \
  --graphLifting 0 --graphLiftingAfterQuantization 1 \
  --huffman 1 \
  --bitlength 1 \
  --threads 2 \
  --subres 2 \
  --elmres 2 \
  --iter 100 \
  --tol 1e-10 \
  --vtk 0 \
  --timing 1 \
  --prefix output/poisson_bddc_compression
```

### Poisson + BDDC + MPI

```bash
mpirun --allow-run-as-root -np 2 ./bddc \
  --mpi 1 \
  --compression 0 \
  --threads 2 \
  --subres 2 \
  --elmres 2 \
  --iter 100 \
  --tol 1e-10 \
  --vtk 0 \
  --timing 1 \
  --prefix output/poisson_bddc_mpi
```

### Poisson + BDDC + compression + MPI

```bash
mpirun --allow-run-as-root -np 2 ./bddc \
  --mpi 1 \
  --compression 1 \
  --compressionBits 16 \
  --graphLifting 0 --graphLiftingAfterQuantization 1 \
  --huffman 1 \
  --bitlength 1 \
  --threads 2 \
  --subres 2 \
  --elmres 2 \
  --iter 100 \
  --tol 1e-10 \
  --vtk 0 \
  --timing 1 \
  --prefix output/poisson_bddc_compression_mpi
```

For fair comparisons, keep `--subres`, `--elmres`, `--iter`, `--tol`,
`--threads`, and the grid configuration fixed while changing only `--mpi`,
`--compression`, and the compression settings.

## Current solver, transform-order and MPI validation

`tests/run_poisson_solver_check.py` targets the geometric `bddc` executable. The
new solver/routing CLI options are not added to the separate algebraic `bddc-add`
driver. The campaign validates the stationary Poisson problem with
CG, one-term FCG and restarted right-preconditioned FGMRES. It exposes the latest
Kaskade destination routing and same-node compression bypass. Existing defaults
and the older frozen-buffer workflow are retained. No EMI time-integration or
physical-step histogram schedule is implied for this stationary example.

The new campaign explicitly uses original-residual CG rho in **all** controls
and compressed cases (`--originalResidual 1`); compare its own controls, not the
iteration count of an older run with different recurrence settings. Constant
splitting (`--constantSplit 1`) is an opt-in prolongation experiment; it does not
prove that a lossy preconditioner is SPD. FCG and FGMRES also require measured
convergence. Split and local-bypass campaigns keep histograms disabled.

Transform orders are:

| Order | Pipeline |
|---|---|
| `before` | Floating graph lifting → quantization → full coding |
| `after` | Quantization → reversible integer graph lifting → full coding |
| `both` | Floating graph lifting → quantization → reversible integer lifting → full coding |

Full coding includes zigzag, bit-length/tail coding and Huffman. `both` explicitly
sets `--graphLiftingBothOrders 1`; merely enabling the two older switches retains
their existing precedence. When both runs are present, the campaign requires
identical exported solutions and convergence traces for `before` and `both`.

The FGMRES kernel and root-owned MPI adapter were ported from the validated EMI
implementation at commit `a850a1f`. Arnoldi vectors are broadcast, local BDDC
corrections are reduced, and the exact assembled operator is applied on rank zero.
The preconditioner uses local matrix-diagonal partition-of-unity weights on input
and output. This is a **validation implementation**, not a scalable distributed
Krylov implementation. Its outer-solver logical communication is reported
separately from BDDC interface bytes. Final errors use an independent distributed
assembly of the true residual. FGMRES stores a global iterate; its zero interface
jump describes that global representation, not a separate audit of every
preconditioner intermediate.

### Build on ZIB

From `poisson_new_scratch`, inside an active Slurm allocation:

```bash
poisson_sif=/data/numerik/people/fchegini/containers/kaskade-robin-20261105/kaskade7-robin.sif
poisson_kaskade=/data/numerik/people/fchegini/kaskade7_new

srun --mpi=none --nodes=1 --ntasks=1 --ntasks-per-node=1 \
  apptainer exec --cleanenv \
  --bind "$poisson_kaskade:/kaskade7" \
  --bind "$(dirname "$PWD"):/BDDC" \
  --pwd /BDDC/poisson_new_scratch "$poisson_sif" \
  /bin/bash tests/build_poisson_cluster.sh
```

The build preserves dependency files, forces compilation of the application, and
runs the transfer, residual, FGMRES and Python tests. Build once before submitting
campaigns; do not replace the executable while any campaign is running.

### Local MPI validation inside the existing Open MPI container

```bash
python3 tests/run_poisson_solver_check.py \
  --ranks 2 --expected-nodes 1 --subres 4 --elmres 8 \
  --bits 16 --solvers cg fcg fgmres --split both \
  --transform-orders before after both \
  --check-routing --check-local-bypass \
  --mpi-args=--allow-run-as-root
```

This runs 61 cases: one direct reference, six uncompressed controls, and 54
compressed cases. Single-node bypass must match the same solver's uncompressed
solution and convergence trace exactly. Destination and legacy routing must also
match exactly, including logical bytes. A failed baseline stops the campaign;
failed compressed cases are recorded and the remaining experiments continue.

For the earlier mesh and additional precisions, use `--subres 4 --elmres 64
--bits 8 12 16 24`. Low-bit failures are experimental outcomes, not silently relaxed
acceptance criteria. Use `--plan-only` to inspect the cases before running them.

### Two, four or eight nodes

Submit from the Poisson repository root after the build; `output` must exist for
Slurm's log file. Each layout uses two MPI ranks per node:

```bash
mkdir -p output
sbatch --nodes=2 --ntasks=4 --ntasks-per-node=2 \
  tests/submit_poisson_solver_check.sbatch --subres 4 --elmres 64
sbatch --nodes=4 --ntasks=8 --ntasks-per-node=2 \
  tests/submit_poisson_solver_check.sbatch --subres 4 --elmres 64
sbatch --nodes=8 --ntasks=16 --ntasks-per-node=2 \
  tests/submit_poisson_solver_check.sbatch --subres 4 --elmres 64
```

Defaults test 16-bit transform-first with split disabled, all three solvers, legacy
versus destination routing, and local bypass. Add `--split both --transform-orders
before after both --bits 12 16` to expand the matrix and request more wall time if
needed. The runner checks actual MPI shared-memory topology, requires inter-node
R/P traffic on multi-node layouts, and bounds each solve (`--timeout 600`).
The script uses `self,vader,tcp` with Vader single-copy disabled, preserving the
configuration that completed the EMI cluster checks. Same-node MPI remains active;
bypass avoids compression, not MPI or all global collectives. Any DOF with owners
sharing a node remains exact on all copies, including remote copies. Therefore
local bypass can change inter-node bytes and the numerical trajectory.

Every output directory contains commands, source/executable hashes, per-case logs,
`summary.json`, `summary.csv`, and `summary.md`. They report iterations, true residuals, raw
relative Euclidean solution errors against direct and matched uncompressed solves,
solver-wrapper/launcher times, and logical same-rank/same-node/inter-node bytes.
These are coefficient norms, not mass-matrix or fractional Sobolev norms.
`PASS` checks completion, residuals and applicable parity; **no scientific
solution-error threshold is imposed**. Multi-node bypass differences are measured,
not required to equal the fully compressed trajectory.

Byte savings use the same solver's uncompressed baseline, including its actual
iteration count. They exclude MPI framing, coarse/setup communication, codebook
training and validation collectives; FGMRES outer communication is separate.
Solver-wrapper time includes setup and validation/output. Application time is
launcher wall time. These checks do not establish a runtime speedup; performance
claims need repeated measurements under matched output and validation settings.

Poisson uses a `1e8` Dirichlet penalty on the outer boundary. In the first real
16-bit, `subres=4`, `elmres=8` split sweep, CG broke down/stalled and FCG broke down
for transform-first/combined splitting, while FGMRES converged. These are recorded
numerical outcomes; splitting is therefore **off by default** in the Poisson
campaign. The interaction of large diagonal weights with the constant split is a
possible explanation, not a verified root cause. Use `--split both` to investigate;
no tolerance or breakdown test is relaxed. An all-local bypass selects the native
uncompressed solver instantiation (as in EMI), preserving bitwise raw parity.

### Validation performed on 8 November 2026

Tests used an isolated source/executable snapshot on `lxdesk-010399`, the existing
`localhost/kaskade7-robin-cluster:20261105` image, and two Open MPI ranks on one
physical node. The working directory's existing binaries were not replaced.

- `subres=4`, `elmres=64`, 12/24 bits, CG/FCG/FGMRES, before/after/both, split off:
  **61/61 cases passed**, including direct/raw controls, legacy/destination parity,
  exact single-node bypass parity and combined/before parity. All iterative cases
  used 10 updates. Fully coded logical payload savings were about 71.3–73.1% at
  12 bits and 53.2–58.1% at 24 bits. FGMRES outer traffic is excluded from these
  interface-payload percentages.
- `subres=4`, `elmres=8`, 16 bits, split off/on with the same full matrix:
  **51/61 cases passed**. Ten failures are five split configurations under both
  routes: CG before/after/both and FCG before/both. CG after reached 1,000 updates;
  the other four broke down. All 18 local-bypass cases matched raw controls
  exactly. FGMRES split cases converged but required 33–50 updates versus six
  without splitting; splitting is not recommended as this campaign's default.
- Six Python tests cover mode validation, reference failures, matrix construction,
  topology/rank/counter rejection and frozen observer checks. Standalone C++ tests
  cover residual assembly, FGMRES numerical failure handling and MPI accounting.
- The existing direct/CG/quantization/capture/before/after/both workflow passed
  **7/7 cases**, including real production frozen packet replay and observer parity.
  The Open MPI adapter test passed on two ranks, including failure synchronization.
  Both `bddc` and `bddc-add` compiled with the current shared options header;
  `bddc-add` was a build-compatibility check, not a new solver-matrix run.

The detailed runs are under
`output/poisson-upgrade-validation-20261108.trb7Np/validation-mesh64` and
`output/poisson-upgrade-validation-20261108.trb7Np/validation-small`; the older
workflow regression is in the same snapshot under `validation-regression`.
These are single-node correctness results, not multi-node Poisson validation or
repeated performance measurements. The two/four/eight-node batch commands above
are the next validation stage. Solution errors remain reported quantities with
no imposed scientific acceptance threshold.

## 3D full-compression validation

The geometric `bddc` driver accepts `--dimension 2` (default) or `--dimension 3`.
3D uses a unit cube split into `subres^3` subdomains, with symmetric tetrahedral
meshes and P1 elements. All three coordinate directions participate in boundary
conditions and shared-DOF matching. The source is the negative Laplacian of
`u = x^2(1-x)y^2(1-y)z^2(1-z)` for `sigma=1`; the existing outer Dirichlet penalty
is retained. The 3D RHS uses degree-eight quadrature for P1 test functions.
The historical 2D RHS quadrature and arithmetic are preserved.

Start with **`--subres 2 --elmres 4`**, not the 2D value `--elmres 64`: 3D mesh
storage grows cubically. The starting mesh has eight subdomains, 12,288 tetrahedra
and 2,969 global potential DOFs. The runner checks the dimension, subdomain count
and global DOF count so an accidentally selected 2D binary cannot pass as 3D.

After rebuilding with `tests/build_poisson_cluster.sh`, submit from this directory:

```bash
mkdir -p output
poisson_job3d=$(sbatch --parsable \
  --nodes=2 --ntasks=4 --ntasks-per-node=2 \
  tests/submit_poisson_solver_check.sbatch \
  --dimension 3 --subres 2 --elmres 4 \
  --bits 12 --solvers cg fcg fgmres --split off \
  --transform-orders before after both \
  --no-check-routing --no-check-local-bypass \
  --iterations 1000)
echo "3D Poisson job: $poisson_job3d"
tail -n 60 -F "output/poisson-check-${poisson_job3d}.log"
```

This selects **full compression on all interfaces**, with destination routing,
including same-node interfaces. Splitting and histograms are off. Zigzag,
bit-length/tail coding and Huffman remain enabled. The two `--no-check-*` flags
override the batch wrapper's broader validation defaults; they remove legacy
routing and local-bypass cases. Direct and solver-matched uncompressed controls
are retained for error and byte comparisons: 13 cases total, of which nine use
12-bit compression. Select `--solvers cg` for a five-case CG-only sweep.

The same `self,vader,tcp` transport with Vader single-copy disabled is used.
FGMRES remains the root-owned validation adapter described above. Reported errors
compare coefficient vectors with direct and uncompressed solutions, not the
manufactured continuum solution; no scientific error threshold is imposed.

3D validation on 8 November 2026: all 13 cases above passed in the existing ZIB
Open MPI container with two ranks on one node. Raw controls used five updates;
12-bit CG/FCG used six and 12-bit FGMRES used five. Combined/before solution and
trace parity passed for all three solvers. All 13 corresponding 2D regression
cases matched the pre-extension solutions and traces exactly. Seven Python tests
and the standalone geometry test passed. This does not yet establish multi-node
3D correctness or performance.

The tested executable and full logs are preserved under
`output/poisson-3d-validation-20261108.nMc1z6/` (`validation-3d` and
`validation-2d-regression`). To run that exact build without replacing the primary
`bddc` executable, add
`--executable "$PWD/output/poisson-3d-validation-20261108.nMc1z6/bddc"`
to the batch command above.

## Final solution visualization

Add `--vtk 1` to `tests/run_poisson_solver_check.py` or its batch wrapper to
write `<case>/solution.vtu` after each converged solve. Default is `--vtk 0`.
This requires rebuilding the executable with the VTK export support; an older
isolated executable will not acquire it just by updating the runner.

For 3D FCG, 12-bit full compression and transform-first, select:

```bash
poisson_vtk_job=$(sbatch --parsable \
  --nodes=2 --ntasks=4 --ntasks-per-node=2 \
  tests/submit_poisson_solver_check.sbatch \
  --dimension 3 --subres 2 --elmres 4 \
  --bits 12 --solvers fcg --split off --transform-orders before \
  --no-check-routing --no-check-local-bypass --vtk 1 --iterations 1000)
tail -n 60 -F "output/poisson-check-${poisson_vtk_job}.log"
```

There are three cases: direct, uncompressed FCG, and compressed FCG. Under the
results directory printed in the log, open
`fcg-before-12-destination/solution.vtu` in ParaView, click Apply, and colour by
`u`. A Slice filter shows the interior of the cube; `subdomain` is also available
as cell data. This is a final solution export, not an iteration animation.

Rank zero writes the entire conforming triangle/tetrahedron mesh using the same
global DOF numbering and Float64 coefficients as `solution.csv`. The runner
checks point/cell counts and exact VTK/CSV coefficient equality. MPI ranks
collectively report write failure. Export does not change the solver state;
reported validation timings include export and should not be used as timing
benchmarks. Historical serial runs without residual monitoring retain their
per-iteration output.

VTK validation on 8 November 2026: direct, uncompressed FCG and compressed FCG
passed with two MPI ranks on one node in both 2D and 3D. All exported coefficient
vectors and solver traces exactly matched the VTK-off runs. Mesh checks verified
connectivity, cell types, nondegenerate cells and total area/volume one. Eight
Python tests passed. The tested VTK-capable executable is
`output/poisson-vtk-validation-20261108.1fEGCC/bddc`; select it with
`--executable "$PWD/output/poisson-vtk-validation-20261108.1fEGCC/bddc"`.
Its `validation-vtk` directory contains the three 3D examples and
`validation-vtk-2d` contains the corresponding 2D checks. This export has not yet
been tested across physical nodes.
