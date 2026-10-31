# Poisson BDDC scratch example

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
default transform configuration. Pre-quantization lifting is rejected when
compression is enabled without integer-after mode: its neighbor-wise transformed
reconstruction does not match the nodal published-residual feedback. To disable
transforms, set both options to zero. In the shell runners, `LIFTING_MODE=1` is
the new default, `0` disables lifting, and legacy `2` is rejected by the solver.

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
| `--graphLifting` | `0` | Legacy pre-quantization transform; rejected without integer-after mode |
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
