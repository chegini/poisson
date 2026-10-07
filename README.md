# Poisson BDDC scratch example

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
make clean KASKADE7=/kaskade7 MPI=1
make -j1 KASKADE7=/kaskade7 MPI=1 MPICXX=mpicxx
```

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
| `--tol` | `1e-10` | BDDC residual tolerance |
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
| `--compressionBits` | `16` | Quantization precision |
| `--graphLifting` | `1` | Enable graph-lifting transform |
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
--graphLifting 1 \
--huffman 1 \
--bitlength 1
```

The transfer pipeline is:

```text
interface vector
-> graph lifting
-> quantization
-> bit-length encoding
-> Huffman coding
-> packed payload
```

The receiver applies the inverse operations. Compression variants can be
tested as follows:

```bash
# Quantization only
--compression 1 --compressionBits 16 \
--graphLifting 0 --huffman 0 --bitlength 0

# Quantization plus Huffman coding
--compression 1 --compressionBits 16 \
--graphLifting 0 --huffman 1 --bitlength 0

# Full compression pipeline
--compression 1 --compressionBits 16 \
--graphLifting 1 --huffman 1 --bitlength 1
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
  --graphLifting 1 \
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
  --graphLifting 1 \
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
  --graphLifting 1 \
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
