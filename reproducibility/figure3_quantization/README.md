# Figure 3: Poisson iteration count versus quantization bits

Legacy source: `Compressed_BDDC/images/itr_vs_bits.png` and
`Compressed_BDDC/code/itr_vs_bits.m`.

The historical setup is a 2-D unit square with 256 subdomains, a fixed global
512 x 512 element mesh, and 16 interface DoFs. The reproduction uses
`SUBRES=16`, `ELMRES=32`, so the current code has the same 256-subdomain and
512 x 512 global-resolution structure.

The first pass is quantization only. This isolates the effect plotted in the
legacy figure:

```text
LIFTING_MODE=0, ZIGZAG=0, HUFFMAN=0, BITLENGTH=0, MIN_INTERFACE_SIZE=0
```

The current executable reports the BDDC/CG iteration count. The old MATLAB
plot contains both Steepest Descent and Conjugate Gradient curves, but the
current driver does not expose the historical SD run as a separate command;
therefore the new data must label the measured solver explicitly.

Run:

```bash
bash submit_sweep.sh
```

