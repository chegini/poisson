# Figure 5: iteration count versus bits for increasing subdomain counts

Legacy source: `Compressed_BDDC/figure5.tex` and
`Compressed_BDDC/code/subdomains_vs_bits.m`.

The sweep uses subdomain counts 16, 64, 256, 1024, and 4096. They correspond
to `SUBRES=4, 8, 16, 32, 64` in the current 2-D driver. To preserve a fixed
global mesh of 512 x 512 elements, the local resolution is selected as
`ELMRES=512/SUBRES`.

As for Figure 3, the first reproduction is quantization only. This keeps the
comparison with the historical figure scientifically interpretable; graph
lifting and entropy coding should be added as a separate, explicitly labelled
production-pipeline study.

Run:

```bash
bash submit_sweep.sh
```

