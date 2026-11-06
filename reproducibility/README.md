# Reproducibility package: Poisson BDDC compression

This directory records the numerical experiments used to reproduce the legacy
Poisson results from `Compressed_BDDC` and to regenerate the corresponding
figures for D3-3.

The source of truth for execution is the current cluster scratch project:

```text
/data/numerik/people/fchegini/project/fatemeh/MicroCard/current_papers/BDDC_petsc_SDC_2023/gatevariables/20Nov2023/BDDC/poisson_new_scratch
```

The scripts in this package submit jobs through the scratch project's
`run_case_mpi.sh`. Every submitted case is recorded in a CSV manifest, and
the application log, Slurm job ID, and extracted compression metrics are kept
under `logs/` and `results/`.

## Experiments

- `figure3_quantization/`: the legacy 2-D Poisson iteration-versus-bits
  study. It uses a fixed global 512 x 512 mesh, 256 subdomains, and the
  quantization-only pipeline. The current executable exposes the BDDC/CG
  iteration count; the historical SD curve is documented as a legacy result
  and is not silently reconstructed.
- `figure5_subdomains/`: the fixed-global-resolution sweep over 16, 64, 256,
  1024, and 4096 subdomains and quantization bits 4, 8, 12, 16, 24, 32, 48,
  and 64. The global mesh is kept at 512 x 512 by choosing
  `ELMRES=512/SUBRES`.

Both experiments deliberately start with quantization only (`LIFTING_MODE=0`,
`HUFFMAN=0`, `BITLENGTH=0`, `ZIGZAG=0`) so that they match the historical
figures. The production pipeline can be tested afterward by changing these
flags in the recorded manifest.

## Cluster workflow

From the cluster login node:

```bash
cd /data/numerik/people/fchegini/project/fatemeh/MicroCard/current_papers/BDDC_petsc_SDC_2023/gatevariables/20Nov2023/BDDC/poisson_new_scratch
build_jid=$(sbatch --parsable --partition=small --nodelist=htc-cmp023 build_poisson.sbatch)
echo "$build_jid"
squeue -j "$build_jid"
sacct -j "$build_jid" --format=JobID,State,ExitCode,Elapsed,NodeList
```

After a successful build, run either sweep script. The scripts submit one
Slurm job per parameter point and append the exact environment and command to
their `run_manifest.csv`.

Set `BUILD_JID` to make every experiment depend on successful completion of
the build job. If the build fails, the experiment jobs remain unstarted.

```bash
export REPORT_REPRO_DIR=/data/numerik/people/fchegini/project/fatemeh/MicroCard/current_papers/BDDC_petsc_SDC_2023/gatevariables/20Nov2023/BDDC/poisson_new_scratch/reproducibility
bash /path/to/reproducibility/figure3_quantization/submit_sweep.sh
bash /path/to/reproducibility/figure5_subdomains/submit_sweep.sh
```

Set `REPORT_REPRO_DIR` to the report checkout when running on a machine where
the report directory is writable. The cluster and the local report checkout
are separate filesystems, so manifests and logs created on the cluster should
be copied back into this report directory after the jobs finish.

The scripts default to the shared cluster path above. For local dry runs or a
different checkout, override `SCRATCH_DIR` and `REPORT_REPRO_DIR`.

## Interpretation

`MIN_INTERFACE_SIZE` is kept at zero in these historical reproductions. The
interface-type threshold experiments are separate sensitivity studies and
must not be mixed into the legacy Figure 3/5 curves. Use the threshold mask
only in a dedicated manifest, e.g. `THRESHOLD_INTERFACE_TYPES=2` for edges,
`4` for faces, `6` for edges plus faces, and `0` for all classified interfaces.

