Use the following workflow to reproduce Figure 7.

### 1. Go to the benchmark directory

cd /data/numerik/people/fchegini/project/fatemeh/MicroCard/current_papers/BDDC_petsc_SDC_2023/gatevariables/20Nov2023/BDDC/poisson_new_scratch

mkdir -p reproducibility/direct_logs

### 2. Run the quantization-bit sweep

This reproduces the measured curve using graph lifting, zigzag, bit-length coding, and Huffman coding:

for bits in 4 5 6 7 8 9 10 12 13 14 15 16 24 27 32 64; do
  echo "===== Running Figure 7: bits=$bits ====="

  srun \
    --mpi=pmi2 \
    --kill-on-bad-exit=0 \
    --partition=big \
    --nodelist=htc-cmp011 \
    --nodes=1 \
    --ntasks=8 \
    --ntasks-per-node=8 \
    --cpus-per-task=1 \
    --time=00:10:00 \
    --output="reproducibility/direct_logs/figure3-b${bits}-%t.out" \
    --error="reproducibility/direct_logs/figure3-b${bits}-%t.err" \
    bash -lc '
      bits="$1"
      eval "exec 3<&${PMI_FD}"

      podman run --rm --pull=never \
        --network=host \
        --ipc=host \
        --pid=host \
        --preserve-fds=1 \
        --env-host \
        -e PMI_FD=3 \
        -e PMI_RANK \
        -e PMI_SIZE \
        -e LD_LIBRARY_PATH=/petsc/lib \
        -v "$PWD:/BDDC" \
        -v /data/numerik/people/fchegini/kaskade7_new:/kaskade7 \
        -v /home/htc/fchegini/opencarp_dependencies/petsc_install:/petsc:ro \
        -w /BDDC \
        localhost/kaskade7-mpi-build-environment:latest \
        ./bddc \
          --mpi 1 \
          --compression 1 \
          --compressionBits "$bits" \
          --compressionMinInterfaceSize 0 \
          --compressionThresholdInterfaceTypes 4 \
          --graphLifting 1 \
          --graphLiftingAfterQuantization 1 \
          --huffman 1 \
          --bitlength 1 \
          --zigzag 1 \
          --threads 1 \
          --subres 16 \
          --elmres 32 \
          --iter 20 \
          --tol 1e-10 \
          --timing 0 \
          --vtk 0 \
          --histogram 0 \
          --sigma 1.0 \
          --interfacetypes 7 \
          --crisscross 1 \
          --bddcExcludeCornersFromCompression 1
    ' _ "$bits"

  rc=$?
  echo "bits=$bits exit_status=$rc"

  if [ "$rc" -ne 0 ]; then
    echo "Figure 7 run failed for bits=$bits" >&2
    break
  fi
done

### 3. Verify all compression reports

for bits in 4 5 6 7 8 9 10 12 13 14 15 16 24 27 32 64; do
  f="reproducibility/direct_logs/figure3-b${bits}-0.out"

  echo "===== bits=$bits ====="

  grep -E \
    'raw sent bytes|transmitted sent bytes|savings:|ratio:' \
    "$f"
done

### 4. Generate the CSV used by LaTeX

csv="reproducibility/results/figure3_compression_data.csv"

mkdir -p reproducibility/results

cat > "$csv" <<'EOF'
bits,savings_percent,compression_ratio,raw_sent_bytes,transmitted_sent_bytes,subres,elmres,global_elements,mpi_ranks,threads,compression,graph_lifting,graph_lifting_after_quantization,huffman,bitlength,zigzag,min_interface_size,threshold_interface_types,exclude_corners,source
EOF

for bits in 4 5 6 7 8 9 10 12 13 14 15 16 24 27 32 64; do
  f="reproducibility/direct_logs/figure3-b${bits}-0.out"

  raw=$(grep -m1 'raw sent bytes:' "$f" | sed -E 's/.*raw sent bytes: *//')
  transmitted=$(grep -m1 'transmitted sent bytes:' "$f" | sed -E 's/.*transmitted sent bytes: *//')
  savings=$(grep -m1 'savings:' "$f" | sed -E 's/.*savings: *//; s/%.*//')
  ratio=$(grep -m1 'ratio:' "$f" | sed -E 's/.*ratio: *//; s/:1.*//')

  printf '%s,%s,%s,%s,%s,16,32,262144,8,1,1,1,1,1,1,1,0,4,1,direct_logs/figure3-b%s-0.out\n' \
    "$bits" \
    "$savings" \
    "$ratio" \
    "$raw" \
    "$transmitted" \
    "$bits" >> "$csv"
done

Check the generated file:

column -s, -t reproducibility/results/figure3_compression_data.csv

### 5. Copy the CSV to the report directory

If the report is on your Mac, run this on the Mac:

scp \
  fchegini@htc-login03:/data/numerik/people/fchegini/project/fatemeh/MicroCard/current_papers/BDDC_petsc_SDC_2023/gatevariables/20Nov2023/BDDC/poisson_new_scratch/reproducibility/results/figure3_compression_data.csv \
  /Users/fchegini/paper/data_compression_report_D3_3/reproducibility/results/figure3_compression_data.csv

### 6. Compile the LaTeX report

On the Mac:

cd /Users/fchegini/paper/data_compression_report_D3_3

make clean
make

The generated PDF is:

/Users/fchegini/paper/data_compression_report_D3_3/D3-3.pdf

The LaTeX figure reads:

reproducibility/results/figure3_compression_data.csv

Therefore, changing the CSV changes Figure 7 automatically when the report is compiled.