#!/usr/bin/env bash

VCF="41-Perognathus.X.multi_individual.filtered_missing.vcf.gz"
POP="pop_assignments.multi_individual.txt"
OUTDIR="X-pixy_results"
CORES=24
STATS=("dxy" "fst" "watterson_theta" "tajima_d")
SIZE="50000"

mkdir -p "$OUTDIR"

for STAT in "${STATS[@]}"; do
  echo "Running pixy for $STAT..."
  pixy --stats "$STAT" \
       --vcf "$VCF" \
       --populations "$POP" \
       --output_folder "$OUTDIR" \
       --n_cores "$CORES" \
       --include_multiallelic_snps \
       --fst_type hudson \
       --bypass_invariant_check \
       --window_size "$SIZE"
done

echo "All pixy runs complete. Outputs in $OUTDIR/"
