#!/usr/bin/env bash
#
# Summarize the MACS3 peaks of the ChIP samples listed in a design file.
#
# Writes three tables to the output directory:
#
#   peak_counts.tsv      one row per ChIP sample: fragments, peaks, FRiP (fraction
#                        of fragments that overlap a peak of that sample)
#   replicate_pairs.tsv  one row per pair of replicates in a group: how many peaks
#                        of each are found in the other, and the Jaccard index
#   WT_peak_support.tsv  one row per merged WT peak: how many WT and sfp1KO
#                        replicates have a peak there, and the nearest gene
#
# Peaks in the rDNA repeat and the mitochondrial genome are flagged in
# WT_peak_support.tsv. They reflect repeats and copy number, not binding.
#
# Usage: bash compare.sh design.csv fasta/yeast.gff results

set -euo pipefail

DESIGN=$1
GFF=$2
OUT=$3

# Name and group of every ChIP sample, in design file order.
CHIP=$(awk -F, '$4 == "ChIP" {print $2 "\t" $3}' "$DESIGN")
GROUPS_LIST=$(echo "$CHIP" | cut -f2 | uniq)

peakfile() { echo "peaks/$1_peaks.narrowPeak"; }

# --- Peak counts and FRiP ---
printf "sample\tgroup\tfragments\tpeaks\tfrip_pct\n" > "$OUT/peak_counts.tsv"
echo "$CHIP" | while IFS=$'\t' read -r s g; do
    frags=$(samtools view -c -f 64 "bam/$s.bam")
    peaks=$(wc -l < "$(peakfile "$s")" | tr -d ' ')
    inpeaks=$(samtools view -c -f 64 -L "$(peakfile "$s")" "bam/$s.bam")
    awk -v s="$s" -v g="$g" -v f="$frags" -v p="$peaks" -v i="$inpeaks" \
        'BEGIN {printf "%s\t%s\t%d\t%d\t%.1f\n", s, g, f, p, 100 * i / f}'
done >> "$OUT/peak_counts.tsv"

# --- Pairwise agreement within each group ---
printf "group\tsample_a\tsample_b\tpeaks_a\tpeaks_b\ta_in_b_pct\tb_in_a_pct\tjaccard\n" > "$OUT/replicate_pairs.tsv"
for g in $GROUPS_LIST; do
    samples=($(echo "$CHIP" | awk -v g="$g" '$2 == g {print $1}'))
    for ((i = 0; i < ${#samples[@]}; i++)); do
        for ((j = i + 1; j < ${#samples[@]}; j++)); do
            a=$(peakfile "${samples[i]}")
            b=$(peakfile "${samples[j]}")
            na=$(wc -l < "$a" | tr -d ' ')
            nb=$(wc -l < "$b" | tr -d ' ')
            ab=$(bedtools intersect -u -a "$a" -b "$b" | wc -l | tr -d ' ')
            ba=$(bedtools intersect -u -a "$b" -b "$a" | wc -l | tr -d ' ')
            jac=$(bedtools jaccard -a <(sort -k1,1 -k2,2n "$a") -b <(sort -k1,1 -k2,2n "$b") | tail -1 | cut -f3)
            awk -v g="$g" -v sa="${samples[i]}" -v sb="${samples[j]}" -v na="$na" -v nb="$nb" \
                -v ab="$ab" -v ba="$ba" -v jac="$jac" \
                'BEGIN {printf "%s\t%s\t%s\t%d\t%d\t%.1f\t%.1f\t%.3f\n", g, sa, sb, na, nb, 100 * ab / na, 100 * ba / nb, jac}'
        done
    done
done >> "$OUT/replicate_pairs.tsv"

# --- Support for every merged WT peak ---
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

# Transcription start sites of the protein-coding genes: chrom, start, end,
# gene name (the systematic name if it has none), strand. Sorted like the peaks.
awk 'BEGIN {FS = OFS = "\t"} $3 == "gene" && $9 ~ /biotype=protein_coding/ {
        id = ""; name = ""
        n = split($9, a, ";")
        for (i = 1; i <= n; i++) {
            if (a[i] ~ /^gene_id=/) id = substr(a[i], 9)
            if (a[i] ~ /^Name=/) name = substr(a[i], 6)
        }
        tss = ($7 == "+") ? $4 - 1 : $5 - 1
        print $1, tss, tss + 1, (name == "" ? id : name), $7
     }' "$GFF" | LC_ALL=C sort -k1,1 -k2,2n > "$WORK/tss.bed"

wt=$(echo "$CHIP" | awk '$2 == "WT" {print $1}')
ko=$(echo "$CHIP" | awk '$2 != "WT" {print $1}')

for s in $wt; do cut -f1-3 "$(peakfile "$s")"; done \
    | LC_ALL=C sort -k1,1 -k2,2n \
    | bedtools merge > "$WORK/merged.bed"

# Number of samples in a list with a peak overlapping each merged region.
count_support() {
    local samples=$1
    cp "$WORK/merged.bed" "$WORK/acc.bed"
    for s in $samples; do
        bedtools intersect -c -a "$WORK/acc.bed" -b "$(peakfile "$s")" \
            | awk 'BEGIN {FS = OFS = "\t"} {if ($NF > 1) $NF = 1; print}' > "$WORK/acc2.bed"
        mv "$WORK/acc2.bed" "$WORK/acc.bed"
    done
    awk 'BEGIN {FS = "\t"} {n = 0; for (i = 4; i <= NF; i++) n += $i; print n}' "$WORK/acc.bed"
}

count_support "$wt" > "$WORK/wt.txt"
count_support "$ko" > "$WORK/ko.txt"

# Nearest transcription start site of each merged region (the first if tied).
bedtools closest -d -t first -a "$WORK/merged.bed" -b "$WORK/tss.bed" \
    | cut -f7,9 > "$WORK/nearest.txt"

printf "chrom\tstart\tend\twt_reps\tko_reps\tnearest_gene\tdistance\tflag\n" > "$OUT/WT_peak_support.tsv"
paste "$WORK/merged.bed" "$WORK/wt.txt" "$WORK/ko.txt" "$WORK/nearest.txt" \
    | awk 'BEGIN {FS = OFS = "\t"} {
          flag = "."
          if ($1 == "Mito") flag = "mitochondrial"
          if ($1 == "XII" && $3 > 451000 && $2 < 469000) flag = "rDNA"
          print $0, flag
      }' >> "$OUT/WT_peak_support.tsv"
