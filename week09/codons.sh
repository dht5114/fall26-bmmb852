#!/usr/bin/env bash
#
# Tabulate the coding SNPs in a snpEff-annotated VCF.
#
# For every PASS SNP inside a CDS, take the first (most severe) snpEff ANN
# entry, find the codon in the reference, swap in the alternate base, and
# translate both with the standard genetic code (RSV uses the standard code).
# All RSV genes are on the + strand of the RefSeq genome, so no reverse
# complement is needed.
#
# snpEff looks at one SNP at a time. When two or three SNPs fall in the same
# codon (and are carried by the same samples) the real codon has all of them,
# so the CODON_SNPS, COMBINED_CODON and COMBINED_AA columns apply every SNP
# in the codon together.
#
# Usage: bash codons.sh vcf/merged.ann.vcf.gz fasta/rsv_a.fa > results/coding_variants.tsv

set -euo pipefail

VCF=$1
FASTA=$2

# Sample names for the genotype columns.
SAMPLES=$(bcftools query -l "$VCF" | paste -sd '\t' -)

printf "POS\tREF\tALT\tGENE\tTRANSCRIPT\tEFFECT\tHGVS_C\tHGVS_P\tREF_CODON\tALT_CODON\tREF_AA\tALT_AA\tCODON_SNPS\tCOMBINED_CODON\tCOMBINED_AA\t%s\n" "$SAMPLES"

# One line per variant: POS REF ALT ANN GT...
bcftools query -i 'FILTER="PASS" && TYPE="snp"' -f '%POS\t%REF\t%ALT\t%INFO/ANN[\t%GT]\n' "$VCF" \
| awk -v fasta="$FASTA" '
# Mark the changed bases of a codon in uppercase, e.g. Aag -> Agg.
function mark(codon, ref,    i, out, c) {
    out = ""
    for (i = 1; i <= 3; i++) {
        c = substr(codon, i, 1)
        out = out (c == substr(ref, i, 1) ? tolower(c) : c)
    }
    return out
}
BEGIN {
    FS = OFS = "\t"

    # Read the (single-sequence) genome.
    while ((getline line < fasta) > 0)
        if (line !~ /^>/) seq = seq toupper(line)

    # Standard genetic code, codons in TCAG order.
    b = "TCAG"
    aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    n = 0
    for (i = 1; i <= 4; i++)
        for (j = 1; j <= 4; j++)
            for (k = 1; k <= 4; k++)
                code[substr(b, i, 1) substr(b, j, 1) substr(b, k, 1)] = substr(aa, ++n, 1)
}
{
    # First ANN entry: Allele|Annotation|Impact|Gene|GeneID|FeatureType|FeatureID|
    # Biotype|Rank|HGVS.c|HGVS.p|cDNA.pos|CDS.pos|AA.pos|Distance|Errors
    split($4, ann, ",")
    split(ann[1], f, "|")
    if (f[13] == "") next                       # not in a CDS

    split(f[13], cds, "/")
    offset = (cds[1] - 1) % 3                   # position of the SNP in its codon
    start = $1 - offset

    gts = ""
    for (i = 5; i <= NF; i++) gts = gts OFS $i

    # Remember the variant; codons are keyed by transcript, start and genotypes.
    m++
    pos[m] = $1; ref[m] = $2; alt[m] = $3; off[m] = offset; gt[m] = gts
    info[m] = f[4] OFS f[7] OFS f[2] OFS f[10] OFS f[11]
    key[m] = f[7] SUBSEP start SUBSEP gts
    st[m] = start

    # Combined codon: start from the reference and apply every SNP in it.
    if (!(key[m] in combined)) combined[key[m]] = substr(seq, start, 3)
    c = combined[key[m]]
    combined[key[m]] = substr(c, 1, offset) $3 substr(c, offset + 2)
    count[key[m]]++
}
END {
    for (i = 1; i <= m; i++) {
        ref_codon = substr(seq, st[i], 3)
        alt_codon = substr(ref_codon, 1, off[i]) alt[i] substr(ref_codon, off[i] + 2)
        comb = combined[key[i]]
        print pos[i], ref[i], alt[i], info[i], \
              mark(ref_codon, alt_codon), mark(alt_codon, ref_codon), \
              code[ref_codon], code[alt_codon], \
              count[key[i]], mark(comb, ref_codon), code[comb] gt[i]
    }
}'
