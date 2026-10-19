# Week 08 - Predict variant effects
#
# Extends the week07 pipeline (reads -> BAM -> VCF) to several RSV-A samples:
# each sample is aligned separately, then variants are called on all BAM
# files at once (bcftools mpileup with every BAM) to make a multi-sample VCF.
# The variants are annotated with snpEff using a database built from the
# RefSeq genome and its GFF, and the coding changes are tabulated with the
# reference and alternate codons.
#
# Usage:
#   make all              # every step: genome, all samples, multi-sample VCF, annotation
#   make samples          # align and call variants for every sample in samples.csv
#   make sample SRR=... SAMPLE_NAME=...   # one sample only (fastq -> bam -> vcf)
#   make merge            # call variants on all BAM files together
#   make annotate         # build the snpEff database and annotate the merged VCF
#   make codons           # table of coding variants with ref/alt codon and amino acid
#   make clean            # remove reads, alignments, variants and annotations
#   make realclean        # also remove the genome, its index and the snpEff database
#
# Requirements: curl, fastq-dump (sra-tools), seqkit, bwa, samtools, bcftools,
#               snpEff, GNU parallel

# URL of the reference genome and its annotation (RSV-A, RefSeq NC_038235.1)
NCBI_URL = https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/815/475/GCF_002815475.1_ASM281547v1/GCF_002815475.1_ASM281547v1
GENOME_URL = $(NCBI_URL)_genomic.fna.gz
GFF_URL = $(NCBI_URL)_genomic.gff.gz

# Name of the reference genome
GENOME_NAME = rsv_a

# Genome FASTA and GFF files
GENOME_FASTA = fasta/$(GENOME_NAME).fa
GFF = fasta/$(GENOME_NAME).gff

# Design file: one line per sample (sample,run)
DESIGN = samples.csv

# Sample names and BAM files from the design file
SAMPLES = $(shell tail -n +2 $(DESIGN) | cut -d, -f1)
BAMS = $(SAMPLES:%=bam/%.bam)

# Default single sample (overridden on the command line or by 'make samples')
SRR = SRR35979456
SAMPLE_NAME = 2024RSV_10

# Target coverage, used by 'make stats' and the filter
COVERAGE = 10

# Number of read pairs (spots) to download per sample
N = 5000

# FASTQ files
FASTQ_1 = fastq/$(SRR)_1.fastq
FASTQ_2 = fastq/$(SRR)_2.fastq

# BAM file derived from the sample name
BAM = bam/$(SAMPLE_NAME).bam

# Read group, so the multi-sample VCF knows which reads belong to which sample
READ_GROUP = @RG\tID:$(SRR)\tSM:$(SAMPLE_NAME)\tPL:ILLUMINA

# BWA index sentinel
INDEX = $(GENOME_FASTA).bwt

# Number of threads for bwa, and number of samples processed in parallel
THREADS = 4
JOBS = 2

# Single-sample VCF
VCF = vcf/$(SAMPLE_NAME).vcf.gz

# Multi-sample VCF (all BAM files called together) and its annotated version
MERGED = vcf/merged.vcf.gz
ANNOTATED = vcf/merged.ann.vcf.gz

# Annotation parameters for mpileup:
#   -d      max reads per position (default 250); raised so deep regions are not downsampled
#   -a      extra per-sample and per-site tags (allele depths, strand-specific depths, strand bias)
ANN = -d 1000 --annotate 'INFO/AD,FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/SP'

# Call parameters for call:
#   --ploidy 1  RSV is a haploid (single-stranded RNA) virus
#   -a GQ       add genotype quality
CALL = --ploidy 1 --annotate 'FORMAT/GQ'

# Filter expression: calls matching it are kept but marked LowQual in FILTER
FILTER = QUAL<30 || INFO/DP<$(COVERAGE)

# snpEff database: name, directory, config file and the built database
SNPEFF_DB = $(GENOME_NAME)
SNPEFF_DIR = snpeff
SNPEFF_CONFIG = $(SNPEFF_DIR)/snpEff.config
SNPEFF_BIN = $(SNPEFF_DIR)/data/$(SNPEFF_DB)/snpEffectPredictor.bin

# snpEff reports
SNPEFF_HTML = results/snpEff_summary.html
SNPEFF_CSV = results/snpEff_summary.csv

# Table of coding variants with codons and amino acid changes
CODONS = results/coding_variants.tsv

# --- All definitions above this line. All code below this line ---

# Run all steps.
all: genome index samples merge annotate codons

# Download the genome and its annotation.
genome: $(GENOME_FASTA) $(GFF)

$(GENOME_FASTA):
	mkdir -p $(dir $@)
	curl -fsSL $(GENOME_URL) | gunzip > $@
	samtools faidx $@
	seqkit stats $@

$(GFF):
	mkdir -p $(dir $@)
	curl -fsSL $(GFF_URL) | gunzip > $@

# Index the genome.
index: $(INDEX)

$(INDEX): $(GENOME_FASTA)
	bwa index $<

# --- Per-sample steps (parametrized by SRR and SAMPLE_NAME) ---

# One sample: download reads, align, call variants on its own.
sample: fastq align vcf

# Download the first N read pairs.
fastq: $(FASTQ_1) $(FASTQ_2)

$(FASTQ_1):
	mkdir -p $(dir $@)
	fastq-dump -X $(N) --split-files -O $(dir $@) $(SRR)
	seqkit stats $(FASTQ_1) $(FASTQ_2)

$(FASTQ_2): $(FASTQ_1) ;

# Align the reads to the genome, then sort and index the BAM file.
align: $(BAM)

$(BAM): $(INDEX) $(FASTQ_1) $(FASTQ_2)
	mkdir -p $(dir $@)
	bwa mem -t $(THREADS) -R '$(READ_GROUP)' $(GENOME_FASTA) $(FASTQ_1) $(FASTQ_2) \
	    | samtools sort -o $@
	samtools index $@

# Alignment statistics for one sample.
stats: $(BAM)
	samtools flagstat $<
	samtools coverage $<
	@echo "Bases covered >= $(COVERAGE)x:"
	@samtools depth -a $< \
	    | awk '{n++; if ($$3 >= $(COVERAGE)) k++} END {printf "  %d of %d (%.1f%%)\n", k, n, 100*k/n}'

# Call variants for one sample (same pipeline as week07).
vcf: $(VCF)

$(VCF): $(BAM) $(GENOME_FASTA)
	mkdir -p $(dir $@)
	bcftools mpileup $(ANN) -O u -f $(GENOME_FASTA) $(BAM) \
	    | bcftools call $(CALL) -mv -O u \
	    | bcftools norm -f $(GENOME_FASTA) -d all -O u \
	    | bcftools filter -s LowQual -e '$(FILTER)' -O u \
	    | bcftools sort --write-index -O z -o $@

# --- All samples ---

# Run the per-sample steps for every sample in the design file.
samples: $(INDEX)
	cat $(DESIGN) | parallel --colsep , --header : -j $(JOBS) \
	    $(MAKE) sample SRR={run} SAMPLE_NAME={sample}

# Coverage summary of every BAM file.
coverage:
	@for b in $(BAMS); do samtools coverage -H $$b | sed "s|^|$$b\t|"; done

# Call variants on all BAM files together, so every sample gets a genotype
# at every variant site (no "holes" as with bcftools merge).
merge: $(MERGED)

$(MERGED): $(BAMS) $(GENOME_FASTA)
	mkdir -p $(dir $@)
	bcftools mpileup $(ANN) -O u -f $(GENOME_FASTA) $(BAMS) \
	    | bcftools call $(CALL) -mv -O u \
	    | bcftools norm -f $(GENOME_FASTA) -d all -O u \
	    | bcftools filter -s LowQual -e '$(FILTER)' -O u \
	    | bcftools sort --write-index -O z -o $@
	bcftools stats $@ | grep '^SN' | cut -f 3-

# --- Variant effect prediction ---

# Build a snpEff database from the RefSeq FASTA and GFF.
#
# The RefSeq GFF attaches the two overlapping M2 ORFs (M2-1 and M2-2) directly
# to the gene, next to an mRNA that has no gene. snpEff then joins both ORFs
# into one bogus transcript. The awk step drops that orphan mRNA (and its exon)
# and gives each M2 ORF its own mRNA (rna-M2-1, rna-M2-2) and exon.
snpeff: $(SNPEFF_BIN)

$(SNPEFF_BIN): $(GENOME_FASTA) $(GFF)
	mkdir -p $(dir $@)
	awk 'BEGIN {FS = OFS = "\t"} \
	    $$9 ~ /^ID=rna-NC_|Parent=rna-NC_/ {next} \
	    $$3 == "CDS" && match($$9, /Parent=gene-[^;]+/) { \
	        gene = substr($$9, RSTART + 7, RLENGTH - 7); rna = "rna-M2-" (++orf); \
	        print $$1, $$2, "mRNA", $$4, $$5, ".", $$7, ".", "ID=" rna ";Parent=" gene ";gbkey=mRNA;gene=M2"; \
	        print $$1, $$2, "exon", $$4, $$5, ".", $$7, ".", "ID=exon-" rna ";Parent=" rna; \
	        sub(/Parent=gene-[^;]+/, "Parent=" rna, $$9) } \
	    {print}' $(GFF) > $(dir $@)genes.gff
	cp $(GENOME_FASTA) $(dir $@)sequences.fa
	printf "data.dir = ./data/\n$(SNPEFF_DB).genome : Human respiratory syncytial virus A\n" > $(SNPEFF_CONFIG)
	snpEff build -c $(SNPEFF_CONFIG) -gff3 -noLog -v $(SNPEFF_DB)

# Annotate the multi-sample VCF.
annotate: $(ANNOTATED)

$(ANNOTATED): $(MERGED) $(SNPEFF_BIN)
	mkdir -p $(dir $(SNPEFF_HTML))
	snpEff ann -noLog -no-upstream -no-downstream -c $(SNPEFF_CONFIG) -s $(SNPEFF_HTML) -csvStats $(SNPEFF_CSV) $(SNPEFF_DB) $< \
	    | bcftools view --write-index -O z -o $@

# Table of the coding variants: the first (most severe) snpEff annotation of
# each PASS SNP, the reference and alternate codon, the codon with every SNP
# in it applied together, and each sample's genotype (see codons.sh).
codons: $(CODONS)

$(CODONS): $(ANNOTATED)
	bash codons.sh $< $(GENOME_FASTA) > $@
	column -t -s $$'\t' $@ | head -40

# Remove the reads, the alignments, the variants and the annotations.
clean:
	rm -rf fastq bam vcf results

# Remove everything that was generated.
realclean: clean
	rm -rf fasta $(SNPEFF_DIR)

.PHONY: all genome index sample fastq align stats vcf samples coverage merge \
        snpeff annotate codons clean realclean
