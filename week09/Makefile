# Week 09 - Create a workflow
#
# Automates the week07/08 variant calling pipeline (reads -> BAM -> VCF) for
# many samples of BioProject PRJNA1358921 (RSV surveillance, Sydney 2024).
# The samples are listed in design.csv; 'make batch' runs the per-sample
# pipeline on the first NSAMPLES samples of one subtype with GNU parallel,
# merges the per-sample VCFs into one VCF and annotates it with snpEff.
#
# Usage:
#   make batch                    # the whole workflow on the first 10 RSV-A samples
#   make batch NSAMPLES=4         # fewer (or more) samples
#   make dryrun                   # print the per-sample commands without running them
#   make all SRR=... SAMPLE=...   # one sample only (fastq -> bam -> vcf)
#   make merge                    # merge the per-sample VCFs into vcf/merged.vcf.gz
#   make joint                    # call variants on all BAM files together (for comparison)
#   make annotate                 # snpEff annotation of the merged VCF
#   make codons                   # table of coding variants with ref/alt codon and amino acid
#   make coverage                 # coverage summary of every BAM file
#   make design                   # rebuild design.csv from the SRA/ENA metadata (needs bio, csvtk)
#   make clean                    # remove reads, alignments, variants and annotations
#   make realclean                # also remove the genome, its index and the snpEff database
#
# Requirements: curl, fastq-dump (sra-tools), seqkit, bwa, samtools, bcftools,
#               snpEff, GNU parallel

# BioProject the samples come from
BIOPROJECT = PRJNA1358921

# URL of the reference genome and its annotation (RSV-A, RefSeq NC_038235.1)
NCBI_URL = https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/815/475/GCF_002815475.1_ASM281547v1/GCF_002815475.1_ASM281547v1
GENOME_URL = $(NCBI_URL)_genomic.fna.gz
GFF_URL = $(NCBI_URL)_genomic.gff.gz

# Name of the reference genome
GENOME_NAME = rsv_a

# Genome FASTA and GFF files
GENOME_FASTA = fasta/$(GENOME_NAME).fa
GFF = fasta/$(GENOME_NAME).gff

# Design file: run_accession,sample,subtype,read_count (one line per SRA run)
DESIGN = design.csv

# Only samples of this subtype are processed (the reference is RSV-A)
SUBTYPE = RSVA

# Number of samples to process
NSAMPLES = 10

# The rows of the design file to process: header + first NSAMPLES of SUBTYPE
SELECT = awk -F, 'NR == 1 || $$3 == "$(SUBTYPE)"' $(DESIGN) | head -n $$(($(NSAMPLES) + 1))

# Sample names and BAM/VCF files of the selected samples
SAMPLES = $(shell $(SELECT) | tail -n +2 | cut -d, -f2)
BAMS = $(SAMPLES:%=bam/%.bam)
VCFS = $(SAMPLES:%=vcf/%.vcf.gz)

# Number of samples processed at the same time, and bwa threads per sample
JOBS = 4
THREADS = 2

# The parallel command that runs the workflow for each row of the design file
PARALLEL = parallel --colsep , --header : --eta --lb -j $(JOBS)

# Default single sample (overridden on the command line or by 'make batch')
SRR = SRR35979456
SAMPLE = RSVA_10

# Number of read pairs (spots) to download per sample (~180x coverage)
N = 10000

# Minimum depth for a call to PASS, and used by 'make stats'
COVERAGE = 10

# FASTQ files
FASTQ_1 = fastq/$(SRR)_1.fastq
FASTQ_2 = fastq/$(SRR)_2.fastq

# BAM file derived from the sample name
BAM = bam/$(SAMPLE).bam

# Read group, so the merged VCF knows which reads belong to which sample
READ_GROUP = @RG\tID:$(SRR)\tSM:$(SAMPLE)\tPL:ILLUMINA

# BWA index sentinel
INDEX = $(GENOME_FASTA).bwt

# Single-sample VCF
VCF = vcf/$(SAMPLE).vcf.gz

# Merged VCF (per-sample VCFs combined), its annotated version,
# and the joint-called VCF (all BAM files called together)
MERGED = vcf/merged.vcf.gz
ANNOTATED = vcf/merged.ann.vcf.gz
JOINT = vcf/joint.vcf.gz

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

# Show the usage.
usage:
	@sed -n '/^# Usage/,/^# Requirements/p' Makefile | sed '$$d'

# --- Per-sample workflow (parametrized by SRR and SAMPLE) ---

# One sample: genome, reads, alignment, variants.
all: genome index fastq align vcf

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

# Call variants for one sample.
vcf: $(VCF)

$(VCF): $(BAM) $(GENOME_FASTA)
	mkdir -p $(dir $@)
	bcftools mpileup $(ANN) -O u -f $(GENOME_FASTA) $(BAM) \
	    | bcftools call $(CALL) -mv -O u \
	    | bcftools norm -f $(GENOME_FASTA) -d all -O u \
	    | bcftools filter -s LowQual -e '$(FILTER)' -O u \
	    | bcftools sort --write-index -O z -o $@

# --- Batch processing using the design file ---

# Rebuild the design file from the BioProject metadata. The sample name
# comes from the library name (RSVA/RSV_10 -> RSVA_10) so that it does not
# start with a number like the SRA alias (2024RSV_10).
design:
	mkdir -p metadata
	bio search -H --csv $(BIOPROJECT) > metadata/$(BIOPROJECT).csv
	csvtk cut -f run_accession,library_name,read_count metadata/$(BIOPROJECT).csv \
	    | awk -F, -v OFS=, 'NR == 1 {print "run_accession,sample,subtype,read_count"; next} \
	        {split($$2, a, "/"); n = a[2]; sub(/RSV_/, "", n); print $$1, a[1] "_" n, a[1], $$3}' \
	    | (read h; echo "$$h"; sort -t, -k3,3 -k2,2V) > $(DESIGN)

# Print the commands that batch would run.
dryrun:
	@$(SELECT) | parallel --colsep , --header : -k --dry-run \
	    $(MAKE) all SRR={run_accession} SAMPLE={sample} N=$(N)

# Run the whole workflow:
#   1. download the genome and build the index once (parallel jobs must not race on it)
#   2. run the per-sample workflow for every selected row of the design file in parallel
#   3. merge the VCFs and annotate the result
batch: genome index
	$(SELECT) | $(PARALLEL) \
	    $(MAKE) all SRR={run_accession} SAMPLE={sample} N=$(N)
	$(MAKE) merge joint annotate codons NSAMPLES=$(NSAMPLES)

# Coverage summary of every BAM file.
coverage:
	@printf "sample\tcoverage_pct\tmean_depth\n"
	@for s in $(SAMPLES); do \
	    samtools coverage -H bam/$$s.bam | awk -v s=$$s -v OFS='\t' '{print s, $$6, $$7}'; done

# Merge the per-sample VCFs into one multi-sample VCF. A sample that had no
# variant at a site gets a missing genotype (.): reference or no coverage.
# -F x: a site is PASS if it passed the filter in at least one sample
# (the default, -F +, marks it LowQual if any sample failed).
merge: $(MERGED)

$(MERGED): $(VCFS)
	mkdir -p $(dir $@)
	bcftools merge -F x -O z -o $@ --write-index $^
	bcftools stats $@ | grep '^SN' | cut -f 3-

# Call variants on all BAM files together, so every sample gets a genotype
# at every variant site (0 = reference, . = no reads).
joint: $(JOINT)

$(JOINT): $(BAMS) $(GENOME_FASTA)
	mkdir -p $(dir $@)
	bcftools mpileup $(ANN) -O u -f $(GENOME_FASTA) $(BAMS) \
	    | bcftools call $(CALL) -mv -O u \
	    | bcftools norm -f $(GENOME_FASTA) -d all -O u \
	    | bcftools filter -s LowQual -e '$(FILTER)' -O u \
	    | bcftools sort --write-index -O z -o $@
	bcftools stats $@ | grep '^SN' | cut -f 3-

# --- Variant effect prediction (from week08) ---

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

# Annotate the merged VCF.
annotate: $(ANNOTATED)

$(ANNOTATED): $(MERGED) $(SNPEFF_BIN)
	mkdir -p $(dir $(SNPEFF_HTML))
	snpEff ann -noLog -no-upstream -no-downstream -c $(SNPEFF_CONFIG) -s $(SNPEFF_HTML) -csvStats $(SNPEFF_CSV) $(SNPEFF_DB) $< \
	    | bcftools view --write-index -O z -o $@

# Table of the coding variants with ref/alt codons and genotypes (see codons.sh).
codons: $(CODONS)

$(CODONS): $(ANNOTATED)
	bash codons.sh $< $(GENOME_FASTA) > $@

# Remove the reads, the alignments, the variants and the annotations.
clean:
	rm -rf fastq bam vcf results

# Remove everything that was generated.
realclean: clean
	rm -rf fasta $(SNPEFF_DIR)

.PHONY: usage all genome index fastq align stats vcf design dryrun batch \
        coverage merge joint snpeff annotate codons clean realclean
