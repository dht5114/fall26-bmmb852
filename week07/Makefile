# Week 07 - Generate a VCF file
#
# Extends the week05 pipeline: aligns Illumina reads from an RSV-A sample to
# the RSV-A reference genome with bwa mem, then calls variants with bcftools
# (mpileup -> call -> norm -> sort), marks low-quality calls, and summarizes
# the VCF with bcftools stats.
#
# Usage:
#   make                  # run all steps: genome fastq index align stats vcf vcfstats
#   make estimate         # show how N was chosen for the target coverage
#   make N=10000          # use a different number of read pairs
#   make SRR=... SAMPLE_NAME=...   # process another paired-end run
#   make clean            # remove the reads, alignments and variants
#   make realclean        # also remove the genome and its index
#
# Requirements: curl, fastq-dump (sra-tools), seqkit, bwa, samtools, bcftools

# URL of the reference genome (RSV-A, RefSeq NC_038235.1)
GENOME_URL = https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/815/475/GCF_002815475.1_ASM281547v1/GCF_002815475.1_ASM281547v1_genomic.fna.gz

# Name of the reference genome
GENOME_NAME = rsv_a

# Genome FASTA file
GENOME_FASTA = fasta/$(GENOME_NAME).fa

# SRR accession number (Illumina MiSeq, paired-end 2x151 bp)
SRR = SRR35979456

# Target coverage, used by 'make estimate' and 'make stats'
COVERAGE = 10

# Number of read pairs (spots) to download (see 'make estimate' and README)
N = 5000

# FASTQ files
FASTQ_1 = fastq/$(SRR)_1.fastq
FASTQ_2 = fastq/$(SRR)_2.fastq

# Sample name
SAMPLE_NAME = 2024RSV_10

# BAM file derived from the sample name
BAM = bam/$(SAMPLE_NAME).bam

# Read group, so IGV and later tools know which sample the reads are from
READ_GROUP = @RG\tID:$(SRR)\tSM:$(SAMPLE_NAME)\tPL:ILLUMINA

# BWA index sentinel
INDEX = $(GENOME_FASTA).bwt

# Number of threads for bwa
THREADS = 4

# VCF file of all variant calls (low-quality calls are marked in FILTER)
VCF = vcf/$(SAMPLE_NAME).vcf.gz

# VCF statistics report
VCF_STATS = vcf/$(SAMPLE_NAME).stats.txt

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

# --- All definitions above this line. All code below this line ---

# Run all steps.
all: genome fastq index align stats vcf vcfstats

# Download the genome.
genome: $(GENOME_FASTA)

$(GENOME_FASTA):
	mkdir -p $(dir $@)
	curl -fsSL $(GENOME_URL) | gunzip > $@
	samtools faidx $@
	seqkit stats $@

# Download the first N read pairs.
fastq: $(FASTQ_1) $(FASTQ_2)

$(FASTQ_1):
	mkdir -p $(dir $@)
	fastq-dump -X $(N) --split-files -O $(dir $@) $(SRR)
	seqkit stats $(FASTQ_1) $(FASTQ_2)

$(FASTQ_2): $(FASTQ_1) ;

# How many read pairs are needed for the target coverage?
#   N = coverage * genome size / (2 * read length)
estimate: $(GENOME_FASTA)
	@SIZE=$$(seqkit fx2tab -n -l $< | awk '{s+=$$NF} END {print s}'); \
	LEN=$$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&rettype=runinfo&retmode=text&id=$(SRR)" \
	       | awk -F, 'NR==2 {print $$7}'); \
	echo "Genome size          $$SIZE bp"; \
	echo "Read pair length     $$LEN bp (from SRA metadata)"; \
	echo "Target coverage      $(COVERAGE)x"; \
	echo "N (read pairs)       $$(( ($(COVERAGE) * SIZE + LEN - 1) / LEN ))"

# Index the genome.
index: $(INDEX)

$(INDEX): $(GENOME_FASTA)
	bwa index $<

# Align the reads to the genome, then sort and index the BAM file.
align: $(BAM)

$(BAM): $(INDEX) $(FASTQ_1) $(FASTQ_2)
	mkdir -p $(dir $@)
	bwa mem -t $(THREADS) -R '$(READ_GROUP)' $(GENOME_FASTA) $(FASTQ_1) $(FASTQ_2) \
	    | samtools sort -o $@
	samtools index $@

# Alignment statistics.
stats: $(BAM)
	samtools flagstat $<
	samtools coverage $<
	@echo "Bases covered >= $(COVERAGE)x:"
	@samtools depth -a $< \
	    | awk '{n++; if ($$3 >= $(COVERAGE)) k++} END {printf "  %d of %d (%.1f%%)\n", k, n, 100*k/n}'

# Call variants: genotype likelihoods -> variant calls -> normalize -> sort,
# then mark (not remove) low-quality calls.
vcf: $(VCF)

$(VCF): $(BAM) $(GENOME_FASTA)
	mkdir -p $(dir $@)
	bcftools mpileup $(ANN) -O u -f $(GENOME_FASTA) $(BAM) \
	    | bcftools call $(CALL) -mv -O u \
	    | bcftools norm -f $(GENOME_FASTA) -d all -O u \
	    | bcftools filter -s LowQual -e '$(FILTER)' -O u \
	    | bcftools sort --write-index -O z -o $@

# Statistics report on the VCF file.
vcfstats: $(VCF_STATS)

$(VCF_STATS): $(VCF)
	bcftools stats -F $(GENOME_FASTA) -s - $< > $@
	@grep '^SN' $@ | cut -f 3-
	@echo "Calls passing the filter: $$(bcftools view -H -f PASS $< | wc -l)"

# Remove the reads, the alignments and the variants.
clean:
	rm -rf fastq bam vcf

# Remove everything that was generated.
realclean: clean
	rm -rf fasta

.PHONY: all genome fastq estimate index align stats vcf vcfstats clean realclean
