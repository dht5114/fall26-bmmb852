# Week 10 - ChIP-seq peak calling
#
# Reuses the week09 workflow (design file + GNU parallel around a single-sample
# Makefile) on the yeast Sfp1 ChIP-seq samples of BioProject PRJNA488473, and
# replaces variant calling with peak calling:
# reads -> BAM -> coverage track, then ChIP BAM + input BAM -> MACS3 peaks.
#
# Usage:
#   make batch                # the whole workflow on every sample in design.csv
#   make dryrun               # print the per-sample commands without running them
#   make background           # process the input (background) samples
#   make chip                 # process the ChIP samples and call their peaks
#   make compare              # peak counts and replicate agreement in results/
#   make all SRR=... SAMPLE=...               # one sample: reads -> BAM -> coverage
#   make peaks SAMPLE=... CONTROL=...         # peaks of one ChIP sample vs its input
#   make clean                # remove reads, alignments, tracks, peaks and results
#   make realclean            # also remove the genome and its index
#
# Requirements: curl, fastq-dump (sra-tools), seqkit, bowtie2, samtools, bedtools,
#               bedGraphToBigWig, macs3, GNU parallel

# URL of the reference genome and its annotation (S. cerevisiae S288C, R64-1-1)
ENSEMBL_URL = https://ftp.ensembl.org/pub/release-116
FASTA_URL = $(ENSEMBL_URL)/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
GFF_URL = $(ENSEMBL_URL)/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.63.gff3.gz

# Name of the reference genome
GENOME_NAME = yeast

# Genome FASTA, its index and the GFF annotation
GENOME_FASTA = fasta/$(GENOME_NAME).fa
GENOME_SIZES = $(GENOME_FASTA).fai
GFF = fasta/$(GENOME_NAME).gff

# Bowtie2 index prefix and sentinel
INDEX = fasta/$(GENOME_NAME)
INDEX_FILE = $(INDEX).1.bt2

# Design file: run,sample,group,type,replicate,control (one line per SRA run)
DESIGN = design.csv

# Number of samples processed at the same time, and threads per sample
JOBS = 2
THREADS = 4

# The parallel command that runs the workflow for each row of the design file
PARALLEL = parallel --colsep , --header : --eta --lb -j $(JOBS)

# Rows of the design file of one type: header + the rows with that type
ROWS = awk -F, -v t=$(1) 'NR == 1 || $$4 == t' $(DESIGN)

# Default single sample (overridden on the command line or by 'make batch')
SRR = SRR7767712
SAMPLE = WT_ChIP_1
CONTROL = WT_input_1

# Number of read pairs (spots) to download per sample (2x25 bp reads)
N = 3000000

# FASTQ files
FASTQ_1 = fastq/$(SRR)_1.fastq
FASTQ_2 = fastq/$(SRR)_2.fastq

# BAM file and the bowtie2 mapping log, derived from the sample name
BAM = bam/$(SAMPLE).bam
BAM_LOG = bam/$(SAMPLE).log

# BAM file of the background sample
CONTROL_BAM = bam/$(CONTROL).bam

# Coverage track, normalized to fragments per million
BW = bw/$(SAMPLE).bw

# Keep concordant pairs that map confidently; drop unmapped, secondary,
# QC-failed, duplicate and supplementary alignments.
MIN_MAPQ = 20
BAM_FILTER = -q $(MIN_MAPQ) -f 2 -F 3844

# MACS3 parameters:
#   -f BAMPE  use the real fragments from the read pairs
#   -g        effective genome size of yeast
#   -q        false discovery rate cutoff
MACS = -f BAMPE -g 1.2e7 -q 0.01

# Narrow peak file (Sfp1 is a transcription factor)
PEAKS = peaks/$(SAMPLE)_peaks.narrowPeak

# Tables written by 'make compare'
COUNTS = results/peak_counts.tsv
PAIRS = results/replicate_pairs.tsv
SUPPORT = results/WT_peak_support.tsv

# --- All definitions above this line. All code below this line ---

# Show the usage.
usage:
	@sed -n '/^# Usage/,/^# Requirements/p' Makefile | sed '$$d'

# --- Per-sample workflow (parametrized by SRR and SAMPLE) ---

# One sample: genome, reads, alignment, coverage track.
all: genome index fastq align bigwig

# Download the genome and its annotation.
genome: $(GENOME_FASTA) $(GFF)

$(GENOME_FASTA):
	mkdir -p $(dir $@)
	curl -fsSL $(FASTA_URL) | gunzip > $@
	samtools faidx $@
	seqkit stats $@

$(GFF):
	mkdir -p $(dir $@)
	curl -fsSL $(GFF_URL) | gunzip > $@

# Index the genome.
index: $(INDEX_FILE)

$(INDEX_FILE): $(GENOME_FASTA)
	bowtie2-build -q $< $(INDEX)

# Download the first N read pairs.
fastq: $(FASTQ_1) $(FASTQ_2)

$(FASTQ_1):
	mkdir -p $(dir $@)
	fastq-dump -X $(N) --split-files -O $(dir $@) $(SRR)
	seqkit stats $(FASTQ_1) $(FASTQ_2)

$(FASTQ_2): $(FASTQ_1) ;

# Align the reads, keep the confident pairs, then sort and index the BAM file.
# bowtie2 instead of bwa mem: at its default score cutoff of 30 bwa mem
# discards 25 bp reads.
align: $(BAM)

$(BAM): $(INDEX_FILE) $(FASTQ_1) $(FASTQ_2)
	mkdir -p $(dir $@)
	bowtie2 -p $(THREADS) -x $(INDEX) -X 1000 --no-mixed --no-discordant \
	    --rg-id $(SRR) --rg SM:$(SAMPLE) -1 $(FASTQ_1) -2 $(FASTQ_2) 2> $(BAM_LOG) \
	    | samtools view -b $(BAM_FILTER) \
	    | samtools sort -o $@
	samtools index $@

# Alignment statistics for one sample.
stats: $(BAM)
	cat $(BAM_LOG)
	samtools flagstat $<

# Coverage track: fragment coverage scaled to fragments per million.
bigwig: $(BW)

$(BW): $(BAM) $(GENOME_FASTA)
	mkdir -p $(dir $@)
	bedtools genomecov -ibam $< -bg -pc \
	    -scale $$(samtools view -c -f 64 $< | awk '{print 1000000 / $$1}') \
	    | LC_ALL=C sort -k1,1 -k2,2n > $(BW:.bw=.bedGraph)
	bedGraphToBigWig $(BW:.bw=.bedGraph) $(GENOME_SIZES) $@
	rm -f $(BW:.bw=.bedGraph)

# Call peaks for a ChIP sample relative to its background sample.
peaks: $(PEAKS)

$(PEAKS): $(BAM) $(CONTROL_BAM)
	mkdir -p $(dir $@)
	macs3 callpeak $(MACS) -t $(BAM) -c $(CONTROL_BAM) -n $(SAMPLE) --outdir $(dir $@) 2> peaks/$(SAMPLE).log

# --- Batch processing using the design file ---

# Print the commands that batch would run.
dryrun:
	@$(call ROWS,input) | parallel --colsep , --header : -k --dry-run \
	    $(MAKE) all SRR={run} SAMPLE={sample} N=$(N)
	@$(call ROWS,ChIP) | parallel --colsep , --header : -k --dry-run \
	    $(MAKE) all peaks SRR={run} SAMPLE={sample} CONTROL={control} N=$(N)

# Process the background samples.
background: genome index
	$(call ROWS,input) | $(PARALLEL) \
	    $(MAKE) all SRR={run} SAMPLE={sample} N=$(N)

# Process the ChIP samples and call their peaks. Needs the background BAM files.
chip: genome index
	$(call ROWS,ChIP) | $(PARALLEL) \
	    $(MAKE) all peaks SRR={run} SAMPLE={sample} CONTROL={control} N=$(N)

# Run the whole workflow. The genome is prepared once and the background runs
# first, so parallel jobs never race on the index or on a control BAM file.
batch: background chip compare

# Peak counts and replicate agreement.
compare:
	mkdir -p results
	bash compare.sh $(DESIGN) $(GFF) results
	column -t -s $$'\t' $(COUNTS)

# Remove the reads, the alignments, the tracks, the peaks and the results.
clean:
	rm -rf fastq bam bw peaks results

# Remove everything that was generated.
realclean: clean
	rm -rf fasta

.PHONY: usage all genome index fastq align stats bigwig peaks dryrun \
        background chip batch compare clean realclean
