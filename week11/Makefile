# RNA-Seq Analysis with the HiSat2 aligner

# The location for the data
DATA_URL=http://data.biostarhandbook.com/data/rnaseq-tutorial-2026.tar.gz

# The reference genome file
FASTA=fasta/chr22.genome.fa

# The index file sentinel
FASTA_IDX = ${FASTA}.1.ht2

# The annotation file
GTF=gtf/chr22.gtf

# The sample name
SAMPLE = HBR_1

# The path to the read file
R1=fastq/${SAMPLE}_R1.fq

# Output alignment file
BAM = bam/${SAMPLE}.bam

# Output coverage file
BW  = bw/${SAMPLE}.bw

# Count files as text
COUNT = counts/${SAMPLE}.txt

# Design file
DESIGN = design.csv

# Current directory
DIR = $(shell pwd)

# Netcat command to send commands to IGV
NC = nc localhost 60151

# The design file
DESIGN = design.csv

# Run two jobs at a time since the alignment is already multithreadedwith 4 threads, thus using 2x4=8 cores
PARALLEL = parallel --eta --colsep , --header : --jobs 2 --lb

# --- All definitions above this line. All code below this line ---

# First target is the help.
help:
	@echo "# "
	@echo "# SAMPLE=${SAMPLE}"
	@echo "# FASTA=${FASTA}"
	@echo "# R1=${R1}"
	@echo "# BAM=${BAM}"
	@echo "# BW=${BW}"
	@echo "# COUNT=${COUNT}"
	@echo "# "
	@echo "# make data"
	@echo "# make index"
	@echo "# make align"
	@echo "# make count"
	@echo "# "
	@echo "# make all"
	@echo "# "


# Download and unpack the data
${FASTA} ${GTF} ${R1} &:
	curl ${DATA_URL} | tar xzvf -

data: ${FASTA} ${GTF} ${R1}
	ls -lh ${FASTA} ${GTF} ${R1}

# The index file
${FASTA_IDX}: ${FASTA}
	hisat2-build ${FASTA} ${FASTA}

# Index the genome.
index: ${FASTA_IDX}
	ls -lh ${FASTA_IDX}

# The rules to generate the BAM file
${BAM}: ${FASTA_IDX} ${R1}
	mkdir -p $(dir ${BAM})
	hisat2 --threads 4 --sensitive -x ${FASTA} -U ${R1} | \
		samtools sort --write-index -o ${BAM}

# The rules to generate the bigwig file
${BW}: ${BAM}
	mkdir -p $(dir ${BW})
	bamCoverage --binSize 1 -b ${BAM} -o ${BW}

align: ${BAM} ${BW}
	ls -lh ${BAM}
	ls -lh ${BW}

${COUNT}: ${BW}
	mkdir -p $(dir ${COUNT})
	featureCounts -a ${GTF} -s 2 -o ${COUNT} ${BAM}

counts: ${COUNT}
	ls -lh ${COUNT}

all: counts

# Visualize with IGV 
igv: 
	echo "new" | ${NC}
	echo "genome ${DIR}/$(FASTA)" | ${NC}
	echo "load ${DIR}/$(GTF)" | ${NC}
	# Loads up the coverage files for each sample.
	cat $(DESIGN) | $(PARALLEL) 'echo "load ${DIR}/bw/{sample}.bw" | ${NC}'
	
# Generate all counts for every sample
batch:
	# Ensure that the index is built before running the alignments.
	make index
	
	# Run the alignments for each sample.
	cat $(DESIGN) | $(PARALLEL) make SAMPLE={sample} all
	

.PHONY: all help download index align bam bw counts

# Additional Makefile settings
SHELL := bash
.DELETE_ON_ERROR:
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
MAKEFLAGS += --warn-undefined-variables --no-print-directory
