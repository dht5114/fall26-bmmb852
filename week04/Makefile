# Week 04 - Get FASTQ from SRA
#
# Surveys the SRA for human respiratory syncytial virus (RSV) sequencing
# data, downloads the first N reads of one SRR run, and runs quality
# control (FastQC -> fastp trimming -> FastQC again).
#
# The Makefile is generic: change ACC to download a run from any
# sequencing platform. The layout (single/paired) and the platform are
# read from the SRA metadata, and the files are renamed after the sample.
#
# Usage:
#   make                        # download + QC the default run (Illumina)
#   make ACC=SRR39304124        # same for another run (Oxford Nanopore)
#   make ACC=SRR39304124 N=5000 # change the number of reads
#   make info                   # show what the metadata says about ACC
#   make stats                  # read statistics before/after trimming
#   make survey                 # summarize all SRA runs for the search
#   make clean                  # remove reads and reports for ACC
#   make realclean              # remove everything that was generated
#
# Requirements: curl, csvtk, fastq-dump (sra-tools), fastqc, fastp, seqkit

# SRR accession of the run to download
ACC ?= SRR35979456

# Number of reads (spots) to download
N ?= 10000

# Search term used on the SRA website
QUERY = human respiratory syncytial virus (RSV)

# NCBI E-utilities (called with plain curl; edirect gave SSL errors)
EUTILS = https://eutils.ncbi.nlm.nih.gov/entrez/eutils

SHELL := bash

# ---------------------------------------------------------------------
# Run metadata
# ---------------------------------------------------------------------

# Metadata (runinfo table) of the run, downloaded once and cached
META = metadata/$(ACC).csv

# The file names below depend on the metadata, so download it while the
# Makefile is being read (skipped for the survey and clean targets).
ifeq ($(filter survey clean realclean,$(MAKECMDGOALS)),)
$(shell mkdir -p metadata; test -s $(META) || \
    curl -s "$(EUTILS)/efetch.fcgi?db=sra&rettype=runinfo&retmode=text&id=$(ACC)" \
    | csvtk grep -f Run -p $(ACC) > $(META) 2>/dev/null)
endif

# Look up one column of the metadata by name
field = $(shell csvtk cut -f $(1) $(META) 2>/dev/null | sed -n 2p)

SAMPLE   := $(call field,SampleName)
PLATFORM := $(call field,Platform)
LAYOUT   := $(call field,LibraryLayout)

# Descriptive file name: sample name + platform, e.g. 2024RSV_10_illumina
# (uses the SRR number if the metadata has no sample name)
NAME ?= $(shell echo "$(or $(SAMPLE),$(ACC))_$$(echo $(PLATFORM) | tr A-Z a-z)" \
          | tr -c 'A-Za-z0-9._\n-' '_')

# ---------------------------------------------------------------------
# File names, organized by data type
# ---------------------------------------------------------------------

ifeq ($(LAYOUT),PAIRED)
  R1 = fastq/$(NAME)_1.fastq
  R2 = fastq/$(NAME)_2.fastq
  T1 = fastq/$(NAME)_1.trimmed.fastq
  T2 = fastq/$(NAME)_2.trimmed.fastq
  FASTP_IO = -i $(R1) -I $(R2) -o $(T1) -O $(T2) --detect_adapter_for_pe
else
  R1 = fastq/$(NAME).fastq
  T1 = fastq/$(NAME).trimmed.fastq
  FASTP_IO = -i $(R1) -o $(T1)
endif

READS   = $(R1) $(R2)
TRIMMED = $(T1) $(T2)

# FastQC reports (one per FASTQ file)
QC_RAW  = $(patsubst fastq/%.fastq,fastqc/%_fastqc.html,$(READS))
QC_TRIM = $(patsubst fastq/%.fastq,fastqc/%_fastqc.html,$(TRIMMED))

# ---------------------------------------------------------------------
# Trimming settings depend on the platform
# ---------------------------------------------------------------------

ifneq ($(filter OXFORD_NANOPORE PACBIO_SMRT,$(PLATFORM)),)
  # Long reads: a sliding window would cut every read at its first
  # low-quality stretch, so only trim the ends and then filter whole
  # reads by average quality and length.
  TRIM_OPTS = --cut_front --cut_tail --cut_mean_quality 10 \
              --average_qual 12 --length_required 200
else
  # Short reads: sliding-window quality trimming from 5' to 3'
  # (like Trimmomatic SLIDINGWINDOW:4:20) plus adapter removal.
  TRIM_OPTS = --cut_right --cut_right_window_size 4 \
              --cut_right_mean_quality 20 --length_required 50
endif

# ---------------------------------------------------------------------
# Targets
# ---------------------------------------------------------------------

all: fastq qc trim qc_trim

# Print what the metadata says about the run
info:
	@echo "Run       $(ACC)"
	@echo "Sample    $(SAMPLE)"
	@echo "Organism  $(call field,ScientificName)"
	@echo "Platform  $(PLATFORM) ($(call field,Model))"
	@echo "Strategy  $(call field,LibraryStrategy)"
	@echo "Layout    $(LAYOUT)"
	@echo "Spots     $(call field,spots)"
	@echo "Files     $(READS)"

# 1. Download the first N reads and rename them after the sample
fastq: $(READS)

$(R1):
	@test -s $(META) || { echo "No SRA metadata found for $(ACC)"; exit 1; }
	mkdir -p fastq
	fastq-dump -X $(N) --split-3 -O fastq $(ACC)
	for f in fastq/$(ACC)*.fastq; do mv $$f $${f/$(ACC)/$(NAME)}; done

# 2. QC report on the raw reads
qc: $(QC_RAW)

# 3. Trim the reads with fastp (also writes an HTML before/after report)
trim: $(TRIMMED)

$(T1): $(READS)
	mkdir -p fastp
	fastp $(FASTP_IO) $(TRIM_OPTS) \
	    -h fastp/$(NAME).html -j fastp/$(NAME).json -R "$(NAME) ($(ACC))"

# The second file of a pair is made together with the first one
ifeq ($(LAYOUT),PAIRED)
$(R2): $(R1) ;
$(T2): $(T1) ;
endif

# 4. QC report on the trimmed reads
qc_trim: $(QC_TRIM)

# FastQC report for any FASTQ file
fastqc/%_fastqc.html: fastq/%.fastq
	mkdir -p fastqc
	fastqc -q -o fastqc $<

# Read statistics before and after trimming
stats: $(READS) $(TRIMMED)
	seqkit stats -a $^

# ---------------------------------------------------------------------
# Survey of all SRA runs that match the search term
# ---------------------------------------------------------------------

SEARCH = metadata/search_runinfo.csv

survey: $(SEARCH)
	@echo "Runs:        $$(csvtk nrow $<)"
	@echo "BioProjects: $$(csvtk cut -f BioProject $< | sed 1d | sort -u | wc -l | tr -d " ")"
	@echo "BioSamples:  $$(csvtk cut -f BioSample $< | sed 1d | sort -u | wc -l | tr -d " ")"
	@echo "Gigabases:   $$(csvtk cut -f bases $< | sed 1d | awk '{s+=$$1} END {printf "%.0f", s/1e9}')"
	@for col in LibraryStrategy Platform Model ScientificName; do \
	    echo; echo "== $$col"; \
	    csvtk cut -f $$col $< | sed 1d | sort | uniq -c | sort -rn | head -10; \
	done
	@echo; echo "== Release year"
	@csvtk cut -f ReleaseDate $< | sed 1d | cut -c1-4 | sort | uniq -c

# Search the SRA, then download the runinfo of every hit in batches
$(SEARCH):
	mkdir -p metadata
	curl -s -G "$(EUTILS)/esearch.fcgi" -d db=sra -d retmax=100000 \
	    --data-urlencode "term=$(QUERY)" \
	    | grep -o '<Id>[0-9]*</Id>' | tr -dc '0-9\n' > metadata/search_ids.txt
	rm -f $@.tmp
	split -l 400 metadata/search_ids.txt metadata/batch_
	for b in metadata/batch_*; do \
	    curl -s -X POST "$(EUTILS)/efetch.fcgi" \
	        -d "db=sra&rettype=runinfo&retmode=text&id=$$(paste -sd, $$b)" \
	        | grep -v '^$$' >> $@.tmp; \
	    sleep 0.5; \
	done
	(head -1 $@.tmp; grep -v '^Run,' $@.tmp) > $@
	rm -f $@.tmp metadata/batch_*

# ---------------------------------------------------------------------

# Remove the reads and reports of the current run
clean:
	rm -f fastq/$(NAME)*.fastq fastqc/$(NAME)*_fastqc.* fastp/$(NAME).*

# Remove everything that was downloaded or generated
realclean:
	rm -rf fastq fastqc fastp metadata

.PHONY: all info fastq qc trim qc_trim stats survey clean realclean
