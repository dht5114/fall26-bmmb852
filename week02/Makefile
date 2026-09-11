# Week 02 - Obtain genomic data
#
# Download the FASTA and GFF files for the human respiratory syncytial
# virus (RSV) reference genome (assembly GCF_002815475.1, ASM281547v1)
# from the NCBI FTP site.
#
# Usage:
#   make          # download everything (same as 'make all')
#   make fasta   # download only the FASTA file
#   make gff     # download only the GFF file
#   make clean   # remove the downloaded data

# Base URL of the assembly on the NCBI FTP site
URL = https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/815/475/GCF_002815475.1_ASM281547v1

# File name prefix used by NCBI for this assembly
NAME = GCF_002815475.1_ASM281547v1

# Local files, organized by data type
FASTA = fasta/$(NAME)_genomic.fna
GFF   = gff/$(NAME)_genomic.gff

all: fasta gff

fasta: $(FASTA)

gff: $(GFF)

# Download and uncompress the genome FASTA file
$(FASTA):
	mkdir -p fasta
	curl -fsSL -o $(FASTA).gz $(URL)/$(NAME)_genomic.fna.gz
	gunzip -f $(FASTA).gz

# Download and uncompress the annotation GFF file
$(GFF):
	mkdir -p gff
	curl -fsSL -o $(GFF).gz $(URL)/$(NAME)_genomic.gff.gz
	gunzip -f $(GFF).gz

clean:
	rm -rf fasta gff

.PHONY: all fasta gff clean