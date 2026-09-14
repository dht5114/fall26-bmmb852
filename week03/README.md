# Week 03 - Peer review of a classmate's repository

## Repository reviewed

- Forked from: Emily Snyder
- Original repository: <https://github.com/EXS5825/appbio-2026>
- My fork: <https://github.com/dht5114/appbio-2026>
- Assignment reviewed: `week02` (obtain genomic data + visualize a genome)
- Pull request: <https://github.com/EXS5825/appbio-2026/pull/2>

## 1. Fork and clone

I forked the repository on GitHub and cloned my fork:

```bash
git clone https://github.com/dht5114/appbio-2026.git
cd appbio-2026/week02
```

## 2. Is the code doing anything dangerous?

No. The code in `week02/` is:

- `Makefile` — downloads the *Genlisea aurea* assembly (GCA_000441915.1)
  with the NCBI `datasets` CLI and extracts the FASTA and GFF with `unzip`.
  The only destructive command is `make clean`, which runs
  `rm -f` on the three files the Makefile itself created.
- `genomebrowser/Makefile` — despite its name, this file is Java source code
  for a Swing genome browser. It only reads the FASTA/GFF files given on the
  command line; it makes no network connections, runs no shell commands and
  writes/deletes no files.

## 3–4. Evaluation of the README

The README is well written and nicely motivates the choice of organism. It
answers every question and includes screenshots from the genome browser.
However, it is not fully clear how to run the code:

- The run instructions say `pixi run makefile`, but the repository has no
  `pixi.toml`, so this command does not work. The Makefile is run with
  `make`.
- The required tools (`datasets`, `unzip`) are not listed.
- The example output shown does not match what the Makefile actually prints
  (it shows a different `unzip` path and omits the GFF step).
- The `genomebrowser/README.md` says to run `make run`, but
  `genomebrowser/Makefile` contains Java code, so `make run` fails with
  `Makefile:1: *** missing separator. Stop.`
- The six reading frames question is answered only with a screenshot, with no
  written description of the frames.

The expected outcomes (file names and sizes) are clearly stated, which made
checking the results easy.

## 5. Are the results reproducible?

Yes, for the main data download. Running `make` inside the `bioinfo`
environment produced exactly the files and sizes reported in the README:

```
-rw-r--r--  21M  GCA_000441915.1.zip
-rw-r--r--  43M  GCA_000441915.1_genomic.fna
-rw-r--r-- 114M  GCA_000441915.1_genomic.gff
```

```bash
awk '$3=="gene"' GCA_000441915.1_genomic.gff | wc -l
# 17685   (matches the README)

grep -v ">" GCA_000441915.1_genomic.fna | tr -d '\n' | wc -c
# 43357795   (10,684 scaffolds)
```

Issues found while checking the numbers:

- The README says the genome is **63.36 Mb**, but the assembled sequence is
  **43.36 Mb**. Using the assembly length, gene density is
  17685 ÷ 43.36 ≈ **408 genes/Mb** (the README reports ≈279 genes/Mb).
- The README says the genome is "over 70 times smaller" than *Arabidopsis*
  (~135 Mb); 135 ÷ 43.4 is only about **3 times** smaller.
- The genome browser could not be reproduced with `make run` (see above). The
  Java code also uses switch expressions, so it needs Java 14 or newer, which
  is not documented.

## 6–7. Comparison with my solution (by the AI agent)

| | Classmate's solution | My solution |
|---|---|---|
| Organism | *Genlisea aurea* (plant, 43 Mb, 10,684 scaffolds) | RSV-A (virus, 15.2 kb, 1 sequence) |
| Download tool | NCBI `datasets` CLI → zip → `unzip` | `curl` from the NCBI FTP site → `gunzip` |
| Dependencies | `datasets` must be installed (checked by the Makefile) | Only standard tools (`curl`, `gunzip`) |
| File layout | All files in the working directory | Organized into `fasta/` and `gff/` |
| Makefile targets | `all`, `clean`, `check_deps` | `all`, `fasta`, `gff`, `clean` |
| Run instructions | `pixi run makefile` (does not work) | `make` (works) |
| Answers | Mostly screenshots, short text | Commands + tables + written explanations |
| Extras | Custom Java genome browser | IGV screenshots and `igv.html` |

**Which is better?** The AI agent judged my solution to be the more
reproducible and better documented one: its run instructions work as
written, it needs no extra tools, the downloaded files are organized by type,
and every answer shows the command used. The classmate's solution has real
strengths — the `check_deps` target gives a helpful error if `datasets` is
missing, the `datasets` approach generalizes easily to any accession, and the
custom genome browser is an ambitious extra — but the incorrect run command,
the broken `make run`, and the genome size error make it harder to reproduce
and verify.

## 8. Summary of findings

The classmate's week 02 repository is safe to run and its core result is
reproducible: the Makefile downloads the *Genlisea aurea* assembly and
produces FASTA and GFF files of exactly the sizes reported, and the gene
count of 17,685 matches. The README is engaging and answers all questions,
and the Makefile is well structured, including a dependency check for the
`datasets` tool.

The main problems are with readability and documentation. The instructions
tell the reader to run `pixi run makefile`, which does not work because the
project has no pixi configuration; the required tools are not listed; and the
example output does not match the Makefile. The genome size used for the gene
density (63.36 Mb) does not match the assembled sequence (43.36 Mb), which
changes the density to about 408 genes/Mb and the *Arabidopsis* comparison to
about 3× rather than 70×. The bonus genome browser is stored as Java source in
a file named `Makefile`, so the documented `make run` fails.

## 9–10. Change made to the fork

I fixed the "How the Makefile should be used" section of `week02/README.md`:
replaced `pixi run makefile` with `make`, listed the required tools, documented
`make clean`, and updated the example output to match what the Makefile
actually prints.

```bash
git add week02/README.md
git commit -m "Fix run instructions in week02 README"
git push origin main
```

## 11. Pull request

Pull request to the original repository: <https://github.com/EXS5825/appbio-2026/pull/2>
