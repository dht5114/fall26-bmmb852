# Week 01 - Setting up the computer and learning the command line

I set up my computer following the Biostar Handbook instructions and created
the `bioinfo` conda environment. For my editor I use VS Code
with the Claude Code extension.

## 1. samtools version in the `bioinfo` environment

```bash
conda activate bioinfo
samtools --version
```

Output (first lines):

```
samtools 1.24
Using htslib 1.24
Copyright (C) 2026 Genome Research Ltd.

Samtools compilation details:
    Features:       build=configure curses=yes 
    CC:             x86_64-apple-darwin13.4.0-clang -std=gnu23
    CPPFLAGS:       -D_FORTIFY_SOURCE=2 -isystem /Users/dorotheatse/micromamba/envs/bioinfo/include -mmacosx-version-min=11.3 -mmacosx-version-min=10.13
    CFLAGS:         -Wall -march=core2 -mtune=haswell -mssse3 -ftree-vectorize -fPIC -fstack-protector-strong -O2 -pipe -isystem /Users/dorotheatse/micromamba/envs/bioinfo/include -fdebug-prefix-map=/opt/mambaforge/envs/bioconda/conda-bld/samtools_1784063242737/work=/usr/local/src/conda/samtools-1.24 -fdebug-prefix-map=/Users/dorotheatse/micromamba/envs/bioinfo=/usr/local/src/conda-prefix
    LDFLAGS:        -Wl,-headerpad_max_install_names -Wl,-dead_strip_dylibs -Wl,-rpath,/Users/dorotheatse/micromamba/envs/bioinfo/lib -L/Users/dorotheatse/micromamba/envs/bioinfo/lib
    HTSDIR:         
    LIBS:           
    CURSES_LIB:     -ltinfow -lncursesw

HTSlib compilation details:
    Features:       build=configure libcurl=yes S3=yes GCS=yes libdeflate=yes lzma=yes bzip2=yes plugins=yes plugin-path=/Users/dorotheatse/micromamba/envs/bioinfo/libexec/htslib htscodecs=1.6.7
    CC:             x86_64-apple-darwin13.4.0-clang -std=gnu23
    CPPFLAGS:       -D_FORTIFY_SOURCE=2 -isystem /Users/dorotheatse/micromamba/envs/bioinfo/include -mmacosx-version-min=11.3 -mmacosx-version-min=10.13
    CFLAGS:         -Wall -march=core2 -mtune=haswell -mssse3 -ftree-vectorize -fPIC -fstack-protector-strong -O2 -pipe -isystem /Users/dorotheatse/micromamba/envs/bioinfo/include -fdebug-prefix-map=/opt/mambaforge/envs/bioconda/conda-bld/htslib_1783614739681/work=/usr/local/src/conda/htslib-1.24 -fdebug-prefix-map=/Users/dorotheatse/micromamba/envs/bioinfo=/usr/local/src/conda-prefix -fvisibility=hidden
    LDFLAGS:        -Wl,-headerpad_max_install_names -Wl,-dead_strip_dylibs -Wl,-rpath,/Users/dorotheatse/micromamba/envs/bioinfo/lib -L/Users/dorotheatse/micromamba/envs/bioinfo/lib -fvisibility=hidden -rdynamic

HTSlib URL scheme handlers present:
    built-in:	 file, preload, data
    Google Cloud Storage:	 gs+http, gs+https, gs
    libcurl:	 gophers, smtp, wss, rtsp, tftp, mqtts, pop3, imaps, pop3s, ws, ftps, ftp, gopher, imap, http, https, sftp, smtps, scp, dict, mqtt, telnet
    Amazon S3:	 s3+https, s3, s3+http
    crypt4gh-needed:	 crypt4gh
    mem:	 mem
```

## 2. Creating a nested directory structure

The `-p` flag makes `mkdir` create all intermediate directories at once:

```bash
mkdir -p week01/project/data/raw
mkdir -p week01/project/results
```

## 3. Creating files in different directories

`echo` with `>` redirects text into a file, creating the file if it does not
exist:

```bash
echo "raw sequencing data goes here" > week01/project/data/raw/reads.txt
echo "analysis results go here" > week01/project/results/summary.txt
```

Verifying the structure with a recursive listing:

```bash
ls -R week01/project
```

```
data    results

week01/project/data:
raw

week01/project/data/raw:
reads.txt

week01/project/results:
summary.txt
```

## 4. Accessing files with relative and absolute paths

A relative path starts from the current directory; an absolute path starts
from `/` and works from anywhere. First I moved into the project directory
and checked where I was:

```bash
cd week01/project
pwd
```

```
/Users/dorotheatse/fall26_bmmb852/fall26-bmmb852/week01/project
```

Reading the same file both ways - relative to `week01/project`:

```bash
cat data/raw/reads.txt
```

```
raw sequencing data goes here
```

and with the absolute path:

```bash
cat /Users/dorotheatse/fall26_bmmb852/fall26-bmmb852/week01/project/data/raw/reads.txt
```

```
raw sequencing data goes here
```

Going back up two levels with the relative path `..` (parent directory):

```bash
cd ../..
```
