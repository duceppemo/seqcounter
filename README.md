# SeqCounter
## Description
Count occurrence of reference sequence(s) from NGS sequencing (fastq) and assembly (fasta) data. Designed to run on Linux. Tested on Ubuntu 18.

## Workflow
1. Single- or paired-end fastq files or fasta assembly are searched for the spoligo spacers using `Seal` from `BBtools`.
2. Reference sequence counts reported using `pandas`.
## Installation
1. Create and activate `seqcounter` virtual environment:
```commandline
# Create virtual environment
conda create -n seqcounter -y -c bioconda bbmap pandas psutil

# Activate virtual environment
conda activate seqcounter
```
2. Clone the `seqcounter` repository and make sure it runs smoothly:
```commandline
# Clone repo
git clone https://github.com/duceppemo/seqcounter

# Go into cloned repo
cd seqcounter

# Test seqcounter
python seqcounter.py -h
```
## Usage
```commandline
usage: python seqcounter.py [-h] [-i /path/to/folder/with/fastq/or/fasta/files] [-r1 /path/to/sample_R1.[fastq|fasta]] [-r2 /path/to/R2/fastq] -ref /path/to/ref.fasta [-k 25] [-m 1] -o /path/to/output/folder/ [-t 16] [-p 2] [-v]

Count occurrence of reference sequence(s) in fastq or fasta.

optional arguments:
  -h, --help            show this help message and exit
  -i /path/to/folder/with/fastq/or/fasta/files, --input /path/to/folder/with/fastq/or/fasta/files
                        Folder containing fastq or fasta file to search. Gzipped or not. Mandatory.
  -r1 /path/to/sample_R1.[fastq|fasta]
                        R1 fastq file from paired-end or single end sequencing data or fasta file. Gzipped or not. Mandatory.
  -r2 /path/to/R2/fastq
                        R2 fastq file from paired-end. Optional.
  -ref /path/to/ref.fasta
                        Fasta file containing reference sequence(s) to find and count occurrence. Optional.
  -k 25, --kmer-size 25
                        Kmer length used for finding reference sequences. reference sequences shorter than k will not be found. Default 25. Optional.
  -m 1, --mismatch 1    Maximum number of mismatches between the reference and the input sequences (Maximum Hamming distance for ref kmers). Default 1. Optional.
  -o /path/to/output/folder/, --output /path/to/output/folder/
                        Folder to hold the result files. Mandatory.
  -t 16, --threads 16   Number of threads. Default is maximum available(16). Optional
  -p 2, --parallel 2    Number of samples to process in parallel. Default is 2. Optional.
  -v, --version         show program's version number and exit

```
## Examples
1. Single-end input fastq file:
```commandline
python seqcounter.py \
    -r1 ERR1744454_1.fastq.gz \
    -o spoligotyper_output/
```
2. Paired-end fastq input file:
```commandline
python seqcounter.py \
    -r1 ERR2747598_R1.fastq.gz \
    -r2 ERR2747598_R2.fastq.gz \
    -o spoligotyper_output/
```
3. Fasta input file (genome assembly):
```commandline
python seqcounter.py \
    -r1 NC_002945.4.fasta \
    -o spoligotyper_output/
```
4. Input folder containing single-end fastq, paired-end fastq and fasta files:
```commandline
python seqcounter.py \
    -i input/folder/ \
    -o spoligotyper_output/
```

