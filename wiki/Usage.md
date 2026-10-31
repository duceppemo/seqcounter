# Usage

```
seqcounter -ref FILE (-r1 FILE | -i FOLDER) -o FOLDER [options]
```

Two things are always needed: the **reference sequences** to count (`-ref`) and the **input** to count them in,
either one sample (`-r1`, with `-r2` for the second mate of a pair) or a folder of samples (`-i`).

```bash
seqcounter -ref refs.fasta -r1 sample_R1.fastq.gz -r2 sample_R2.fastq.gz -o results/  # Paired-end reads
seqcounter -ref refs.fasta -r1 sample.fastq.gz -o results/                            # Single-end reads
seqcounter -ref refs.fasta -r1 assembly.fasta -o results/                             # An assembly
seqcounter -ref refs.fasta -i samples/ -o results/ -p 4                               # A folder, 4 at a time
```

## The reference file
A fasta file, gzipped or not, holding the sequences to count:

```
>spacer01 Kamerbeek 1997
ATAGAGGGTCGCCGGTTCTGGATCA
>IS6110
AGATGCACCGATGCAGGCATAGTCG
```

* Each sequence gives **one row** of the counts table, named after the **first word** of its `>` line
  (`spacer01`, `IS6110`). The rest of the line is kept in the JSON report as the description.
* Two sequences whose first word is the same are **rejected**: Seal would add their counts together in one row.
* A sequence **shorter than the k-mer size** (25 by default) is rejected too: Seal cannot find it at all. Either
  lower `-k` to at most the length of the shortest reference, or remove it.
* There is no limit on how many sequences the file holds, but Seal keeps every reference k-mer in memory: see
  `--memory` below.

## Input
### One sample
`-r1` takes a fastq or fasta file, gzipped or not; `-r2` the second mate of a paired-end run. Both mates are
counted, once each. The sample name is the file name without its extension and without the read suffix
(`S1_R1.fastq.gz` → `S1`); `-s` overrides it.

Seal recognizes fastq and fasta by file extension, so a fastq file named `.fasta` (or a file with no known
extension, such as a Galaxy `dataset_1.dat`) would be misread. seqcounter looks at the first character of the
content and presents the file to Seal under a matching name, so these files are counted correctly.

### A folder of samples
`-i` searches a folder and its subfolders for fasta and fastq files. Hidden files and folders are skipped, symbolic
links to folders are followed (without looping), and the output folder is skipped, so running twice into a
subfolder of the input does not count the reports as samples.

Each sample is one fasta file, one single-end fastq file, or a pair of R1/R2 fastq files. Pairs are recognized from
the read suffix: `_R1`/`_R2`, `_1`/`_2`, and the Illumina form `_R1_001`/`_R2_001`. Everything before the suffix is
the sample name.

If two files give the same sample name (the same name in two subfolders, a fasta and a fastq file of the same name,
two R1 files...), the run stops with the list of the files concerned rather than merging or overwriting them.

## Options
| Option | Default | What it does |
| --- | --- | --- |
| `-ref`, `--ref FILE` | required | Fasta file of the reference sequences to count |
| `-r1`, `--r1 FILE` | | One sample: reads (single-end or R1) or an assembly |
| `-r2`, `--r2 FILE` | | R2 reads, for paired-end fastq files |
| `-i`, `--input FOLDER` | | A folder of samples, searched recursively |
| `-s`, `--sample NAME` | the file name | With `-r1`: the sample name used in the reports |
| `-o`, `--output FOLDER` | required | Folder to hold the reports; created if needed |
| `-k`, `--kmer-size N` | 25 | K-mer size used to find the references |
| `-m`, `--mismatch N` | 1 | Mismatches allowed per k-mer (Hamming distance) |
| `-t`, `--threads N` | all available | Threads, shared between the samples counted at the same time |
| `-p`, `-j`, `--parallel`, `--jobs N` | 1 | With `-i`: samples counted at the same time |
| `--memory SIZE` | from the reference size | Java heap given to Seal, per job, e.g. `2g` |
| `--verbose` | | Debug messages, including the exact Seal command |
| `-v`, `--version` | | Print the version |
| `-h`, `--help` | | Print the help |

### `-k` and `-m`: how exact the match has to be
A read is counted for a reference when it contains a k-mer of that reference, on either strand, with at most `-m`
mismatches. Raising `-k` or lowering `-m` makes the match stricter; `-m 0` requires an exact k-mer.

`-k` must be at most the length of the shortest reference sequence, and `-m` smaller than `-k`. A short `-k` (below
about 20) starts matching unrelated sequence by chance.

### `-t` and `-p`: threads and samples
With `-i`, `-p` samples are counted at the same time and the `-t` threads are shared between them (`-t 16 -p 4`
gives each sample 4 threads). Counting a few large fastq files is usually fastest with `-p 1` and all the threads;
many small samples go faster with a higher `-p`.

### `--memory`: the Java heap
Seal holds every reference k-mer in memory. Without `--memory`, seqcounter sets the heap from the total length of
the reference sequences (1 GB at least, about 40 MB per megabase of reference), capped at 85% of the machine's
memory divided by the number of jobs. This is what keeps `-p` usable: `seal.sh` left to itself reserves most of the
free memory **per job**. Set `--memory` yourself for a very large reference set, or when counting on a shared
machine.

## Exit codes
| Code | Meaning |
| --- | --- |
| 0 | Every sample was counted |
| 1 | The run could not start (bad reference file, no input, no `seal.sh`...), or at least one sample failed |
| 2 | The command line is wrong |
| 130 | Interrupted (Ctrl-C) |

A sample that fails does not stop the others: its column is left empty in the counts table, its error is written in
the per-sample report and the JSON, and the exit code is 1.

## Coming from the first version
Every flag of the first version (`-i`, `-r1`, `-r2`, `-ref`, `-k`, `-m`, `-o`, `-t`, `-p`, `-v`) still means the
same thing, `-v` still prints the version, and `refseqcounts.tsv` is still the counts table, with the command line
on its first line. `python seqcounter.py` still works from a clone. See the [Changelog](Changelog) for what changed
around that.
