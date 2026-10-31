# How it works

## What a count is
For each sample, seqcounter runs [Seal](https://sourceforge.net/projects/bbmap/) once over the input, with the
reference file as the reference set. Seal breaks every reference sequence into k-mers, then reads the input and
counts, for each reference, **how many reads contain one of its k-mers**.

So a count is a number of **reads**, or of **contigs** for a fasta file — not a number of copies:

* A reference present twice in one contig of an assembly is counted **once**: one contig contains it.
* In reads, the same reference is covered by about twice as many reads as a reference present once, so the counts
  are roughly proportional to the number of copies, with the depth and the read length as the scale.
* A read containing two different references is counted once for **each** of them (`ambiguous=all`,
  `clearzone=999999`): seqcounter counts every reference a read matches, not only its best match.
* Both mates of a paired-end sample are read, and each counts separately.

A reference is found on **either strand** (`rcomp=t`), and with up to `-m` mismatches per k-mer (a Hamming
distance, so substitutions, not insertions or deletions). A reference **shorter than `-k`** has no k-mer and cannot
be found at all, which is why seqcounter refuses to run on one.

## The Seal command
With `--verbose`, the exact command is printed. It looks like this:

```
seal.sh -Xmx1024m in=sample_R1.fastq.gz in2=sample_R2.fastq.gz ref=refs.fasta k=25 rcomp=t hdist=1 \
        maskmiddle=f clearzone=999999 ambiguous=all nzo=f qin=33 ow=t stats=/tmp/seqcounter_xxx/stats.tsv threads=8
```

| Setting | Why |
| --- | --- |
| `k=25`, `hdist=1` | `-k` and `-m`: the k-mer size and the mismatches allowed per k-mer |
| `rcomp=t` | Count a match on either strand |
| `maskmiddle=f` | Do not treat the middle base of a k-mer as a wildcard, which is Seal's default and would allow one more mismatch than asked for |
| `clearzone=999999`, `ambiguous=all` | Count a read for every reference it matches, instead of assigning it to the best one. **Both are needed**: see below |
| `nzo=f` | Report references with no match, so every reference has a row |
| `qin=33` | Force the quality encoding: autodetection fails on some low quality nanopore reads |
| `int=f` (single input file) | Never guess that one fastq file is interleaved, which would halve the counts |
| `-Xmx` | The Java heap: see `--memory` in [Usage](Usage) |

Seal writes its stats to a temporary file, which seqcounter parses and removes; nothing but the reports is left in
the output folder.

### Why `clearzone` is not optional
Seal assigns each read to the reference it shares the most k-mers with, and `ambiguous=` only decides what happens
to reads it considers *ambiguous*. With the default `clearzone=0`, a read is ambiguous only when two references tie
exactly, so `ambiguous=all` on its own is not enough: a read holding an 80 bp reference and a 30 bp one is not a
tie, and the short reference gets nothing. Checked on one read containing both, with BBTools 39.01:

| Settings | 80 bp reference | 30 bp reference |
| --- | --- | --- |
| `ambiguous=all clearzone=999999` (seqcounter) | 1 | **1** |
| `ambiguous=all` alone (`clearzone=0`) | 1 | **0** |

`clearzone=999999` makes every reference sharing a k-mer with the read fall inside the ambiguity zone of the best
match, and `ambiguous=all` then counts the read for all of them. Lower either one and references that share a read
with a better match lose their counts silently. Seal's own default is `ambiguous=random`, which would scatter tied
reads between references at random.

## Before Seal runs
* **The reference file** is read and checked: it must be fasta, no sequence may be empty or shorter than `-k`, and
  no two sequences may share their first word (Seal reports one row per name, so their counts would be added
  together).
* **The samples** are resolved: one sample from `-r1`/`-r2`, or a recursive search of `-i` with R1/R2 pairing. Two
  files giving the same sample name stop the run.
* **The file formats** are checked against the content, because Seal recognizes fastq and fasta by extension only.
  A fastq file named `.fasta`, or a file with no known extension, is linked under a name matching what is actually
  in it.
* **The file names** are made safe. `seal.sh` builds its Java command line and runs it through the shell, and it
  splits `ref=` lists on commas; a path containing a space, a quote, `$`, a backquote, `;`, `&`, parentheses or a
  comma would be re-interpreted, and in BBTools 40 any argument containing `xmx` or `xms` is taken for a memory
  setting. Such paths are symlinked into a temporary folder under a neutral name (`in0.fastq.gz`), and the
  temporary folder itself is chosen so that its own path is safe.

## In parallel
With `-i` and `-p` greater than 1, `-p` samples are counted at the same time in threads of one process, each
running its own `seal.sh`. The `-t` threads are divided between them, and each job gets the same `--memory`, which
is why the automatic heap size is divided by the number of jobs. Results are collected in the order of the samples,
so the counts table does not depend on which sample finished first.

A sample whose Seal run fails keeps its error in its result instead of raising: the other samples are counted and
reported, and the exit code is 1.

## Why Seal
Counting k-mer matches over a fastq file is exactly what Seal is written for, and it is fast at it (a few seconds
for a bacterial run). Doing the same in Python would be tens of times slower for no gain in accuracy. The cost is
the dependency on BBTools and Java, which is why seqcounter checks for `seal.sh` before anything else and reports
the most informative line of Seal's output when it fails.
