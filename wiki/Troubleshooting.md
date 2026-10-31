# Troubleshooting

## "seal.sh" was not found in your PATH
BBTools is not installed, or not in the environment you are running in:

```bash
conda install -c conda-forge -c bioconda bbmap
seal.sh --version
```

seqcounter also looks for a `seal.sh` next to the Python running it, so `/path/to/env/bin/seqcounter` works without
activating that environment. See [Installation](Installation).

## Every count is 0
In order of likelihood:

1. **The references are not there.** Check with a positive control: count the reference file in itself,
   `seqcounter -ref refs.fasta -r1 refs.fasta -o test_out`. Every reference should be found once. If that works,
   seqcounter and Seal are fine and the sequences really are absent from your samples.
2. **`-k` is too high for a divergent target.** A 25-mer has to match (with at most `-m` mismatches) end to end. On
   a reference that differs from your samples by a few percent, few 25-mers survive: try `-k 19 -m 1`.
3. **`-m 0` with real data.** One sequencing error in a k-mer is enough to lose it. The default `-m 1` is there for
   that reason.
4. **The wrong file.** The per-sample report gives `Reads` and `Bases`: if they are not what you expect, the input
   is not what you think it is.

## A count is far higher than expected
* A short `-k` (below about 20) matches unrelated sequence by chance, and a repeated element inside a reference
  (an IS element, a tandem repeat) matches everything that carries it elsewhere in the genome.
* Remember that a read matching two references is counted for both, so columns do not add up to `Matched`.
* With `--verbose`, check the Seal command: `hdist` and `k` are the two settings that change what matches.

## A sample failed
The error is in `refseqcounts_samples.tsv` and in the JSON report, and the run exits with 1 while still reporting
the samples that worked. The most common causes:

| Message | Cause |
| --- | --- |
| `java.io.IOException: Unexpected end of ZLIB input stream` | A truncated or corrupt gzipped file |
| `Exception ... OutOfMemoryError` | Not enough Java heap: raise `--memory`, or lower `-p` |
| `java.io.FileNotFoundException` | The file disappeared, or is not readable |
| `No reference sequences` | The reference file is empty or not fasta |

Run with `--verbose` to see the whole Seal command and output.

## Seal runs out of memory
Seal keeps every reference k-mer in memory. seqcounter sizes the heap from the length of the reference sequences
(about 40 MB per megabase, 1 GB at least), divided by the number of parallel jobs. For a very large reference set,
set it yourself: `--memory 16g -p 1`.

If seqcounter warns that the references need more memory than is available, the machine is too small for that
reference set and that `-p`: lower `-p`, or split the reference file and merge the resulting tables.

## "Several files give the same sample name"
Two input files reduce to the same sample name, so one would overwrite the other in the table. The message lists
them. Rename, move, or point `-i` at a narrower folder. Typical cases: the same sample in two subfolders, a fasta
and a fastq file of the same name, or an R1 file present twice.

## "reference sequence(s) ... are shorter than the k-mer size"
Seal cannot find a sequence shorter than `-k`. Lower `-k` to at most the length of the shortest reference, or take
that sequence out of the reference file. With a primer or probe set of 18 to 22 bp, `-k 18` is the usual answer.

## "Several reference sequences ... have the same name"
Seal reports one row per reference name, taken up to the first space, so two sequences whose first word is the same
would have their counts added together. Rename them in the reference file.

## The temporary folder contains a space, a comma, "xmx" or "xms"
`seal.sh` cannot be given such a path. Point `TMPDIR` at a folder whose path is made of plain characters:
`TMPDIR=/var/tmp seqcounter ...`.

## Nothing is printed, or the table is mixed with the messages
The counts table goes to standard output, the messages to standard error. Redirect them separately:

```bash
seqcounter -ref refs.fasta -i samples/ -o results/ > counts.tsv 2> seqcounter.log
```

## Still stuck
[Open an issue](https://github.com/duceppemo/seqcounter/issues/new/choose) with the version
(`seqcounter --version`), the command, and the output of a run with `--verbose`.
