# FAQ

### Does a count tell me how many copies of the sequence are there?
No. A count is the number of **reads**, or contigs, containing the reference sequence. In an assembly, a sequence
present twice in one contig is counted once. In reads, the count grows with the number of copies, but it also grows
with the sequencing depth and with the length of the reference, so copy number has to be worked out from the depth,
not read off the table. See [How it works](How-it-works).

### Can I count nanopore or PacBio reads?
Yes, long reads are counted like any other fastq. Two things to keep in mind: a long read counts **once** however
many copies it carries, so counts are much lower than with short reads at the same depth; and on noisy reads a
25-mer often carries more than one error, which loses the match. `-k 19 -m 1` recovers most of it. For
high-accuracy data, or when you want to count per contig, counting the assembly is usually cleaner.

### Can I use it for presence and absence?
Yes, but the threshold is yours to choose. Read the count against `Reads` in `refseqcounts_samples.tsv`: on a deep
run a single matching read is as likely to be an error or a contaminant as a real hit, while on a shallow one it
may be all you get. seqcounter deliberately does not call presence for you.

### Why is `Matched` not the sum of a column?
`Matched` counts reads matching **at least one** reference; the table counts a read once for **each** reference it
matches. A read containing two references adds 1 to `Matched` and 1 to two rows.

### How many reference sequences can I give?
As many as fit in Seal's memory: a few thousand short sequences are nothing, a whole genome needs a bigger heap.
seqcounter sizes the heap from the total length of the references (about 40 MB per megabase), and `--memory` sets it
by hand. A short reference set is also easier to read: the table has one row per sequence.

### Can two references share sequence?
Yes. A read matching both is counted for both, so overlapping references give correlated counts. If you want each
read attributed to one reference only, seqcounter is not the right tool: that is a classification problem, not a
counting one.

### Does an interleaved fastq file work?
Yes. seqcounter tells Seal explicitly not to treat a single file as interleaved, so each read in it is counted
once, which is what you want for counting. Paired files are given as `-r1`/`-r2`.

### Does it need a reference genome, an index or a database?
No. The only input besides your samples is the fasta file of sequences to count, and nothing is cached between
runs.

### Are the results reproducible?
Yes. The same inputs and options give the same counts, whatever `-p` and `-t` are: rows are sorted by reference
name, columns follow the samples, and no step depends on which sample finished first. Only the run metadata
differs between runs: the first line of `refseqcounts.tsv` holds the command, and the JSON report the timestamps.

### Can I count protein sequences?
No, nucleotide sequences only. Both strands are searched, which only makes sense for DNA.

### How is this different from `bbduk.sh` or `kmercountexact.sh`?
Seal, which seqcounter drives, is the BBTools program for counting reads per reference sequence; `bbduk.sh` filters
or trims reads, and `kmercountexact.sh` counts k-mers of the input, not reads per reference. seqcounter adds sample
discovery and pairing, the checks on the reference file, the per-sample totals, the three reports and the parallel
runs, and it pins the Seal settings that make a count mean what this wiki says it means.

### Is it tied to spoligotyping?
No, although that is where it came from: counting the 43 spacers of the direct repeat locus is one reference file
among others. [spoligotyper](https://github.com/duceppemo/spoligotyper) does the whole spoligotyping job, including
reading the pattern, the species and the lineage.

### Can I call it from Python?
Yes, although the command line is the supported interface. The pieces are small:

```python
from seqcounter.counts import count_sample, write_reports
from seqcounter.references import load_references
from seqcounter.samples import single_sample

references = load_references('refs.fasta', 25)
result = count_sample(single_sample('sample_R1.fastq.gz', 'sample_R2.fastq.gz'), 'refs.fasta', threads=4)
print(result.counts, result.reads, result.matched)
write_reports([result], references, 'results/')
```
