# Output files

Every run writes three files to the folder given by `-o`, and prints the counts table to standard output.

| File | What it holds |
| --- | --- |
| `refseqcounts.tsv` | The counts: one row per reference sequence, one column per sample |
| `refseqcounts_samples.tsv` | One row per sample: how much was read, how much matched, and any error |
| `refseqcounts.json` | The whole run: parameters, software versions, references and counts |

The files are overwritten by the next run into the same folder.

## refseqcounts.tsv
```
#seqcounter -ref refs.fasta -i samples/ -o results/ -p 2
Reference	assembly	no_match	paired	single_end
absent	0	0	0	0
present	1	0	86	34
twice	1	0	154	65
```

* The **first line** is the command that produced the file, prefixed with `#`, for the record. Skip it when reading
  the file (`pandas.read_csv(path, sep='\t', comment='#', index_col=0)`, `grep -v '^#'`, ...).
* One row per reference sequence, named after the first word of its `>` line, **sorted by name**. A reference found
  in no sample has a row of zeros: nothing is dropped.
* One column per sample, in the order the samples were counted (alphabetical for `-i`).
* A cell is the number of reads, or contigs, containing that reference: see [How it works](How-it-works).
* A sample that **failed** has empty cells, so the samples that worked stay readable as numbers.

The same table is printed to standard output, without the `#` line, so it can be piped or redirected:
`seqcounter ... -o results/ > counts.tsv`. Messages go to standard error.

## refseqcounts_samples.tsv
One row per sample, which is what a count has to be read against: 34 matching reads mean something different in a
thousand reads and in a million.

| Column | Meaning |
| --- | --- |
| `Sample` | Sample name |
| `Type` | `fastq` (reads) or `fasta` (assembly) |
| `Reads` | Reads in the input, both mates together; contigs for a fasta file |
| `Bases` | Bases in the input |
| `Matched` | Reads matching at least one reference sequence (a read matching two is counted once here) |
| `MatchedPct` | `Matched` as a percentage of `Reads` |
| `Files` | The input file(s), separated by `;` |
| `Error` | Empty, or why the sample failed |

`Matched` is not the sum of the column in `refseqcounts.tsv`: a read containing two different references is counted
once in `Matched` and once for each reference in the counts table.

## refseqcounts.json
For pipelines, and for the record: everything the run knew.

```json
{
  "tool": "seqcounter",
  "version": "1.0.0",
  "references": [
    {"name": "present", "description": "example reference", "length": 60}
  ],
  "samples": [
    {"sample": "assembly", "type": "fasta", "files": ["/data/samples/assembly.fasta"],
     "reads": 1, "bases": 4180, "matched": 1, "matched_percent": 100.0,
     "counts": {"absent": 0, "present": 1, "twice": 1}, "error": null}
  ],
  "command": "seqcounter -ref refs.fasta -i samples/ -o results/ -p 2",
  "started": "2026-10-30T12:00:00-04:00",
  "finished": "2026-10-30T12:00:03-04:00",
  "software": {"seqcounter": "1.0.0", "Python": "3.12.3", "BBTools": "39.13", "Java": "openjdk version \"17.0.9\""},
  "parameters": {"Reference": "/data/refs.fasta", "Input": "/data/samples", "Output folder": "/data/results",
                 "K-mer size": 25, "Mismatches": 1, "Threads": 16, "Parallel samples": 2, "Seal memory": "1024m"}
}
```

`description` is the part of the `>` line after the first word. `error` is `null` for a sample that was counted, and
the message otherwise.

## Reading the counts
A count is a number of reads, which depends on how much was sequenced. To compare samples, normalize by `Reads`
from the per-sample report (counts per million reads, say), and remember that a longer reference is covered by more
reads than a short one at the same depth. For presence/absence, pick a minimum count that matches the depth of your
data: on a 30× bacterial run, a single matching read is as likely to be an error or a contaminant as a real hit.
