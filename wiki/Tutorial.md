# Tutorial

Count one reference sequence in three public samples, and read the result properly. About ten minutes, most of it
downloading, and nothing is needed but seqcounter, BBTools and `curl`.

The sequence is **IS6110**, the insertion element of the *Mycobacterium tuberculosis* complex. It is a good first
reference because every isolate has several copies, and not the same number: counting it says something real about
the sample.

## 1. Install

```bash
conda create -n seqcounter -c conda-forge -c bioconda bbmap python">=3.10"
conda activate seqcounter
pip install seqcounter
seqcounter --version
```

See [Installation](Installation) for other routes.

## 2. Make the reference file

Any fasta file of sequences works. Take IS6110 from GenBank:

```bash
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=X17348.1&rettype=fasta&retmode=text" \
  -o refs.fasta
head -1 refs.fasta
```
```
>X17348.1 Mycobacterium tuberculosis IS6110 IS-like element
```

The row in the counts table is named after the **first word** of the `>` line, so this would give a row called
`X17348.1`. Rename it to something you want to read in a table:

```bash
sed -i '1s/^>X17348.1/>IS6110 X17348.1/' refs.fasta
```

## 3. Get some samples

Three *M. tuberculosis* runs from the ENA, about 150 MB in all:

```bash
mkdir -p reads
for spec in ERR140/010/ERR14057010 ERR140/033/ERR14057033 ERR140/006/ERR14057206; do
  run=${spec##*/}
  for mate in 1 2; do
    curl -s "https://ftp.sra.ebi.ac.uk/vol1/fastq/${spec}/${run}_${mate}.fastq.gz" -o "reads/${run}_R${mate}.fastq.gz"
  done
done
ls reads/
```
```
ERR14057010_R1.fastq.gz  ERR14057033_R1.fastq.gz  ERR14057206_R1.fastq.gz
ERR14057010_R2.fastq.gz  ERR14057033_R2.fastq.gz  ERR14057206_R2.fastq.gz
```

The `_R1`/`_R2` names are all seqcounter needs to pair the files; three samples will be found, not six.

## 4. Count

```bash
seqcounter -ref refs.fasta -i reads -o results -p 3
```
```
3 sample(s) found in reads
3 sample(s), 3 at a time, 21 thread(s) each
ERR14057033: 3,550 of 586,752 read(s) matched
ERR14057010: 5,442 of 599,360 read(s) matched
ERR14057206: 10,994 of 1,055,736 read(s) matched
Reference	ERR14057010	ERR14057033	ERR14057206
IS6110	5442	3550	10994
Report saved in results/refseqcounts.tsv
Done in 3s
```

Three seconds for 2.2 million reads. The table went to standard output as well as to `results/`, so
`seqcounter ... > counts.tsv` works, and the messages went to standard error.

## 5. Read the reports

`results/refseqcounts.tsv` is the table above, with the command that produced it on its first line.
`results/refseqcounts_samples.tsv` is what the counts have to be read against:

```
Sample       Type   Reads    Bases      Matched  MatchedPct
ERR14057010  fastq  599360   88796451   5442     0.9080
ERR14057033  fastq  586752   87002595   3550     0.6050
ERR14057206  fastq  1055736  156986138  10994    1.0414
```

ERR14057206 has the largest count, 10,994 — and also nearly twice as many reads as the others. A raw count is not
comparable between samples of different size. Divide by the reads:

| sample | IS6110 | total reads | per million reads |
| --- | --- | --- | --- |
| ERR14057010 | 5,442 | 599,360 | 9,080 |
| ERR14057033 | 3,550 | 586,752 | 6,050 |
| ERR14057206 | 10,994 | 1,055,736 | 10,414 |

Now the order is different from the raw counts, and the real result appears: ERR14057033 carries clearly less
IS6110 than the other two, per read sequenced. That is a property of the isolates — IS6110 copy number varies
between strains, which is why it was the standard fingerprint for tuberculosis for decades. The
[Validation](Validation) page turns these same numbers into copy numbers (14.7, 22.8 and 24.9) by adding
single-copy references to the file and calibrating against them.

`results/refseqcounts.json` holds all of it again, with the parameters, the software versions and the input files,
for a pipeline or for the record.

## 6. The same reference in an assembly

```bash
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000962.3&rettype=fasta&retmode=text" \
  -o H37Rv.fasta
seqcounter -ref refs.fasta -r1 H37Rv.fasta -o assembly_results
```
```
H37Rv: 1 of 1 read(s) matched
Reference	H37Rv
IS6110	1
```

**1**, not 16 — and H37Rv does hold 16 copies of IS6110. A fasta file is counted the same way as reads: a count is
the number of sequences containing the reference, and this genome is a single contig, so the answer is one.
Assemblies tell you presence or absence; the copy number is in the reads. This is the one thing about seqcounter
worth remembering, and [How it works](How-it-works) explains it in full.

## Where to go next

* [Usage](Usage) — every option, and what `-k` and `-m` actually change
* [Output files](Output-files) — the three reports, column by column
* [Validation](Validation) — what the counts were checked against, with the numbers
* [Troubleshooting](Troubleshooting) — when everything comes back 0, or a count looks too high
