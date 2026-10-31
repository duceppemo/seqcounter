<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/duceppemo/seqcounter/main/assets/logo_dark.svg">
    <img src="https://raw.githubusercontent.com/duceppemo/seqcounter/main/assets/logo.svg" width="460" alt="seqcounter: count reference sequences in reads or assemblies">
  </picture>
</p>

Count how many sequencing reads (fastq) or contigs (fasta) contain each sequence of a reference fasta file, using
[Seal](https://sourceforge.net/projects/bbmap/) from BBTools.

```
seqcounter -ref refs.fasta -i samples/ -o results/
```

```
Reference   assembly  no_match  paired  single_end
absent      0         0         0       0
present     1         0         86      34
twice       1         0         154     65
```

Any set of short sequences works: spoligotype spacers, marker or resistance genes, primer or probe sequences,
plasmid or phage fragments, barcodes, spike-ins. seqcounter does not interpret them; it counts them, the same way
for reads and for assemblies, and writes the counts as a table, a per-sample summary and a JSON report.

## Where to start
* **[Installation](Installation)** — conda, pip, or from a clone
* **[Tutorial](Tutorial)** — count IS6110 in three public samples, in ten minutes
* **[Usage](Usage)** — the input, every option, and what to expect
* **[Output files](Output-files)** — the three reports, column by column
* **[How it works](How-it-works)** — what a count is, exactly, and how Seal is called
* **[Validation](Validation)** — the counts checked against known answers, in an assembly, in simulated reads and
  in public read sets
* **[Troubleshooting](Troubleshooting)** — when a count is 0, a sample fails, or Seal runs out of memory
* **[FAQ](FAQ)**
* **[Changelog](Changelog)**

## What it is not
seqcounter reports counts; it does not call a genotype, a species or a presence/absence threshold. A count has to be
read in the context of how much was sequenced, which the per-sample summary gives. Tools built on top of these
counts decide what they mean: [spoligotyper](https://github.com/duceppemo/spoligotyper), for example, counts the 43
spacers of the direct repeat locus and reads them as a spoligotype.

## Requirements
Python 3.10 or newer, and BBTools (`seal.sh`, which needs Java). seqcounter itself uses nothing but the Python
standard library.

## Citing
Please cite seqcounter and BBTools:

> Duceppe M-O. seqcounter: count the occurrences of reference sequences in sequencing reads or genome assemblies.
> Zenodo. https://doi.org/10.5281/zenodo.23826933

That DOI resolves to the latest version. The "Cite this repository" button on the repository page gives the same
citation in APA and BibTeX, from [CITATION.cff](https://github.com/duceppemo/seqcounter/blob/main/CITATION.cff).
