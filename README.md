<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/duceppemo/seqcounter/main/assets/logo_dark.svg">
    <img src="https://raw.githubusercontent.com/duceppemo/seqcounter/main/assets/logo.svg" width="520" alt="seqcounter: count reference sequences in reads or assemblies">
  </picture>
</p>

<p align="center">
  <a href="https://github.com/duceppemo/seqcounter/actions/workflows/tests.yml"><img src="https://github.com/duceppemo/seqcounter/actions/workflows/tests.yml/badge.svg" alt="Tests"></a>
  <a href="https://codecov.io/gh/duceppemo/seqcounter"><img src="https://codecov.io/gh/duceppemo/seqcounter/graph/badge.svg" alt="Coverage"></a>
  <a href="https://github.com/duceppemo/seqcounter/releases/latest"><img src="https://img.shields.io/github/v/release/duceppemo/seqcounter?cacheSeconds=3600" alt="Release"></a>
  <a href="https://pypi.org/project/seqcounter/"><img src="https://img.shields.io/pypi/v/seqcounter?cacheSeconds=3600" alt="PyPI"></a>
  <img src="https://img.shields.io/badge/python-3.10%E2%80%933.14-blue" alt="Python 3.10–3.14">
  <a href="https://github.com/duceppemo/seqcounter/blob/main/LICENSE"><img src="https://img.shields.io/github/license/duceppemo/seqcounter" alt="License"></a>
  <a href="https://github.com/duceppemo/seqcounter/wiki"><img src="https://img.shields.io/badge/docs-wiki-informational" alt="Documentation"></a>
  <a href="https://doi.org/10.5281/zenodo.23826933"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.23826933.svg" alt="DOI"></a>
</p>

Count how many sequencing reads (fastq) or contigs (fasta) contain each sequence of a reference file, with
[Seal](https://sourceforge.net/projects/bbmap/) from BBTools. Point it at a reference fasta file and a sample, or a
folder of samples, and it writes one table of counts per reference sequence and per sample.

Any set of short sequences works: spoligotype spacers, marker or resistance genes, primer or probe sequences,
plasmid or phage fragments, barcodes, spike-ins. seqcounter does not interpret them; it counts them, reproducibly
and in the same way for reads and for assemblies.

```
$ seqcounter -ref refs.fasta -i samples/ -o results/
Reference   assembly  no_match  paired  single_end
absent      0         0         0       0
present     1         0         86      34
twice       1         0         154     65
```

## Features
* **Reads or assemblies**: single-end or paired-end fastq, or fasta, gzipped or not. Both mates of a pair are
  counted, once each.
* **Batch mode**: point it at a folder; fastq and fasta files are found recursively and R1/R2 files paired
  automatically. Several samples can be counted at the same time.
* **Both strands, and mismatches**: a reference is found on either strand, and with up to `-m` mismatches per k-mer
  (one by default).
* **Three reports**: the counts table (`refseqcounts.tsv`, the table the first version of seqcounter wrote, with the
  command that produced it on its first line), a per-sample summary of how much was read and how much matched, and a
  JSON report with the parameters, the software versions and the counts.
* **Nothing is hidden**: references Seal cannot find (shorter than the k-mer size), names that would collide, and
  samples that failed are reported, not silently dropped. A failed sample leaves the others readable.
* **Nothing to configure**: no database, no index, no reference genome. Only the Python standard library, plus
  BBTools.
* **Tested**: 105 tests, 100% coverage, end to end on synthetic references, genomes and reads with known counts.
* **Validated** on real data: counts identical to an independent k-mer scan, read for read; 16.1 copies recovered
  for a 16-copy element in simulated reads; a negative control at zero throughout
  ([Validation](https://github.com/duceppemo/seqcounter/wiki/Validation)).

## Installation
```
conda create -n seqcounter -c conda-forge -c bioconda bbmap python">=3.10"
conda activate seqcounter
pip install seqcounter
```
BBTools (which provides `seal.sh`) is not on PyPI, so it comes from bioconda; seqcounter itself needs nothing else.
See [Installation](https://github.com/duceppemo/seqcounter/wiki/Installation) for other options, including running
it from a clone with `python seqcounter.py`, as the first version did.

## Quick start
```
seqcounter -ref refs.fasta -r1 sample_R1.fastq.gz -r2 sample_R2.fastq.gz -o results/  # Paired-end reads
seqcounter -ref refs.fasta -r1 sample.fastq.gz -o results/                            # Single-end reads
seqcounter -ref refs.fasta -r1 assembly.fasta -o results/                             # An assembly
seqcounter -ref refs.fasta -i samples/ -o results/ -p 4                               # A folder, 4 at a time
```
`-ref` is a fasta file of the sequences to count: each one gives a row, named after the first word of its `>` line.
The counts are printed and saved in `results/`. To see it work on data with known counts, run
`bash examples/run_example.sh`, which generates its own references and samples.

**What a count is**: the number of reads, or contigs, that contain the reference sequence, not the number of copies
it has in them. A reference present twice in a single-contig assembly is counted once; in reads, about twice as many
reads match it. See [How it works](https://github.com/duceppemo/seqcounter/wiki/How-it-works).

## Documentation
The **[wiki](https://github.com/duceppemo/seqcounter/wiki)** covers
[installation](https://github.com/duceppemo/seqcounter/wiki/Installation), a
[tutorial](https://github.com/duceppemo/seqcounter/wiki/Tutorial) that counts IS6110 in three public samples in ten
minutes, [usage and options](https://github.com/duceppemo/seqcounter/wiki/Usage),
[output files](https://github.com/duceppemo/seqcounter/wiki/Output-files),
[how it works](https://github.com/duceppemo/seqcounter/wiki/How-it-works),
[validation](https://github.com/duceppemo/seqcounter/wiki/Validation),
[troubleshooting](https://github.com/duceppemo/seqcounter/wiki/Troubleshooting) and the
[FAQ](https://github.com/duceppemo/seqcounter/wiki/FAQ).

## Citing
If you use seqcounter, please cite it and BBTools, which does the counting:

> Duceppe M-O. seqcounter: count the occurrences of reference sequences in sequencing reads or genome assemblies.
> Zenodo. https://doi.org/10.5281/zenodo.23826933

This DOI always points to the latest version; each release also has its own DOI, listed on
[Zenodo](https://doi.org/10.5281/zenodo.23826933).

> Bushnell B. BBTools. https://sourceforge.net/projects/bbmap/

GitHub's **"Cite this repository"** button gives the citation in APA and BibTeX formats, from
[`CITATION.cff`](https://github.com/duceppemo/seqcounter/blob/main/CITATION.cff).

## Related
* [spoligotyper](https://github.com/duceppemo/spoligotyper) — *in silico* spoligotyping of the *M. tuberculosis*
  complex. It counts the 43 spacers the way seqcounter counts any reference, and then reads them as a spoligotype,
  with the species, the lineage and a PDF report.
* [brucetyper](https://github.com/duceppemo/brucetyper) — *in silico* Bruce-ladder, Suis-ladder, AMOS and RB51
  junction PCR for *Brucella*, from assemblies.

## Contributing
Bug reports, questions and pull requests are welcome: see [CONTRIBUTING.md](https://github.com/duceppemo/seqcounter/blob/main/CONTRIBUTING.md).

## Author
Marc-Olivier Duceppe, Canadian Food Inspection Agency (CFIA): marc-olivier.duceppe@inspection.gc.ca

## License
[MIT](https://github.com/duceppemo/seqcounter/blob/main/LICENSE)
