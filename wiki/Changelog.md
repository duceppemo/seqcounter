# Changelog

## Unreleased
* A [Tutorial](Tutorial): counting IS6110 in three public *M. tuberculosis* read sets, from installing to reading
  the reports and normalizing the counts, with the real output of every command.
* A [Validation](Validation) page, and `validation/run_validation.sh`, which reproduces it from public data in a
  few minutes. What it establishes: the counts are identical, read for read, to an independent exact k-mer scan
  over 200,000 reads and 77 references; a reference is counted once per sequence of a fasta file, even at 16
  copies; with even coverage the counts follow the model to a median of 1.00 and recover a 16-copy element as 16.1
  copies; a sequence absent from the genome stays at zero everywhere; and `-m 1` adds reads from near-matches
  elsewhere in the genome for 22 of 77 references, by nothing at the median window and up to +43%. Two figures,
  drawn from the results by `validation/make_figures.py`, show the counts against the model and the spread of
  real coverage.
* The logo is drawn as outlines rather than text, so it renders without the font it was designed in, and the wiki
  front page carries it.

## 1.0.1 (2026-10-30)
* The default branch is now `main`, and the links that named the branch follow it: the logo in the README, which
  `raw.githubusercontent.com` serves per branch and does not redirect after a rename, and the `blob/` links from
  this wiki.
* The Zenodo DOI of the archived releases is recorded in `CITATION.cff` and in the README badge:
  [10.5281/zenodo.23826933](https://doi.org/10.5281/zenodo.23826933) resolves to the latest version.
* No change to the code: 1.0.1 and 1.0.0 count identically.

## 1.0.0 (2026-10-30)
seqcounter becomes an installable, tested package. Everything the first version did, it still does, with the same
flags and the same `refseqcounts.tsv`; what changed is around that.

**Packaging and interface**
* Installable with `pip install seqcounter` (or from bioconda): a `seqcounter` command, and `python -m seqcounter`.
  `python seqcounter.py` still works from a clone, with the same arguments, so existing scripts keep working.
* **No Python dependency at all**: pandas and psutil are gone (the table is written directly, and the memory check
  no longer needs psutil). Only BBTools and Java are needed, as before.
* The counts table is also **printed to standard output**, with the messages on standard error, so
  `seqcounter ... > counts.tsv` works.
* New options: `-s`/`--sample` (the sample name used in the reports), `--memory` (Seal's Java heap),
  `--verbose` (debug messages, including the Seal command), `-j`/`--jobs` as another name for `-p`/`--parallel`.
  `-v` still prints the version.
* Two new reports beside `refseqcounts.tsv`: `refseqcounts_samples.tsv` (reads, bases, matched reads and the error
  per sample) and `refseqcounts.json` (parameters, software versions, references and counts).
  See [Output files](Output-files).
* Messages are logged with a timestamp and a level, and an expected error is reported without a traceback.
  Exit codes: 0, 1 (failed), 2 (bad command line), 130 (interrupted).

**Counts that can be trusted**
* **A single fastq file is never read as interleaved** (`int=f`). Seal guesses interleaving from the read names,
  and a file it guessed wrong on had its pairs counted as one read each: the counts were halved.
* **The quality encoding is forced** (`qin=33`), instead of being autodetected, which fails on some low quality
  nanopore reads.
* **The reference file is checked** before anything runs: a sequence shorter than `-k` (which Seal cannot find at
  all), an empty sequence, two sequences whose first word is the same (Seal would add their counts together in one
  row), a file that is not fasta. Each is reported with what to do about it.
* **Sample name collisions are reported** instead of being resolved silently. The first version put the files of a
  folder in a dictionary keyed by sample name, so two files reducing to the same name (the same sample in two
  subfolders, a fasta and a fastq of the same name) were mixed up or overwritten. R1/R2 pairing is now explicit,
  and an R2 without its R1 is counted as single-end with a warning.
* **File names Seal cannot be given are handled.** `seal.sh` runs its command line through the shell and splits
  `ref=` lists on commas, so a path with a space, a quote, `$`, `;`, `&`, parentheses or a comma was
  re-interpreted — in the best case the run failed, in the worst a file name could run a command. Such paths, and
  paths containing `xmx` or `xms` (which BBTools 40 takes for a memory setting), are now linked under a neutral
  name.
* **Formats are recognized by content**, not only by extension: a fastq file named `.fasta`, or a file with no
  known extension such as a Galaxy `dataset_1.dat`, is presented to Seal under a matching name instead of being
  misread.
* **The column order is deterministic** (the samples in order, alphabetical for `-i`), where it used to depend on
  which sample finished first.
* **The stats file** Seal writes goes to a private temporary folder, which is always removed, instead of being
  written into the output folder under the sample name and deleted afterwards.
* `maskmiddle=f` is kept from the first version: without it Seal ignores the middle base of each k-mer, which
  allows one more mismatch than `-m` asks for.

**Running it**
* **Seal's memory is bounded**: the Java heap is set from the total length of the reference sequences (1 GB at
  least), capped at 85% of the machine's memory divided by the number of parallel jobs. The first version gave
  `seal.sh` no `-Xmx`, so each job reserved most of the free memory and `-p 2` could put the machine into swap.
  `--memory` overrides it.
* **One sample failing no longer stops the run**: its column is left empty, its error is reported in the per-sample
  and JSON reports, and the exit code is 1.
* `-t` respects the CPUs the process is actually allowed to use (`sched_getaffinity`), which matters inside a
  SLURM job or a container.
* The output folder is skipped when it is inside the input folder, so running twice into `samples/results/` does
  not count the reports as samples.

**Development**
* 105 tests, **100% coverage** of the package, end to end on synthetic references, genomes and reads with known
  counts; the tests that need Seal are skipped when BBTools is absent.
* GitHub Actions: tests on Python 3.10, 3.12 and 3.14 with ruff and Codecov, PyPI publishing on release, and this
  wiki published from `wiki/`.
* A [wiki](Home), a [bioconda recipe](https://github.com/duceppemo/seqcounter/tree/main/recipe), a
  [`CITATION.cff`](https://github.com/duceppemo/seqcounter/blob/main/CITATION.cff), and
  `examples/run_example.sh`, which generates its own data and counts it.

## 0.1 (2023-05-04, tagged 2026-10-30)
The first version: `seqcounter.py` and `seqcounter_methods.py`, run from a clone, counting reference sequences in
fastq or fasta with Seal and writing `refseqcounts.tsv` with pandas. Tagged as
[v0.1](https://github.com/duceppemo/seqcounter/releases/tag/v0.1) so that work based on the clone-and-run layout has
a fixed point to refer to.
