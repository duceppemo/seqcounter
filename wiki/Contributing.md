# Contributing

Bug reports, questions, ideas and pull requests are all welcome. The full guide is in
[CONTRIBUTING.md](https://github.com/duceppemo/seqcounter/blob/master/CONTRIBUTING.md); this page is the short
version, plus what the code looks like from the inside.

## Reporting something
[Open an issue](https://github.com/duceppemo/seqcounter/issues/new/choose) with the version
(`seqcounter --version`), the command you ran and the output of a run with `--verbose`. For a count that looks
wrong, add the row of `refseqcounts.tsv`, the length of the reference sequence and the `-k`/`-m` you used.
[Troubleshooting](Troubleshooting) covers the usual causes.

## Development setup
```bash
git clone https://github.com/duceppemo/seqcounter
cd seqcounter
conda env create -f environment.yml
conda activate seqcounter
pip install --no-deps -e .
pytest --cov
ruff check .
```
The tests that drive Seal are skipped when `seal.sh` is not installed; the rest of the suite runs anywhere.

## The layout
```
seqcounter/
  __init__.py      the version
  __main__.py      python -m seqcounter
  errors.py        SeqCounterError: the errors reported without a traceback
  references.py    reading and checking the reference fasta file
  samples.py       finding samples in a folder, pairing R1/R2, naming them
  seal.py          everything about Seal: the command, safe paths, the stats file, the heap size
  counts.py        counting the samples, and the three reports
  cli.py           the arguments, and the run from end to end
seqcounter.py      a shim at the root, so python seqcounter.py still works
tests/             pytest, 100% coverage; conftest.py builds the synthetic data
wiki/              this wiki, published by .github/workflows/wiki.yml
recipe/            the bioconda recipe
examples/          run_example.sh: generates references and samples, and counts them
```

One quirk to know about: the root `seqcounter.py` shim has the same name as the package, so setuptools resolves
`seqcounter` to that file when it reads the version. That is why `pyproject.toml` points at
`seqcounter.__init__.__version__`, and why the shim imports nothing until it is run as a script. Renaming either
one breaks `python seqcounter.py` or the build.

Two rules worth knowing before changing anything:

* **Everything Seal needs is in `seal.py`**, including the settings that define what a count means
  (`clearzone`, `ambiguous`, `rcomp`, `maskmiddle`, `nzo`, `int`, `qin`) and the path handling that keeps
  `seal.sh`'s shell expansion harmless. Change a setting there, and say so in [How it works](How-it-works).
* **`refseqcounts.tsv` is a compatibility surface.** Its name, its `#command` first line and its shape go back to
  the first version and are what other people's scripts read. Add columns and files rather than changing it.

## Pull requests
* Branch from `master`, one change per pull request.
* Add or update tests; the suite is at 100% coverage and should stay there.
* Update `wiki/` if users are affected, and add a line to [Changelog](Changelog) under "Unreleased".
