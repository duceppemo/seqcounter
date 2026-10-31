# Contributing to seqcounter

Thanks for your interest! Bug reports, questions, ideas and pull requests are all welcome.
Please follow the [code of conduct](CODE_OF_CONDUCT.md).

## Reporting a bug or asking a question
[Open an issue](https://github.com/duceppemo/seqcounter/issues/new/choose) and pick the matching template.
For bugs, please include the version (`seqcounter --version`), the command you ran and the messages it printed
(ideally from a run with `--verbose`). For a count that looks wrong, include the row of `refseqcounts.tsv`, the
length of the reference sequence and the k-mer size you used. Check the
[troubleshooting page](https://github.com/duceppemo/seqcounter/wiki/Troubleshooting) and the
[FAQ](https://github.com/duceppemo/seqcounter/wiki/FAQ) first.

## Development setup
```
git clone https://github.com/duceppemo/seqcounter
cd seqcounter
conda env create -f environment.yml
conda activate seqcounter
pip install --no-deps -e .
pytest --cov
```
seqcounter itself needs nothing but the Python standard library; the environment is for BBTools and the test tools.
The end-to-end tests in `tests/test_integration.py` run Seal on small synthetic references, genomes and reads, and
are skipped if `seal.sh` is not installed. The suite covers 100% of the package, so please keep it there.

Lint with [ruff](https://docs.astral.sh/ruff/): `ruff check .`, or install the git hook with `pre-commit install`.
`ruff format` is **not** used: the code is wrapped by hand to read as prose, and reformatting it would rewrite
most of the files for nothing.
The tests run on GitHub Actions for Python 3.10, 3.12 and 3.14.

To see it work on data you can look at, run `bash examples/run_example.sh`: it builds a small reference file and a
few samples with known counts, runs seqcounter on them and prints the reports.

## Pull requests
* Branch from `main` and keep each pull request focused on one change.
* Add or update tests for any change in behaviour.
* Match the style of the surrounding code.
* Update the documentation in `wiki/` if the change affects users (see below), and add a line to
  `wiki/Changelog.md` under "Unreleased" (the section is created by the first change after a release).

## Documentation
The [wiki](https://github.com/duceppemo/seqcounter/wiki) is maintained in the [`wiki/`](wiki/) folder of this
repository and published automatically when changes reach `main`. **Do not edit the wiki on GitHub directly**:
those changes would be overwritten. Edit the files in `wiki/` in your pull request instead.

## Releases (maintainers)
1. Update the version in `seqcounter/__init__.py` and `CITATION.cff` (`version` and `date-released`).
2. Rename "Unreleased" in `wiki/Changelog.md` to the version and date.
3. Commit, tag (`git tag -a vX.Y.Z`), push the commit and the tag, then create the GitHub release once the tests pass.
4. Publishing the GitHub release automatically uploads the package to
   [PyPI](https://pypi.org/project/seqcounter/) (`.github/workflows/publish.yml`, trusted publishing: no token
   needed; the release tag must match the package version), and Zenodo archives it and mints a DOI. Add that DOI to
   `CITATION.cff` and to the badge in `README.md`.
5. Bioconda: update `recipe/meta.yaml` (version and sha256 of the new GitHub tarball). Bioconda's bot usually opens
   the pull request to `bioconda-recipes` on its own for new releases; otherwise copy the recipe and open it by hand.
