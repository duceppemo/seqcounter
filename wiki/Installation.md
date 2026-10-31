# Installation

seqcounter needs **Python 3.10 or newer** and **BBTools**, which provides `seal.sh` and needs Java. seqcounter
itself has no Python dependency: it uses only the standard library.

## conda (recommended)
BBTools is not on PyPI, so the simplest route is a conda environment holding both:

```bash
conda create -n seqcounter -c conda-forge -c bioconda bbmap python">=3.10"
conda activate seqcounter
pip install seqcounter
seqcounter --version
```

`mamba` or `micromamba` work the same way and are faster.

## pip, with BBTools already installed
```bash
pip install seqcounter
```
or, without touching your environment:
```bash
pipx install seqcounter
```
`seal.sh` must be in your `PATH`. seqcounter also finds a `seal.sh` installed next to the Python that runs it, so
`/path/to/env/bin/seqcounter` works without activating that environment.

To install BBTools by hand instead, download it from [SourceForge](https://sourceforge.net/projects/bbmap/), unpack
it and add its folder to your `PATH`:
```bash
tar -xzf BBMap_39.13.tar.gz
export PATH="$PWD/bbmap:$PATH"
seal.sh --version
```

## From a clone
Nothing has to be installed to run seqcounter from a clone, as the first version did:

```bash
conda create -n seqcounter -y -c conda-forge -c bioconda bbmap python">=3.10"
conda activate seqcounter
git clone https://github.com/duceppemo/seqcounter
cd seqcounter
python seqcounter.py -h
```

`python seqcounter.py` takes exactly the same arguments as the installed `seqcounter` command, so scripts written
for the first version keep working. `pip install -e .` from the clone adds the `seqcounter` command as well.

## Checking the installation
```bash
bash examples/run_example.sh
```
from a clone generates a reference file and four samples whose counts are known, counts them and prints the
reports. Without a clone:
```bash
printf '>ref1\nTTCGAACGCTGACGTCAGCATCGATTGCAGCTAGCTTAGCGATCGGATCGA\n' > refs.fasta
seqcounter -ref refs.fasta -r1 refs.fasta -o test_out   # ref1 is found once, in its own file
```

## Which version is installed
```bash
seqcounter --version
seal.sh --version          # BBTools and Java
```
The JSON report of every run also records the versions of seqcounter, Python, BBTools and Java that produced it.
