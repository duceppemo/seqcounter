#!/usr/bin/env bash
# Check seqcounter's counts against counts that are known, on the Mycobacterium tuberculosis H37Rv genome:
#
#   1. an assembly, where every reference present must be counted once, whatever its copy number
#   2. simulated reads, where coverage is even and the expected count follows from the copy number
#   3. public read sets, where coverage is not even
#
# Everything is downloaded (about 150 MB) and nothing is kept outside validation/data and validation/results.
# Takes a few minutes. Results are on the wiki: https://github.com/duceppemo/seqcounter/wiki/Validation
#
#   bash validation/run_validation.sh
#
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"
data=data
results=results
mkdir -p "$data" "$results"

command -v seal.sh >/dev/null || { echo "seal.sh (BBTools) is not installed" >&2; exit 1; }
seqcounter="seqcounter"
command -v seqcounter >/dev/null || seqcounter="python3 $PWD/../seqcounter.py"

eutils=https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi
fetch() {  # fetch <url> <file>
  [ -s "$2" ] || { echo "downloading $(basename "$2")"; curl -sS --max-time 900 "$1" -o "$2"; }
}

echo "== reference genome and insertion sequence =="
fetch "$eutils?db=nuccore&id=NC_000962.3&rettype=fasta&retmode=text" "$data/H37Rv.fasta"
fetch "$eutils?db=nuccore&id=X17348.1&rettype=fasta&retmode=text" "$data/IS6110.fasta"

echo
echo "== reference set: IS6110, single-copy windows, and a sequence that is not in the genome =="
python3 validate.py references "$data/H37Rv.fasta" "$data/IS6110.fasta" "$data/refs.fasta"

echo
echo "== simulated reads, even coverage =="
[ -s "$data/sim.fastq.gz" ] || python3 validate.py simulate "$data/H37Rv.fasta" "$data/sim.fastq.gz"

echo
echo "== public read sets of M. tuberculosis (ENA) =="
mkdir -p "$data/reads"
for spec in ERR140/010/ERR14057010 ERR140/033/ERR14057033 ERR140/006/ERR14057206; do
  run=${spec##*/}
  for mate in 1 2; do
    fetch "https://ftp.sra.ebi.ac.uk/vol1/fastq/${spec}/${run}_${mate}.fastq.gz" "$data/reads/${run}_R${mate}.fastq.gz"
  done
done

echo
echo "== counting =="
$seqcounter -ref "$data/refs.fasta" -r1 "$data/H37Rv.fasta" -s assembly -o "$results/assembly" >/dev/null
for mismatch in 0 1; do
  $seqcounter -ref "$data/refs.fasta" -r1 "$data/sim.fastq.gz" -s "simulated_m$mismatch" \
              -o "$results/sim_m$mismatch" -m "$mismatch" >/dev/null
done
$seqcounter -ref "$data/refs.fasta" -i "$data/reads" -o "$results/reads" -p 3 >/dev/null

echo
echo "== the assembly: one contig, so every reference present is counted once =="
grep -v '^#' "$results/assembly/refseqcounts.tsv" | awk -F'\t' '
  NR == 1 { next } { count[$2]++ } END { for (c in count) printf "  %d reference(s) counted %s\n", count[c], c }'

echo
echo "== the same counts, asked twice: seqcounter -m 0 against an exact k-mer scan in Python =="
[ -s "$data/sim_subset.fastq.gz" ] || python3 validate.py subset "$data/sim.fastq.gz" "$data/sim_subset.fastq.gz"
$seqcounter -ref "$data/refs.fasta" -r1 "$data/sim_subset.fastq.gz" -s subset -o "$results/subset" -m 0 >/dev/null
python3 validate.py exact "$data/refs.fasta" "$data/sim_subset.fastq.gz" "$results/subset/refseqcounts.json"

echo
echo "== simulated reads, even coverage: the model should hold =="
python3 validate.py report "$data/refs.fasta" "$results"/sim_m*/refseqcounts.json

echo
echo "== what -m 1 adds over -m 0, on the same simulated reads =="
python3 validate.py mismatch "$results/sim_m0/refseqcounts.json" "$results/sim_m1/refseqcounts.json"

echo
echo "== public read sets: coverage is not even =="
python3 validate.py report "$data/refs.fasta" "$results/reads/refseqcounts.json"
