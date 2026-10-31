#!/usr/bin/env bash
# Count three reference sequences in four samples whose counts are known, and show the reports.
# Everything is generated here: no download, and nothing but seqcounter, Seal and Python is needed.
#
#   bash examples/run_example.sh
#
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
data="$here/data"
results="$here/results"
rm -rf "$data" "$results"
mkdir -p "$data/samples"

command -v seal.sh >/dev/null || { echo "seal.sh (BBTools) is not installed: conda install -c bioconda bbmap" >&2;
                                   exit 1; }
seqcounter="seqcounter"
command -v seqcounter >/dev/null || seqcounter="python3 $here/../seqcounter.py"

# Three 60 bp reference sequences, and a 4 kb genome holding "present" once and "twice" twice
python3 - "$data" <<'PYEOF'
import gzip
import random
import sys

folder = sys.argv[1]
rng = random.Random(42)
random_seq = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
refs = {name: random_seq(60) for name in ('absent', 'present', 'twice')}
with open(folder + '/refs.fasta', 'w') as f:
    for name, seq in refs.items():
        f.write('>{} example reference\n{}\n'.format(name, seq))

genome = (random_seq(1000) + refs['present'] + random_seq(1000) + refs['twice'] + random_seq(1000)
          + refs['twice'] + random_seq(1000))
with open(folder + '/samples/assembly.fasta', 'w') as f:
    f.write('>contig1\n{}\n'.format(genome))
with open(folder + '/samples/no_match.fasta', 'w') as f:
    f.write('>contig1\n{}\n'.format(random_seq(4000)))

revcomp = lambda s: s[::-1].translate(str.maketrans('ACGT', 'TGCA'))
def write_reads(path, depth, mate=None):
    with gzip.open(path, 'wt') as f:
        for i in range(int(depth * len(genome) / 100)):
            start = rng.randint(0, len(genome) - 100)
            read = genome[start:start + 100]
            if (rng.random() < 0.5) != (mate == 2):
                read = revcomp(read)
            f.write('@read{}/{}\n{}\n+\n{}\n'.format(i, mate or 1, read, 'I' * 100))

write_reads(folder + '/samples/single_end.fastq.gz', 30)
write_reads(folder + '/samples/paired_R1.fastq.gz', 15, mate=1)
write_reads(folder + '/samples/paired_R2.fastq.gz', 15, mate=2)
print('Reference sequences and samples written to ' + folder)
PYEOF

echo
echo "### $seqcounter -ref $data/refs.fasta -i $data/samples -o $results -p 2"
$seqcounter -ref "$data/refs.fasta" -i "$data/samples" -o "$results" -p 2

echo
echo "### $results/refseqcounts.tsv (without its first line, which holds the command that was run)"
grep -v '^#' "$results/refseqcounts.tsv" | column -t -s $'\t'
echo
echo "### $results/refseqcounts_samples.tsv"
cut -f 1-6 "$results/refseqcounts_samples.tsv" | column -t -s $'\t'
echo
echo "What to expect: 'absent' is in no sample, 'present' and 'twice' are in every sample but no_match, and about"
echo "twice as many reads match 'twice' as 'present' because the genome holds it twice. In the assembly, which is"
echo "one contig, both are counted once: the count is the number of sequences containing the reference, not the"
echo "number of copies. The full run, with the software versions and parameters, is in refseqcounts.json."
