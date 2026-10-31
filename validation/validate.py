#!/usr/bin/env python3
"""
Build the reference set for the validation, and report what seqcounter counted against what is expected.

    python3 validate.py references H37Rv.fasta IS6110.fasta refs.fasta
    python3 validate.py simulate H37Rv.fasta sim.fastq.gz --reads 1000000 --length 150 --error 0.002
    python3 validate.py subset sim.fastq.gz sim_subset.fastq.gz --reads 200000
    python3 validate.py exact refs.fasta sim_subset.fastq.gz out_subset/refseqcounts.json
    python3 validate.py report refs.fasta out_assembly/refseqcounts.json out_*/refseqcounts.json

The expected counts come from a uniform-coverage model: a read of length R matches a reference of length L when it
overlaps one of its copies by at least k bases, which happens for L + R - 2k + 1 of the G start positions in the
genome, so a run of N reads is expected to give N * copies * (L + R - 2k + 1) / G. See the wiki's Validation page.
"""
import argparse
import gzip
import json
import random
import statistics
import sys

K = 25  # The k-mer size used throughout, matching seqcounter's default


def read_fasta(path):
    records, name, parts = {}, None, []
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                if name:
                    records[name] = ''.join(parts)
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line.strip().upper())
    records[name] = ''.join(parts)
    return records


def reverse_complement(seq):
    return seq[::-1].translate(str.maketrans('ACGT', 'TGCA'))


def counts_of(report_path):
    """The counts of the first (or only) sample of a seqcounter JSON report."""
    with open(report_path) as f:
        return json.load(f)['samples'][0]['counts']


def samples_of(report_path):
    with open(report_path) as f:
        return json.load(f)['samples']


def occurrences(genome, sequences):
    """{name: [occurrences of each k-mer in the genome, both strands]}, in one pass over the genome."""
    wanted = {}
    for seq in sequences.values():
        for i in range(len(seq) - K + 1):
            kmer = seq[i:i + K]
            wanted[kmer] = 0
            wanted[reverse_complement(kmer)] = 0
    for i in range(len(genome) - K + 1):
        kmer = genome[i:i + K]
        if kmer in wanted:
            wanted[kmer] += 1
    return {name: [wanted[seq[i:i + K]] + wanted[reverse_complement(seq[i:i + K])]
                   for i in range(len(seq) - K + 1)] for name, seq in sequences.items()}


def build_references(genome_path, element_path, out_path, windows=60, lengths=(600, 300, 150, 60)):
    """
    The reference set: the insertion sequence, single-copy windows of the genome, and a sequence that is not in it.

    A window is kept only when every one of its k-mers occurs exactly once in the genome, on either strand, so its
    copy number is one by construction rather than by assumption.
    """
    genome = next(iter(read_fasta(genome_path).values()))
    element = next(iter(read_fasta(element_path).values()))
    size = len(genome)

    candidates = {'w1000_{}'.format(i): (40_000 + i * (size - 80_000) // windows, 1000) for i in range(windows)}
    for i, length in enumerate(lengths):
        for j in range(6):
            candidates['w{}_{}'.format(length, j)] = (60_000 + (i * 6 + j) * (size - 120_000) // (6 * len(lengths)),
                                                      length)
    sequences = {name: genome[start:start + length] for name, (start, length) in candidates.items()}
    sequences['IS6110'] = element
    counts = occurrences(genome, sequences)

    single = sorted(name for name, c in counts.items() if name != 'IS6110' and min(c) == 1 and max(c) == 1)
    element_copies = int(statistics.median(counts['IS6110']))
    rng = random.Random(1)
    with open(out_path, 'w') as f:
        f.write('>IS6110 GenBank X17348.1, {} copies in this genome\n{}\n'.format(element_copies, element))
        for name in single:
            start, length = candidates[name]
            f.write('>{} {}-{} of the genome, single copy\n{}\n'.format(
                name, start + 1, start + length, sequences[name]))
        f.write('>absent_300bp random sequence, not in the genome\n{}\n'.format(
            ''.join(rng.choice('ACGT') for _ in range(300))))

    print('genome {:,} bp'.format(size))
    print('IS6110 {} bp: k-mer occurrences min {} median {} max {} -> {} copies'.format(
        len(element), min(counts['IS6110']), element_copies, max(counts['IS6110']), element_copies))
    print('{} of {} windows are strictly single copy; {} references written to {}'.format(
        len(single), len(candidates), len(single) + 2, out_path))
    return 0


def simulate(genome_path, out_path, reads, length, error, seed=42):
    """Reads drawn uniformly from the genome, on both strands: coverage is even by construction."""
    genome = next(iter(read_fasta(genome_path).values()))
    rng = random.Random(seed)
    with gzip.open(out_path, 'wt', compresslevel=1) as f:
        for i in range(reads):
            start = rng.randrange(len(genome) - length)
            read = genome[start:start + length]
            if error:
                read = ''.join(rng.choice('ACGT') if rng.random() < error else base for base in read)
            if rng.random() < 0.5:
                read = reverse_complement(read)
            f.write('@r{}\n{}\n+\n{}\n'.format(i, read, 'I' * length))
    print('{:,} reads of {} bp with {:.1%} errors -> {} ({:.1f}x)'.format(
        reads, length, error, out_path, reads * length / len(genome)))
    return 0


def subset(fastq_path, out_path, reads):
    """The first reads of a fastq file, so that the Python scan and seqcounter can be given the same input."""
    written = 0
    with gzip.open(fastq_path, 'rt') as f, gzip.open(out_path, 'wt', compresslevel=1) as out:
        for i, line in enumerate(f):
            out.write(line)
            if i % 4 == 3:
                written += 1
                if written == reads:
                    break
    print('{:,} reads -> {}'.format(written, out_path))
    return 0


def exact(refs_path, fastq_path, report_path):
    """
    Count, in Python, the reads sharing an exact k-mer with each reference, and compare with what seqcounter
    counted with -m 0. The two must agree read for read: it is the same question asked twice.
    """
    kmers = {}
    for name, seq in read_fasta(refs_path).items():
        own = {seq[i:i + K] for i in range(len(seq) - K + 1)}
        kmers[name] = own | {reverse_complement(k) for k in own}
    counted = dict.fromkeys(kmers, 0)
    reads = 0
    with gzip.open(fastq_path, 'rt') as f:
        for i, line in enumerate(f):
            if i % 4 != 1:
                continue
            read = line.strip()
            reads += 1
            read_kmers = {read[j:j + K] for j in range(len(read) - K + 1)}
            for name, reference_kmers in kmers.items():
                if read_kmers & reference_kmers:
                    counted[name] += 1

    seqcounter_counts = counts_of(report_path)
    differences = {name: (counted[name], seqcounter_counts[name])
                   for name in kmers if counted[name] != seqcounter_counts[name]}
    print('{:,} reads scanned in Python against {} references'.format(reads, len(kmers)))
    if differences:
        print('DIFFERENT for {} reference(s):'.format(len(differences)))
        for name, (mine, theirs) in sorted(differences.items())[:10]:
            print('  {:<16} python {:>7}  seqcounter {:>7}'.format(name, mine, theirs))
        return 1
    print('identical for every reference, including IS6110 ({:,} reads) and the control ({})'.format(
        counted['IS6110'], counted['absent_300bp']))
    return 0


def mismatch(exact_report, mismatch_report):
    """What allowing one mismatch per k-mer adds: the same sample counted with -m 0 and with -m 1."""
    a = counts_of(exact_report)
    b = counts_of(mismatch_report)
    gained = {name: (a[name], b[name]) for name in a if b[name] != a[name]}
    windows = [n for n in a if n.startswith('w')]
    extra = [(b[n] - a[n]) / a[n] for n in windows if a[n]]
    print('{} of {} references counted more reads with -m 1 than with -m 0'.format(len(gained), len(a)))
    print('single-copy windows: median {:+.1%}, worst {:+.1%}'.format(
        statistics.median(extra), max(extra)))
    worst = sorted(gained.items(), key=lambda kv: (kv[1][1] - kv[1][0]) / max(kv[1][0], 1), reverse=True)[:5]
    for name, (exact_count, mismatch_count) in worst:
        print('  {:<16} {:>7} -> {:<7} {:+.0%}'.format(
            name, exact_count, mismatch_count, (mismatch_count - exact_count) / max(exact_count, 1)))
    print('  the control stays at {}'.format(b['absent_300bp']))
    return 0


def reference_lengths(refs_path):
    return {name: len(seq) for name, seq in read_fasta(refs_path).items()}


def report(refs_path, genome_size, reports):
    """Compare the counts of each run with the uniform-coverage model."""
    lengths = reference_lengths(refs_path)
    windows = [name for name in lengths if name.startswith('w')]
    print('{} single-copy windows, {} bp to {} bp\n'.format(
        len(windows), min(lengths[w] for w in windows), max(lengths[w] for w in windows)))
    header = '{:<16} {:>9} {:>7} {:>8} {:>13} {:>11} {:>8} {:>8}'
    print(header.format('sample', 'reads', 'depth', 'median', 'quartiles', 'IS6110', 'copies', 'control'))
    for path in reports:
        for sample in samples_of(path):
            counts, n, bases = sample['counts'], sample['reads'], sample['bases']
            read_length = bases / n
            rates = [counts[w] / (n * (lengths[w] + read_length - 2 * K + 1) / genome_size) for w in windows]
            median = statistics.median(rates)
            quartiles = statistics.quantiles(rates, n=4)
            per_copy = n * (lengths['IS6110'] + read_length - 2 * K + 1) / genome_size
            print(header.format(
                sample['sample'], '{:,}'.format(n), '{:.1f}x'.format(bases / genome_size),
                '{:.2f}'.format(median), '{:.2f}-{:.2f}'.format(quartiles[0], quartiles[2]),
                '{:,}'.format(counts['IS6110']), '{:.1f}'.format(counts['IS6110'] / per_copy / median),
                counts['absent_300bp']))
    print('\nmedian and quartiles: counted reads divided by the reads the model expects, over the single-copy')
    print('windows. copies: IS6110 counts divided by the rate of one copy in the same sample. control: the')
    print('reference that is not in the genome, which must be 0 everywhere.')
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest='command', required=True)
    build = sub.add_parser('references')
    build.add_argument('genome')
    build.add_argument('element')
    build.add_argument('out')
    sim = sub.add_parser('simulate')
    sim.add_argument('genome')
    sim.add_argument('out')
    sim.add_argument('--reads', type=int, default=1_000_000)
    sim.add_argument('--length', type=int, default=150)
    sim.add_argument('--error', type=float, default=0.002)
    ex = sub.add_parser('exact')
    ex.add_argument('refs')
    ex.add_argument('fastq')
    ex.add_argument('report')
    mm = sub.add_parser('mismatch')
    mm.add_argument('exact_report')
    mm.add_argument('mismatch_report')
    sub_parser = sub.add_parser('subset')
    sub_parser.add_argument('fastq')
    sub_parser.add_argument('out')
    sub_parser.add_argument('--reads', type=int, default=200_000)
    rep = sub.add_parser('report')
    rep.add_argument('refs')
    rep.add_argument('--genome-size', type=int, default=4_411_532)
    rep.add_argument('reports', nargs='+')
    args = parser.parse_args()
    if args.command == 'references':
        return build_references(args.genome, args.element, args.out)
    if args.command == 'simulate':
        return simulate(args.genome, args.out, args.reads, args.length, args.error)
    if args.command == 'mismatch':
        return mismatch(args.exact_report, args.mismatch_report)
    if args.command == 'subset':
        return subset(args.fastq, args.out, args.reads)
    if args.command == 'exact':
        return exact(args.refs, args.fastq, args.report)
    return report(args.refs, args.genome_size, args.reports)


if __name__ == '__main__':
    sys.exit(main())
