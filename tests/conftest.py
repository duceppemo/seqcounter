"""Synthetic reference sequences, genomes and reads with known counts."""

import gzip
import random
import shutil

import pytest

REF_LENGTH = 60  # Longer than the default k-mer size of 25


def random_seq(rng, length):
    return ''.join(rng.choice('ACGT') for _ in range(length))


def reverse_complement(seq):
    return seq[::-1].translate(str.maketrans('ACGT', 'TGCA'))


def write_fasta(path, records):
    with open(path, 'w') as f:
        for name, seq in records.items():
            f.write('>{}\n{}\n'.format(name, seq))


def reads(seq, rng, depth, length=100):
    """Random reads over seq, on both strands."""
    for _ in range(int(depth * len(seq) / length)):
        p = rng.randint(0, len(seq) - length)
        read = seq[p:p + length]
        yield reverse_complement(read) if rng.random() < 0.5 else read


def write_reads(path, seq, rng, depth=20, length=100, mate=None):
    """Gzipped fastq."""
    with gzip.open(path, 'wt') as f:
        for i, read in enumerate(reads(seq, rng, depth, length)):
            if mate == 2:
                read = reverse_complement(read)
            f.write('@r{}/{}\n{}\n+\n{}\n'.format(i, mate or 1, read, 'I' * length))


@pytest.fixture(scope='session')
def references():
    """Three reference sequences: "present" and "twice" are in the genome, "absent" is not."""
    rng = random.Random(7)
    return {name: random_seq(rng, REF_LENGTH) for name in ('absent', 'present', 'twice')}


@pytest.fixture(scope='session')
def data(tmp_path_factory, references):
    """
    The reference files, and a "samples" folder whose counts are known: an assembly, single-end reads,
    paired-end reads, and a sample containing none of the reference sequences.
    """
    rng = random.Random(11)
    folder = tmp_path_factory.mktemp('data')
    samples = folder / 'samples'
    samples.mkdir()
    write_fasta(folder / 'refs.fasta', references)
    with gzip.open(folder / 'refs.fasta.gz', 'wt') as f:
        f.write((folder / 'refs.fasta').read_text())

    # The genome holds "present" once and "twice" twice, with "absent" nowhere
    genome = ''.join([random_seq(rng, 1000), references['present'], random_seq(rng, 1000), references['twice'],
                      random_seq(rng, 1000), references['twice'], random_seq(rng, 1000)])
    empty_genome = random_seq(rng, 4000)
    write_fasta(samples / 'assembly.fasta', {'contig1': genome})
    write_fasta(samples / 'nothing.fasta', {'contig1': empty_genome})
    write_reads(samples / 'single.fastq.gz', genome, rng, depth=30)
    write_reads(samples / 'paired_R1.fastq.gz', genome, rng, depth=15, mate=1)
    write_reads(samples / 'paired_R2.fastq.gz', genome, rng, depth=15, mate=2)
    return folder


@pytest.fixture(scope='session')
def seal_or_skip():
    """Skip the test when seal.sh is not installed: the unit tests do not need BBTools."""
    if shutil.which('seal.sh') is None:
        pytest.skip('seal.sh (BBTools) is not installed')
    return shutil.which('seal.sh')
