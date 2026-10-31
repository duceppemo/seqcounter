import gzip

import pytest

from seqcounter.errors import SeqCounterError
from seqcounter.references import load_references, read_fasta

from .conftest import write_fasta


def test_load_references(tmp_path, references):
    path = tmp_path / 'refs.fasta'
    write_fasta(path, references)
    loaded = load_references(path, 25)
    assert [r.name for r in loaded] == ['absent', 'present', 'twice']  # File order is kept
    assert all(r.length == 60 for r in loaded)
    assert all(r.description == '' for r in loaded)


def test_gzipped_and_wrapped(tmp_path):
    """A gzipped reference file, with the sequence wrapped over several lines."""
    path = tmp_path / 'refs.fasta.gz'
    with gzip.open(path, 'wt') as f:
        f.write('>ref1 the first one\nACGTACGTAC\nGTACGTACGT\n\n>ref2\nACGT\n')
    loaded = load_references(path, 4)
    assert [(r.name, r.description, r.length) for r in loaded] == [('ref1', 'the first one', 20), ('ref2', '', 4)]


def test_read_fasta_ignores_leading_junk(tmp_path):
    path = tmp_path / 'refs.fasta'
    path.write_text('\nACGT\n>ref1\nACGT\n')  # Sequence before any ">" line
    assert read_fasta(path) == [('ref1', 4)]


def test_missing_file(tmp_path):
    with pytest.raises(SeqCounterError, match='Reference file not found'):
        load_references(tmp_path / 'missing.fasta', 25)


def test_not_fasta(tmp_path):
    path = tmp_path / 'refs.txt'
    path.write_text('ACGT\nACGT\n')
    with pytest.raises(SeqCounterError, match='No sequence found'):
        load_references(path, 25)


def test_unreadable_file(tmp_path):
    """A file that claims to be gzipped but is not."""
    path = tmp_path / 'refs.fasta.gz'
    path.write_bytes(b'\x1f\x8bnot gzip')
    with pytest.raises(SeqCounterError, match='cannot be read'):
        load_references(path, 25)


def test_empty_sequence(tmp_path):
    path = tmp_path / 'refs.fasta'
    path.write_text('>ref1\n>ref2\nACGTACGT\n')
    with pytest.raises(SeqCounterError, match='is empty'):
        load_references(path, 4)


def test_sequence_shorter_than_kmer(tmp_path):
    path = tmp_path / 'refs.fasta'
    path.write_text('>long\n{}\n>short\nACGTACGT\n'.format('A' * 30))
    with pytest.raises(SeqCounterError, match=r'shorter than the k-mer size \(25 bp\): "short" \(8 bp\)'):
        load_references(path, 25)
    assert len(load_references(path, 8)) == 2  # Exactly k is long enough


def test_duplicate_names(tmp_path):
    """Two sequences whose first word is the same: Seal would add their counts together."""
    path = tmp_path / 'refs.fasta'
    path.write_text('>ref1 copy one\nACGTACGT\n>ref1 copy two\nACGTACGT\n')
    with pytest.raises(SeqCounterError, match='have the same name: "ref1"'):
        load_references(path, 8)


def test_name_with_tabs(tmp_path):
    path = tmp_path / 'refs.fasta'
    path.write_text('>\nACGTACGT\n>a\tb\nACGTACGT\n')
    assert [r.name for r in load_references(path, 8)] == ['unnamed', 'a_b']
