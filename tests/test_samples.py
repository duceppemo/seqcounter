import pytest

from seqcounter.errors import SeqCounterError
from seqcounter.samples import find_samples, sample_name, single_sample, split_extension


def touch(folder, *names):
    for name in names:
        path = folder / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text('x')


def test_split_extension():
    assert split_extension('S1_R1.fastq.gz') == ('S1_R1', 'fastq')
    assert split_extension('a.b.FNA') == ('a.b', 'fasta')
    assert split_extension('notes.txt') is None
    assert split_extension('.fastq') is None


def test_sample_name():
    assert sample_name('/data/S1_R1.fastq.gz') == 'S1'
    assert sample_name('S1_S7_L001_R2_001.fastq') == 'S1_S7_L001'
    assert sample_name('Iso_R1.fasta') == 'Iso_R1'  # No read suffix on a fasta file
    assert sample_name('E.coli.v2.fasta') == 'E.coli.v2'
    assert sample_name('reads.dat') == 'reads.dat'  # An unknown extension is kept: Seal reads the content
    assert sample_name('_R1.fastq') == '_R1'  # Stripping the suffix would leave nothing


def test_find_samples(tmp_path):
    touch(tmp_path, 'run1/S1_R1.fastq.gz', 'run1/S1_R2.fastq.gz',  # Paired
          'run1/S2_S7_L001_R1_001.fastq.gz', 'run1/S2_S7_L001_R2_001.fastq.gz',  # Illumina names
          'run2/SRR1_1.fq', 'run2/SRR1_2.fq',  # SRA names
          'run2/single.fastq',  # Single-end
          'assemblies/Iso_R1.fasta', 'assemblies/Iso_R10.fna.gz',  # Fasta: no read suffix
          'notes.txt', '.hidden.fasta', '.snapshots/old.fasta')
    samples = {s.name: s for s in find_samples(tmp_path)}
    assert sorted(samples) == ['Iso_R1', 'Iso_R10', 'S1', 'S2_S7_L001', 'SRR1', 'single']
    assert [p.split('/')[-1] for p in samples['S1'].files] == ['S1_R1.fastq.gz', 'S1_R2.fastq.gz']
    assert samples['S1'].paired and samples['SRR1'].paired and samples['S2_S7_L001'].paired
    assert not samples['single'].paired
    assert samples['Iso_R1'].file_type == 'fasta'


def test_r2_without_r1(tmp_path, caplog):
    touch(tmp_path, 'lonely_R2.fastq')
    samples = find_samples(tmp_path)
    assert [s.name for s in samples] == ['lonely'] and not samples[0].paired
    assert 'R2 file without R1' in caplog.text


def test_output_folder_excluded(tmp_path):
    touch(tmp_path, 'S1.fasta', 'results/S2.fasta')
    assert [s.name for s in find_samples(tmp_path, exclude=[tmp_path / 'results'])] == ['S1']


@pytest.mark.parametrize('names', [
    ['a/S1.fasta', 'b/S1.fasta'],  # Same name in two folders
    ['S1.fasta', 'S1_R1.fastq', 'S1_R2.fastq'],  # Fasta and fastq
    ['a/S1_R1.fastq', 'b/S1_R1.fastq'],  # Two R1
    ['S1.fastq', 'S1_R1.fastq'],  # Single-end and R1
])
def test_ambiguous_names(tmp_path, names):
    touch(tmp_path, *names)
    with pytest.raises(SeqCounterError, match='same sample name'):
        find_samples(tmp_path)


def test_no_samples(tmp_path):
    touch(tmp_path, 'notes.txt')
    with pytest.raises(SeqCounterError, match='No fasta or fastq'):
        find_samples(tmp_path)
    with pytest.raises(SeqCounterError, match='not found'):
        find_samples(tmp_path / 'missing')


def test_symlinked_folders(tmp_path):
    touch(tmp_path, 'real/S1.fasta', 'input/S2.fasta')
    (tmp_path / 'input' / 'linked').symlink_to(tmp_path / 'real')
    (tmp_path / 'input' / 'loop').symlink_to(tmp_path / 'input')  # Must not loop forever
    assert [s.name for s in find_samples(tmp_path / 'input')] == ['S1', 'S2']


def test_names_without_tabs(tmp_path):
    """Tabs and line breaks in file names would break the table."""
    assert sample_name('my\tsample\n1.fasta') == 'my_sample_1'
    (tmp_path / 'my\tsample.fasta').write_text('>a\nACGT\n')
    assert [s.name for s in find_samples(tmp_path)] == ['my_sample']
    (tmp_path / 'my_sample.fasta').write_text('>a\nACGT\n')
    with pytest.raises(SeqCounterError, match='same sample name'):
        find_samples(tmp_path)


def test_single_sample(tmp_path):
    touch(tmp_path, 'S1_R1.fastq.gz', 'S1_R2.fastq.gz', 'asm.fasta', 'reads.dat')
    sample = single_sample(tmp_path / 'S1_R1.fastq.gz', tmp_path / 'S1_R2.fastq.gz')
    assert sample.name == 'S1' and sample.paired and sample.file_type == 'fastq'
    assert all(p.startswith('/') for p in sample.files)  # Absolute paths
    assert single_sample(tmp_path / 'asm.fasta').file_type == 'fasta'
    assert single_sample(tmp_path / 'asm.fasta', name='mine').name == 'mine'
    assert single_sample(tmp_path / 'reads.dat').file_type == 'fastq'  # Unknown extension: Seal reads the content


def test_single_sample_errors(tmp_path):
    touch(tmp_path, 'S1_R1.fastq.gz', 'asm.fasta')
    with pytest.raises(SeqCounterError, match='Input file not found'):
        single_sample(tmp_path / 'missing.fastq.gz')
    with pytest.raises(SeqCounterError, match='Input file not found'):
        single_sample(tmp_path / 'S1_R1.fastq.gz', tmp_path / 'missing.fastq.gz')
    with pytest.raises(SeqCounterError, match='-r2 is for paired-end fastq'):
        single_sample(tmp_path / 'asm.fasta', tmp_path / 'S1_R1.fastq.gz')
    with pytest.raises(SeqCounterError, match='are the same file'):
        single_sample(tmp_path / 'S1_R1.fastq.gz', tmp_path / 'S1_R1.fastq.gz')
