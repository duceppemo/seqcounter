"""End to end runs with the real Seal, from the command line entry point. Skipped without BBTools installed."""

import gzip
import json
import subprocess
import sys
from pathlib import Path

import pytest

from seqcounter.cli import main

pytestmark = pytest.mark.usefixtures('seal_or_skip')

REPO = Path(__file__).resolve().parent.parent


def run(*argv):
    """Run the command line as the installed "seqcounter" does, and return its exit code."""
    try:
        main([str(a) for a in argv])
    except SystemExit as exit_info:
        return exit_info.code or 0
    return 0


def read_counts(folder):
    """{reference name: {sample: count}} from refseqcounts.tsv, and its "#command" line."""
    lines = (Path(folder) / 'refseqcounts.tsv').read_text().splitlines()
    command = lines.pop(0)
    header = lines.pop(0).split('\t')[1:]
    table = {}
    for line in lines:
        fields = line.split('\t')
        table[fields[0]] = dict(zip(header, fields[1:], strict=True))
    return table, command


def test_assembly(data, tmp_path, capsys):
    """One assembly: the count is the number of contigs containing the reference, and 0 for the absent one."""
    out = tmp_path / 'out'
    assert run('-ref', data / 'refs.fasta', '-r1', data / 'samples' / 'assembly.fasta', '-o', out) == 0
    table, command = read_counts(out)
    assert table == {'absent': {'assembly': '0'}, 'present': {'assembly': '1'}, 'twice': {'assembly': '1'}}
    assert command.startswith('#seqcounter -ref ')
    assert 'assembly' in capsys.readouterr().out  # The table is printed too

    samples = [line.split('\t') for line in (out / 'refseqcounts_samples.tsv').read_text().splitlines()]
    assert samples[1][:2] == ['assembly', 'fasta'] and samples[1][2] == '1'  # One contig
    report = json.loads((out / 'refseqcounts.json').read_text())
    assert report['samples'][0]['counts'] == {'absent': 0, 'present': 1, 'twice': 1}
    assert report['samples'][0]['matched'] == 1 and report['samples'][0]['error'] is None
    assert report['software']['BBTools'] != 'unknown'
    assert report['parameters']['K-mer size'] == 25


def test_single_end_reads(data, tmp_path):
    out = tmp_path / 'out'
    assert run('-ref', data / 'refs.fasta', '-r1', data / 'samples' / 'single.fastq.gz', '-o', out, '-s', 'mine') == 0
    table, _ = read_counts(out)
    assert set(table) == {'absent', 'present', 'twice'}
    assert table['absent'] == {'mine': '0'}
    # "twice" is in the genome twice, so about twice as many reads match it as "present"
    assert int(table['present']['mine']) > 5
    assert int(table['twice']['mine']) > int(table['present']['mine'])


def test_paired_end_reads_are_counted_once_each(data, tmp_path):
    """Both mates are read: the counts are those of the two files together, not of R1 twice."""
    out = tmp_path / 'out'
    assert run('-ref', data / 'refs.fasta', '-r1', data / 'samples' / 'paired_R1.fastq.gz', '-r2', data / 'samples' / 'paired_R2.fastq.gz',
               '-o', out) == 0
    table, _ = read_counts(out)
    assert list(table['present']) == ['paired']  # R1 and R2 are one sample
    report = json.loads((out / 'refseqcounts.json').read_text())
    assert len(report['samples'][0]['files']) == 2
    reads = report['samples'][0]['reads']

    single = tmp_path / 'single_out'
    assert run('-ref', data / 'refs.fasta', '-r1', data / 'samples' / 'paired_R1.fastq.gz', '-o', single) == 0
    assert json.loads((single / 'refseqcounts.json').read_text())['samples'][0]['reads'] == reads / 2


def test_folder_of_samples(data, tmp_path):
    out = tmp_path / 'out'
    assert run('-ref', data / 'refs.fasta.gz', '-i', data / 'samples', '-o', out, '-p', '2', '-t', '4') == 0
    table, _ = read_counts(out)
    assert sorted(table['present']) == ['assembly', 'nothing', 'paired', 'single']
    assert table['present']['nothing'] == '0' and table['twice']['nothing'] == '0'
    assert table['present']['assembly'] == '1'
    assert int(table['twice']['single']) > 0


def test_output_folder_inside_the_input_folder(data, tmp_path):
    """The reports of an earlier run are not taken for samples: the output folder is skipped."""
    folder = tmp_path / 'reads'
    folder.mkdir()
    (folder / 'S1.fasta').write_text((data / 'samples' / 'assembly.fasta').read_text())
    out = folder / 'results'
    assert run('-ref', data / 'refs.fasta', '-i', folder, '-o', out) == 0
    assert run('-ref', data / 'refs.fasta', '-i', folder, '-o', out) == 0
    table, _ = read_counts(out)
    assert list(table['present']) == ['S1']


def test_mismatches_and_kmer_size(data, tmp_path, references):
    """A reference with one mismatch is found with -m 1 (the default), not with -m 0."""
    ref = tmp_path / 'mutated.fasta'
    mutated = references['present']
    mutated = mutated[:30] + ('T' if mutated[30] != 'T' else 'A') + mutated[31:]
    ref.write_text('>mutated\n{}\n'.format(mutated))

    with_mismatch = tmp_path / 'with'
    assert run('-ref', ref, '-r1', data / 'samples' / 'assembly.fasta', '-o', with_mismatch, '-m', '1') == 0
    assert read_counts(with_mismatch)[0]['mutated'] == {'assembly': '1'}

    without = tmp_path / 'without'
    assert run('-ref', ref, '-r1', data / 'samples' / 'assembly.fasta', '-o', without, '-m', '0', '-k', '31') == 0
    assert read_counts(without)[0]['mutated'] == {'assembly': '0'}  # No 31-mer of it is exact


def test_unsafe_file_names(data, tmp_path):
    """A file name Seal cannot be given safely (spaces, quotes, "xmx") is still counted."""
    folder = tmp_path / 'in put'
    folder.mkdir()
    name = "it's xmx 1.fasta"
    (folder / name).write_text((data / 'samples' / 'assembly.fasta').read_text())
    refs = folder / 'my refs.fasta'
    refs.write_text((data / 'refs.fasta').read_text())
    out = tmp_path / 'out'
    assert run('-ref', refs, '-r1', folder / name, '-o', out) == 0
    assert read_counts(out)[0]['present'] == {"it's xmx 1": '1'}


def test_fastq_named_fasta(data, tmp_path):
    """Seal recognizes formats by extension only: a fastq file named .fasta is read as fastq anyway."""
    misnamed = tmp_path / 'misnamed.fasta'
    with gzip.open(data / 'samples' / 'single.fastq.gz', 'rt') as f:
        misnamed.write_text(f.read())
    out = tmp_path / 'out'
    assert run('-ref', data / 'refs.fasta', '-r1', misnamed, '-o', out) == 0
    report = json.loads((out / 'refseqcounts.json').read_text())
    assert report['samples'][0]['reads'] > 100  # Read as reads, not as one long contig


def test_one_sample_fails(data, tmp_path, caplog):
    """A broken file does not stop the run: the other samples are reported, and the exit code is 1."""
    folder = tmp_path / 'reads'
    folder.mkdir()
    (folder / 'good.fasta').write_text((data / 'samples' / 'assembly.fasta').read_text())
    (folder / 'broken.fastq.gz').write_bytes(b'\x1f\x8b' + b'not gzip at all' * 10)
    out = tmp_path / 'out'
    assert run('-ref', data / 'refs.fasta', '-i', folder, '-o', out) == 1
    table, _ = read_counts(out)
    assert table['present'] == {'broken': '', 'good': '1'}
    assert '1 of 2 sample(s) failed: broken' in caplog.text
    samples = {line.split('\t')[0]: line.split('\t') for line in
               (out / 'refseqcounts_samples.tsv').read_text().splitlines()}
    assert samples['broken'][-1]  # The error is reported in the samples table
    assert json.loads((out / 'refseqcounts.json').read_text())['samples'][0]['error']


def test_reference_shorter_than_the_kmer_size(data, tmp_path, caplog):
    out = tmp_path / 'out'
    assert run('-ref', data / 'refs.fasta', '-r1', data / 'samples' / 'assembly.fasta', '-o', out, '-k', '61') == 1
    assert 'shorter than the k-mer size' in caplog.text
    assert not out.exists()  # Nothing is written when the run cannot start


def test_entry_points(data, tmp_path):
    """The module and the seqcounter.py shim of a clone run the same command line."""
    shim = ['seqcounter.py'] if (REPO / 'seqcounter.py').is_file() else []  # Absent from the distribution
    for argv in (['-m', 'seqcounter'], *([shim] if shim else [])):
        out = tmp_path / ('out' + argv[-1])
        proc = subprocess.run([sys.executable, *argv, '-ref', str(data / 'refs.fasta'),
                               '-r1', str(data / 'samples' / 'assembly.fasta'), '-o', str(out)],
                              capture_output=True, text=True, cwd=REPO)
        assert proc.returncode == 0, proc.stderr
        assert 'present\t1' in proc.stdout
        assert read_counts(out)[0]['absent'] == {'assembly': '0'}
