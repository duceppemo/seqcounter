import json
from datetime import datetime

import pytest

from seqcounter import counts
from seqcounter.counts import (
    SAMPLES_HEADER,
    Result,
    count_sample,
    count_samples,
    counts_table,
    run_info,
    write_reports,
)
from seqcounter.references import Reference
from seqcounter.samples import Sample
from seqcounter.seal import SealError, SealStats

REFERENCES = [Reference('present', 'in the genome', 60), Reference('absent', '', 60)]


def result(sample, counts_dict, **kwargs):
    return Result(sample, 'fastq', ['{}.fastq.gz'.format(sample)], counts=counts_dict, **kwargs)


def test_matched_percent():
    assert result('S1', {}, reads=200, matched=50).matched_percent == 25.0
    assert result('S1', {}, reads=0, matched=0).matched_percent is None  # An empty file: no division by zero
    assert result('S1', {}, reads=200).matched_percent is None


def test_counts_table():
    results = [result('S1', {'present': 7, 'absent': 0}), result('S2', {'present': 3, 'absent': 1})]
    assert counts_table(results, REFERENCES) == [['Reference', 'S1', 'S2'],
                                                 ['absent', '0', '1'],  # Rows sorted by reference name
                                                 ['present', '7', '3']]


def test_counts_table_with_a_failed_sample():
    """A sample that failed has empty cells: the others stay readable as numbers."""
    results = [result('S1', {'present': 7, 'absent': 0}), result('S2', {}, error='Seal failed')]
    assert counts_table(results, REFERENCES) == [['Reference', 'S1', 'S2'], ['absent', '0', ''], ['present', '7', '']]


def test_counts_table_reports_an_unknown_reference():
    """A name Seal reports that is not in the reference file still gets a row, rather than being dropped."""
    rows = counts_table([result('S1', {'present': 1, 'absent': 0, 'surprise': 2})], REFERENCES)
    assert [row[0] for row in rows[1:]] == ['absent', 'present', 'surprise']


def test_count_sample(monkeypatch):
    monkeypatch.setattr(counts, 'run_seal',
                        lambda *a, **k: SealStats({'present': 5, 'absent': 0}, reads=100, bases=10000, matched=5))
    counted = count_sample(Sample('S1', 'fastq', ['S1.fastq.gz']), 'refs.fasta')
    assert counted.counts == {'present': 5, 'absent': 0}
    assert (counted.reads, counted.bases, counted.matched, counted.error) == (100, 10000, 5, None)


def test_count_sample_keeps_an_error(monkeypatch, caplog):
    """One sample failing does not stop a run over a folder."""
    def fail(*a, **k):
        raise SealError('Seal failed (exit code 1): out of memory')

    monkeypatch.setattr(counts, 'run_seal', fail)
    counted = count_sample(Sample('S1', 'fastq', ['S1.fastq.gz']), 'refs.fasta')
    assert counted.counts == {} and 'out of memory' in counted.error
    assert 'out of memory' in caplog.text


def test_count_sample_without_totals(monkeypatch):
    monkeypatch.setattr(counts, 'run_seal', lambda *a, **k: SealStats({'present': 5}))
    assert count_sample(Sample('S1', 'fastq', ['S1.fastq.gz']), 'refs.fasta').reads is None


def test_count_samples_shares_the_threads(monkeypatch):
    seen = []

    def fake_run_seal(inputs, refs, threads=1, memory='1g', k=25, hdist=1):
        seen.append((inputs[0], threads, k, hdist, memory))
        return SealStats({'present': 1})

    monkeypatch.setattr(counts, 'run_seal', fake_run_seal)
    samples = [Sample('S{}'.format(i), 'fastq', ['S{}.fastq.gz'.format(i)]) for i in range(4)]
    results = count_samples(samples, 'refs.fasta', jobs=2, threads=9, kmer_size=31, mismatch=0, memory='2g')
    assert [r.sample for r in results] == ['S0', 'S1', 'S2', 'S3']  # In the order given
    assert {t[1] for t in seen} == {4}  # 9 threads over 2 jobs
    assert {(t[2], t[3], t[4]) for t in seen} == {(31, 0, '2g')}


def test_count_samples_never_asks_for_more_jobs_than_samples(monkeypatch):
    seen = []
    monkeypatch.setattr(counts, 'run_seal', lambda inputs, refs, threads=1, **k: seen.append(threads)
                        or SealStats({}))
    count_samples([Sample('S1', 'fastq', ['S1.fastq.gz'])], 'refs.fasta', jobs=8, threads=4)
    assert seen == [4]  # One sample: all the threads, not 4 // 8 = 0


def test_write_reports(tmp_path, monkeypatch):
    monkeypatch.setattr(counts, 'versions', lambda: {'BBTools': '39.01', 'Java': 'openjdk 17'})
    results = [result('S1', {'present': 7, 'absent': 0}, reads=200, bases=20000, matched=7),
               result('S2', {}, error='Seal failed (exit code 1): boom\nReads:\t100\tTime:\t1.2s')]
    run = run_info('seqcounter -ref refs.fasta -i reads/ -o out/', {'K-mer size': 25},
                   datetime(2026, 10, 30, 12, 0).astimezone())
    rows, reports = write_reports(results, REFERENCES, tmp_path / 'out', command=run['command'], run=run)
    assert [p.name for p in reports] == ['refseqcounts.tsv', 'refseqcounts_samples.tsv', 'refseqcounts.json']
    assert rows[0] == ['Reference', 'S1', 'S2']

    # The counts table, with the command line of the run as its first line, as the first version of seqcounter wrote it
    lines = reports[0].read_text().splitlines()
    assert lines[0] == '#seqcounter -ref refs.fasta -i reads/ -o out/'
    assert lines[1:] == ['Reference\tS1\tS2', 'absent\t0\t', 'present\t7\t']

    samples = [line.split('\t') for line in reports[1].read_text().splitlines()]
    assert samples[0] == ['Sample', 'Type', 'Reads', 'Bases', 'Matched', 'MatchedPct', 'Files', 'Error']
    assert samples[1] == ['S1', 'fastq', '200', '20000', '7', '3.5000', 'S1.fastq.gz', '']
    assert samples[2][:6] == ['S2', 'fastq', '', '', '', '']
    # Seal's output is tab aligned and several lines long: it must not add columns to a tab separated file
    assert len(samples[2]) == len(SAMPLES_HEADER)
    assert samples[2][7] == 'Seal failed (exit code 1): boom Reads: 100 Time: 1.2s'

    report = json.loads(reports[2].read_text())
    assert report['tool'] == 'seqcounter' and report['command'] == run['command']
    assert report['software']['BBTools'] == '39.01' and report['parameters'] == {'K-mer size': 25}
    assert report['references'][0] == {'name': 'present', 'description': 'in the genome', 'length': 60}
    assert report['samples'][0]['counts'] == {'present': 7, 'absent': 0}
    assert report['samples'][1]['error'].startswith('Seal failed')
    assert report['started'] < report['finished']


def test_write_reports_without_a_command(tmp_path):
    rows, reports = write_reports([result('S1', {'present': 1, 'absent': 0})], REFERENCES, tmp_path / 'out')
    assert rows[0] == ['Reference', 'S1']
    assert reports[0].read_text().splitlines()[0] == 'Reference\tS1'  # No "#command" line
    assert 'command' not in json.loads(reports[2].read_text())


def test_print_table(capsys):
    counts.print_table([['Reference', 'S1'], ['present', '7']])
    assert capsys.readouterr().out == 'Reference\tS1\npresent\t7\n'


@pytest.mark.parametrize('memory', ['1g', '512m'])
def test_run_info(monkeypatch, memory):
    monkeypatch.setattr(counts, 'versions', lambda: {'BBTools': '39.01', 'Java': 'openjdk 17'})
    info = run_info('seqcounter -h', {'Seal memory': memory}, datetime(2026, 1, 1, 9, 0).astimezone())
    assert info['parameters']['Seal memory'] == memory
    assert info['software']['seqcounter'] and info['software']['Python']
    assert info['started'].startswith('2026-01-01T09:00:00')
