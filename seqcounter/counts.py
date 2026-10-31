"""Count the reference sequences in each sample, and write the reports."""

import json
import logging
import platform
import sys
from concurrent import futures
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

from . import __version__
from .errors import SeqCounterError
from .seal import DEFAULT_KMER_SIZE, DEFAULT_MISMATCHES, run_seal, versions

log = logging.getLogger(__name__)

COUNTS_FILE = 'refseqcounts.tsv'  # The name the first version of seqcounter wrote, and the same table
SAMPLES_FILE = 'refseqcounts_samples.tsv'
JSON_FILE = 'refseqcounts.json'


@dataclass
class Result:
    """What one sample gave: a count per reference sequence, or an error."""
    sample: str
    file_type: str
    files: list
    counts: dict = field(default_factory=dict)
    reads: int | None = None  # Reads, or contigs for a fasta file
    bases: int | None = None
    matched: int | None = None  # Reads matching at least one reference sequence
    error: str | None = None

    @property
    def matched_percent(self):
        if self.matched is None or not self.reads:
            return None
        return 100.0 * self.matched / self.reads


def count_sample(sample, ref, kmer_size=DEFAULT_KMER_SIZE, mismatch=DEFAULT_MISMATCHES, threads=1, memory='1g'):
    """
    Count the reads (or contigs) of one sample matching each reference sequence.

    An error on one sample is kept in its Result rather than raised, so that a run over a folder reports the
    samples that worked.
    """
    log.info('Counting %s (%s, %d file(s))', sample.name, sample.file_type, len(sample.files))
    try:
        stats = run_seal(sample.files, ref, threads=threads, memory=memory, k=kmer_size, hdist=mismatch)
    except SeqCounterError as e:
        log.error('%s: %s', sample.name, e)
        return Result(sample.name, sample.file_type, list(sample.files), error=str(e))
    result = Result(sample.name, sample.file_type, list(sample.files), counts=stats.counts,
                    reads=stats.reads, bases=stats.bases, matched=stats.matched)
    log.info('%s: %s of %s read(s) matched', sample.name,
             '{:,}'.format(stats.matched) if stats.matched is not None else '?',
             '{:,}'.format(stats.reads) if stats.reads is not None else '?')
    return result


def count_samples(samples, ref, jobs=1, threads=1, **options):
    """
    Count the reference sequences in several samples, up to jobs samples at a time.

    The threads are shared between the jobs, as Seal uses them within one sample.

    :return: list of Result, in the order of samples
    """
    jobs = max(1, min(jobs, len(samples)))
    per_job = max(1, threads // jobs)
    if jobs > 1:
        log.info('%d sample(s), %d at a time, %d thread(s) each', len(samples), jobs, per_job)
    with futures.ThreadPoolExecutor(max_workers=jobs) as executor:
        return list(executor.map(lambda s: count_sample(s, ref, threads=per_job, **options), samples))


def counts_table(results, references):
    """
    The table written to refseqcounts.tsv and printed: one row per reference sequence, one column per sample.

    Rows are sorted by reference name, as the first version of seqcounter sorted them. A sample that failed has
    empty cells, so that the samples that worked are still readable as numbers.

    :return: list of rows, the first one being the header
    """
    names = sorted({r.name for r in references} | {name for result in results for name in result.counts})
    rows = [['Reference'] + [result.sample for result in results]]
    for name in names:
        row = [name]
        for result in results:
            row.append('' if result.error else str(result.counts.get(name, 0)))
        rows.append(row)
    return rows


def write_counts(results, references, path, command=None):
    """Write the counts table, keeping the "#command line" first line of the first version of seqcounter."""
    rows = counts_table(results, references)
    with open(path, 'w') as f:
        if command:
            f.write('#{}\n'.format(command))
        for row in rows:
            f.write('\t'.join(row) + '\n')
    return rows


SAMPLES_HEADER = ['Sample', 'Type', 'Reads', 'Bases', 'Matched', 'MatchedPct', 'Files', 'Error']


def write_samples(results, path):
    """One row per sample: how much was read, how much matched, and the error if it failed."""
    with open(path, 'w') as f:
        f.write('\t'.join(SAMPLES_HEADER) + '\n')
        for r in results:
            percent = r.matched_percent
            f.write('\t'.join([r.sample, r.file_type,
                               '' if r.reads is None else str(r.reads),
                               '' if r.bases is None else str(r.bases),
                               '' if r.matched is None else str(r.matched),
                               '' if percent is None else '{:.4f}'.format(percent),
                               ';'.join(r.files),
                               (r.error or '').replace('\n', ' ')]) + '\n')


def write_json(results, references, path, run=None):
    """Everything the run produced, for a pipeline: parameters, versions, references and counts."""
    report = {
        'tool': 'seqcounter',
        'version': __version__,
        'references': [{'name': r.name, 'description': r.description, 'length': r.length} for r in references],
        'samples': [{'sample': r.sample, 'type': r.file_type, 'files': r.files, 'reads': r.reads, 'bases': r.bases,
                     'matched': r.matched, 'matched_percent': r.matched_percent, 'counts': r.counts,
                     'error': r.error} for r in results],
    }
    report.update(run or {})
    with open(path, 'w') as f:
        json.dump(report, f, indent=2, sort_keys=False)
        f.write('\n')


def run_info(command, parameters, started, finished=None):
    """The run metadata of the JSON report: what was run, with which versions and parameters."""
    software = {'seqcounter': __version__, 'Python': platform.python_version()}
    software.update(versions())
    return {'command': command, 'started': started.isoformat(timespec='seconds'),
            'finished': (finished or datetime.now().astimezone()).isoformat(timespec='seconds'),
            'software': software, 'parameters': parameters}


def print_table(rows, stream=None):
    """Print the counts table, tab separated: "seqcounter ... > counts.tsv" works."""
    stream = stream or sys.stdout
    for row in rows:
        print('\t'.join(row), file=stream)
    stream.flush()


def write_reports(results, references, output, command=None, run=None):
    """Write the three reports in the output folder and return their paths, the counts table first."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    counts, samples, report = output / COUNTS_FILE, output / SAMPLES_FILE, output / JSON_FILE
    rows = write_counts(results, references, counts, command=command)
    write_samples(results, samples)
    write_json(results, references, report, run=run)
    return rows, [counts, samples, report]
