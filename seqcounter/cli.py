"""Command line entry point."""

import logging
import os
import re
import shlex
import sys
import time
from argparse import ArgumentParser, ArgumentTypeError, RawDescriptionHelpFormatter
from datetime import datetime
from pathlib import Path

from . import __version__
from .counts import count_samples, print_table, run_info, write_reports
from .errors import SeqCounterError
from .references import load_references
from .samples import find_samples, single_sample
from .seal import DEFAULT_KMER_SIZE, DEFAULT_MISMATCHES, auto_memory, check_seal

log = logging.getLogger('seqcounter')

# Expected errors are reported without a traceback
USER_ERRORS = (SeqCounterError, OSError)

EPILOG = '''examples:
  seqcounter -ref refs.fasta -r1 sample_R1.fastq.gz -r2 sample_R2.fastq.gz -o results/
  seqcounter -ref refs.fasta -r1 sample.fastq.gz -o results/
  seqcounter -ref refs.fasta -r1 assembly.fasta -o results/
  seqcounter -ref refs.fasta -i folder_of_fastq_and_fasta/ -o results/

documentation: https://github.com/duceppemo/seqcounter/wiki'''


def setup_logging(verbose=False):
    """Configure the package logger only, so importing seqcounter never alters the root logger."""
    for handler in list(log.handlers):
        log.removeHandler(handler)
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)-7s %(message)s', datefmt='%H:%M:%S'))
    log.addHandler(handler)
    log.setLevel(logging.DEBUG if verbose else logging.INFO)


def positive_int(value):
    try:
        value = int(value)
    except ValueError:
        raise ArgumentTypeError('"{}" is not an integer'.format(value)) from None
    if value < 1:
        raise ArgumentTypeError('must be >= 1')
    return value


def java_memory(value):
    if not value[:-1].isdigit() or value[-1].lower() not in 'mg':
        raise ArgumentTypeError('"{}" is not a Java memory size such as 500m or 2g'.format(value))
    return value


def available_cpus():
    try:
        return len(os.sched_getaffinity(0))  # Respects cgroup/taskset limits on Linux (e.g. SLURM jobs)
    except AttributeError:
        return os.cpu_count() or 1


def format_elapsed(seconds):
    """A readable duration, as the first version of seqcounter printed it: "1h2m3s"."""
    minutes, seconds = divmod(int(round(seconds)), 60)
    hours, minutes = divmod(minutes, 60)
    days, hours = divmod(hours, 24)
    periods = [('d', days), ('h', hours), ('m', minutes), ('s', seconds)]
    return ''.join('{}{}'.format(value, name) for name, value in periods if value) or '0s'


def build_parser():
    max_cpu = available_cpus()
    parser = ArgumentParser(prog='seqcounter', formatter_class=RawDescriptionHelpFormatter, epilog=EPILOG,
                            description='Count the occurrences of reference sequences in sequencing reads (fastq)\n'
                                        'or genome assemblies (fasta), with Seal from BBTools.')
    inputs = parser.add_argument_group('input (one sample with -r1, or a folder of samples with -i)')
    source = inputs.add_mutually_exclusive_group(required=True)
    source.add_argument('-r1', '--r1', metavar='FILE',
                        help='Reads (single-end or R1 of paired-end) or assembly. fastq or fasta, gzipped or not.')
    source.add_argument('-i', '--input', metavar='FOLDER',
                        help='Folder of fastq and fasta files, searched recursively. Each fasta file, single-end '
                             'fastq file or pair of R1/R2 fastq files is a sample.')
    inputs.add_argument('-r2', '--r2', metavar='FILE',
                        help='R2 reads, for paired-end fastq files.')
    inputs.add_argument('-ref', '--ref', metavar='FILE', required=True,
                        help='Fasta file of the reference sequences to count, gzipped or not. One row per sequence '
                             'is reported, named after the first word of its ">" line.')
    inputs.add_argument('-s', '--sample', metavar='NAME',
                        help='With -r1: sample name used in the reports. Default: the input file name, without '
                             'extension and read suffix (_R1, _1, ...).')

    output = parser.add_argument_group('output')
    output.add_argument('-o', '--output', metavar='FOLDER', required=True,
                        help='Folder to hold the reports. Created if needed.')

    counting = parser.add_argument_group('counting')
    counting.add_argument('-k', '--kmer-size', metavar='N', type=positive_int, default=DEFAULT_KMER_SIZE,
                          help='K-mer size used to find the reference sequences. A reference shorter than k cannot '
                               'be found. Default: {}.'.format(DEFAULT_KMER_SIZE))
    counting.add_argument('-m', '--mismatch', metavar='N', type=int, default=DEFAULT_MISMATCHES,
                          help='Mismatches allowed between a reference k-mer and the input (maximum Hamming '
                               'distance). Default: {}.'.format(DEFAULT_MISMATCHES))

    other = parser.add_argument_group('performance and other options')
    other.add_argument('-t', '--threads', metavar='N', type=positive_int, default=max_cpu,
                       help='Number of threads. Default: all available ({}).'.format(max_cpu))
    other.add_argument('-p', '-j', '--parallel', '--jobs', metavar='N', type=positive_int, default=1, dest='parallel',
                       help='With -i: number of samples counted at the same time, sharing the threads. Default: 1.')
    other.add_argument('--memory', metavar='SIZE', type=java_memory,
                       help='Memory for Seal (Java heap size), per job, e.g. 2g. Default: from the size of the '
                            'reference sequences, within 85%% of the machine.')
    other.add_argument('--verbose', action='store_true', help='Show debug messages, including the Seal command.')
    other.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)
    return parser


def check_arguments(parser, args):
    if args.input and (args.r2 or args.sample):
        parser.error('-r2 and --sample are only used with -r1')
    if args.r1 and args.parallel > 1:
        parser.error('--parallel is only used with -i')
    if args.mismatch < 0:
        parser.error('--mismatch cannot be negative')
    if args.mismatch >= args.kmer_size:
        parser.error('--mismatch ({}) must be smaller than --kmer-size ({})'.format(args.mismatch, args.kmer_size))
    if args.sample is not None and (not args.sample.strip() or '/' in args.sample or args.sample in ('.', '..')):
        parser.error('--sample must be a name, not a path: "{}"'.format(args.sample))
    if args.sample is not None and re.search(r'[\t\r\n]', args.sample):
        parser.error('--sample cannot contain tabs or line breaks (they would break the table)')


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    check_arguments(parser, args)
    setup_logging(args.verbose)
    command = shlex.join([parser.prog, *(sys.argv[1:] if argv is None else argv)])
    log.debug('seqcounter %s: %s', __version__, command)
    started = datetime.now().astimezone()
    start = time.monotonic()
    output = Path(args.output).expanduser()
    try:
        check_seal()
        references = load_references(args.ref, args.kmer_size)
        if args.input:
            samples = find_samples(args.input, exclude=[output])
            log.info('%d sample(s) found in %s', len(samples), args.input)
        else:
            samples = [single_sample(args.r1, args.r2, name=args.sample)]
        if any(s.file_type == 'fasta' for s in samples):
            log.info('fasta input: the count of a reference sequence is the number of sequences (contigs) that '
                     'contain it, not the number of copies it has in them.')
        jobs = min(args.parallel, len(samples))
        memory = args.memory or auto_memory(sum(r.length for r in references), jobs=jobs)
        results = count_samples(samples, args.ref, jobs=jobs, threads=args.threads, kmer_size=args.kmer_size,
                                mismatch=args.mismatch, memory=memory)
        run = run_info(command, {'Reference': os.path.abspath(args.ref),
                                 'Input': ' '.join(os.path.abspath(f) for f in (args.input, args.r1, args.r2) if f),
                                 'Output folder': str(output.resolve()),
                                 'K-mer size': args.kmer_size, 'Mismatches': args.mismatch,
                                 'Threads': args.threads, 'Parallel samples': jobs, 'Seal memory': memory}, started)
        rows, reports = write_reports(results, references, output, command=command, run=run)
    except USER_ERRORS as e:
        log.error('%s', e)
        log.debug('Traceback:', exc_info=True)  # Shown with --verbose
        sys.exit(1)
    except KeyboardInterrupt:
        log.error('Interrupted')
        sys.exit(130)

    # The counts go to stdout, the messages to stderr: "seqcounter ... > counts.tsv" works
    print_table(rows)
    for report in reports:
        log.info('Report saved in %s', report)
    log.info('Done in %s', format_elapsed(time.monotonic() - start))
    failed = [r.sample for r in results if r.error]
    if failed:
        log.error('%d of %d sample(s) failed: %s', len(failed), len(results), ', '.join(failed))
        sys.exit(1)
