"""Count the reads (or contigs) matching each reference sequence with Seal, from BBTools."""

import gzip
import logging
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from functools import cache
from pathlib import Path

from .errors import SeqCounterError

log = logging.getLogger(__name__)

DEFAULT_KMER_SIZE = 25  # A reference sequence shorter than k cannot be found
DEFAULT_MISMATCHES = 1  # Maximum Hamming distance between a reference k-mer and the input


class SealError(SeqCounterError):
    pass


@cache
def executable():
    """
    Path of seal.sh: the one in the PATH, otherwise the one installed next to the running Python.
    The fallback makes "/path/to/env/bin/seqcounter" work without activating the conda environment.

    :return: path to seal.sh, or None if it cannot be found
    """
    found = shutil.which('seal.sh')
    if found:
        return found
    # Do not resolve symlinks: a virtualenv's python is a symlink to the base interpreter
    candidate = Path(sys.executable).parent / 'seal.sh'
    if candidate.is_file() and os.access(candidate, os.X_OK):
        return str(candidate)
    return None


def check_seal():
    """Make sure seal.sh is installed and return its path."""
    path = executable()
    if path is None:
        raise SealError('"seal.sh" was not found in your PATH or in "{}". Install BBTools in the environment '
                        'with "conda install -c bioconda bbmap".'.format(Path(sys.executable).parent))
    return path


def seal_command(inputs, refs, stats_file, threads, memory, k=DEFAULT_KMER_SIZE, hdist=DEFAULT_MISMATCHES):
    """
    :param inputs: one fasta/fastq file, or two paired-end fastq files
    :param refs: reference fasta file(s): a path or a list of paths
    :param k: k-mer size. A reference sequence shorter than k cannot be found
    :param hdist: number of mismatches allowed in a k-mer
    """
    refs = [refs] if isinstance(refs, (str, os.PathLike)) else list(refs)
    cmd = [check_seal(), '-Xmx{}'.format(memory), 'in={}'.format(inputs[0])]
    if len(inputs) > 1:
        cmd.append('in2={}'.format(inputs[1]))
    else:
        cmd.append('int=f')  # Never guess that a single fastq file is interleaved: that would halve the counts
    cmd += ['ref={}'.format(','.join(str(r) for r in refs)),
            'k={}'.format(k),
            'rcomp=t',  # Count a match on either strand
            'hdist={}'.format(hdist),
            'maskmiddle=f',  # Do not treat the middle base of a k-mer as a wildcard
            'clearzone=999999',
            'ambiguous=all',  # Count a read for every reference sequence it matches
            'nzo=f',  # Also report reference sequences with no match
            'qin=33',  # Force the quality encoding: autodetection fails on some low quality nanopore reads
            'ow=t',
            'stats={}'.format(stats_file),
            'threads={}'.format(threads)]
    return cmd


# Paths that seal.sh cannot pass to Seal safely: it runs its command line with eval, so any shell character (spaces,
# quotes, $, `, ;, &, |, parentheses...) is re-interpreted, it splits ref= lists on commas, and BBTools 40 takes any
# argument containing "xmx" or "xms" (in any case, anywhere) for a Java memory setting. Only plain characters are
# passed as they are: the other paths are linked under a neutral name.
UNSAFE = re.compile(r'[^A-Za-z0-9._/+:@-]|xm[xs]', re.IGNORECASE)


def extensions(name):
    """The extensions Seal uses to recognize the format, e.g. ".fastq.gz" for "S1_R1.fastq.gz"."""
    suffixes = [s for s in Path(name).suffixes[-2:] if re.fullmatch(r'\.[A-Za-z0-9]+', s) and not UNSAFE.search(s)]
    return ''.join(suffixes)


SEQUENCE_EXTENSION = re.compile(r'\.(fastq|fq|fasta|fa|fna|fas|fsa)(\.gz)?$', re.IGNORECASE)


def content_extension(path):
    """
    Extension matching the content of a sequence file, e.g. ".fastq.gz". Seal recognizes formats and compression by
    extension only, and Galaxy, for example, names its files "dataset_1.dat".
    """
    with open(path, 'rb') as f:
        gzipped = f.read(2) == b'\x1f\x8b'
    try:
        with (gzip.open(path, 'rt') if gzipped else open(path)) as f:
            first = next((line for line in f if line.strip()), '>')[0]
    except (OSError, EOFError, UnicodeDecodeError):
        first = '>'  # Unreadable: let Seal report the error
    return ('.fastq' if first == '@' else '.fasta') + ('.gz' if gzipped else '')


def safe_paths(paths, folder, prefix, check_extension=False):
    """
    Paths that Seal can read: unsafe ones, and with check_extension those whose extension does not match their
    content, are linked into folder under a neutral name (e.g. "in0.fastq.gz").
    """
    safe = []
    for i, path in enumerate(paths):
        path = str(path)
        extension = extensions(path)
        if check_extension:
            detected = content_extension(path)
            match = SEQUENCE_EXTENSION.search(Path(path).name)
            fastq = match and match.group(1).lower() in ('fastq', 'fq')
            if not match or fastq != detected.startswith('.fastq') or bool(match.group(2)) != detected.endswith('.gz'):
                extension = detected
        if UNSAFE.search(path) or extension != extensions(path):
            link = Path(folder) / '{}{}{}'.format(prefix, i, extension)
            link.symlink_to(Path(path).resolve())
            path = str(link)
        safe.append(path)
    return safe


def safe_temporary_folder(attempts=20):
    """A new temporary folder whose path Seal can read: random names can contain "xmx" or "xms", so retry."""
    for _ in range(attempts):
        folder = tempfile.mkdtemp(prefix='seqcounter_')
        if not UNSAFE.search(folder):
            return folder
        os.rmdir(folder)
    raise SealError('The temporary folder "{}" contains a space, a comma, "xmx" or "xms", which Seal cannot handle. '
                    'Set TMPDIR to another folder.'.format(tempfile.gettempdir()))


GENERIC_ERRORS = re.compile(r'Exception in thread "main"|terminated in an error state')


def failure_reason(lines):
    """
    The most informative line of Seal's error output: the first error that is not Seal's generic final message.
    Java exceptions often end with ":" and give their message on the next line, which is then added.
    """
    candidates = [i for i, line in enumerate(lines) if re.search(r'Exception|Error|error', line)
                  and not line.strip().startswith('at ')]
    specific = [i for i in candidates if not GENERIC_ERRORS.search(lines[i])]
    if not candidates:
        return lines[-1].strip()
    i = (specific or candidates)[0]
    reason = lines[i].strip()
    if reason.endswith(':'):
        following = [line.strip() for line in lines[i + 1:] if line.strip() and not line.strip().startswith('at ')]
        if following:
            reason += ' ' + following[0]
    return reason


@dataclass
class SealStats:
    counts: dict = field(default_factory=dict)  # {reference name: number of matching reads}
    reads: int | None = None  # Total number of reads (or contigs) in the input
    bases: int | None = None
    matched: int | None = None  # Reads matching at least one reference sequence


def parse_stats(stats_file):
    """
    Read Seal's stats file.

    #File	sample_R1.fastq.gz
    #Total	822714	123407100
    #Matched	799	0.09712%
    #Name	Reads	ReadsPct
    spacer25	62	0.00754%
    """
    stats = SealStats()
    with open(stats_file) as f:
        for line in f:
            fields = line.rstrip('\n').split('\t')
            if not line.strip():
                continue
            try:
                if fields[0] == '#Total':
                    stats.reads = int(fields[1])
                    stats.bases = int(fields[2]) if len(fields) > 2 else None
                elif fields[0] == '#Matched':
                    stats.matched = int(fields[1])
                elif not line.startswith('#'):
                    stats.counts[fields[0].split()[0]] = int(fields[1])  # Name without the description
            except (IndexError, ValueError):
                raise SealError('Unexpected line in Seal stats file {}: {}'.format(stats_file, line.strip())) from None
    return stats


def run_seal(inputs, refs, threads=1, memory='1g', k=DEFAULT_KMER_SIZE, hdist=DEFAULT_MISMATCHES):
    """
    Count the reads (or contigs) that contain each reference sequence.

    :param inputs: one fasta/fastq file, or two paired-end fastq files
    :param refs: reference fasta file(s)
    :param memory: Java heap size given to Seal, e.g. "1g"
    :return: SealStats
    """
    refs = [refs] if isinstance(refs, (str, os.PathLike)) else list(refs)
    tmp = safe_temporary_folder()
    try:
        stats_file = Path(tmp) / 'stats.tsv'
        cmd = seal_command(safe_paths(inputs, tmp, 'in', check_extension=True), safe_paths(refs, tmp, 'ref'),
                           stats_file, threads, memory, k=k, hdist=hdist)
        log.debug('Running: %s', shlex.join(cmd))
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0 or not stats_file.exists():
            detail = (proc.stderr or proc.stdout).strip().splitlines() or ['no output']
            raise SealError('Seal failed (exit code {}): {}\nCommand: {}\n{}'.format(
                proc.returncode, failure_reason(detail), shlex.join(cmd), '\n'.join(detail[-10:])))
        log.debug('Seal output:\n%s', proc.stderr.strip())
        return parse_stats(stats_file)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def versions():
    """{"BBTools": version, "Java": version} of the Seal installation, "unknown" when they cannot be read."""
    found = {'BBTools': 'unknown', 'Java': 'unknown'}
    path = executable()
    if path is None:
        return found
    try:  # A small heap: without -Xmx, seal.sh reserves most of the free memory, even for --version
        out = subprocess.run([path, '-Xmx64m', '--version'], capture_output=True, text=True, timeout=60)
        match = re.search(r'(?:BBMap|BBTools) version (\S+)', out.stdout + out.stderr)
        if match:
            found['BBTools'] = match.group(1)
        java = Path(path).parent / 'java'
        java = str(java) if java.exists() else shutil.which('java')
        if java:
            out = subprocess.run([java, '-version'], capture_output=True, text=True, timeout=60)
            first = (out.stderr or out.stdout).strip().splitlines()
            if first:
                found['Java'] = first[0]
    except (OSError, subprocess.SubprocessError):
        pass
    return found


def system_memory_mb():
    """Total memory of the machine in MB, or None when it cannot be read."""
    try:
        return os.sysconf('SC_PHYS_PAGES') * os.sysconf('SC_PAGE_SIZE') // (1024 * 1024)
    except (AttributeError, ValueError, OSError):
        return None


def auto_memory(ref_bases, jobs=1):
    """
    Java heap size for Seal, from the size of the reference sequences: Seal holds every reference k-mer in memory,
    at a few tens of bytes each. Never more than 85% of the machine's memory, shared between the jobs: without
    -Xmx, seal.sh would reserve most of the free memory for each job at once.

    :return: a size for -Xmx, e.g. "1024m"
    """
    needed = 512 + ref_bases * 40 / 1e6  # MB
    megabytes = max(1024, -(-int(needed) // 256) * 256)  # Rounded up to 256 MB
    total = system_memory_mb()
    if total:
        allowed = max(512, int(total * 0.85) // max(1, jobs))
        if megabytes > allowed:
            log.warning('The reference sequences need about %d MB of memory for Seal, but only %d MB is available '
                        '(85%% of the machine, shared between %d job(s)). Seal may run out of memory: use fewer '
                        'jobs, or --memory to set the size yourself.', megabytes, allowed, jobs)
            megabytes = allowed
    return '{}m'.format(megabytes)
