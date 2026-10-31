"""Read the reference sequences to count, and check that Seal can count them."""

import gzip
import logging
import re
from dataclasses import dataclass
from pathlib import Path

from .errors import SeqCounterError

log = logging.getLogger(__name__)


@dataclass
class Reference:
    name: str  # The name Seal reports: the first word of the fasta header
    description: str  # The rest of the header, or ""
    length: int


def open_text(path):
    """Open a text file, gzipped or not."""
    with open(path, 'rb') as f:
        gzipped = f.read(2) == b'\x1f\x8b'
    return gzip.open(path, 'rt') if gzipped else open(path)


def read_fasta(path):
    """[(header, sequence length)] of a fasta file, gzipped or not, in file order."""
    records, header, length = [], None, 0
    with open_text(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header is not None:
                    records.append((header, length))
                header, length = line[1:].strip(), 0
            elif line and header is not None:
                length += len(re.sub(r'\s', '', line))
    if header is not None:
        records.append((header, length))
    return records


def load_references(path, kmer_size):
    """
    The reference sequences of a fasta file, checked against how Seal counts them.

    Rejected rather than counted wrongly: an empty file, a sequence shorter than the k-mer size (Seal cannot
    find it at all), and two sequences whose first word is the same (their counts would collide in one row).

    :return: list of Reference, in file order
    """
    if not Path(path).is_file():
        raise SeqCounterError('Reference file not found: {}'.format(path))
    try:
        records = read_fasta(path)
    except (OSError, EOFError, UnicodeDecodeError) as e:
        raise SeqCounterError('Reference file {} cannot be read: {}'.format(path, e)) from None
    if not records:
        raise SeqCounterError('No sequence found in the reference file {}. It must be in fasta format, with a '
                              '">name" line before each sequence.'.format(path))

    references, seen, duplicates, short = [], {}, [], []
    for header, length in records:
        name, _, description = header.partition(' ')
        name = re.sub(r'[\t\r\n]', '_', name) or 'unnamed'
        if length == 0:
            raise SeqCounterError('The reference sequence "{}" of {} is empty.'.format(name, path))
        if length < kmer_size:
            short.append('"{}" ({} bp)'.format(name, length))
        if name in seen:
            duplicates.append('"{}"'.format(name))
        seen[name] = True
        references.append(Reference(name, description.strip(), length))

    if duplicates:
        raise SeqCounterError('Several reference sequences of {} have the same name: {}. Seal reports one row per '
                              'name, up to the first space, so their counts would be added together. Rename '
                              'them.'.format(path, ', '.join(sorted(set(duplicates)))))
    if short:
        raise SeqCounterError('{} reference sequence(s) of {} are shorter than the k-mer size ({} bp): {}. Seal '
                              'cannot find them: lower -k to at most the shortest reference, or remove '
                              'them.'.format(len(short), path, kmer_size, ', '.join(short)))
    log.info('%d reference sequence(s) read from %s (%d-%d bp)', len(references), path,
             min(r.length for r in references), max(r.length for r in references))
    return references
