import gzip
import os
import subprocess

import pytest

from seqcounter import seal
from seqcounter.seal import SealError


@pytest.fixture(autouse=True)
def clear_cache():
    """executable() is cached: the tests change the PATH."""
    seal.executable.cache_clear()
    yield
    seal.executable.cache_clear()


def test_executable_in_path(tmp_path, monkeypatch):
    fake = tmp_path / 'seal.sh'
    fake.write_text('#!/bin/sh\n')
    fake.chmod(0o755)
    monkeypatch.setenv('PATH', str(tmp_path))
    assert seal.executable() == str(fake)
    assert seal.check_seal() == str(fake)


def test_executable_next_to_python(tmp_path, monkeypatch):
    """Found in the folder of the running Python: "env/bin/seqcounter" works without activating the environment."""
    bin_folder = tmp_path / 'bin'
    bin_folder.mkdir()
    fake = bin_folder / 'seal.sh'
    fake.write_text('#!/bin/sh\n')
    fake.chmod(0o755)
    monkeypatch.setenv('PATH', str(tmp_path / 'empty'))
    monkeypatch.setattr(seal.sys, 'executable', str(bin_folder / 'python'))
    assert seal.executable() == str(fake)


def test_executable_not_found(tmp_path, monkeypatch):
    monkeypatch.setenv('PATH', str(tmp_path / 'empty'))
    monkeypatch.setattr(seal.sys, 'executable', str(tmp_path / 'bin' / 'python'))
    assert seal.executable() is None
    with pytest.raises(SealError, match='was not found in your PATH'):
        seal.check_seal()


def test_not_executable(tmp_path, monkeypatch):
    """A seal.sh without the execute bit is not usable."""
    (tmp_path / 'bin').mkdir()
    (tmp_path / 'bin' / 'seal.sh').write_text('#!/bin/sh\n')
    monkeypatch.setenv('PATH', str(tmp_path / 'empty'))
    monkeypatch.setattr(seal.sys, 'executable', str(tmp_path / 'bin' / 'python'))
    assert seal.executable() is None


def test_seal_command(tmp_path, monkeypatch):
    monkeypatch.setattr(seal, 'check_seal', lambda: 'seal.sh')
    single = seal.seal_command(['in.fastq.gz'], 'refs.fasta', tmp_path / 'stats.tsv', 4, '2g')
    assert 'in=in.fastq.gz' in single and 'int=f' in single and 'in2=' not in ' '.join(single)
    assert 'ref=refs.fasta' in single and 'k=25' in single and 'hdist=1' in single and 'threads=4' in single
    assert single[1] == '-Xmx2g'
    paired = seal.seal_command(['r1.fq', 'r2.fq'], ['a.fasta', 'b.fasta'], tmp_path / 'stats.tsv', 1, '1g',
                               k=31, hdist=0)
    assert 'in2=r2.fq' in paired and 'int=f' not in paired
    assert 'ref=a.fasta,b.fasta' in paired and 'k=31' in paired and 'hdist=0' in paired


def test_extensions():
    assert seal.extensions('S1_R1.fastq.gz') == '.fastq.gz'
    assert seal.extensions('S1.fasta') == '.fasta'
    assert seal.extensions('/path/to/dataset_1.dat') == '.dat'
    assert seal.extensions('weird name.fastq') == '.fastq'  # The space is handled by safe_paths
    assert seal.extensions('S1.xmx.fastq') == '.fastq'  # A suffix Seal would read as a memory setting


def test_content_extension(tmp_path):
    (tmp_path / 'a').write_text('@r1\nACGT\n+\nIIII\n')
    (tmp_path / 'b').write_text('>c1\nACGT\n')
    with gzip.open(tmp_path / 'c', 'wt') as f:
        f.write('@r1\nACGT\n+\nIIII\n')
    (tmp_path / 'd').write_bytes(b'\x1f\x8bbroken')
    (tmp_path / 'e').write_text('\n\n@r1\nACGT\n+\nIIII\n')  # Blank lines first
    assert seal.content_extension(tmp_path / 'a') == '.fastq'
    assert seal.content_extension(tmp_path / 'b') == '.fasta'
    assert seal.content_extension(tmp_path / 'c') == '.fastq.gz'
    assert seal.content_extension(tmp_path / 'd') == '.fasta.gz'  # Unreadable: let Seal report the error
    assert seal.content_extension(tmp_path / 'e') == '.fastq'


def test_safe_paths_leaves_plain_paths(tmp_path):
    (tmp_path / 'S1.fasta').write_text('>c\nACGT\n')
    paths = seal.safe_paths([tmp_path / 'S1.fasta'], tmp_path, 'ref')
    assert paths == [str(tmp_path / 'S1.fasta')]


@pytest.mark.parametrize('name', ['my sample.fastq', "quote'.fastq", 'xmx.fastq', 'a,b.fastq', 'd$t.fastq'])
def test_safe_paths_links_unsafe_names(tmp_path, name):
    """A name Seal cannot be given safely is linked under a neutral one."""
    source = tmp_path / name
    source.write_text('@r\nACGT\n+\nIIII\n')
    folder = tmp_path / 'tmp'
    folder.mkdir()
    linked = seal.safe_paths([source], folder, 'in')[0]
    assert linked == str(folder / 'in0.fastq')
    assert os.path.realpath(linked) == str(source.resolve())


def test_safe_paths_fixes_wrong_extension(tmp_path):
    """A fastq file named .fasta (or without a known extension) is linked under the extension of its content."""
    (tmp_path / 'reads.fasta').write_text('@r\nACGT\n+\nIIII\n')
    (tmp_path / 'reads.dat').write_text('@r\nACGT\n+\nIIII\n')
    with gzip.open(tmp_path / 'reads.fastq', 'wt') as f:  # Gzipped, but not named .gz
        f.write('@r\nACGT\n+\nIIII\n')
    folder = tmp_path / 'tmp'
    folder.mkdir()
    linked = seal.safe_paths([tmp_path / 'reads.fasta', tmp_path / 'reads.dat', tmp_path / 'reads.fastq'],
                             folder, 'in', check_extension=True)
    assert [os.path.basename(p) for p in linked] == ['in0.fastq', 'in1.fastq', 'in2.fastq.gz']


def test_safe_paths_keeps_right_extension(tmp_path):
    (tmp_path / 'reads.fastq').write_text('@r\nACGT\n+\nIIII\n')
    assert seal.safe_paths([tmp_path / 'reads.fastq'], tmp_path, 'in', check_extension=True) == \
        [str(tmp_path / 'reads.fastq')]


def test_safe_temporary_folder(monkeypatch):
    folder = seal.safe_temporary_folder()
    assert os.path.isdir(folder) and not seal.UNSAFE.search(folder)
    os.rmdir(folder)

    made = []

    def unsafe_mkdtemp(prefix):
        folder = '/tmp/{}xmx{}'.format(prefix, len(made))
        os.mkdir(folder)
        made.append(folder)
        return folder

    monkeypatch.setattr(seal.tempfile, 'mkdtemp', unsafe_mkdtemp)
    with pytest.raises(SealError, match='Set TMPDIR'):
        seal.safe_temporary_folder(attempts=3)
    assert len(made) == 3 and not any(os.path.exists(f) for f in made)  # Every attempt was removed


def test_failure_reason():
    assert seal.failure_reason(['Loading', 'Exception in thread "main" java.lang.RuntimeException: Nope',
                                '\tat jgi.Seal.main(Seal.java:1)']) == \
        'Exception in thread "main" java.lang.RuntimeException: Nope'
    # A specific error is preferred to Seal's generic final message
    assert seal.failure_reason(['java.io.FileNotFoundException: refs.fasta', 'Exception in thread "main" ...',
                                'This program terminated in an error state.']) == \
        'java.io.FileNotFoundException: refs.fasta'
    # A message on the next line is added
    assert seal.failure_reason(['java.lang.RuntimeException:', '\tat jgi.Seal.main', 'No reference sequence']) == \
        'java.lang.RuntimeException: No reference sequence'
    assert seal.failure_reason(['java.lang.RuntimeException:']) == 'java.lang.RuntimeException:'
    assert seal.failure_reason(['Loading', 'nothing useful']) == 'nothing useful'  # No error line at all


def test_parse_stats(tmp_path):
    stats = tmp_path / 'stats.tsv'
    stats.write_text('#File\tS1.fastq.gz\n#Total\t822714\t123407100\n#Matched\t799\t0.09712%\n'
                     '#Name\tReads\tReadsPct\nref1 some description\t62\t0.00754%\nref2\t0\t0.00000%\n\n')
    parsed = seal.parse_stats(stats)
    assert parsed.counts == {'ref1': 62, 'ref2': 0}
    assert (parsed.reads, parsed.bases, parsed.matched) == (822714, 123407100, 799)

    stats.write_text('#Total\t10\n#Name\tReads\nref1\t1\n')  # Older Seal: no bases column
    assert seal.parse_stats(stats).bases is None


def test_parse_stats_unexpected_line(tmp_path):
    stats = tmp_path / 'stats.tsv'
    stats.write_text('#Name\tReads\nref1\tnot a number\n')
    with pytest.raises(SealError, match='Unexpected line'):
        seal.parse_stats(stats)


def test_run_seal_reports_failure(tmp_path, monkeypatch):
    """A Seal that fails gives its command and the most informative line of its output."""
    (tmp_path / 'in.fastq').write_text('@r\nACGT\n+\nIIII\n')
    monkeypatch.setattr(seal, 'check_seal', lambda: '/bin/false')
    monkeypatch.setattr(seal.subprocess, 'run', lambda *a, **k: subprocess.CompletedProcess(
        a[0], 1, '', 'java.io.FileNotFoundException: refs.fasta\nThis program terminated in an error state.'))
    with pytest.raises(SealError, match='FileNotFoundException'):
        seal.run_seal([tmp_path / 'in.fastq'], 'refs.fasta')


def test_run_seal_reports_silent_failure(tmp_path, monkeypatch):
    (tmp_path / 'in.fastq').write_text('@r\nACGT\n+\nIIII\n')
    monkeypatch.setattr(seal, 'check_seal', lambda: '/bin/false')
    monkeypatch.setattr(seal.subprocess, 'run',
                        lambda *a, **k: subprocess.CompletedProcess(a[0], 0, '', ''))  # No stats file written
    with pytest.raises(SealError, match='no output'):
        seal.run_seal([tmp_path / 'in.fastq'], 'refs.fasta')


def test_auto_memory(monkeypatch):
    monkeypatch.setattr(seal, 'system_memory_mb', lambda: 64000)
    assert seal.auto_memory(1000) == '1024m'  # Small references: the floor
    assert seal.auto_memory(100_000_000) == '4608m'  # 100 Mb of reference
    monkeypatch.setattr(seal, 'system_memory_mb', lambda: None)  # Unknown system memory: no cap
    assert seal.auto_memory(100_000_000) == '4608m'


def test_auto_memory_capped(monkeypatch, caplog):
    monkeypatch.setattr(seal, 'system_memory_mb', lambda: 2000)
    assert seal.auto_memory(100_000_000, jobs=2) == '850m'  # 85% of 2000, shared between 2 jobs
    assert 'may run out of memory' in caplog.text


def test_system_memory_mb():
    memory = seal.system_memory_mb()
    assert memory is None or memory > 0


def test_system_memory_unknown(monkeypatch):
    monkeypatch.setattr(seal.os, 'sysconf', lambda name: (_ for _ in ()).throw(ValueError()))
    assert seal.system_memory_mb() is None


def test_versions_without_seal(tmp_path, monkeypatch):
    monkeypatch.setenv('PATH', str(tmp_path / 'empty'))
    monkeypatch.setattr(seal.sys, 'executable', str(tmp_path / 'bin' / 'python'))
    assert seal.versions() == {'BBTools': 'unknown', 'Java': 'unknown'}


def test_versions_handles_a_broken_seal(tmp_path, monkeypatch):
    fake = tmp_path / 'seal.sh'
    fake.write_text('#!/bin/sh\nexit 1\n')
    fake.chmod(0o755)
    monkeypatch.setenv('PATH', str(tmp_path))
    monkeypatch.setattr(seal.shutil, 'which', lambda name: str(fake) if name == 'seal.sh' else None)
    assert seal.versions()['BBTools'] == 'unknown'


def test_versions(seal_or_skip):
    found = seal.versions()
    assert found['BBTools'] != 'unknown' and found['Java'] != 'unknown'


def test_versions_when_the_commands_say_nothing(tmp_path, monkeypatch):
    """A Seal and a Java that print nothing leave the versions unknown rather than failing."""
    fake = tmp_path / 'seal.sh'
    fake.write_text('#!/bin/sh\n')
    fake.chmod(0o755)
    (tmp_path / 'java').write_text('#!/bin/sh\n')
    (tmp_path / 'java').chmod(0o755)
    monkeypatch.setenv('PATH', str(tmp_path))
    monkeypatch.setattr(seal.subprocess, 'run',
                        lambda *a, **k: subprocess.CompletedProcess(a[0], 0, '', ''))
    assert seal.versions() == {'BBTools': 'unknown', 'Java': 'unknown'}


def test_versions_when_the_commands_cannot_run(tmp_path, monkeypatch):
    fake = tmp_path / 'seal.sh'
    fake.write_text('#!/bin/sh\n')
    fake.chmod(0o755)
    monkeypatch.setenv('PATH', str(tmp_path))
    monkeypatch.setattr(seal.subprocess, 'run', lambda *a, **k: (_ for _ in ()).throw(OSError('no fork')))
    assert seal.versions() == {'BBTools': 'unknown', 'Java': 'unknown'}
