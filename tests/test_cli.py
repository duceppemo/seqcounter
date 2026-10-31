import runpy

import pytest

from seqcounter import cli
from seqcounter.errors import SeqCounterError


def parse(*argv):
    parser = cli.build_parser()
    args = parser.parse_args(argv)
    cli.check_arguments(parser, args)
    return args


def test_defaults():
    args = parse('-ref', 'refs.fasta', '-r1', 'S1.fastq.gz', '-o', 'out')
    assert (args.kmer_size, args.mismatch, args.parallel, args.memory) == (25, 1, 1, None)
    assert args.threads >= 1


def test_memory_size():
    assert parse('-ref', 'r.fasta', '-r1', 'S1.fastq', '-o', 'out', '--memory', '2g').memory == '2g'
    assert parse('-ref', 'r.fasta', '-r1', 'S1.fastq', '-o', 'out', '--memory', '500M').memory == '500M'


def test_the_first_version_flags_still_work():
    """Every flag of the first version of seqcounter, which is how scripts call it."""
    args = parse('-i', 'reads/', '-ref', 'refs.fasta', '-k', '31', '-m', '0', '-o', 'out', '-t', '8', '-p', '4')
    assert (args.input, args.ref, args.kmer_size, args.mismatch, args.output) == ('reads/', 'refs.fasta', 31, 0, 'out')
    assert (args.threads, args.parallel) == (8, 4)
    long_form = parse('--input', 'reads/', '--ref', 'refs.fasta', '--kmer-size', '31', '--mismatch', '0',
                      '--output', 'out', '--threads', '8', '--jobs', '4')
    assert (long_form.threads, long_form.parallel) == (8, 4)


def test_version_is_minus_v(capsys):
    """-v printed the version in the first version of seqcounter: --verbose is the long one."""
    with pytest.raises(SystemExit) as exit_info:
        parse('-v')
    assert exit_info.value.code == 0
    assert capsys.readouterr().out.startswith('seqcounter ')


@pytest.mark.parametrize('argv, message', [
    (['-o', 'out', '-ref', 'r.fasta'], 'one of the arguments'),  # No input
    (['-i', 'reads/', '-r1', 'S1.fastq', '-ref', 'r.fasta', '-o', 'out'], 'not allowed with'),  # Both inputs
    (['-i', 'reads/', '-ref', 'r.fasta', '-o', 'out', '-r2', 'S2.fastq'], '-r2 and --sample are only used'),
    (['-i', 'reads/', '-ref', 'r.fasta', '-o', 'out', '-s', 'name'], '-r2 and --sample are only used'),
    (['-r1', 'S1.fastq', '-ref', 'r.fasta', '-o', 'out', '-p', '2'], '--parallel is only used with -i'),
    (['-r1', 'S1.fastq', '-ref', 'r.fasta', '-o', 'out', '-m', '-1'], 'cannot be negative'),
    (['-r1', 'S1.fastq', '-ref', 'r.fasta', '-o', 'out', '-m', '25'], 'must be smaller than --kmer-size'),
    (['-r1', 'S1.fastq', '-ref', 'r.fasta', '-o', 'out', '-s', 'a/b'], 'must be a name, not a path'),
    (['-r1', 'S1.fastq', '-ref', 'r.fasta', '-o', 'out', '-s', ' '], 'must be a name, not a path'),
    (['-r1', 'S1.fastq', '-ref', 'r.fasta', '-o', 'out', '-s', 'a\tb'], 'cannot contain tabs'),
    (['-r1', 'S1.fastq', '-o', 'out'], 'required'),  # No reference
    (['-r1', 'S1.fastq', '-ref', 'r.fasta'], 'required'),  # No output folder
    (['-r1', 'S1.fastq', '-ref', 'r.fasta', '-o', 'out', '-t', '0'], 'must be >= 1'),
    (['-r1', 'S1.fastq', '-ref', 'r.fasta', '-o', 'out', '-t', 'many'], 'not an integer'),
    (['-r1', 'S1.fastq', '-ref', 'r.fasta', '-o', 'out', '--memory', '2tb'], 'not a Java memory size'),
])
def test_argument_errors(capsys, argv, message):
    with pytest.raises(SystemExit) as exit_info:
        parse(*argv)
    assert exit_info.value.code == 2
    assert message in capsys.readouterr().err


def test_format_elapsed():
    assert cli.format_elapsed(0) == '0s'
    assert cli.format_elapsed(45.4) == '45s'
    assert cli.format_elapsed(3723) == '1h2m3s'
    assert cli.format_elapsed(90061) == '1d1h1m1s'


def test_available_cpus(monkeypatch):
    assert cli.available_cpus() >= 1
    monkeypatch.delattr(cli.os, 'sched_getaffinity', raising=False)
    monkeypatch.setattr(cli.os, 'cpu_count', lambda: None)
    assert cli.available_cpus() == 1


def test_setup_logging_adds_one_handler_only():
    cli.setup_logging()
    cli.setup_logging(verbose=True)
    assert len(cli.log.handlers) == 1 and cli.log.level == 10


def test_main_reports_a_missing_seal(tmp_path, monkeypatch, caplog):
    monkeypatch.setattr(cli, 'check_seal', lambda: (_ for _ in ()).throw(SeqCounterError('no seal.sh here')))
    with pytest.raises(SystemExit) as exit_info:
        cli.main(['-r1', str(tmp_path / 'S1.fastq'), '-ref', str(tmp_path / 'r.fasta'), '-o', str(tmp_path / 'out')])
    assert exit_info.value.code == 1
    assert 'no seal.sh here' in caplog.text


def test_main_reports_an_interruption(tmp_path, monkeypatch, caplog):
    monkeypatch.setattr(cli, 'check_seal', lambda: (_ for _ in ()).throw(KeyboardInterrupt()))
    with pytest.raises(SystemExit) as exit_info:
        cli.main(['-r1', str(tmp_path / 'S1.fastq'), '-ref', str(tmp_path / 'r.fasta'), '-o', str(tmp_path / 'out')])
    assert exit_info.value.code == 130
    assert 'Interrupted' in caplog.text


def test_module_entry_point(monkeypatch, capsys):
    """"python -m seqcounter" runs the same command line."""
    monkeypatch.setattr(cli.sys, 'argv', ['seqcounter', '--version'])
    with pytest.raises(SystemExit) as exit_info:
        runpy.run_module('seqcounter', run_name='__main__')
    assert exit_info.value.code == 0
    assert capsys.readouterr().out.startswith('seqcounter ')


def test_module_is_importable():
    """Importing seqcounter.__main__ must not run anything."""
    import seqcounter.__main__ as entry
    assert entry.main is cli.main
