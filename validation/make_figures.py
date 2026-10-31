#!/usr/bin/env python3
"""
Draw the two figures on the wiki's Validation page from the results of run_validation.sh.

    python3 make_figures.py        # needs matplotlib; run run_validation.sh first

It reads results/*/refseqcounts.json and writes, in assets/, a light and a dark version of each figure:

    validation_model.svg      counted against expected, on simulated reads (section 3 of the page)
    validation_coverage.svg   counted over expected per sample, simulated and real (section 5)

Nothing is drawn that is not in the reports: the expected counts come from the same uniform-coverage model as
`validate.py report`, and the text that rides the marks is computed from the counts, never typed in. The text is
written as outlines (svg.fonttype='path'), so the files render identically without DejaVu Sans installed, and the
background is left transparent so the dark version sits on a dark page.
"""
import pathlib
import statistics
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

from validate import K, read_fasta, samples_of  # noqa: E402  (after sys.path)

GENOME_SIZE = 4_411_532
ELEMENT = 'IS6110'
CONTROL = 'absent_300bp'
COPIES = 16  # Measured, not quoted: the median IS6110 k-mer occurs 16 times in the genome

# Two series only, so two hues: the windows take the project's teal, snapped to the step that clears the
# chroma floor and the lightness band of each mode, and IS6110 the warm hue furthest from it.
LIGHT = {
    'windows': '#0086a5', 'element': '#eb6834', 'surface': '#ffffff',
    'ink': '#243b45', 'muted': '#5d6a72', 'guide': '#aab4bb', 'grid': '#e4e8ea',
}
DARK = {
    'windows': '#009cc0', 'element': '#d95926', 'surface': '#0d1117',
    'ink': '#eef4f6', 'muted': '#a9b6bd', 'guide': '#5d6a72', 'grid': '#232a31',
}


def expected(length, read_length, reads):
    """The reads a single copy of a reference of this length is expected to get, under even coverage."""
    return reads * (length + read_length - 2 * K + 1) / GENOME_SIZE


def measurements(report_path, lengths):
    """Per sample: the window ratios, IS6110's count and implied copy number, and the control."""
    out = []
    for sample in samples_of(report_path):
        counts, reads = sample['counts'], sample['reads']
        read_length = sample['bases'] / reads
        windows = sorted(name for name in counts if name.startswith('w'))
        ratios = {w: counts[w] / expected(lengths[w], read_length, reads) for w in windows}
        median = statistics.median(ratios.values())
        per_copy = expected(lengths[ELEMENT], read_length, reads)
        out.append({
            'sample': sample['sample'],
            'reads': reads,
            'read_length': read_length,
            'depth': sample['bases'] / GENOME_SIZE,
            'expected': {w: expected(lengths[w], read_length, reads) for w in windows},
            'counts': counts,
            'ratios': ratios,
            'median': median,
            'quartiles': statistics.quantiles(ratios.values(), n=4),
            'element_expected': per_copy,
            'element_copies': counts[ELEMENT] / per_copy / median,
            'control': counts[CONTROL],
        })
    return out


def style(plt, theme):
    plt.rcParams.update({
        'svg.fonttype': 'path', 'font.family': 'DejaVu Sans', 'font.size': 10,
        'figure.facecolor': 'none', 'axes.facecolor': 'none',
        'text.color': theme['ink'], 'axes.labelcolor': theme['muted'],
        'xtick.color': theme['muted'], 'ytick.color': theme['muted'],
        'xtick.labelsize': 9, 'ytick.labelsize': 9, 'axes.labelsize': 9.5,
        'xtick.major.size': 0, 'ytick.major.size': 0, 'xtick.minor.size': 0, 'ytick.minor.size': 0,
    })


def recede(ax, theme, spines=('bottom',), axis='both'):
    for side, spine in ax.spines.items():
        spine.set_visible(side in spines)
        spine.set_color(theme['grid'])
        spine.set_linewidth(1)
    ax.grid(True, axis=axis, color=theme['grid'], linewidth=1, linestyle='-')
    ax.set_axisbelow(True)


def figure_model(plt, theme, sim, out_path):
    """Counted against expected, log-log, on the simulated reads: the windows on the 1-copy line, IS6110 on 16."""
    fig, ax = plt.subplots(figsize=(7.4, 4.3))
    recede(ax, theme)
    xs = [sim['expected'][w] for w in sim['ratios']]
    ys = [sim['counts'][w] for w in sim['ratios']]
    limits, top = (25, 600), 11000
    for copies, dash in ((1, None), (COPIES, (4, 3))):
        line, = ax.plot(limits, [copies * g for g in limits], color=theme['guide'], linewidth=1.2, zorder=1)
        if dash:
            line.set_dashes(dash)
        at = min(limits[1], 0.55 * top / copies)  # Label the line where it is still inside the view
        ax.annotate('{} cop{}'.format(copies, 'y' if copies == 1 else 'ies'),
                    (at, copies * at), textcoords='offset points', xytext=(0, 6),
                    ha='right', fontsize=9, color=theme['muted'])
    ax.scatter(xs, ys, s=34, color=theme['windows'], edgecolor=theme['surface'], linewidth=1,
               alpha=0.85, zorder=3, label='{} single-copy windows, 60-1,000 bp'.format(len(xs)))
    ax.scatter([sim['element_expected']], [sim['counts'][ELEMENT]], s=95, marker='D',
               color=theme['element'], edgecolor=theme['surface'], linewidth=1.6, zorder=4,
               label='IS6110, 16 copies in the genome')
    ax.annotate('{:,} reads = {:.1f} copies'.format(sim['counts'][ELEMENT], sim['element_copies']),
                (sim['element_expected'], sim['counts'][ELEMENT]), textcoords='offset points',
                xytext=(-13, 1), fontsize=9.5, color=theme['ink'], ha='right', va='center')
    ax.set(xscale='log', yscale='log', xlim=limits, ylim=(25, top))
    ax.set_xticks([30, 50, 100, 200, 500])
    ax.set_yticks([30, 100, 300, 1000, 3000])
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_formatter(lambda value, _: '{:,.0f}'.format(value))
        axis.set_minor_formatter(lambda value, _: '')
    ax.set_xlabel('reads expected for one copy   N \u00d7 (L + R \u2212 2k + 1) / G')
    ax.set_ylabel('reads counted')
    ax.set_title('Every reference lands on the copy number it has',
                 loc='left', fontsize=12, fontweight='bold', color=theme['ink'], pad=30)
    ax.annotate('{:,} simulated reads of {:.0f} bp, {:.0f}\u00d7, counted with -m 0'.format(
        sim['reads'], sim['read_length'], sim['depth']),
        (0, 1), xycoords='axes fraction', xytext=(0, 9), textcoords='offset points',
        fontsize=9.5, color=theme['muted'])
    ax.annotate('the negative control is absent from the genome:\n0 expected, 0 counted, so off the log axes',
                (0.985, 0.03), xycoords='axes fraction', ha='right', fontsize=9, color=theme['muted'])
    legend = ax.legend(loc='upper left', frameon=False, fontsize=9.5, handletextpad=0.4,
                       borderpad=0, labelspacing=0.5)
    for text in legend.get_texts():
        text.set_color(theme['ink'])
    fig.tight_layout()
    fig.savefig(out_path, format='svg', transparent=True)
    plt.close(fig)


def figure_coverage(plt, theme, samples, out_path):
    """One strip per sample: where the single-copy windows fall against the model, and what IS6110 implies."""
    import numpy as np

    fig, ax = plt.subplots(figsize=(7.4, 4.1))
    recede(ax, theme, spines=(), axis='x')
    rng = np.random.default_rng(7)
    rows = list(enumerate(samples))
    for index, sample in rows:
        y = len(samples) - index
        values = list(sample['ratios'].values())
        jitter = rng.uniform(-0.17, 0.17, len(values))
        ax.scatter(values, [y + j for j in jitter], s=26, color=theme['windows'], alpha=0.8,
                   edgecolor=theme['surface'], linewidth=0.8, zorder=3,
                   label='one single-copy window' if index == 0 else None)
        low, high = sample['quartiles'][0], sample['quartiles'][2]
        ax.plot([low, high], [y - 0.34, y - 0.34], color=theme['muted'], linewidth=2,
                solid_capstyle='butt', zorder=4)
        ax.plot([sample['median']] * 2, [y - 0.46, y - 0.22], color=theme['ink'], linewidth=2,
                solid_capstyle='butt', zorder=5,
                label='median and quartiles' if index == 0 else None)
        ax.annotate('{:.2f}'.format(sample['median']), (sample['median'], y - 0.52),
                    ha='center', va='top', fontsize=9, color=theme['ink'])
        ax.annotate('{:.1f}'.format(sample['element_copies']), (1.04, y), xycoords=('axes fraction', 'data'),
                    va='center', fontsize=10, color=theme['ink'], fontweight='bold')
    ax.axvline(1, color=theme['guide'], linewidth=1.2, zorder=2)
    ax.annotate('the model', (1, len(samples) + 0.55), ha='center', fontsize=9, color=theme['muted'])
    ax.annotate('IS6110\ncopies implied', (1.04, len(samples) + 0.5), xycoords=('axes fraction', 'data'),
                va='center', fontsize=9, color=theme['muted'])
    ax.set_yticks([len(samples) - index for index, _ in rows])
    ax.set_yticklabels(['{}\n{:.0f}\u00d7'.format(sample['sample'], sample['depth']) for _, sample in rows],
                       fontsize=9.5, color=theme['ink'])
    ax.tick_params(axis='y', pad=8)
    ax.set_ylim(0.3, len(samples) + 0.9)
    ax.set_xlim(0, 2.6)
    ax.set_xticks([0, 0.5, 1, 1.5, 2])
    ax.set_xlabel('reads counted / reads expected')
    ax.set_title('Simulated coverage lands on the model; real libraries spread',
                 loc='left', fontsize=12, fontweight='bold', color=theme['ink'], pad=30)
    ax.annotate('each single-copy window counted, divided by what even coverage predicts',
                (0, 1), xycoords='axes fraction', xytext=(0, 9), textcoords='offset points',
                fontsize=9.5, color=theme['muted'])
    legend = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.17), ncol=2, frameon=False,
                       fontsize=9.5, handletextpad=0.4, borderpad=0, columnspacing=1.6)
    for text in legend.get_texts():
        text.set_color(theme['ink'])
    fig.tight_layout()
    fig.savefig(out_path, format='svg', transparent=True)
    plt.close(fig)


def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    here = pathlib.Path(__file__).resolve().parent
    assets = here.parent / 'assets'
    lengths = {name: len(seq) for name, seq in read_fasta(here / 'data' / 'refs.fasta').items()}
    simulated = measurements(here / 'results' / 'sim_m0' / 'refseqcounts.json', lengths)[0]
    reads = measurements(here / 'results' / 'reads' / 'refseqcounts.json', lengths)
    simulated['sample'] = 'simulated'
    for theme, suffix in ((LIGHT, ''), (DARK, '_dark')):
        style(plt, theme)
        figure_model(plt, theme, simulated, assets / 'validation_model{}.svg'.format(suffix))
        figure_coverage(plt, theme, [simulated] + reads, assets / 'validation_coverage{}.svg'.format(suffix))
    print('wrote {}/validation_model[_dark].svg and validation_coverage[_dark].svg'.format(assets))
    return 0


if __name__ == '__main__':
    sys.exit(main())
