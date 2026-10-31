#!/usr/bin/env python3
"""
Draw the three figures on the wiki's Validation page from the results of run_validation.sh.

    python3 make_figures.py        # needs matplotlib; run run_validation.sh first

It reads results/*/refseqcounts.json, the genome and its feature table, and writes, in assets/, a light and a
dark version of each figure:

    validation_model.svg      counted against expected, on simulated reads (section 3 of the page)
    validation_coverage.svg   counted over expected per sample, simulated and real (section 5)
    validation_chromosome.svg the chromosome, its IS6110 copies and its repeated kilobases (the data section)

Nothing is drawn that is not in the data: the expected counts come from the same uniform-coverage model as
`validate.py report`, the copy numbers and the IS6110 sites are counted in the genome here, the named loci are at
the coordinates NCBI gives for them, and the text that rides the marks is computed, never typed in. The text is
written as outlines (svg.fonttype='path'), so the files render identically without DejaVu Sans installed, and the
background is left transparent so the dark version sits on a dark page.
"""
import pathlib
import re
import statistics
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

from validate import K, read_fasta, samples_of  # noqa: E402  (after sys.path)

GENOME_SIZE = 4_411_532
ELEMENT = 'IS6110'
CONTROL = 'absent_300bp'
COPIES = 16  # Measured, not quoted: the median IS6110 k-mer occurs 16 times in the genome
BIN = 1_000  # The copy-number ring of the chromosome figure, in bases per bin

# The few loci the chromosome figure names, by the gene the feature table calls them. Everything here is drawn at
# the coordinates NCBI gives; nothing is placed by hand. The direct-repeat locus has no gene, so it is taken from
# the repeat_region that spans it.
LOCI = [
    ('dnaA', 'oriC (dnaA)'),
    ('gyrA', 'gyrA'),
    ('rpoB', 'rpoB'),
    ('rrs', 'rrn operon'),
    ('katG', 'katG'),
    ('pncA', 'pncA'),
    ('esxA', 'esx-1 (RD1)'),
]
DR_NOTE = 'direct repeat region'  # The spoligotyping locus, as the feature table describes it

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


def kmer_codes(sequence):
    """The canonical 2-bit code of every k-mer of a sequence, and which of them are free of ambiguous bases."""
    import numpy as np

    table = np.full(256, 255, dtype=np.uint8)
    for value, base in enumerate('ACGT'):
        table[ord(base)] = value
    bases = table[np.frombuffer(sequence.encode(), dtype=np.uint8)]
    clean = bases < 4
    values = np.where(clean, bases, 0).astype(np.uint64)
    n = len(values) - K + 1
    forward = np.zeros(n, dtype=np.uint64)
    reverse = np.zeros(n, dtype=np.uint64)
    for j in range(K):
        forward = (forward << np.uint64(2)) | values[j:j + n]
        reverse = (reverse << np.uint64(2)) | (np.uint64(3) - values[K - 1 - j:K - 1 - j + n])
    usable = np.lib.stride_tricks.sliding_window_view(clean, K).all(axis=1)
    return np.minimum(forward, reverse), usable


def copy_number_profile(genome, element):
    """
    How many times each k-mer of the genome occurs in it, binned, and where the element sits.

    The same question the counts answer, asked of the genome alone: a k-mer that occurs once belongs to a
    single-copy locus, and a k-mer that occurs sixteen times belongs to something there are sixteen of. The
    element's positions are every place a k-mer of the element occurs, grouped into sites.
    """
    import numpy as np

    codes, usable = kmer_codes(genome)
    unique, inverse, counts = np.unique(codes[usable], return_inverse=True, return_counts=True)
    occurrence = np.zeros(len(codes), dtype=np.int32)
    occurrence[usable] = counts[inverse]
    bins = len(occurrence) // BIN
    binned = np.median(occurrence[:bins * BIN].reshape(bins, BIN), axis=1)
    element_codes, element_usable = kmer_codes(element)
    # The element comes from its own GenBank record, which is not identical to this genome's copies of it: some of
    # its k-mers are not in the genome at all. They count as zero, as validate.py counts them, rather than being
    # read off whichever k-mer sorts next.
    wanted = element_codes[element_usable]
    found = np.searchsorted(unique, wanted)
    present = np.take(unique, found, mode='clip') == wanted
    median_copies = int(np.median(np.where(present, counts[np.clip(found, 0, len(counts) - 1)], 0)))
    positions = np.flatnonzero(np.isin(codes, wanted))
    sites = []
    for run in np.split(positions, np.flatnonzero(np.diff(positions) > 2000) + 1):
        span = run[-1] + K - run[0]
        sites.append({'start': int(run[0]), 'end': int(run[-1] + K), 'copies': max(1, round(span / len(element)))})
    return {'binned': binned, 'edges': np.arange(bins) * BIN + BIN / 2, 'sites': sites,
            'median_copies': median_copies, 'total_copies': sum(site['copies'] for site in sites)}


def annotated_loci(features_path):
    """The midpoint of each named locus, from the feature table NCBI serves for the record."""
    wanted = dict(LOCI)
    found, span = {}, None
    for line in features_path.read_text().splitlines():
        fields = line.split('\t')
        if len(fields) >= 3 and fields[0] and fields[2] in ('gene', 'repeat_region'):
            try:
                span = (int(fields[0].lstrip('<>')), int(fields[1].lstrip('<>')), fields[2])
            except ValueError:
                span = None
        elif len(fields) >= 5 and span:
            start, end, kind = span
            middle = (min(start, end) + max(start, end)) / 2
            if kind == 'gene' and fields[3] == 'gene' and fields[4] in wanted:
                found.setdefault(wanted[fields[4]], middle)
            # The locus itself, not one of the 36 bp repeats inside it
            elif (kind == 'repeat_region' and fields[3] == 'note' and DR_NOTE in fields[4]
                    and abs(end - start) > 1000):
                found.setdefault('DR locus', middle)
    return sorted(found.items(), key=lambda item: item[1])


def window_positions(report_path):
    """Where the single-copy windows came from, as their descriptions in the report record it."""
    import json

    with open(report_path) as f:
        references = json.load(f)['references']
    starts = []
    for reference in references:
        match = re.match(r'(\d+)-(\d+) of the genome', reference.get('description', ''))
        if match:
            starts.append((int(match.group(1)) + int(match.group(2))) / 2)
    return starts


def spread(angles, minimum):
    """Nudge label angles apart, keeping their order, so that neighbours do not overprint."""
    out = list(angles)
    for _ in range(200):
        moved = False
        for i in range(len(out) - 1):
            gap = out[i + 1] - out[i]
            if gap < minimum:
                shift = (minimum - gap) / 2
                out[i] -= shift
                out[i + 1] += shift
                moved = True
        if not moved:
            break
    return out


def style(plt, theme):
    plt.rcParams.update({
        'svg.fonttype': 'path', 'font.family': 'DejaVu Sans', 'font.size': 10,
        'svg.hashsalt': 'seqcounter',  # Fixed, so that redrawing unchanged data rewrites the same bytes
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
    fig.savefig(out_path, format='svg', transparent=True, metadata={'Date': None})
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
    fig.savefig(out_path, format='svg', transparent=True, metadata={'Date': None})
    plt.close(fig)


def figure_chromosome(plt, theme, profile, loci, windows, out_path):
    """The chromosome as a circle: where the references sit, and how many copies of each the genome holds."""
    import numpy as np
    from matplotlib.lines import Line2D

    fig = plt.figure(figsize=(7.6, 7.6))
    ax = fig.add_subplot(projection='polar')
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_axis_off()
    ax.set_ylim(0, 1)
    ring = np.linspace(0, 2 * np.pi, 721)
    ceiling = 20  # The copy-number ring is logarithmic, and nothing in this genome reaches 20 per kilobase

    def angle(position):
        return 2 * np.pi * np.asarray(position) / GENOME_SIZE

    def radius(copies):
        return 0.44 + 0.22 * np.log10(copies) / np.log10(ceiling)

    # The chromosome, with its ruler reading inwards so that the locus names have the outside to themselves
    ax.fill_between(ring, 0.90, 0.95, color=theme['grid'], linewidth=0, zorder=2)
    for position in range(0, GENOME_SIZE, 500_000):
        ax.plot([angle(position)] * 2, [0.885, 0.90], color=theme['guide'], linewidth=1, zorder=3)
        ax.text(angle(position), 0.835, '{:.1f}'.format(position / 1e6) if position else '0 Mb',
                ha='center', va='center', fontsize=8, color=theme['muted'])

    # Every place a 25-mer occurs more than once, at the height of how often it occurs
    ax.plot(ring, [radius(1)] * len(ring), color=theme['windows'], linewidth=1.6, zorder=3)
    repeated = np.flatnonzero(profile['binned'] > 1)
    for index in repeated:
        at = angle(profile['edges'][index])
        ax.plot([at] * 2, [radius(1), radius(profile['binned'][index])], color=theme['windows'],
                linewidth=1.8, solid_capstyle='butt', zorder=4)
    for copies in (COPIES,):
        line, = ax.plot(ring, [radius(copies)] * len(ring), color=theme['guide'], linewidth=1, zorder=2)
        line.set_dashes((3, 3))
    for copies, label in ((1, '1 copy'), (COPIES, '{} copies'.format(COPIES))):
        ax.text(angle(260_000), radius(copies), ' ' + label, ha='left', va='center',
                fontsize=8.5, color=theme['muted'], zorder=6)

    # Where the single-copy references were taken from, and IS6110, which crosses the chromosome
    for position in windows:
        ax.plot([angle(position)] * 2, [0.70, 0.745], color=theme['guide'], linewidth=1, zorder=3)
    for site in profile['sites']:
        middle = angle((site['start'] + site['end']) / 2)
        ax.plot([middle] * 2, [0.875, 0.975], color=theme['element'], linewidth=2.4,
                solid_capstyle='butt', zorder=6)
        ax.plot([middle], [0.87], marker='o', markersize=5, color=theme['element'],
                markeredgecolor=theme['surface'], markeredgewidth=1, zorder=6)
        if site['copies'] > 1:
            ax.text(middle, 1.0, '\u00d7{}'.format(site['copies']), ha='center', va='bottom',
                    fontsize=9, color=theme['ink'], zorder=7, clip_on=False)

    # The named loci, outside, pulled apart just enough not to overprint. The spreading starts just after the
    # widest gap between two of them, so that the wrap at the origin is not treated as a crowd.
    places = [float(angle(position)) for _, position in loci]
    turn = 2 * np.pi
    first = int(np.argmax([(places[(i + 1) % len(places)] - place) % turn for i, place in enumerate(places)])) + 1
    ordered = [loci[(first + i) % len(loci)] for i in range(len(loci))]
    rotated = [places[(first + i) % len(places)] for i in range(len(places))]
    angles = spread([rotated[0] + (place - rotated[0]) % turn for place in rotated], np.radians(14))
    for (name, position), label_angle in zip(ordered, angles, strict=True):
        at = float(angle(position))
        ax.plot([at, at, label_angle, label_angle], [0.95, 0.99, 1.03, 1.05],
                color=theme['guide'], linewidth=1, clip_on=False, zorder=4)
        degrees = np.degrees(label_angle) % 360
        flip = 180 < degrees < 360
        ax.text(label_angle, 1.07, name, ha='right' if flip else 'left', va='center', fontsize=9.5,
                color=theme['ink'], rotation=(90 - degrees) + (180 if flip else 0),
                rotation_mode='anchor', clip_on=False, zorder=4)

    ax.text(0, 0.16, 'M. tuberculosis', ha='center', va='center', fontsize=10.5,
            fontstyle='italic', color=theme['ink'])
    ax.text(0, 0.02, 'H37Rv\nNC_000962.3\n{:,} bp'.format(GENOME_SIZE), ha='center', va='center',
            fontsize=9.5, color=theme['muted'], linespacing=1.6)
    fig.suptitle('Sixteen copies of IS6110, and 75 references with exactly one', x=0.5, y=0.985,
                 fontsize=13, fontweight='bold', color=theme['ink'])
    fig.text(0.5, 0.945, 'the chromosome the validation counts against, and the few places it repeats',
             ha='center', fontsize=10, color=theme['muted'])
    handles = [
        Line2D([], [], color=theme['element'], linewidth=2.4,
               label='IS6110, {} copies at {} sites'.format(profile['total_copies'], len(profile['sites']))),
        Line2D([], [], color=theme['windows'], linewidth=1.8,
               label='25-mer copies per kilobase'),
        Line2D([], [], color=theme['guide'], linewidth=1,
               label='{} single-copy references'.format(len(windows))),
    ]
    legend = fig.legend(handles=handles, loc='upper left', ncol=1, frameon=False, fontsize=9.5,
                        handletextpad=0.7, labelspacing=0.7, bbox_to_anchor=(0.015, 0.90))
    for text in legend.get_texts():
        text.set_color(theme['ink'])
    fig.subplots_adjust(left=0.12, right=0.88, top=0.79, bottom=0.13)
    fig.savefig(out_path, format='svg', transparent=True, metadata={'Date': None})
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
    genome = next(iter(read_fasta(here / 'data' / 'H37Rv.fasta').values()))
    element = next(iter(read_fasta(here / 'data' / 'IS6110.fasta').values()))
    profile = copy_number_profile(genome, element)
    loci = annotated_loci(here / 'data' / 'H37Rv_features.txt')
    windows = window_positions(here / 'results' / 'sim_m0' / 'refseqcounts.json')
    print('{} copies of {} at {} sites, median k-mer occurrence {}; {} loci named'.format(
        profile['total_copies'], ELEMENT, len(profile['sites']), profile['median_copies'], len(loci)))
    for theme, suffix in ((LIGHT, ''), (DARK, '_dark')):
        style(plt, theme)
        figure_model(plt, theme, simulated, assets / 'validation_model{}.svg'.format(suffix))
        figure_coverage(plt, theme, [simulated] + reads, assets / 'validation_coverage{}.svg'.format(suffix))
        figure_chromosome(plt, theme, profile, loci, windows,
                          assets / 'validation_chromosome{}.svg'.format(suffix))
    print('wrote {}/validation_[model|coverage|chromosome][_dark].svg'.format(assets))
    return 0


if __name__ == '__main__':
    sys.exit(main())
