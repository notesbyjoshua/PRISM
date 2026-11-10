import hashlib
import json
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, FancyBboxPatch, Polygon
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
TABLES = ROOT / 'publication/final_research_120'
SYNTHETIC = ROOT / 'publication/synthetic_consistency'
FIGURES = ROOT / 'manuscript/figures'
COLORS = ['#27647B', '#D1843E']
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.titlesize': 10,
                     'axes.labelsize': 9, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
                     'svg.fonttype': 'none', 'svg.hashsalt': 'PRISM120_FINAL_RESEARCH_20261109',
                     'axes.spines.top': False, 'axes.spines.right': False})
INVENTORY = []


def short(syndrome):
    name = syndrome.split('__')[0].split(';')[0].replace(' syndrome', '').replace(' SYNDROME', '').title()
    return name.replace('Kbg', 'KBG').replace('Cornelia De Lange', 'Cornelia de Lange')


def save(fig, name, sources, caption, alt):
    FIGURES.mkdir(parents=True, exist_ok=True)
    for ext in ['svg', 'png']:
        fig.savefig(FIGURES / f'{name}.{ext}', dpi=300, bbox_inches='tight',
                    metadata={'Date': None} if ext == 'svg' else {'Software': 'PRISM aggregate plotting script'})
        if ext == 'svg':
            path = FIGURES / f'{name}.{ext}'
            path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines()) + '\n')
    plt.close(fig)
    INVENTORY.append({'name': name, 'sources': {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in sources},
                      'caption': caption, 'alt_text': alt,
                      'outputs': {ext: hashlib.sha256((FIGURES / f'{name}.{ext}').read_bytes()).hexdigest() for ext in ['svg', 'png']}})


def workflow():
    fig, ax = plt.subplots(figsize=(7.1, 7.1))
    ax.set(xlim=(0, 10), ylim=(0, 10))
    ax.axis('off')
    def box(x, y, w, h, text, color):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=.12', facecolor=color, edgecolor='#66727D', linewidth=.7))
        ax.text(x + w / 2, y + h / 2, text, ha='center', va='center', fontsize=9)
    def arrow(a, b):
        ax.annotate('', xy=b, xytext=a, arrowprops={'arrowstyle': '-|>', 'color': '#5F6B76'})
    ax.text(.2, 9.75, 'A  Inputs and FaceKit measurement platform', weight='bold', fontsize=11)
    box(.4, 8.5, 9.2, .7, 'GMDB photographs     •     FairFace margin 0.25     •     PDIDB synthetic images', '#F3F5F6')
    arrow((5, 8.35), (5, 7.95))
    box(.4, 6.3, 9.2, 1.5, 'FaceKit (Chen et al., 2026; pinned source)\n478 landmarks + facial transformation matrix\nCanonical-depth pose correction\n120 feature-specific normalized measurements', '#E0EEF1')
    arrow((9.85, 6.2), (9.85, 5.5))
    ax.text(.2, 5.85, 'B  PRISM cohort design and statistical investigation', weight='bold', fontsize=11)
    box(.4, 4.3, 9.2, 1.1, 'QC: one unambiguous face; 120 finite outputs; yaw/pitch/roll ≤ 15°/15°/10°\nGMDB: age 3 ≤ age < 10; select one image per recorded patient\nFairFace: source age label 3–9; expanded balanced reference', '#EDF2E5')
    arrow((2.6, 4.15), (2.6, 3.65))
    arrow((7.4, 4.15), (7.4, 3.65))
    box(.4, 1.45, 4.4, 2.05, 'Real-patient study\n1,845 selected GMDB patients\n3,409 FairFace reference paths\n48 × 120 = 5,760 primary tests\nDemographics; within-GMDB\nCovariates; photo/reference checks', '#EDF2E5')
    box(5.2, 1.45, 4.4, 2.05, 'Separate conditional fidelity\n1,563 synthetic images pass QC\n9 matched labels; 1,424 images\nAt least 10 rows per source\nDistribution and profile summaries\nModel tags; heavy-tail checks', '#FBEDD8')
    ax.text(5, .85, 'Synthetic rows contribute zero observations to real-patient inference.', ha='center', weight='bold', fontsize=9)
    ax.text(5, .35, 'Original schematic. No real or generated facial photographs are displayed.', ha='center', color='#53606B', fontsize=8)
    save(fig, 'fig5_facekit_prism_workflow', [ROOT / 'ACTIVE_ANALYSIS.json', SYNTHETIC / 'cohort_and_mapping_audit.csv'],
         'FaceKit measurement platform and PRISM study design. FaceKit provides landmark detection, pose correction and 120 measurements; PRISM supplies cohort selection and the statistical branches. Synthetic fidelity is separate and conditional on a GMDB-trained generator (mentor statement; checkpoint-specific provenance unverified). FaceKit: arXiv:2610.09288v1, §§2.1–2.2. This original schematic contains no patient images.',
         'Flow diagram separates the FaceKit measurement steps from PRISM real-patient testing and the synthetic-only fidelity branch.')


def face(ax):
    ax.add_patch(Ellipse((0, 0), 1.7, 2.3, facecolor='#F7F3ED', edgecolor='#7E8B92', lw=1))
    for x in [-.35, .35]:
        ax.add_patch(Ellipse((x, .27), .42, .13, facecolor='white', edgecolor='#7E8B92', lw=.8))
        ax.add_patch(Ellipse((x, .27), .055, .055, color='#7E8B92'))
        ax.plot([x - .2, x, x + .2], [.55, .6, .55], color='#7E8B92', lw=1.5)
    ax.plot([-.08, 0, .08], [.08, -.18, .08], color='#7E8B92', lw=.8)
    ax.plot([-.25, 0, .25], [-.5, -.46, -.5], color='#7E8B92', lw=1)
    ax.set(xlim=(-1.1, 1.1), ylim=(-1.2, 1.25), aspect='equal')
    ax.axis('off')


def measure(ax, a, b, color, label=None):
    ax.annotate('', xy=b, xytext=a, arrowprops={'arrowstyle': '<->', 'color': color, 'lw': 1.2,
                                               'shrinkA': 0, 'shrinkB': 0, 'mutation_scale': 6})
    if label:
        ax.text((a[0] + b[0]) / 2 + .04, (a[1] + b[1]) / 2, label, color=color, fontsize=9)


def anatomy():
    fig, axes = plt.subplots(2, 3, figsize=(7.1, 6.8))
    titles = ['A  Philtrum / midface height', 'B  Interpupillary /\nbizygomatic width', 'C  Outer-canthal /\nbizygomatic width',
              'D  Eyebrow polygon area', 'E  Malar geometric index']
    formulas = ['distance(2, 0) / midface height\nmidface height = |y168 − y2|',
                'distance(468, 473)\n/ bizygomatic width',
                'distance(33, 263)\n/ bizygomatic width',
                'mean(star-polygon area of each brow)\n/ bizygomatic width²',
                'Mean signed horizontal cheek-to-oval\ninterpolation distance / bizygomatic width']
    for ax, title, formula in zip(axes.flat, titles, formulas):
        face(ax)
        ax.set_title(title, fontsize=9, loc='left', pad=8)
        ax.text(.5, -.08, formula, transform=ax.transAxes, ha='center', va='top', fontsize=7, wrap=True)
    a, b, c, d, e, unused = axes.flat
    measure(a, (.12, -.18), (.12, -.46), COLORS[0], 'N')
    measure(a, (-.12, .45), (-.12, -.18), COLORS[1], 'D')
    for ax in [b, c, d, e]:
        measure(ax, (-.8, -.12), (.8, -.12), COLORS[1])
    measure(b, (-.35, .27), (.35, .27), COLORS[0])
    measure(c, (-.56, .27), (.56, .27), COLORS[0])
    for x in [-.35, .35]:
        d.add_patch(Polygon([(x-.21,.54),(x,.63),(x+.21,.54),(x+.15,.5),(x,.56),(x-.15,.5)],
                            facecolor=COLORS[0], alpha=.65, edgecolor=COLORS[0]))
    e.plot([.77, .68], [-.1, -.6], '--', color='#7E8B92')
    e.scatter([.47], [-.36], color=COLORS[0], s=12)
    measure(e, (.47, -.36), (.723, -.36), COLORS[0])
    unused.axis('off')
    unused.text(.02, .92, 'Schematic only', transform=unused.transAxes, fontsize=10, weight='bold')
    unused.text(.02, .8, 'Blue: numerator geometry\nOrange: normalization length\n\nN = numerator; D = denominator\n\nNot a patient, detector output,\nor calibrated anthropometry.\nSee exact definitions in S6.',
                transform=unused.transAxes, fontsize=9, va='top', linespacing=1.6)
    fig.subplots_adjust(wspace=.3, hspace=.65, bottom=.1)
    save(fig, 'fig6_anatomical_definitions', [ROOT / 'publication/feature_definitions/current_120_feature_formulas.csv'],
         'Original anatomical schematic for five FaceKit indices. A: subnasale-to-upper-lip distance divided by vertical nasion-to-subnasale midface height. B: iris-center distance divided by bizygomatic width. C: outer-canthal distance divided by bizygomatic width. D: bilateral mean star-polygon eyebrow area divided by squared bizygomatic width; this measures geometry, not hair density. E: bilateral signed horizontal cheek-to-interpolated-face-oval distance divided by bizygomatic width, not three-dimensional malar projection. Blue indicates numerator geometry, orange the length scale. Points are illustrative artwork; exact landmark sets, signs and branches are defined in S6 and the pinned FaceKit source.',
         'Five original face-outline panels illustrate philtrum, iris and canthal distances, eyebrow polygons and cheek-to-oval geometry with their normalization lengths.')


def regional():
    path = TABLES / 'regional_discovery_source.csv'
    frame = pd.read_csv(path).sort_values('fraction_significant')
    fig, ax = plt.subplots(figsize=(7.1, 4.2))
    positions = np.arange(len(frame))
    ax.barh(positions, 100 * frame.fraction_significant, color=COLORS[0], height=.64)
    ax.set_yticks(positions, [textwrap.fill(label, 25) for label in frame.anatomical_region])
    for y, row in zip(positions, frame.itertuples()):
        ax.text(row.fraction_significant * 100 + .6, y, f'{row.significant}/{row.tested} ({row.fraction_significant:.1%})', va='center', fontsize=8)
    ax.set(xlim=(0, 51), xlabel='Primary hypotheses with global BH q < 0.05 (%)')
    ax.set_title('Regional discovery rates across 48 syndrome cohorts', loc='left', pad=12)
    ax.grid(axis='x', alpha=.18)
    ax.set_axisbelow(True)
    fig.tight_layout()
    save(fig, 'fig7_regional_discovery_rates', [path],
         'Significant/tested syndrome-feature counts and rates in nine anatomical regions. All tests belong to the single 5,760-hypothesis primary global-BH family (q < .05). Region denominators equal 48 syndromes times the number of features in the region and therefore differ; correlated features are not independent anatomical discoveries. Source: regional_discovery_source.csv, independently cross-checked against S1 and regional_summary_120.csv.',
         'Horizontal bars show discovery rates for all nine regions, with counts and percentages. Chin/jaw and forehead have the highest rates; denominators differ.')


def distributions():
    path = SYNTHETIC / 'anchor_distribution_source.csv'
    frame = pd.read_csv(path)
    fig, axes = plt.subplots(2, 3, figsize=(7.1, 5.6))
    for index, (ax, row) in enumerate(zip(axes.flat, frame.itertuples())):
        for x, source, color in [(0, 'real', COLORS[0]), (1, 'synthetic', COLORS[1])]:
            lower, median, upper = [getattr(row, f'{source}_{v}') for v in ['q1', 'median', 'q3']]
            ax.add_patch(plt.Rectangle((x - .22, lower), .44, upper - lower, facecolor=color, alpha=.65))
            ax.plot([x - .22, x + .22], [median, median], color='#182A35', lw=1.3)
        ax.set_xticks([0, 1], [f'GMDB\nn={row.real_n}', f'PDIDB\nn={row.synthetic_n}'])
        ax.set_xlim(-.6, 1.6)
        ax.set_title(f'{chr(65 + index)}  {short(row.syndrome)}\n{row.feature.replace("_", " ")}', fontsize=9)
        ax.set_ylabel('Normalized feature units')
        ax.ticklabel_format(axis='y', style='plain', useOffset=False)
        ax.grid(axis='y', alpha=.15)
    axes.flat[-1].axis('off')
    axes.flat[-1].text(.05, .9, 'All five predeclared anchors\n\nBoxes: Q1–Q3\nLine: median\nNo individual points or extrema\n\nModel tags pooled as the\nobserved PDIDB release mixture.\nAges unavailable for PDIDB.', va='top', fontsize=9, linespacing=1.5)
    fig.tight_layout(h_pad=2, w_pad=1.2)
    save(fig, 'fig8_synthetic_anchor_distributions', [path],
         'Real and synthetic interquartile distributions for all five predeclared syndrome-feature anchors. Boxes span Q1–Q3 and horizontal lines show medians in actual normalized feature units; no whiskers, extrema, individual observations or photographs are displayed. GMDB contributes one selected age-eligible image per recorded patient. PDIDB model tags are pooled as an observed release mixture; all synthetic ages are missing. These are within-source distribution summaries, not clinical replication.',
         'Five paired interquartile boxes compare GMDB and PDIDB for the predeclared philtrum, pupil-spacing, malar, eyebrow-area and canthal-width anchors.')


def effects():
    path = SYNTHETIC / 'effect_scatter_source.csv'
    summary = SYNTHETIC / 'syndrome_fidelity_summary.csv'
    frame = pd.read_csv(path)
    fidelity = pd.read_csv(summary)
    fidelity = fidelity.loc[fidelity.model_tag.eq('pooled') & fidelity.sensitivity.eq('unmodified')].set_index('syndrome')
    fig, axes = plt.subplots(3, 3, figsize=(7.1, 7.1), sharex=True, sharey=True)
    limit = np.ceil(np.abs(frame[['real_g', 'synthetic_g']]).max().max() * 2) / 2
    for index, (ax, (syndrome, group)) in enumerate(zip(axes.flat, frame.groupby('syndrome', sort=True))):
        row = fidelity.loc[syndrome]
        ax.plot([-limit, limit], [-limit, limit], '--', color='#929BA1', lw=.8)
        same = group.real_g * group.synthetic_g >= 0
        ax.scatter(group.real_g, group.synthetic_g, c=np.where(same, COLORS[0], COLORS[1]), s=10, alpha=.65, linewidths=0)
        ax.axhline(0, color='#CBD1D6', lw=.5)
        ax.axvline(0, color='#CBD1D6', lw=.5)
        ax.set(xlim=(-limit, limit), ylim=(-limit, limit))
        ax.set_title(f'{chr(65 + index)}  {short(syndrome)}\nr={row.pearson_g:.2f}; sign={row.direction_agreement:.0%}', fontsize=9)
    fig.supxlabel("Real GMDB vs fixed FairFace: Hedges’ g", fontsize=10)
    fig.supylabel("Synthetic PDIDB vs the SAME fixed FairFace: Hedges’ g", fontsize=10)
    fig.tight_layout(rect=(.02, .02, 1, 1), h_pad=1.5)
    save(fig, 'fig9_synthetic_effect_profiles', [path, summary],
         "Real-versus-synthetic effect profiles against the identical frozen 3,409-path FairFace reference. Each point is a syndrome-feature aggregate (120 per panel), not an individual; dashed lines indicate equality. Orange marks strictly opposite effect directions. Titles report within-syndrome Pearson correlation and direction agreement. Hedges’ g uses syndrome-minus-reference mean difference and pooled SD with the stated small-sample correction. Correlated features, a shared reference and GMDB-trained generation induce dependence; no independent confirmatory p-values or confidence intervals are inferred.",
         'Nine effect-profile scatter panels show real and synthetic syndrome-minus-reference Hedges g, an identity line and per-syndrome correlation and sign agreement.')


def fidelity_plot():
    path = SYNTHETIC / 'syndrome_fidelity_summary.csv'
    frame = pd.read_csv(path)
    frame = frame.loc[frame.model_tag.eq('pooled') & frame.sensitivity.eq('unmodified')].sort_values('syndrome')
    fig, axes = plt.subplots(1, 3, figsize=(7.1, 4.4), sharey=True)
    metrics = [('median_absolute_standardized_shift', 'A  Absolute mean shift\n(real SD units)', (.0, .35)),
               ('median_variance_ratio', 'B  Synthetic / real\nvariance ratio', (0, 1.3)),
               ('direction_agreement', 'C  Effect-direction\nagreement', (0, 1.05))]
    for ax, (metric, title, limits) in zip(axes, metrics):
        y = np.arange(len(frame))
        ax.barh(y, frame[metric], color=COLORS[0], height=.6)
        ax.set(xlim=limits, title=title)
        ax.set_yticks(y, [short(s) for s in frame.syndrome])
        ax.grid(axis='x', alpha=.15)
        if metric == 'median_variance_ratio':
            ax.axvline(1, color=COLORS[1], ls='--', lw=1)
    axes[0].invert_yaxis()
    fig.tight_layout(w_pad=1)
    save(fig, 'fig10_synthetic_fidelity', [path],
         'Per-syndrome descriptive fidelity across 120 correlated features. A: median absolute synthetic-minus-real mean shift divided by the real sample SD. B: median sample variance ratio; dashed line is equal dispersion. C: agreement of the signs of real/synthetic Hedges’ g against the same fixed FairFace reference. Model tags are pooled as the empirical release mixture; stratified and 1%/99% winsorized results are supplied separately. No bootstrap intervals or biological independence across features are assumed.',
         'Three aligned bar panels show standardized mean shift, variance preservation and effect-direction agreement for nine matched syndromes.')


def main():
    for function in [workflow, anatomy, regional, distributions, effects, fidelity_plot]:
        function()
    manifest = {'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), 'figures': INVENTORY,
                'format': 'SVG with editable vector text and 300 dpi PNG; 7.1-inch manuscript width; no image-derived art'}
    (TABLES / 'FIGURE_SOURCE_MANIFEST.json').write_text(json.dumps(manifest, indent=2) + '\n')
    text = '# Final research figure captions and alt text\n\n'
    for item in INVENTORY:
        text += f"## {item['name']}\n\n{item['caption']}\n\nAlt text: {item['alt_text']}\n\n"
    (TABLES / 'FIGURE_CAPTIONS.md').write_text(text.rstrip() + '\n')
    print('Generated six original aggregate figures as SVG and PNG')


if __name__ == '__main__':
    main()
