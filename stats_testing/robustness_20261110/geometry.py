import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from .common import (DEFAULT_OUT, SEED, canonicalize, case_rows, corrected, effect, features, load_private,
                     offline, private_output, protect_baseline, require, u_test, write_aggregate)
from .stability import patient_samples, fixed_reference_effect


def alternative_ratios(frame, feature):
    require(frame.face_aspect_ratio.gt(0).all(), 'Invalid face-aspect normalizer')
    if feature == 'outer_canthal_distance':
        require(frame.midface_width.gt(0).all(), 'Invalid midface-width normalizer')
        return {'original_bizygomatic': frame[feature].to_numpy(),
                'height_normalized': (frame[feature] * frame.face_aspect_ratio).to_numpy(),
                'midface_width_normalized': (frame[feature] / frame.midface_width).to_numpy()}
    require(feature == 'philtrum_width' and frame.mouth_width.gt(0).all(), 'Unsupported clinical-geometry case or mouth normalizer')
    return {'original_bizygomatic': frame[feature].to_numpy(),
            'height_normalized': (frame[feature] * frame.face_aspect_ratio).to_numpy(),
            'mouth_width_normalized': (frame[feature] / frame.mouth_width).to_numpy()}


def numerical_geometry_simulation():
    rows = []
    for numerator_change in [-.1, 0, .1]:
        for denominator_change in np.linspace(-.15, .15, 31):
            baseline_numerator, baseline_denominator = 60., 140.
            numerator = baseline_numerator * (1 + numerator_change)
            denominator = baseline_denominator * (1 + denominator_change)
            original = baseline_numerator / baseline_denominator
            ratio = numerator / denominator
            rows.append(dict(simulated_numerator_change=numerator_change, simulated_denominator_change=denominator_change,
                             relative_ratio_change=ratio / original - 1,
                             interpretation='numerical image-space geometry only; no patient measurement or millimeter inference'))
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    offline()
    out = private_output(args.output_dir)
    protect_baseline(out)
    real, reference, raw, alternate, paths = load_private()
    definitions = features().set_index('feature')
    require('bizygomatic_width /' in definitions.loc['face_aspect_ratio', 'formula'], 'Face-aspect ratio direction differs from the verified formula')
    require('UPPER_LIP_PEAK_R' in definitions.loc['philtrum_width', 'formula']
            and 'UPPER_LIP_PEAK_L' in definitions.loc['philtrum_width', 'formula'], 'Philtrum-width construct differs from the verified upper-lip peaks')
    cases = [('WILLIAMS-BEUREN SYNDROME; WBS__93e4a782', 'outer_canthal_distance'),
             ('Cornelia de Lange syndrome__fb4a389d', 'philtrum_width')]
    rows = []
    for case_index, (syndrome, feature) in enumerate(cases):
        target = real.loc[real.disease.eq(syndrome)]
        for comparator, other in [('FairFace', reference), ('other_GMDB', real.loc[real.disease.ne(syndrome)])]:
            x, y = alternative_ratios(target, feature), alternative_ratios(other, feature)
            for definition_index, definition in enumerate(x):
                a, b = x[definition], y[definition]
                g, p = effect(a, b), u_test(a, b)[1]
                samples = patient_samples(a[:, None], target.patient_id, 200,
                    np.random.default_rng([SEED, 40, case_index, definition_index]))
                draws = fixed_reference_effect(samples, b[:, None])[:, 0]
                rows.append(dict(syndrome=syndrome, feature=feature, definition=definition, comparator=comparator,
                                 n_target=len(a), n_comparator=len(b), mean_target=a.mean(), mean_comparator=b.mean(),
                                 median_difference=np.median(a) - np.median(b), hedges_g=float(g), p_value=float(p),
                                 supported=min(len(a), len(b)) >= 10, units='dimensionless derived ratio',
                                 case_resampling_q025=np.quantile(draws, .025), case_resampling_q975=np.quantile(draws, .975),
                                 range_role='conditional 200-case-patient bootstrap quantile range; not a calibrated clinical CI',
                                 scope='new exploratory normalization comparison; original primary estimates/q remain unchanged'))
    table = canonicalize(pd.DataFrame(rows))
    table = corrected(table, 'p_value', table.supported, 'two predeclared geometry cases x three definitions x two comparators; 12 planned')
    write_aggregate(table, out / 'aggregate_candidates/P1_normalization_cases.csv')
    write_aggregate(numerical_geometry_simulation(), out / 'aggregate_candidates/P1_ratio_geometry_simulation.csv')
    matrix = pd.read_csv(out / 'evidence/P0_exhaustive_evidence_matrix.csv', low_memory=False)
    write_aggregate(case_rows(matrix), out / 'aggregate_candidates/P1_geometry_case_original_evidence.csv')
    text = ('# P1 clinical geometry case studies\n\n'
        'Raw component experiment: BLOCKED_PRIVATE_INPUT. Compatible saved landmark coordinates and the separate scale denominators are not present in the current measurement exports. '
        'Normalized columns cannot identify the original numerator and denominator independently. Real-image transformations: BLOCKED_IMAGE_REPROCESSING_APPROVAL; '
        'no pixel access, new extraction, photo export, cropping or brightness manipulation was performed.\n\n'
        'The following algebraic alternative normalizations are identifiable from existing current features. Face aspect is bizygomatic_width/(face_height+epsilon), '
        'so multiplying a width-normalized feature by face_aspect_ratio produces a height-normalized distance, NOT dividing by that ratio. '
        'The same-width-normalized midface_width and mouth_width permit ratio-of-ratios definitions. These are new exploratory estimands, not replacements of the original feature. '
        'The 12-cell new BH family includes both cases, their original and two alternate definitions, and both comparators; it does not alter any frozen primary q.\n\n'
        'Williams–Beuren outer_canthal_distance is outer-corner separation divided by bizygomatic width. Its original positive FairFace effect and within-GMDB association must be interpreted '
        'with different comparator distributions and its full 5,880-model age/sex screen q, not the separate four-model q (which did not target this feature). '
        'Absolute anthropometric reports use a different scale and cannot by themselves contradict this normalized ratio.\n\n'
        'Cornelia de Lange philtrum_width uses upper-lip-peak separation (landmarks 37 and 267), not necessarily the two clinically defined philtral columns. '
        'The long-philtrum description instead concerns length; philtrum_length is a positive related-trait control, not validation of the lip-peak width construct. '
        'No absence of a term in a clinical summary is taken as evidence for novelty. No absolute millimeter dimensions can be inferred.\n\n'
        'P1_ratio_geometry_simulation.csv is a tested numerical demonstration: hold numerator fixed and change the denominator, or change both. '
        'A ratio changes under normalization alone; this does not establish that such a change caused the observed patient association. '
        'Component distributions, real-photo perturbation effects and clinical correspondence of the exact ratios remain unvalidated. '
        'Positive controls retain Cornelia de Lange length, Noonan pupil spacing and Williams–Beuren eyebrow polygon area; eyebrow area is not hair density.\n\n')
    for row in table.itertuples():
        text += (f'{row.syndrome.split("__")[0]}, {row.feature}, {row.definition}, {row.comparator}: '
                 f'g={row.hedges_g:+.6g}, exploratory 12-cell q={row.exploratory_q:.6g}, n={row.n_target}/{row.n_comparator}.\n\n')
    (out / 'P1_clinical_geometry_case_studies.md').write_text(text)
    protect_baseline(out)
    print('Existing-feature ratio sensitivities and numerical geometry simulation complete; raw/image experiments explicitly blocked')


if __name__ == '__main__':
    main()
