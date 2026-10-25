#!/usr/bin/env python3
import json
import re
from pathlib import Path
import sys

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import followup as f


def main():
    raw, _, primary, _, _, _, _, _ = f.load()
    manifest = json.loads((f.FROZEN/'publication_manifest.json').read_text())
    rows = []
    hashes = manifest.get('files_sha256', {})
    if not hashes:
        candidates = [value for value in manifest.values() if isinstance(value, dict) and value and
                      all(isinstance(v, str) and len(v) == 64 for v in value.values())]
        if len(candidates) != 1:
            raise ValueError('Unrecognized frozen manifest schema')
        hashes = candidates[0]
    for name, expected in hashes.items():
        path = f.FROZEN/name
        actual = f.sha(path)
        rows.append({'file': name, 'expected_sha256': expected, 'actual_sha256': actual, 'passed': expected == actual})
    f.check(all(r['passed'] for r in rows), 'Frozen input integrity failed')
    f.write(pd.DataFrame(rows), 'frozen_manifest_validation.csv')
    metadata = pd.read_csv(f.OUT/'metadata_availability.csv')
    fields = ['pose_yaw', 'pose_pitch', 'pose_roll', 'frontal_ok', 'expression', 'occlusion', 'image_width',
              'image_height', 'resolution', 'gaze_annotation', 'landmarks']
    imaging = metadata.loc[metadata.field.isin(fields)].copy()
    imaging['sensitivity_status'] = imaging.field.map(lambda x: 'stored estimate only; selection/audit available' if x.startswith('pose_') or x == 'frontal_ok'
                                                     else 'blocked: no verified measured field in supplied master')
    f.write(imaging, 'imaging_covariate_inventory.csv')
    f.write(pd.DataFrame([{'analysis': family, 'status': 'not computed', 'blocker':
                          'No validated expression/gaze/occlusion/resolution fields or authorized historical landmark/image source supplied',
                          'required_next_step': 'data-owner authorization, measured flags and blinded validation before analysis'}
                         for family in ['iris_offset_y', 'forehead_hairline', 'cheek_malar']]), 'imaging_covariate_sensitivity.csv')
    records = []
    for source, frame in [('GMDB_original', raw), ('GMDB_selected', primary)]:
        for field in ['gene_names', 'gene_entrez_ids', 'hgvs', 'test_type', 'hpo_terms', 'present_features', 'absent_features', 'genetic_confirmation']:
            if field in frame:
                text = frame[field].astype('string').fillna('').str.strip()
                usable = ~text.str.lower().isin(['', 'nan', 'none', 'null', 'unknown', '[]', '{}', 'na', 'n/a'])
                n = int(usable.sum())
                codes = text.str.contains(r'HP[:_]\d{7}', regex=True).sum() if field in ['hpo_terms', 'present_features', 'absent_features'] else 0
            else:
                n, codes = 0, 0
            records.append({'source': source, 'field': field, 'column_present': field in frame,
                            'n_rows': len(frame), 'nonempty_n': n, 'missing_or_empty_n': len(frame)-n,
                            'rows_with_explicit_HPO_code_pattern': int(codes),
                            'confirmation_status': 'not inferred from field population'})
    f.write(pd.DataFrame(records), 'genetics_hpo_feasibility.csv')
    (f.OUT/'component_analysis_report.md').write_text('''# Numerator/denominator feasibility

**Status: not computable from the normalized 125-feature master alone.** No verified historical raw-landmark or unnormalized component table was supplied for this task. No patient images or landmark tracks were opened or processed. A feature dictionary and its dependency graph establish formulas, not per-patient numerator/denominator values. A ratio does not uniquely identify its two components.

`highlighted_components.csv` and `feature_dictionary.csv` preserve actual historical definitions. CdLS philtrum length is the 2-0 landmark distance divided by aligned midface height |y168-y2|. Noonan interpupillary distance is the 468-473 distance divided by the 234-454 bizygomatic distance. Eyebrow thickness is a convex-hull area normalized by squared bizygomatic width, with the mean averaging bilateral ratios; it is not a single eyebrow-width distance or hair-density measurement. Mean malar bulge averages the two historical contour-based descriptors, which linearly interpolate the two same-side face-oval points nearest the cheek landmark in aligned y (or average x when their y values coincide), then divide the signed horizontal contour-to-cheek difference by bizygomatic width. Refer to the full expressions and dependencies, rather than replacing them with a simplified cheek-width arrow.

Next research, after an authorized historical geometry source is identified: audit detector/model and image-scale provenance; estimate numerator and denominator sensitivity separately; retain within-face scale normalization or calibrated image scale; check correspondence and missingness. Uncalibrated raw pixels from different photographs cannot be compared as anatomical sizes. Do not invert normalized columns to fabricate component magnitudes. No component-analysis p-values or biological numerator conclusions are produced.
''')
    (f.OUT/'imaging_covariate_sensitivity_report.md').write_text('''# Imaging-covariate inventory and blocked analyses

**Verified existing fields:** stored frontal flag and detector-derived yaw/pitch/roll estimates, plus historical geometric iris/gaze proxies. The primary cohort requires stored frontal pass. Pose ranking chooses one eligible photograph; it is not a new quality flag. `imaging_covariate_inventory.csv` reports presence and missingness for raw and selected/reference cohorts.

**Unavailable measured fields:** the masters do not contain explicit width/height/resolution, expression, occlusion, validated clinical gaze, or raw-landmark columns. The presence of iris-offset/gaze-asymmetry features does not establish true gaze labels or classifier accuracy. No historical raw geometry source has been verified and authorized here. Existing path/ID columns are not evidence that a requested image pilot is authorized.

**Not run:** expression/gaze/occlusion/resolution flags, blinded manual validation, and flagged/unflagged sensitivity for vertical iris offset, forehead/hairline, or cheek/malar families. `imaging_covariate_sensitivity.csv` records these blockers, not placeholder findings. No model-estimated flags are mislabeled as observed annotations.

**Future protocol:** obtain data-owner approval and an explicitly bounded image/landmark source, prespecify objective flags and a blinded manual sample with assessors unaware of feature effects, record missingness and abstention, quantify agreement and error with intervals, then define a retained cohort and correction families before comparing effects. Report both patient and feature-level attrition; account for within-patient repeats. This would be a measurement-validation study and cannot establish elimination of dataset-source confounding. No patient photograph processing is performed in this follow-up.
''')
    (f.OUT/'HPO_GENETICS_FEASIBILITY.md').write_text('''# HPO/genetics feasibility

`genetics_hpo_feasibility.csv` gives original/selected nonempty counts for available gene, variant/test and HPO fields without exporting their contents. A string with a nonempty gene/variant or test-type field is not a confirmed diagnosis. No explicit genetic-confirmation field exists; unknown, clinical and genetic confirmation categories cannot be assigned from syndrome names or guessed from metadata. Explicit HPO-code patterns are counted only as syntactic presence, not validated clinical annotation or completeness.

`candidate_hpo_mapping.csv` proposes eight directional correspondences for highlighted feature families. All identifiers/labels are checked against the official Human Phenotype Ontology OBO source; its version and SHA-256 are recorded in `hpo_source_provenance.json`. The download contains public ontology data only; no patient or study metadata is transmitted. Identifier verification is distinct from clinical mapping approval. Multiple feature definitions can correspond to one candidate term, and one measure can have multiple candidate directional correspondences. All rows remain uncertain and require clinician approval. No HPO recovery rate is calculated.

In particular, interpupillary/bizygomatic ratio is not an age-normalized diagnostic criterion for hypertelorism, eyebrow landmark-hull area is not hair density, and a 2D malar contour measure is not a clinical assessment of bone projection. Clinicians must decide whether an intended term, direction and threshold apply, accounting for age, facial normalization and acquisition. This package does not assign HPO terms to individual patients or assert that the four syndromes have those phenotypes.

Public ontology source: https://github.com/obophenotype/human-phenotype-ontology ; exact downloaded URL/version/checksum are in the provenance JSON.
''')
    print(f'Frozen integrity: {len(rows)} hashes passed; component/imaging/genetics feasibility reports written')


if __name__ == '__main__':
    main()
