import argparse
import csv
import hashlib
import json
import shutil
from pathlib import Path
from xml.sax.saxutils import escape

from reportlab import Version
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image, PageBreak

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TITLE = 'PRISM_FaceKit_120_vs_125_Comparison_Report'


def read_csv(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('output', type=Path)
    out = parser.parse_args().output.resolve()
    result = json.loads((out/'comparison_metrics.json').read_text())
    result['direction_agreement_all_finite'] = {'status': 'not_run', 'numerator': None, 'denominator': None, 'fraction': None, 'wilson_95': None}
    result['new_discovery_historical_agreement'] = {'status': 'not_run', 'numerator': None, 'denominator': None, 'fraction': None}
    result['effect_distribution'] = {'status': 'not_run', 'pearson': None, 'spearman': None, 'median_absolute_change': None, 'p90_absolute_change': None}
    result['effect_tolerance_sensitivity'] = {'status': 'not_run', 'stricter': {'absolute': .10, 'relative': .15},
                                           'main': {'absolute': .20, 'relative': .25}, 'looser': {'absolute': .30, 'relative': .35},
                                           'policy': 'Package planning rule; no external preregistration or observed fractions'}
    (out/'comparison_metrics.json').write_text(json.dumps(result, indent=2)+'\n')
    manifest = json.loads((out/'reproducibility_manifest.json').read_text())
    model = json.loads((out/'model_asset_audit.json').read_text())
    baseline = result['baseline']
    mapping = result['mapping']
    highlights = read_csv(out/'clinical_highlights_comparison.csv')
    access = read_csv(out/'input_access_summary.csv')
    flow = read_csv(out/'cohort_flow.csv')
    old = manifest['historical_facekit']
    current = manifest['current_facekit_audited']
    application = manifest['application_facekit']
    source = HERE/'local_only/sources'/f'FaceKit-{current}'
    for name, path in [('current_dependency_manifest.toml', source/'pyproject.toml'),
                       ('historical_source_LICENSE.txt', ROOT/'publication/historical_facekit/LICENSE'),
                       ('current_source_LICENSE.txt', source/'LICENSE'),
                       ('canonical_depth_source.txt', source/'src/facekit/core/geometric/canonical_depth.py'),
                       ('landmark_loader_source.txt', source/'src/facekit/core/morph/landmarks.py'),
                       ('batch_loader_source.txt', source/'src/facekit/core/geometric/batch.py')]:
        shutil.copyfile(path, out/name)
    json_checked = 0
    json_hits = 0
    for folder in [ROOT/'data', ROOT/'local_review']:
        for path in folder.rglob('*.json'):
            if any(part in {'node_modules', 'test_deps', '__pycache__'} for part in path.parts) or path.stat().st_size > 15000000:
                continue
            json_checked += 1
            try:
                value = json.loads(path.read_text())
            except (ValueError, UnicodeDecodeError):
                continue
            records = value if isinstance(value, list) else [value]
            json_hits += any(isinstance(item, dict) and any(key in item for key in ['landmarks', 'landmarks_2d_px', 'landmarks_3d']) for item in records)
    (out/'supplemental_input_checks.json').write_text(json.dumps({'json_documents_inspected': json_checked,
        'top_level_landmark_key_candidates': json_hits, 'size_limit_bytes': 15000000,
        'scope': 'data and local_review JSON documents; not an exhaustive search of the computer'}, indent=2))
    changes = [
        ['All 120 outputs', 'Roll alignment -> full yaw/pitch/roll correction with fixed depth; fallback remains roll-only.'],
        ['Brow position', 'Upper brow border -> lower border, with midface-height denominator retained.'],
        ['Lateral brow thickness', 'Three-point area / width squared -> vertical upper-to-lower brow gap / width.'],
        ['Lid-to-iris measurements', 'Fissure-height denominator -> bizygomatic width; bilateral derivatives inherit the change.'],
        ['Eye area; nose/face area', 'Polar-angle reordering -> ordered contour area where a cyclic boundary exists.'],
        ['Nose length and tip width', 'Length: face height -> width; tip: landmarks 49/279 -> 45/275. Bridge rename retains 196/419.'],
        ['Philtrum and lip family', 'Philtrum width and most lip distances: mouth width -> bizygomatic width. Philtrum length retains midface height.'],
        ['Chin and forehead height', 'Chin upper anchor: 17 -> 18; heights change from face height to bizygomatic width.'],
        ['Whole-face proportions', 'Forehead width anchor in two proportions: 103/332 -> 54/284.'],
        ['Area-helper rename', 'Historical _hull_area -> _star_polygon_area retains polar-sort geometry; neither computes a convex hull.'],
    ]
    missing = [
        ['120 primary syndrome vs FairFace', 'NOT RUN', 'No authorized image-level measurements; planned 5,880 only if all 49 groups remain eligible.'],
        ['120 demographic A / B / C', 'NOT RUN', 'No 120-feature cohort. Frozen mapping and analysis-specific BH rules must be reused.'],
        ['120 photograph selection / covariates', 'NOT RUN', 'Historical selection audit only; no current measurements or estimable adjusted models.'],
        ['120 GMDB: all other syndromes', 'NOT RUN', 'Permission for this comparison does not supply images or authorize re-extraction.'],
        ['120 GMDB: restricted eligible 49', 'NOT RUN', 'Requires a separate 120-feature BH family; no carried-over historical totals.'],
        ['120 regional/top associations/effect distributions', 'UNAVAILABLE', 'No current statistical results.'],
        ['125 vs 120 scatter, signed changes, transition matrix, syndrome agreement', 'UNAVAILABLE', 'No paired effects or significance classifications; no placeholder observations plotted.'],
    ]
    claims = [
        {'claim': 'Historical cohort: 1,850 patients; 942 reference images; 379 labels; 49 eligible', 'source': 'baseline_audit.json; cohort_flow.csv', 'verification': 'Frozen selection code rerun in memory'},
        {'claim': 'Historical 6,125 tests; 811 global-BH significant', 'source': 'historical_125_point_estimate_checks.csv.gz', 'verification': 'All g, U, p, q and sample counts reconstructed'},
        {'claim': 'Max g/q error 2.22e-16; max p error 1.11e-16', 'source': 'baseline_audit.json', 'verification': 'Comparison with frozen aggregate results'},
        {'claim': 'Independent rank calculations verify 6,125 tests', 'source': 'historical_125_point_estimate_checks.csv.gz', 'verification': 'Rank sums, tie variance, normal tail, manual BH'},
        {'claim': '94 alternative photographs', 'source': 'baseline_audit.json', 'verification': 'Image comparison aligned by patient key'},
        {'claim': '120 finite outputs in synthetic schema probe', 'source': 'current_feature_order.csv; comparison_metrics.json', 'verification': 'Nine base methods; no detector or plugins'},
        {'claim': '121 changed mappings; 4 removed; 5 absent old names; no added names', 'source': 'feature_crosswalk.csv; comparison_metrics.json', 'verification': 'Includes hairline -> forehead merge'},
        {'claim': 'Strict verified mapping: 0/120 new; 0/125 old', 'source': 'feature_crosswalk.csv', 'verification': 'Full default definition includes changed alignment'},
        {'claim': 'Looser anatomical relation: 121/125 old; 120/120 new', 'source': 'feature_crosswalk.csv', 'verification': 'No equivalence inferred from names'},
        {'claim': '78/120 agree on ten synthetic roll-only probes', 'source': 'synthetic_fallback_diagnostics.csv', 'verification': 'Diagnostic only; not a clinical or equivalence metric'},
        {'claim': '0/1,850 GMDB and 0/942 FairFace references resolve', 'source': 'input_access_summary.csv', 'verification': 'Recorded paths and repository-relative candidates checked'},
        {'claim': f'{json_checked} JSON documents; {json_hits} landmark-key candidates', 'source': 'supplemental_input_checks.json', 'verification': 'Audited data/local_review; scope and size limit recorded'},
        {'claim': 'No photo/array/JSONL candidates; six tracked data CSVs flagged', 'source': 'input_inventory_summary.json; tracked_data_privacy_audit.csv', 'verification': 'Scoped inventory; values never exported to PDF'},
        {'claim': 'No missing historical feature cells or gender entries in selected cohorts', 'source': 'historical_missingness.csv', 'verification': 'Raw missing-cell counts; unknown labels not assessed'},
        {'claim': 'All three 90% fractions: N/A; current extraction gate blocked', 'source': 'comparison_metrics.json', 'verification': 'Numerators, denominators and fractions explicitly null'},
    ]
    write_csv(out/'claim_to_source.csv', claims)
    gate_rows = [['Direction agreement, |g125| >= 0.2', 'N/A / N/A', 'Not estimable; no Wilson interval'],
                 ['Effect similarity, meaningful old effects', 'N/A / N/A', 'Not estimable'],
                 ['Retention of old matched BH discoveries', 'N/A / N/A', 'Not estimable'],
                 ['Strict new-feature mapping coverage', '0 / 120 (0%)', 'Below the 90% certification gate'],
                 ['Extraction / missingness quality', 'BLOCKED', 'No authorized current extraction'],
                 ['Clinical-example gate', 'MANUAL REVIEW', 'All five full definitions change'],
                 ['Historical point estimates', 'PASS', '6,125 tests reconstructed']]
    sections = [
        ('Executive verdict', [
            ('p', 'INDETERMINATE - retain the frozen 125-feature primary analysis.'),
            ('p', 'As of October 28, 2026, the required authorized original photographs or compatible landmark inputs are unavailable in the audited repository locations. The author also confirms that image reprocessing permission remains unconfirmed. No 120-feature patient measurements or statistical results were generated.'),
            ('table', (['Decision condition', 'Observed', 'Assessment'], gate_rows)),
            ('p', 'The historical cohort and all point estimates reproduce. This establishes that the stored historical measurements remain usable for their original analysis; it does not reconstruct original image extraction or validate the current pipeline.'),
            ('p', 'The 90% threshold is an author-chosen operational rule. There is no measured 125-versus-120 agreement fraction, and this report makes no claim of clinical interchangeability.'),
            ('p', f'Historical source: {old}. Current candidate audited: {current}. Application checkout: {application}. The audited current extractor is byte-identical to the application extractor; the paper-to-commit link remains unverified.'),
        ]),
        ('Cohort reconstruction and input access', [
            ('table', (['Historical stage', 'n'], [[r['stage'], r['n']] for r in flow])),
            ('p', 'Selection is unchanged: 3 <= year + month/12 < 10; invalid, unresolved and zero-zero ages excluded; stored quality-pass images only; one image per patient chosen by the lowest max(|yaw|/15, |pitch|/15, |roll|/10), then lexical image ID. Identical reference paths are deduplicated. All selected FairFace images have stored age category 3-9; this does not establish exact-age balance.'),
            ('p', 'The eligible syndrome count is 49 under the historical valid-n and missingness policy. Recomputed sample counts match every frozen hypothesis. The selected cohort has no missing historical feature cells or raw gender entries. This does not verify annotation accuracy, exact age/sex balance, or independence of undocumented duplicate identities.'),
            ('table', (['Selected image references', 'Selected n', 'Resolvable n'], [[r['source'], r['selected_n'], r['referenced_files_resolvable_n']] for r in access])),
            ('p', f'No original-photo candidates or saved array/JSONL landmarks were found in data/, local_review/, FaceKit/ or the pinned source. A separate scan inspected {json_checked} JSON documents in data/local_review and found {json_hits} top-level landmark-key candidates. Existing PNG graphics are report/test assets, not verified original study inputs. Other folders on the computer were not searched.'),
            ('p', 'Alternative historical selection changes 94 photographs when compared by patient key. This is a reconstructed selection count, not a rerun of the 120-feature photograph sensitivity analysis. The original cohort-selection and reference rules must remain frozen even when new quality flags are recorded.'),
        ]),
        ('Feature compatibility and strict coverage', [
            ('p', 'The source-only probe confirms the exact ordered historical schema of 125 columns and a current base schema of 120 finite outputs on synthetic coordinates. Plugins, metadata and pose-quality columns are excluded from these feature totals. These probes test code behavior without processing a patient photograph.'),
            ('image', 'feature_mapping_summary.png'),
            ('p', 'Crosswalk status counts are 121 changed-definition mappings and four removed features without a direct replacement. Five old column names disappear: cupid_bow_peak_asym, hairline_height, lip_outer_to_inner_area, mouth_upper_fit_residual and nostril_region_area. hairline_height is represented by an anatomical many-to-one mapping to forehead_height, which explains the difference between five absent names and four unmapped removals. No new column names are added.'),
            ('p', 'Strict verified full-pipeline coverage is 0/120 new features and 0/125 old features. Every current default measurement uses a changed canonicalization path when a transformation matrix is present. No unconditional equivalence was certified across that change. This is a certification result, not evidence that the anatomy is unrelated or that every numerical result must differ.'),
            ('p', 'The looser anatomical mapping covers 121/125 old features (96.8%) and 120/120 current features (100%). On ten synthetic roll-only probes, 78/120 common columns agree numerically. Those restricted probes neither prove mathematical equivalence nor establish patient-level similarity; they are excluded from the strict coverage and all 90% agreement denominators.'),
            ('p', 'Full expressions, support assignments, constant indices, helper definitions and line locations are preserved in feature_crosswalk.csv, formula_details.csv, the complete old/current extractor excerpts, and canonical_depth_source.txt. The separate formula_difference_report.md provides a feature-by-feature audit.'),
        ]),
        ('Formula, landmark and preprocessing changes', [
            ('table', (['Measurement family', 'Source-grounded change'], changes)),
            ('p', 'The historical source centers at the inner-canthal midpoint and corrects roll. The current source uses a nearest pure rotation obtained by SVD and a fixed canonical depth vector scaled to observed bizygomatic width, then recenters at the same midpoint. Per-image MediaPipe z is ignored. Without a usable matrix it falls back to roll-only and reports derotated=False.'),
            ('p', 'The current frontal_ok field reflects the pose threshold even when rejection is disabled; the historical source could mark any matrix-bearing image as frontal when the check was disabled. Reapplying current flags must not silently reselect the historical cohort. Store new failures and analyze retention explicitly.'),
            ('p', 'The direct image extractor reads the supplied image; it does not establish that the original GMDB and FairFace crops, restoration, resizing, color handling or detector assets were identical. Recovered landmarks need pixel coordinates or original image dimensions to convert normalized coordinates, plus the transformation matrix for full correction. A matrix alone does not recover landmarks from normalized feature columns.'),
        ]),
        ('Agreement rules and highlighted syndromes', [
            ('p', 'All three primary 90% fractions are unavailable. Their numerators and denominators are null, rather than zero. No finite paired effect universe exists, so direction agreement for all effects, meaningful-effect agreement, Wilson intervals, effect correlations, difference quantiles and significance transitions cannot be computed.'),
            ('table', (['Syndrome / feature', 'Historical g', 'Historical global q', 'Current result'],
                       [[r['syndrome']+' / '+r['old_feature'], f'{float(r["g_125"]):.3f}', f'{float(r["q_125_global"]):.3g}', 'NOT RUN; review'] for r in highlights])),
            ('p', 'All five highlighted measurements require manual review under the complete pipeline definition: Cornelia de Lange philtrum length; Noonan interpupillary distance; Williams-Beuren eyebrow thickness and outer-canthal distance; Rubinstein-Taybi malar bulge. Their local formulas retain the same definitions on common canonical coordinates, including the renamed but equivalent brow area helper, while the upstream coordinate correction changes. No clinical contradiction or direction reversal has been measured.'),
            ('p', 'When a valid comparable set is established, use finite effects in both versions and |g125| >= 0.2 for the primary direction and effect rules. Main effect tolerance is |g120-g125| <= max(0.20, 0.25|g125|). Proposed follow-up sensitivities, fixed in this package before current-data extraction, use stricter max(0.10, 0.15|g125|) and looser max(0.30, 0.35|g125|); all are currently NOT RUN. Also report all-finite direction agreement.'),
            ('p', 'For significance retention, use original globally significant matched hypotheses as the denominator and recompute each version\'s global BH family independently. Also report the fraction of new discoveries shared historically and the complete 2x2 transition table. Newly introduced, invalid, missing or unmappable measurements must be listed separately.'),
        ]),
        ('Missing results, quality and privacy', [
            ('table', (['Requested output', 'Status', 'Reason / retained requirement'], missing)),
            ('p', 'Extraction success, new missingness, retained populations, detector artifacts and current age/sex missingness are unmeasured. The current paper describes 120 measurements and a fixed-depth pose correction, but the manuscript version does not identify an exact software commit or asset checksum. Its cohort results must not be substituted for PRISM. Source: https://arxiv.org/html/2610.09288v1, Methods 2.1-2.3.'),
            ('p', 'Six already tracked data CSVs contain identifiers, locators, demographics, genetic/clinical metadata or row-level measurements and are flagged for data-owner review in tracked_data_privacy_audit.csv. Their presence in Git does not establish redistribution permission. No source data or history was deleted or rewritten.'),
            ('p', 'All generated outputs are local and Git-ignored. The PDF contains aggregate quantities, feature names and source-code provenance only. Patient keys, selected-image locators and private candidate inventories remain outside the PDF. The author reports mentor permission for within-GMDB comparisons and derived measurements; image reprocessing and release of row-level files or linking metadata still require separate confirmation.'),
            ('p', 'GMDB-versus-FairFace source confounding remains unresolved by software agreement. Within-GMDB comparisons address a different question and require their own correction families. Clinician review of HPO correspondence, detector quality and demographic support is still needed before interpreting a geometric result clinically.'),
        ]),
        ('Prioritized next steps and restart contract', [
            ('p', '1. Recover the exact selected original GMDB and FairFace photographs, or compatible saved 478-point landmarks, together with image dimensions, transformation matrices, crop/resize history and extraction provenance. Match them against the private frozen selection manifest and raw hashes without changing patient, reference or syndrome membership.'),
            ('p', '2. Obtain explicit data-owner/mentor confirmation that reprocessing is permitted. Confirm the paper-specific FaceKit code revision and detector/model assets with the maintainer; the audited upstream commit is a candidate, not a verified paper execution revision.'),
            ('p', '3. Review the changed normalizers, landmark anchors, polygon areas, removed features and pose correction. Separate exact formula equivalence on a common coordinate system from equivalence of the full measurement process. Document many-to-one relations and clinical mappings. Reconsider whether a 90% strict-coverage gate is scientifically appropriate if intentional changes affect the whole catalogue; changing that rule requires an explicit protocol amendment, not relabeling changed definitions as unchanged.'),
            ('p', '4. Reconstruct the historical features from the same recovered inputs as an upstream control, then extract current features without silent record dropping. Report failures, flags, derotation paths and image quality by source and syndrome. If a detector, preprocessing or asset differs, describe a pipeline comparison rather than a software-only comparison.'),
            ('p', '5. Rerun the original two-sided Mann-Whitney U tests, syndrome-minus-reference Hedges g, marginal bootstrap intervals and planned BH families. Preserve valid-n >= 10 and maximum missing fraction 0.30. The planned primary family is 49 x 120 = 5,880 only if the same 49 groups remain eligible; document every invalid test and its treatment. Rerun demographic A/B/C, both within-GMDB comparator designs and supported photograph/covariate sensitivities.'),
            ('p', '6. Calculate all agreement fractions with exclusions, denominators, Wilson intervals and stricter/looser tolerance sensitivities. Add effect scatterplots, signed changes, Bland-Altman differences, significance transitions, regional summaries and syndrome-wise agreement from actual paired results. Review discordant key examples and targeted age/sex effects with the mentor and clinician.'),
            ('p', '7. Only after those steps, apply the original decision rule: propose a 120-feature manuscript package if every gate passes; retain 125 and present both if a valid comparison fails agreement; remain INDETERMINATE if mapping or execution is still unverified. Publication, live-document edits and patient-level release require author review. No replacement manuscript text is justified now.'),
        ]),
        ('Reproducibility and numerical audit', [
            ('p', f'PRISM commit: {manifest["prism_commit"]}. Historical FaceKit: {old}. Current candidate and updated application submodule: {current}. Current source archive SHA-256: {manifest["source_acquisition"]["sha256"]}. At the author\'s request, the application was updated from 63a35875bd7e633d67a495af5f7248cc74b82667 before the final audit; its previous commit is retained on a local backup branch.'),
            ('p', f'Historical extractor SHA-256: {sha(ROOT/"publication/historical_facekit/extractor.py")}. Current extractor SHA-256: {sha(source/"src/facekit/core/geometric/extractor.py")}. Canonical depth SHA-256: {sha(source/"src/facekit/core/geometric/canonical_depth.py")}. Cached detector SHA-256: {model["sha256"]}. The cached detector was not executed; use of this asset in the historical extraction is unverified.'),
            ('p', 'The declared detector asset is MediaPipe face_landmarker float16, version 1. Dependency ranges for geometric extraction include NumPy >=2,<3, OpenCV-contrib >=4.10, MediaPipe >=0.10.30, pandas >=2 and oaklib >=0.5. These ranges are not a paper execution lockfile. The full historical dependency environment and original detector checksums remain unresolved.'),
            ('p', 'Analysis environment: Python '+manifest['python']+'; '+', '.join(f'{k} {v}' for k,v in manifest['versions'].items())+'. Synthetic probe seed: 20261028; fallback seeds: 0-9. Frozen bootstrap settings are 500 draws, rank bootstrap 200, seed 2026, but intervals and ranking bootstraps were not rerun in this feasibility audit.'),
            ('p', 'The numerical control reconstructed all 6,125 historical g, U, p and q values and sample counts; maximum absolute g/q error was 2.22e-16 and p error 1.11e-16. A separate rank-sum and tie-corrected normal-tail calculation agreed for every test, and a manual BH routine was checked against statsmodels. These checks verify the stored-data calculation, not source-level causation or clinical validity.'),
            ('p', 'Twelve audit tests passed, including frozen publication-manifest hashes, unchanged source hashes, schema and merge accounting, null current-result safeguards, independent ranks/BH, photo alignment and Git-ignore controls. Relevant updated-FaceKit and PRISM tests separately report 38 passed and three skipped: two unavailable optional image fixtures and the absent full HPO mapping fixture. Dependency consistency and CLI checks pass. These are software checks, not patient-image validation. PDF text extraction and rendered-page inspection are additional artifact checks.'),
            ('table', (['Headline claims', 'Aggregate source'], [[r['claim'], r['source']] for r in claims[:6]])),
        ]),
        ('Claim-to-source appendix', [
            ('table', (['Headline claim', 'Aggregate source / check'], [[r['claim'], r['source']+'; '+r['verification']] for r in claims[6:]])),
            ('p', 'The machine-readable claim_to_source.csv includes every headline claim in this report. baseline_audit.json and cohort_flow.csv also provide the input counts (10,942 GMDB rows; 1,000 reference rows). historical_missingness.csv records 231,250 GMDB and 117,750 reference feature cells, with zero missing cells. This report does not treat those stored measurements as newly extracted current data.'),
            ('p', 'Source trail: frozen exact_feature_names.csv and stats_each_disease_vs_healthy.csv.gz; historical feature formulas and licensed extractor snapshot; stats_testing/prepare_cohorts.py and statistical_analysis.py; ancestry_analysis.py and ancestry_disease_interaction.py; reviewer_followup_2026 source for the two comparator designs and covariate/photograph extensions. Source hashes identify each inspected file in reproducibility_manifest.json.'),
            ('p', 'Deliverables include the crosswalk, complete formula audit, aggregate historical checks, clinical highlights, explicit missing-results records, decision report, feasibility blockers, source-linked mapping graphic and reproducibility manifests. The 120 analysis directory contains a NOT RUN status only. No connected Google Doc, frozen result, original photograph or external publication was edited.'),
        ]),
    ]
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name='ReportTitle', fontName='Helvetica-Bold', fontSize=25, leading=29, textColor=colors.HexColor('#123a4c'), spaceAfter=18))
    styles.add(ParagraphStyle(name='ReportHeading', fontName='Helvetica-Bold', fontSize=18, leading=22, textColor=colors.HexColor('#123a4c'), spaceAfter=15))
    styles.add(ParagraphStyle(name='ReportBody', fontName='Helvetica', fontSize=9.4, leading=13.3, spaceAfter=10, splitLongWords=True))
    styles.add(ParagraphStyle(name='Cell', fontName='Helvetica', fontSize=8.1, leading=10.7, alignment=TA_LEFT, splitLongWords=True))
    story, markdown = [], ['# PRISM FaceKit 120 vs 125 Comparison Report', '', 'October 28, 2026 | Aggregate-only feasibility and source audit', '']
    for index, (heading, blocks) in enumerate(sections):
        if index:
            story.append(PageBreak())
        else:
            story += [Paragraph('PRISM', styles['ReportTitle']), Paragraph('FaceKit 120 vs 125', styles['ReportTitle']),
                      Paragraph('October 28, 2026 | Aggregate-only feasibility audit', styles['ReportBody']), Spacer(1, 12)]
        story.append(Paragraph(heading, styles['ReportHeading']))
        markdown += ['## '+heading, '']
        for kind, content in blocks:
            if kind == 'p':
                story.append(Paragraph(escape(content), styles['ReportBody']))
                markdown += [content, '']
            elif kind == 'image':
                story.append(Image(str(out/content), width=6.55*inch, height=2.08*inch))
                markdown += [f'![Feature mapping audit]({content})', '', 'Source: mapping_summary_data.csv; feature_crosswalk.csv.', '']
            else:
                headers, rows = content
                data = [[Paragraph(escape(str(cell)), styles['Cell']) for cell in row] for row in [headers]+rows]
                weights = {2: [2.45, 4.2], 3: [2.65, 1.0, 3.0], 4: [3.1, .8, 1.0, 1.75]}[len(headers)]
                table = Table(data, colWidths=[weight*inch for weight in weights], repeatRows=1, hAlign='LEFT')
                table.setStyle(TableStyle([('BACKGROUND', (0,0), (-1,0), colors.HexColor('#dcecef')),
                    ('ROWBACKGROUNDS', (0,1), (-1,-1), [colors.HexColor('#f3f6f8'), colors.white]),
                    ('VALIGN', (0,0), (-1,-1), 'TOP'), ('LEFTPADDING', (0,0), (-1,-1), 7),
                    ('RIGHTPADDING', (0,0), (-1,-1), 7), ('TOPPADDING', (0,0), (-1,-1), 7),
                    ('BOTTOMPADDING', (0,0), (-1,-1), 7), ('LINEBELOW', (0,0), (-1,0), .6, colors.HexColor('#247c87'))]))
                story += [table, Spacer(1, 10)]
                markdown += ['| '+' | '.join(headers)+' |', '| '+' | '.join(['---']*len(headers))+' |']
                markdown += ['| '+' | '.join(str(cell).replace('|','/') for cell in row)+' |' for row in rows]+['']
    def footer(canvas, doc):
        canvas.setStrokeColor(colors.HexColor('#247c87'))
        canvas.line(42, 40, 570, 40)
        canvas.setFont('Helvetica', 8)
        canvas.setFillColor(colors.HexColor('#506875'))
        canvas.drawString(42, 27, 'PRISM | INDETERMINATE | Local aggregate review | 2026-10-28')
        canvas.drawRightString(570, 27, str(doc.page))
    doc = SimpleDocTemplate(str(out/(TITLE+'.pdf')), pagesize=(612, 792), rightMargin=42, leftMargin=42,
                            topMargin=43, bottomMargin=53, title='PRISM FaceKit 120 vs 125: Feasibility Report', author='PRISM research audit')
    doc.build(story, onFirstPage=footer, onLaterPages=footer)
    (out/(TITLE+'.md')).write_text('\n'.join(markdown))
    (out/'decision_report.md').write_text('# Decision: INDETERMINATE\n\n'+ '\n'.join(markdown[4:])+'\n')
    (out/'feasibility_blockers_and_next_steps.md').write_text('# Feasibility blockers and next steps\n\n'+
        '\n\n'.join(content for heading, blocks in sections if heading in ['Executive verdict', 'Prioritized next steps and restart contract'] for kind, content in blocks if kind == 'p')+'\n')
    detail = ['# Detailed formula-difference report', '',
        'Status is based on the complete measurement definition, including canonicalization. Anatomical relations are not certified equivalences. No 120-feature patient values were computed.', '',
        'All local source expressions are accompanied by the complete licensed extractor sources and helper/constant definitions. Names such as face_height and mouth_width denote canonical-space scales; expressions retain the exact source assignment locations.', '',
        '| Family | Change |', '| --- | --- |']
    detail += ['| '+' | '.join(row)+' |' for row in changes]
    detail += ['', '## Full crosswalk formulas', '']
    for row in read_csv(out/'formula_details.csv'):
        detail += ['### '+row['old_feature']+' -> '+(row['new_feature'] or '(removed)'), '',
                   'Historical source: '+row['old_source'], '', '```python', row['old_formula'], '```', '',
                   'Current source: '+(row['new_source'] or 'No verified direct match.'), '', '```python', row['new_formula'], '```', '']
    (out/'formula_difference_report.md').write_text('\n'.join(detail))
    (out/'aggregate_audit_report.md').write_text('# Aggregate numerical audit\n\n'+
        json.dumps(baseline, indent=2)+'\n\nTwelve automated tests passed. Independent checks cover every historical point estimate; original detector execution and bootstrap intervals were not reproduced. See claim_to_source.csv.\n')
    for path in HERE.glob('*.py'):
        manifest['source_sha256'][str(path.relative_to(ROOT))] = sha(path)
    manifest['report_environment'] = {'reportlab': Version, 'python_executable_role': 'bundled Codex document runtime; separate from analysis environment'}
    (out/'reproducibility_manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    output_hashes = {p.name: sha(p) for p in out.iterdir() if p.is_file() and p.name != 'output_hashes.json'}
    (out/'output_hashes.json').write_text(json.dumps(output_hashes, indent=2)+'\n')
    print('Created aggregate PDF, Markdown, decision, formula audit, blockers and claim-to-source table.')


if __name__ == '__main__':
    main()
