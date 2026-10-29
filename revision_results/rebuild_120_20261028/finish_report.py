import csv
import gzip
import hashlib
import json
import os
import platform
import shutil
import subprocess
from pathlib import Path
from xml.sax.saxutils import escape

from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image, PageBreak

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
RUN = HERE/'local_only/run_20261028'
EXTRACTION = ROOT/'local_review/feature_extraction_120/20261028_current120'
os.umask(0o077)


def read_json(path):
    return json.loads(Path(path).read_text())


def rows(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', newline='') as handle:
        return list(csv.DictReader(handle))


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8*1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def write(name, text):
    (RUN/name).write_text(text.strip()+'\n')


def documentation():
    write('SOURCE_AND_PERMISSION_AUDIT.md', '''# Source and permission audit

Status: completed private research extraction and numerical validation; active 120-feature analysis promoted. No public release or live Google Doc edit was performed.

## GMDB

Source: `GestaltGANwglabshare/gmdb_enhanced_v1.1.0_new`; cached build `a94fb856b30eaeeea670b9f06db903afc44e29db`. The two Arrow streams contain 8,888 and 2,092 rows (10,980 total), 8,346 recorded patient IDs and 581 syndrome IDs. `image_restored` holds embedded PNG bytes; all images decoded as 224 x 224 RGB. No encoded-byte duplicates were found. These are restored images; undocumented earlier transformations and byte identity to historical extraction inputs remain unknown.

The user's handoff conveys mentor authorization for internal research reprocessing under relevant IRB paperwork. This run relies on that attestation; the agent did not independently inspect the IRB protocol or expand the permitted study scope. Original photographs must not appear in papers, supplements or public repositories. Raw images, patient keys, measurements and linkage files remain private.

Original cache: `~/.cache/huggingface/datasets/GestaltGANwglabshare___gmdb_enhanced_v1.1.0_new/default/0.0.0/a94fb856b30eaeeea670b9f06db903afc44e29db/`. Independent private backup: `revision_results/facekit_120_comparison_20261028/local_only/hf_cache_audit_20261028/private_backup/`; its parent contains hashes, sizes and timestamps. Reprocessing read that backup, with no in-place cache edits. All 10,980 metadata rows and all 10,942 historical measurement keys matched; 38 cache rows had no historical measurement. All 1,850 historical chosen keys were available.

## FairFace

Official source: https://github.com/joojs/fairface. The user supplied the official margin-0.25 archive and train/validation label links and explicitly authorized their download. Files are private under `local_review/feature_extraction_120/20261028_current120/inputs/`. ZIP size 578,232,678 bytes; CRC passes. All 97,698 unique label paths correspond to archive entries, with no duplicate/conflicting paths or cross-split path collisions. There are 11,764 age-category 3-9 candidates. The source repository states CC BY 4.0; no photographs were redistributed. Dataset race, gender and age-bin labels are annotations. FairFace is a reference dataset, not verified clinically healthy or verified unique people.

## PDIDB

Source: https://huggingface.co/datasets/GestaltGANwglabshare/PDIDB, revision `d0907a56525dc3b25712d360f9ce3031150eb535`, manually gated. The existing authenticated Hugging Face session successfully accessed the official provider; no gate was bypassed. The SDK used `token=True`; credentials were not printed, copied into outputs or committed. An initial metadata serialization KeyError was corrected and was not a provider access denial.

One train Parquet file, 140,243,455 bytes, contains 1,727 rows in 18 row groups. The image struct holds bytes/path. Other fields are image_id, disease, disease_abbr, model, verification, age_range, gender, race, facial_features and prompt. All age_range/gender/race values are missing. Model tags v2/v1 have 912/815 records; verification tags low/clinician have 1,673/54. These tags do not establish generator identity, clinical validity or source independence. The user identifies the source as synthetic. No patient-lineage field is supplied, and generation method/real-source linkage remain unresolved. Normalized exact disease-label matching is only a taxonomy crosswalk, not a patient link.

The card provides no license or substantive generation description. Access does not establish illustration redistribution rights. No PDIDB image was inserted in publication/report figures. There are no exact encoded-byte duplicates or byte matches to GMDB; this check cannot establish privacy or independence. All synthetic observations are excluded from primary and real-patient sensitivity inference. Optional synthetic-source agreement remains unperformed because linkage/design are unresolved.

## Existing tracked-data concern

Historical per-row files, including `data/master_dataset_disease.csv`, `data/metadata_with_images.csv` and phenotype CSVs, remain in the repository. They were preserved unchanged as requested, not copied into new aggregate deliverables. PI review is needed for redistribution/history remediation. No history rewrite, commit, push, image upload or cloud patient-data logging was performed. The private run and extraction folders are gitignored; new private files use owner-only access. Historical archives are read-only.
''')
    schema = RUN/'FEATURE_SCHEMA_120.csv'
    shutil.copy2(ROOT/'publication/feature_definitions/current_120_feature_formulas.csv', schema)
    write('FEATURE_SCHEMA_AND_EXTRACTION_AUDIT.md', '''# Feature schema and extraction audit

Exactly 120 ordered numeric morphology features are defined in `FEATURE_SCHEMA_120.csv`, including exact formulas, units, current source line pointers and anatomical grouping. This table is a derived transcription; the pinned source is authoritative. Quality/pose and metadata columns are excluded from the feature count. No feature values, p-values or q-values were inferred from historical columns.

FaceKit source pin: `5764bdd4af9360b55a046a16579fb5a890a8182f`; extractor SHA-256 `8657c0fa5e9aad74d300c8c9bbf023e35f6c5f0a0f74559b2febb43f57fe0a0e`. This is the verified current checkout used by every source. The exact paper execution commit is unverified; arXiv 2610.09288v1 alone does not prove it. Historical extractor `e474628f3cb17a83edad59ac6e1d9e5438753032` is preserved, with zero certified unconditional full-pipeline equivalences to the new definition.

Model: local MediaPipe Face Landmarker float16 v1, 478 landmarks; `~/.cache/facekit/face_landmarker.task`, SHA-256 `64184e229b263107bc2b804c6625db1341ff2bb731874b0bcc2fe6544e0bc9ff`. Decoding uses Pillow RGB, without new resizing, cropping, restoration or EXIF transpose. An up-to-two-face detector screens ambiguity, followed by the default one-face detector for full floating-point pixel landmarks and pose matrix. Current FaceKit full matrix canonicalization and normalizations produce measurements with `frontal_check=False`, followed by explicit absolute yaw/pitch/roll gates 15/15/10 degrees. Acceptance also requires derotation and all 120 finite features. Images, embeddings and landmark arrays are not output. Extraction blocks networking and cannot silently fetch assets.

| Source | Attempted | 120 finite vectors | QC accepted | Pose failure | No face | Ambiguous |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| GMDB | 10,980 | 10,936 | 9,404 | 1,532 | 39 | 5 |
| FairFace 3-9 | 11,764 | 10,594 | 5,864 | 4,730 | 1,145 | 25 |
| PDIDB synthetic | 1,727 | 1,726 | 1,563 | 163 | 0 | 1 |

Pose-failed finite vectors remain explicitly marked in all-record catalogs. Missing detection vectors remain missing. Every primary/reference role table has zero missing feature cells. The full PDIDB catalog has 120 missing cells from one ambiguous record.

GMDB primary selection uses frozen conditional-age handling (3 <= age < 10), recorded patient IDs, original syndrome/annotation mappings, current QC, one image per recorded patient, normalized absolute pose score and lexical image-key tie break. There are 1,845 selected recorded patients across 379 labels, including 872 patients in the 48 main eligible syndromes. Compared with historical selection: six patients lost, one newly measurable patient added, one selected photo changes among shared patients. FBXW7 has nine patients and is below the ten-observation main threshold.

FairFace: full candidate pool screened first. A is up to 500 per seven source race labels (3,409; Middle Eastern 409 and others 500). B uses all 5,864 accepted paths. C uses 151 in each of 14 race-by-Female/Male cells (2,114). Sampling is seeded, without replacement. C misses the nominal 250-per-cell goal by 99 per cell; A's shortfall is 91 in one category. Source labels/unknowns are preserved. Exact person counts cannot be established from unique paths. Historical 942 paths yield 942 finite vectors but only 941 pass current QC.

The extraction validator passed 28 independent completed-run checks; 59 replayed pilot statuses and 55 finite vectors reproduce within 1e-10. Thirteen extraction unit tests pass. PDIDB uses the identical pipeline and saved pilot replay diagnostics. All candidates, exclusions, model/code/input hashes, checkpoints and per-source provenance are private. See `UNIFIED_120_MANIFEST.json`, `CATALOG/`, `PDIDB/`, the earlier extraction `PACKAGE_SHA256.json`, and `VALIDATION_FULL.json` for executed evidence.
''')
    write('MANUSCRIPT_HANDOFF.md', '''# Manuscript handoff: validated current 120-feature analysis

The active `manuscript/main.tex`, `publication/` and `stats_results/` now use the validated new analysis. The live Formal Paper Google Doc was not edited. The new TeX draft compiled successfully in the desktop compiler. Aggregate figures contain no patient photographs.

| Section / asset | Replacement and source |
| --- | --- |
| Abstract | 120 features; 1,845 selected patients; 3,409 FairFace reference paths; 48 eligible syndromes; 5,760 primary hypotheses; 955 global-BH discoveries, 482 also abs(g)>=0.8. `REPORT/REPORT_NUMBERS.json` and main table. |
| Data / ethics | Restored GMDB, margin-0.25 FairFace and separate synthetic PDIDB. User-conveyed mentor research authorization; author must verify exact IRB scope/identifier. `SOURCE_AND_PERMISSION_AUDIT.md`. |
| Extraction methods | Current pinned FaceKit, one cached detector/model, 478 landmarks, RGB, explicit QC, exact schema and formulas. `FEATURE_SCHEMA_AND_EXTRACTION_AUDIT.md`, `FEATURE_SCHEMA_120.csv`. |
| Cohort flow | Six historical patients lost, one enters; 48 main syndromes after FBXW7 falls to nine. References A/B/C 3,409/5,864/2,114. `ANALYSIS_PLAN.json`, private cohort eligibility tables. |
| Statistical methods | Historical test signs, eligibility and correction scopes; main marginal bootstrap 500 and rank bootstrap 200, seed 2026; separate demographic and extension families. Main table is recomputed, never translated from old q-values. |
| Main results / Figure 1 | 48-syndrome x 12 selected-feature descriptive heatmap. Full 5,760-test table remains S1; the heatmap is not a separate discovery family. `fig1_cross_syndrome_heatmap.png`. |
| Highlighted examples / Figure 2 | Current g/95% marginal CI: Cornelia philtrum 1.163 [0.886,1.412]; Noonan interpupillary 1.244 [0.757,1.684]; Williams eyebrow area 1.266 [0.987,1.587], outer-canthal 1.039 [0.780,1.307]; RTS malar -1.009 [-1.530,-0.514]. These are normalized geometric indices, not diagnoses. `clinical_highlights_120.csv`. |
| Demographics / Figure 3 | Broad and detailed A/B/C tables and global correction scopes. Say recorded demographic categories, not verified ancestry or causal effects. Sparse interaction cells require caution. |
| Sensitivities / Figure 4 | Alternative photo 996 discoveries, 918 shared/37 primary-only/78 alternative-only; 93 patient selections change. Reference B/C/fixed941: 971/945/919. Overlapping reference specifications are not independent replications. |
| Within-GMDB / extensions | 222 discoveries against all other age-eligible syndromes; 259 with current-eligible comparators. Fixed49 adjusted/pairwise/age families include unsupported cells. Full tables in `RESULTS_120/extensions/`. No diagnostic specificity claim. |
| Historical comparison | Original 49x125/6,125/811 preserved only as frozen history. Common-record benchmarks use 1,844 GMDB and 941 reference paths: old125 798/6,000 vs new120 919/5,760. Matched names are not identical definitions. |
| Limitations | Dataset/source confounding, restored-image provenance, changed normalizations, pose/landmark/expression effects, unknown FairFace identity and health, sparse models, correlated measurements and biased concordance resampling. |
| PDIDB / data availability | 1,727 synthetic records, 1,563 accepted, zero pooled primary rows. Generation lineage, source independence and image-publication rights unresolved; no image illustrations used. |

Concordance resampling must not be presented as clean uncertainty confirmation. For all-feature comparison against other GMDB, observed direction agreement is 0.713 and Pearson r=0.688, while frozen percentile intervals are [0.731,0.775] and [0.725,0.773]. Both exclude their observed point, and bootstrap bias exceeds estimated SE. Bias flags and full draws summaries are retained. Use descriptive points and explicitly disclose the failure; evaluate an improved uncertainty design separately rather than hiding these intervals.

Author review still needed: exact permitted study scope and IRB identifier/document; authorship, funding and contributions; clinical interpretation and proxy measurement wording; provider confirmation of PDIDB generation/lineage and redistribution terms; privacy remediation of historical tracked per-row files; intended venue and final release statement. No new clinical-performance claim, image-publication permission or submission readiness is implied.

Primary citations in the draft: FaceKit arXiv 2610.09288v1 (exact paper execution pin unresolved), GestaltMatcher (2022), FairFace (2021). Use the existing reference bibliographic metadata and official source links; verify venue-specific citation format before submission.
''')
    write('CHANGELOG_120_MIGRATION.md', '''# 120-feature migration changelog

Decision: **PROMOTED / VALIDATED_ACTIVE**, analysis `PRISM_120_20261028`. The user explicitly allowed a validated new 120-feature primary analysis and then requested replacing the active 125-feature analysis. This supersedes the earlier 90% interchangeability promotion rule. Numerical validation and source-role checks passed; the new definitions are not asserted interchangeable with old values.

## Active replacements

- `stats_results/significance_testing/`: new 120-feature primary tables, bootstraps, robust tests, ranking, reliability and score outputs. `ethnicity_analysis/`, `extensions_120/` and real 120 sensitivities are regenerated. Old125 benchmark is private comparison evidence only.
- `publication/`: current schema/formulas, fresh aggregate plots and supplementary tables. Stale figures/audits removed from active exports. Historical licensed extractor snapshot remains explicitly labeled. Generator entry points export the verified current artifacts and check actual 120/5,760/955 counts.
- `manuscript/`: full current draft, fresh figures/supplements and reviewer handoff; compiled successfully. Old draft moved to archive.
- `pds_results_v2/`: old standalone score exports/calculator retired; the export entry point and README use the current score in main tables.
- `ACTIVE_ANALYSIS.json`, `stats_testing/analysis_config.py`: default inputs and schema resolve to private validated current cohorts, with no historical fallback. Main/A/B/C defaults, cohort audit, data-distribution plots and prepare_cohorts updated. Main requires exactly the ordered 120 schema.
- `backend/app/main.py`: dashboard measurement distributions use the current private 120 cohorts. SQLite imported from current results; old database backed up privately. Frontend feature counts/text use 120 and reference terminology.
- `README.md`: current counts, paths and reproduction commands replace old default-input recipes.
- Bootstrap optimization preserves literal historical draw ordering and numerical results. Original executed code and regression helpers were frozen before active routing changes. A checkpoint routing migration points one helper import to its frozen copy; original script and statuses remain in `EXECUTED_CODE/`.

## Preserved history

`archive/historical_125/` is an independently checksummed read-only copy of original measurements, source, full historical outputs and manuscript. `archive/retired_active_exports/` preserves the moved active folders and prepromotion SQLite backup. Modes are owner-read-only for files and owner read/execute for folders. The original Hugging Face shards are untouched and have an independently verified private backup. Git history was not rewritten. Historical `data/*.csv` remain unchanged and require a PI privacy-remediation decision.

## Executed verification

- Historical original 6,125 tests / 811 discoveries reproduced with U/g/p/q numerical error <=2.22e-16.
- Full current numerical validator: PASS; main U exact, p error <=1.11e-16, g error <=4.44e-16; all six comparator specifications, demographic scopes, adjusted and pairwise models, selected literal CI replays and concordance arithmetic checked.
- System Python repository tests: 41 passed, one API test skipped because system FastAPI is absent. Virtualenv targeted API/bootstrap/publication/reanalysis tests: 19 passed, including live in-process API distribution checks.
- Extraction completed-run validator: 28 checks passed; 13 extraction unit tests included in the system suite.
- Publication current audit: PASS, 120 features / 5,760 hypotheses / 955 discoveries. Backend import: completed.
- New manuscript: desktop compiler success.
- Frontend type/build validation: blocked by missing node_modules. Offline `npm ci --offline --ignore-scripts --no-audit --no-fund` returned ENOTCACHED for which-typed-array-1.1.22.tgz. No registry download was attempted to override the offline constraint. TSX edits were inspected; do not claim a successful frontend build.

See `VALIDATION_FULL.json`, `FINAL_DELIVERY_AUDIT.json`, private logs and `REPRODUCIBILITY_MANIFEST.json` for exact evidence. The report was rendered and visually reviewed after generation. No commit, push, deployment, patient-data upload or live document edit was performed.
''')
    write('README.md', '''# Validated 120-feature rebuild

Active analysis: `PRISM_120_20261028`. Main inputs: 1,845 GMDB recorded patients and 3,409 FairFace reference paths; 48 eligible syndromes x 120 features = 5,760 primary tests, 955 global-BH discoveries. PDIDB is separate synthetic extraction only. `ACTIVE_ANALYSIS.json` at repository root selects the active private data/schema.

Read the PDF report, `SOURCE_AND_PERMISSION_AUDIT.md`, `FEATURE_SCHEMA_AND_EXTRACTION_AUDIT.md`, `MANUSCRIPT_HANDOFF.md`, `CHANGELOG_120_MIGRATION.md` and `VALIDATION_FULL.json`. Exact schema: `FEATURE_SCHEMA_120.csv`. Role tables and catalog: `CATALOG/`, `UNIFIED_120_MANIFEST.json`. Full source/error/linkage records remain private. The unified catalog is not an inferential pool.

## Reproduction and checkpoints

Run commands from `/Users/joshua/Documents/PRISM`. Extraction and Arrow work use the existing system Python with PyArrow/MediaPipe; statistics use `.venv/bin/python`. The PDF generator uses bundled Python with ReportLab. Environments and hashes are in `REPRODUCIBILITY_MANIFEST.json` and `ENVIRONMENT/`.

```sh
python3 revision_results/rebuild_120_20261028/prepare.py --dry-run
.venv/bin/python revision_results/rebuild_120_20261028/run_sensitivities.py --dry-run
.venv/bin/python revision_results/rebuild_120_20261028/extensions.py --dry-run
PYTHONPATH="$PWD/FaceKit/src" python3 -m pytest tests -q
.venv/bin/python -m pytest tests/test_active_120.py tests/test_bootstrap_pair_effects.py tests/test_publication_source.py tests/test_reanalysis.py -q
.venv/bin/python revision_results/rebuild_120_20261028/validate.py --scope full
.venv/bin/python publication/manuscript_results_audit.py
```

To reproduce the main statistical analysis without touching validated results:

```sh
.venv/bin/python stats_testing/statistical_analysis.py --output-dir local_review/NEW_120_REPRODUCTION
```

The primary defaults retain 500 effect/200 ranking bootstraps and seed 2026. For the exact six historical/current benchmark specifications, `run_sensitivities.py` uses the frozen `EXECUTED_CODE/stats_testing/` and hashes all inputs. Completed jobs verify output hashes and skip work; interrupted jobs rerun that specification. `extensions.py` similarly verifies completed stage signatures/hashes, or reruns an interrupted stage. Changed numerical code/inputs require a new run; do not overwrite checkpoints. `EXECUTED_CODE/ROUTING_MIGRATION.json` documents the source-location-only adjustment after retirement of old active code.

Image extraction resumes with `--resume` and exact input/model/code/environment signatures. See repository `data/dataset_creation/CURRENT_120_EXTRACTION.md` for FairFace/GMDB pilot/full commands. PDIDB resumes with `python3 revision_results/rebuild_120_20261028/pdidb_extract.py --resume`; it reads the authenticated download's verified private backup, not credentials. No network fetch occurs in extraction. No landmark arrays or image thumbnails are exported.

Prepared frozen cohorts and completed catalog must not be regenerated with a different schema in this run. `prepare.py` skips existing cohort preparation, verifies the archive and maintains an inventory. `promote.py` is a one-time operation and refuses to replace its existing retirement archive; do not rerun it. Current publication wrappers copy the verified aggregate artifacts. Full original generation code is versioned beside this README, with all inputs available only to an authorized local user.

Historical 125 outputs are in `archive/historical_125/` and `archive/retired_active_exports/`. Do not remove them or rewrite history. Source images and per-person tables are restricted even when upstream is hosted online. Publication release, image illustration rights and IRB identifiers need author review; current promotion is a local research workflow decision.
''')
    (HERE/'README.md').write_text('''# Current PRISM 120-feature rebuild

The validated new analysis is active. It contains 120 features, 1,845 selected
GMDB patients, 3,409 FairFace reference paths, 48 eligible syndromes, 5,760 main
hypotheses and 955 global-BH-significant associations. Synthetic PDIDB is kept
separate. See repository `ACTIVE_ANALYSIS.json` for the active configuration.

All protected inputs, row-level results, read-only historical archives and the
full PDF/report package are under ignored `local_only/run_20261028/`. Its README
documents provenance, commands, dry runs, checkpoint resumption and validation.
The historical 125-feature archive is preserved; active exports have been replaced.
No live Google Doc edit, deployment, public upload, commit or push was performed.

Scripts here prepare authorized data, extract synthetic PDIDB, build the separate
source catalog, rerun statistics/extensions, independently validate results,
generate manuscript materials, promote once and produce the private report.
Use the private run README before executing them. The exact paper execution pin
and image-publication permissions remain distinct from technical validation.
''')


def build_pdf():
    number = read_json(RUN/'REPORT/REPORT_NUMBERS.json')
    old = read_json(RUN/'REPORT/historical_named_feature_comparison.json')
    clinical = rows(ROOT/'publication/supplementary/clinical_highlights_120.csv')
    sensitivity = rows(RUN/'REPORT/reference_and_photo_sensitivity.csv')
    paired = rows(RUN/'REPORT/paired_old_new_common_record_results.csv.gz')
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name='Body120', fontName='Helvetica', fontSize=9.4, leading=13.5, spaceAfter=8))
    styles.add(ParagraphStyle(name='Small120', fontName='Helvetica', fontSize=7.5, leading=10, spaceAfter=5))
    styles.add(ParagraphStyle(name='Title120', fontName='Helvetica-Bold', fontSize=25, leading=29, textColor=colors.HexColor('#183b50'), spaceAfter=17))
    styles.add(ParagraphStyle(name='Heading120', fontName='Helvetica-Bold', fontSize=17, leading=21, textColor=colors.HexColor('#183b50'), spaceAfter=13))
    styles.add(ParagraphStyle(name='Cell120', fontName='Helvetica', fontSize=8, leading=10.5, alignment=TA_LEFT))
    story = []
    markdown = ['# PRISM 120-feature reanalysis and source comparison', '', 'Validated active research analysis | 28 October 2026', '']

    def title(text, first=False):
        story.append(Paragraph(escape(text), styles['Title120' if first else 'Heading120']))
        markdown.extend(['', '## '+text, ''])

    def p(text, small=False):
        story.append(Paragraph(escape(text), styles['Small120' if small else 'Body120']))
        markdown.extend([text, ''])

    def table(header, data, widths):
        wrapped = [[Paragraph(escape(str(v)), styles['Cell120']) for v in row] for row in [header]+data]
        item = Table(wrapped, colWidths=[v*inch for v in widths], repeatRows=1, hAlign='LEFT')
        item.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e5eff3')),
            ('VALIGN',(0,0),(-1,-1),'TOP'),('LINEBELOW',(0,0),(-1,0),.7,colors.HexColor('#527386')),
            ('LINEBELOW',(0,1),(-1,-1),.3,colors.HexColor('#d7dfe3')),
            ('LEFTPADDING',(0,0),(-1,-1),6),('RIGHTPADDING',(0,0),(-1,-1),6),
            ('TOPPADDING',(0,0),(-1,-1),6),('BOTTOMPADDING',(0,0),(-1,-1),6)]))
        story.extend([item, Spacer(1, 10)])
        markdown.append('| '+' | '.join(map(str,header))+' |')
        markdown.append('| '+' | '.join(['---']*len(header))+' |')
        markdown.extend('| '+' | '.join(map(str,row))+' |' for row in data)
        markdown.append('')

    def figure(path, max_height=3.1):
        item = Image(str(path))
        scale = min(6.6*inch/item.imageWidth, max_height*inch/item.imageHeight)
        item.drawWidth = item.imageWidth*scale
        item.drawHeight = item.imageHeight*scale
        story.extend([item, Spacer(1,8)])
        markdown.extend(['![Aggregate statistical figure]('+str(path.relative_to(RUN))+')', ''])

    def next_page():
        story.append(PageBreak())

    title('PRISM 120-feature reanalysis and source comparison', True)
    p('Validated active research analysis | 28 October 2026 | Private author review')
    p('The current 120-feature analysis replaces the historical 125-feature analysis in active PRISM results, manuscript materials and dashboard inputs. All measurements were extracted from images with one pinned current pipeline, and all reported tests were recomputed. The old research record remains in a checksummed read-only private archive.')
    table(['Primary analysis', 'Verified value'], [
        ['Selected real GMDB recorded patients', '1,845 across 379 recorded syndrome labels'],
        ['FairFace reference', '3,409 unique image paths, age category 3-9'],
        ['Eligible main syndromes / morphology features', '48 / 120'],
        ['Two-sided primary hypotheses', '5,760; one global Benjamini-Hochberg family'],
        ['Significant q < 0.05', '955 (16.58%); 557 positive and 398 negative g'],
        ['Significant with absolute Hedges g >= 0.8', '482'],
        ['Synthetic observations in primary tests', '0']], [3.25,3.35])
    p('The 48 eligible syndromes contain 872 selected patients. All 1,845 patients contribute to the pooled disease analysis and eligible within-GMDB comparator designs; the remaining sparse labels are not forced into the main syndrome-by-feature family. The alternative-photo analysis and reference variants are sensitivities, not additional independent replications.')
    p('Promotion follows the user\'s explicit decision to allow a validated new 120-feature primary analysis. The earlier 90% compatibility rule is superseded. This decision does not certify interchangeability of old and new measurements, publication rights for photographs, or clinical performance.')
    p('Key artifacts: ACTIVE_ANALYSIS.json; RESULTS_120/; VALIDATION_FULL.json; UNIFIED_120_MANIFEST.json; MANUSCRIPT_HANDOFF.md; CHANGELOG_120_MIGRATION.md. The active TeX draft compiled successfully. No public release, live Google Doc edit or patient-image illustration was performed.', True)

    next_page(); title('Sources, access and private storage')
    table(['Source', 'Available data', 'Role and constraints'], [
        ['Restored real GMDB', 'Two readable Arrow streams: 8,888 + 2,092 = 10,980 images; 8,346 recorded patients; 581 syndrome IDs', 'Restricted research disease source. Mentor/IRB authorization conveyed by user; actual protocol not independently reviewed. No original images in publication.'],
        ['Official FairFace margin 0.25', '97,698 label/image paths; 11,764 candidates labeled age 3-9', 'Reference source with recorded age/race/gender annotations. Clinical health and person uniqueness are not verified.'],
        ['Gated synthetic PDIDB', '1,727 rows, one train Parquet, 18 row groups', 'Authenticated official access succeeded. Synthetic-only catalog; generation lineage and illustration redistribution rights unresolved.']], [1.4,2.2,3.0])
    p('GMDB image_restored contains embedded PNG bytes, not just external filenames. Every cached image decoded as 224 x 224 RGB. All historical metadata keys and all 10,942 historical measurement keys matched the cache, leaving 38 cache-only records; all 1,850 historical primary selected keys were available. Restored bytes and matching keys do not establish identity to the image bytes used in the old extraction.')
    p('Before reprocessing, both shards were independently copied to private storage and verified by SHA-256, size and timestamps. The original cache remained unchanged. Old measurements, statistical code, results, figures and manuscript were frozen in archive/historical_125; retired active folders and the prepromotion SQLite database were separately preserved. New row-level data and checkpoints are gitignored and owner-only.')
    p('Existing tracked historical data/*.csv files were preserved unchanged. Their patient-level metadata and redistribution exposure require a PI remediation decision; no history rewrite was performed. No patient data were uploaded. Hugging Face access used the existing SDK session without printing, copying or committing credentials.')
    p('Source references: github.com/joojs/fairface; huggingface.co/datasets/GestaltGANwglabshare/gmdb_enhanced_v1.1.0_new; huggingface.co/datasets/GestaltGANwglabshare/PDIDB. Dataset access and research use do not automatically authorize photograph publication. See SOURCE_AND_PERMISSION_AUDIT.md for local paths and exact restrictions.', True)

    next_page(); title('One extraction definition and observed attrition')
    p('FaceKit is pinned to 5764bdd4af9360b55a046a16579fb5a890a8182f. Every source uses the same cached MediaPipe Face Landmarker float16 v1, 478 landmarks, full matrix pose correction and current geometric normalizations. The exact paper execution commit is unverified; the source pin and exact ordered formulas are recorded rather than inferred from the paper\'s feature count.')
    p('Pillow decodes RGB; no additional crop, resize, restoration or EXIF transpose is applied. A two-face ambiguity screen precedes default one-face detection. Acceptance requires successful derotation, all 120 finite features, and absolute yaw/pitch/roll within 15/15/10 degrees. Model, code and input signatures prevent silent changes on resume. Images, landmarks and embeddings are not exported.')
    table(['Source', 'Attempted', 'QC pass', 'Pose fail', 'No face', 'Ambiguous'], [
        ['GMDB', '10,980', '9,404', '1,532', '39', '5'],
        ['FairFace 3-9', '11,764', '5,864', '4,730', '1,145', '25'],
        ['PDIDB synthetic', '1,727', '1,563', '163', '0', '1']], [1.6,1,1,1,1,1])
    p('Finite vectors also exist for pose-failed records: GMDB 10,936 total and PDIDB 1,726. Failed detections/ambiguity retain missing values rather than proxies from old columns. The unified catalog contains 24,471 records with source and synthetic flags. Primary disease and all reference role tables have zero missing morphology cells.')
    p('GMDB retains the original conditional-age rule, 3 <= age < 10, recorded patient grouping and syndrome/annotation mappings. Current QC and deterministic normalized-pose/lexical-key ranking select one photograph per recorded patient. Six historical patients are lost, one newly measurable patient enters, and one photo changes among shared patients, yielding 1,845. FBXW7 has nine patients and is excluded by the unchanged n >= 10 main threshold.')
    table(['Reference specification', 'Sampling and size'], [
        ['A: primary balanced by seven source race labels', '500 each except Middle Eastern 409; total 3,409. Shortfall 91.'],
        ['B: all eligible', 'All 5,864 accepted age-category paths.'],
        ['C: category-by-reported-sex balance', '151 in each of 14 Female/Male cells; total 2,114.'],
        ['Historical fixed-key sensitivity', '942 finite vectors; 941 accepted by current QC.']], [2.75,3.85])
    p('ZIP CRC, label schema, split/path uniqueness and full path correspondence pass. Sampling is seeded without replacement, after full screening. These reference specifications overlap. The same people cannot be ruled out across different FairFace paths.', True)

    next_page(); title('Primary results and highlighted measurements')
    p('Two-sided Mann-Whitney U tests compare each eligible syndrome with reference images. Hedges g is syndrome minus reference, using the historical small-sample correction; positive values indicate a larger current geometric index. Eligibility requires at least ten finite values per side and at most 30% missingness. All 5,760 main cells are supported. One global BH correction yields 955 discoveries.')
    figure(RUN/'ACTIVE_PUBLICATION/figures/fig2_clinical_indices.png',2.65)
    labels = ['Cornelia / philtrum', 'Noonan / interpupillary', 'Williams / eyebrow area', 'Williams / outer-canthal', 'RTS / malar']
    table(['Predefined example', 'Disease n', 'g [95% CI]', 'Global BH q'], [
        [label,row['n_disease'],f"{float(row['hedges_g']):.3f} [{float(row['ci_low']):.3f}, {float(row['ci_high']):.3f}]",f"{float(row['q_global']):.3g}"]
        for label,row in zip(labels,clinical)], [2.45,.65,2.25,1.25])
    p('All five examples use 3,409 reference paths. Intervals are marginal percentile bootstraps with 500 resamples, seed 2026; they are not multiplicity-adjusted confidence intervals. Rank stability uses 200 resamples. The optimized Cliff delta bootstrap preserves the literal historical interleaved random-draw stream, including subsequent draws.')
    p('Philtrum is normalized by midface height; interpupillary and outer-canthal widths by bizygomatic width. Eyebrow thickness is a geometric polygon-area index, not hair thickness. Malar geometry is not a calibrated clinical cheek phenotype. The exact source formulas are supplied in FEATURE_SCHEMA_120.csv. These associations do not establish diagnosis, calibrated anthropometry or clinical prediction accuracy.', True)

    next_page(); title('Reference and photograph sensitivities')
    table(['Specification', 'Reference n', 'Hypotheses', 'BH discoveries'], [
        [r['specification'].replace('_',' '),r['reference_n'],r['hypotheses'],r['significant']]
        for r in sensitivity[:5]], [3.4,1,1.1,1.1])
    figure(RUN/'ACTIVE_PUBLICATION/figures/fig4_photo_selection.png',2.8)
    p('The deterministic alternative photograph changes the selected record for 93 of 1,845 patients. It yields 996 discoveries, with 918 shared, 37 primary-only and 78 alternative-only. These changes are reported without selecting the more favorable result. The reference sensitivities retain the same new disease cohort and pipeline; they show consequences of reference composition/size, not independent validation.')
    p('One hundred random-photo replicates also preserve one photo per recorded patient and apply the full 5,760-test family each time. Per-hypothesis effect ranges, significance frequencies and Monte Carlo standard errors remain in random_photo_sensitivity.csv. The procedure measures image-selection sensitivity conditional on this source and cohort; it is not an external replication.')
    p('Full source: REPORT/reference_and_photo_sensitivity.csv, RESULTS_120/alternative_photo, reference_all_eligible, reference_race_sex, reference_fixed_historical and extensions/random_photo_*.csv. All point tests and BH corrections were independently checked.', True)

    next_page(); title('Historical comparison and changed definitions')
    p('The frozen original analysis reproduces 49 eligible syndromes x 125 features = 6,125 tests and 811 BH discoveries. Current cohort attrition, reference expansion and changed measurement definitions prevent treating 955 versus 811 as an isolated improvement due to software.')
    table(['Benchmark', 'Disease / reference', 'Test family / discoveries'], [
        ['Frozen original historical', '1,850 / 942', '49 x 125 = 6,125 / 811'],
        ['Common selected keys, old values', '1,844 / 941', '48 x 125 = 6,000 / 798'],
        ['Common selected keys, new extraction', '1,844 / 941', '48 x 120 = 5,760 / 919'],
        ['New primary, expanded reference', '1,845 / 3,409', '48 x 120 = 5,760 / 955']], [2.8,1.5,2.3])
    p(f"Across {old['paired_hypotheses']:,} matched anatomical names, effect correlation is Pearson r={old['pearson']:.3f} and Spearman rho={old['spearman']:.3f}. Among {old['meaningful_old_effects']:,} old effects with absolute g >= 0.2, direction agreement is {100*old['direction_agreement_fraction']:.1f}% and {100*old['effect_tolerance_fraction']:.1f}% fall within max(0.2, 25% of old absolute g). Only {old['retained_significant_n']} of {old['old_matched_significant_n']} matched old discoveries remain significant ({100*old['retained_fraction']:.1f}%).")
    comparisons=[]
    for label,item in zip(labels,clinical):
        match=[r for r in paired if r['feature']==item['feature'] and r['disease'].split('__')[0]==item['syndrome']]
        assert len(match)==1
        row=match[0]
        comparisons.append([label,f"{float(row['hedges_g_old']):.3f}",f"{float(row['hedges_g_new']):.3f}",f"{float(row['mannwhitney_q_fdr_global_old']):.3g} / {float(row['mannwhitney_q_fdr_global_new']):.3g}"])
    table(['Common-key example', 'Old g', 'New g', 'Old / new q'],comparisons,[2.8,.75,.75,2.3])
    p('Matching a feature name does not certify the same construct. Full pose correction, scale denominators and local formulas changed; zero unconditional full-pipeline equivalences are certified. The old/new tables are descriptive compatibility evidence, not converted measurements. Matched keys also do not prove original pixel-byte identity, so the common-key contrast cannot perfectly isolate code from historical image preprocessing.')
    p('A new 120-feature primary analysis was expressly approved despite these incompatibilities. No statement of 90% measurement identity or retained significance is supported. Full paired effects, including disagreements, remain in REPORT/paired_old_new_common_record_results.csv.gz.', True)

    next_page(); title('Demographic and approved extension families')
    p('A is reference-only Kruskal-Wallis; B compares syndrome/reference within recorded groups; C fits HC3 disease-by-group interactions. Broad and detailed annotation maps and correction scopes are frozen from the historical workflow. Main A uses 120 features per layer; pairwise A corrections are separate within each reference-group pair. B and C each use a separate global BH family per layer.')
    table(['Layer / analysis', 'Hypotheses', 'BH discoveries'], [
        [r['layer'].replace('layer1_broad_ethnicity','Broad').replace('layer2_subcategory','Detailed')+' / '+r['analysis'].replace('_',' '),f"{r['hypotheses']:,}",f"{r['significant']:,}"]
        for r in number['demographic_families']], [4.15,1.25,1.2])
    p('C retains the original computational support minimum of two disease and twenty reference observations in required cells. This can yield unstable sparse estimates; HC3 does not establish biological ancestry, causal effects or reliable small-sample inference. B uses 1,000 effect bootstraps with seed 42. All subgroup tables and failures remain available.')
    table(['Extension family', 'Planned', 'Supported', 'Discoveries'], [
        [label,f"{r['planned']:,}",f"{r['supported']:,}",f"{r['significant']:,}"]
        for label,r in zip(['Within GMDB: all other labels','Within GMDB: other eligible labels','Fixed49 age/sex-adjusted screen','186 fixed anchor pairs x 120','Two age strata x fixed49 x 120'],number['extension_families'])], [3.45,1,1.05,1.1])
    p('Unsupported fixed-family cells retain correction p=1, missing effect/q and an explicit status, preserving the planned denominator. Four predefined age/sex-adjusted within-GMDB examples and three reference models with sex, broad category and yaw/pitch/roll retain separate four- and three-model families. Their coefficients are adjusted mean differences in normalized units, not Hedges g.')
    p('Correlation clustering uses current balanced-reference Spearman correlation, average linkage on 1 - absolute rho and cut 0.3. It is descriptive; no cluster-based reduction of the discovery denominator is made. Family definitions and complete negative/unsupported rows are retained in RESULTS_120/extensions/FAMILY_DEFINITIONS.json and tables.', True)

    next_page(); title('Source confounding, concordance and synthetic data')
    p('Against all other age-eligible GMDB patients, 222 associations pass BH, of which 173 overlap the 955 FairFace discoveries. There are 782 FairFace-only and 49 within-GMDB-only discoveries. The observed effect direction agreement across all 5,760 cells is 0.713; Pearson r is 0.688. Within-GMDB effects are smaller in absolute value in 63.4% of cells, with median absolute attenuation 0.109. This is a comparison of two estimands, not a diagnostic specificity test.')
    p('The frozen joint 1,000-resample concordance method has a material bias problem: its direction-agreement percentile interval [0.731, 0.775] and Pearson interval [0.725, 0.773] exclude their observed points. Bias exceeds the estimated bootstrap SE. These failures are explicitly flagged in the complete summary. The intervals must not be used as clean uncertainty confirmation. A revised resampling/measurement-error design is a separate proposed experiment; no failed interval was silently repaired or omitted.')
    table(['Separate PDIDB extraction', 'Observed result'], [
        ['Source revision', 'd0907a56525dc3b25712d360f9ce3031150eb535'],
        ['All / accepted / finite-vector records', '1,727 / 1,563 / 1,726'],
        ['Model tags', 'v2: 912; v1: 815; generator identity unknown'],
        ['Verification tags', 'low: 1,673; clinician: 54; clinical validity not established'],
        ['Missing age_range / gender / race', '1,727 for each field'],
        ['Source patient lineage / independent subjects', 'Not supplied / unknown'],
        ['Encoded-byte duplicates / GMDB byte overlap', '0 / 0; neither result proves privacy or independence'],
        ['Primary or real-patient inferential contribution', '0; no pooled synthetic disease observations']], [3.0,3.6])
    p('The provider card lacks a license and substantive generation description. Synthetic illustration use is therefore not established, and no image illustrations were produced. A disease-label crosswalk is not source linkage. Optional synthetic-source agreement remains blocked by missing linkage and a justified independent design; extraction completion does not resolve those scientific or rights questions.')
    p('Disease status is confounded with GMDB-versus-FairFace source. Restored images, residual pose, expression, occlusion and landmark errors can affect results. Normalized indices may change through denominators and are not calibrated 3D anthropometry. FairFace has binned ages and unverified health/identity. Correlated features and overlapping comparators do not establish realized FDR.', True)

    next_page(); title('Validation, reproducibility and author handoff')
    table(['Executed check', 'Outcome'], [
        ['Full independent numerical validator', 'PASS: primary U exact, p error <=1.11e-16 and g <=4.44e-16; all six comparator specifications and BH checked.'],
        ['Demographics / regression / intervals', 'A/B/C scopes; independent HC3 matrix checks; pairwise tests; seven predefined models; actual literal CI replays checked.'],
        ['Repository / integration tests', 'System: 41 pass, 1 API skip (no system FastAPI). Virtualenv targeted: 19 pass, including API current-reference check.'],
        ['Extraction / checkpoints', '28 completed-run checks; deterministic pilot replay; six sensitivity checkpoints reverified.'],
        ['Publication / manuscript', 'Current audit passes 120 / 5,760 / 955. New TeX draft compiles successfully.'],
        ['Frontend build', 'Unverified: node_modules absent; offline npm cache missing which-typed-array-1.1.22.tgz (ENOTCACHED).']], [2.1,4.5])
    p('Exact role tables, schema, pinned model/source hashes, output hashes, environment versions and validation logs are available in the private run. The statistical default reads ACTIVE_ANALYSIS.json and rejects a non-120 schema. Reproduce to a new private output directory; do not overwrite historical inputs or validated checkpoints. Interrupted extraction resumes with --resume after signature checks. Statistical comparator/extension stages verify completed output hashes and rerun interrupted stages.')
    p('The active replacement includes results, current publication generators, full draft, supplementary tables, score routing and dashboard database/distribution inputs. Historical original outputs and executed source remain read-only. Old per-row tracked files remain unchanged for PI remediation review. No live document edit, commit, push, deployment or public upload was performed.')
    p('Before submission, the authors must verify the exact IRB scope/identifier, authorship, funding and contributions; review clinical wording; decide historical tracked-data remediation; obtain PDIDB generation/lineage and redistribution terms; and approve aggregate data/code release. A scientifically useful next study would use a held-out real disease/reference source with compatible imaging, identity-aware splitting and blinded geometric validation. Concordance uncertainty needs a bias-aware design before inferential use.')
    p('Reading map: SOURCE_AND_PERMISSION_AUDIT.md documents rights; FEATURE_SCHEMA_AND_EXTRACTION_AUDIT.md and FEATURE_SCHEMA_120.csv define measurements; MANUSCRIPT_HANDOFF.md maps sections/figures; CHANGELOG_120_MIGRATION.md lists replacements; README.md gives commands; REPRODUCIBILITY_MANIFEST.json and FINAL_DELIVERY_AUDIT.json identify final files and checks.', True)
    p('References: FaceKit, arXiv:2610.09288v1 (exact paper execution commit unresolved); Hsieh et al., GestaltMatcher, Nature Genetics (2022); Karkkainen and Joo, FairFace, WACV (2021). Local numerical evidence and source revisions govern the reported results. No patient images, identifiers, landmark arrays or access credentials are included in this report.', True)

    def footer(canvas, doc):
        canvas.saveState()
        canvas.setStrokeColor(colors.HexColor('#c4d3db'))
        canvas.line(.65*inch,.53*inch,7.95*inch,.53*inch)
        canvas.setFont('Helvetica',7)
        canvas.setFillColor(colors.HexColor('#526a76'))
        canvas.drawString(.65*inch,.35*inch,'PRISM | Current 120-feature analysis | Private author review')
        canvas.drawRightString(7.95*inch,.35*inch,str(doc.page))
        canvas.restoreState()

    pdf=RUN/'PRISM_120_Reanalysis_and_Source_Comparison_Report.pdf'
    doc=SimpleDocTemplate(str(pdf),pagesize=(8.5*inch,11*inch),rightMargin=.65*inch,leftMargin=.65*inch,
        topMargin=.65*inch,bottomMargin=.72*inch,title='PRISM 120-feature reanalysis and source comparison',author='PRISM research workflow')
    doc.build(story,onFirstPage=footer,onLaterPages=footer)
    write('PRISM_120_Reanalysis_and_Source_Comparison_Report.md','\n'.join(markdown))


def provenance():
    env=RUN/'ENVIRONMENT';env.mkdir(exist_ok=True)
    interpreters={'statistics':ROOT/'.venv/bin/python','extraction':Path('/Library/Frameworks/Python.framework/Versions/3.11/bin/python3')}
    environments={}
    for label,python in interpreters.items():
        result=subprocess.run([str(python),'-m','pip','freeze'],capture_output=True,text=True,check=True)
        (env/(label+'_pip_freeze.txt')).write_text(result.stdout)
        environments[label]={'python':str(python),'version':subprocess.check_output([str(python),'--version'],text=True).strip(),
            'pip_freeze_sha256':sha(env/(label+'_pip_freeze.txt'))}
    files=[]
    for folder in [RUN/'RESULTS_120',RUN/'COHORTS',RUN/'CATALOG',RUN/'PDIDB',RUN/'EXECUTED_CODE',ROOT/'publication',ROOT/'manuscript',ROOT/'stats_testing',ROOT/'backend/app',ROOT/'frontend/components']:
        for path in sorted(folder.rglob('*')):
            if path.is_file() and path.suffix not in {'.pyc','.sqlite3'} and '__pycache__' not in path.parts and not path.name.endswith(('-wal','-shm')):
                files.append({'path':str(path.relative_to(ROOT)),'size_bytes':path.stat().st_size,'sha256':sha(path)})
    for path in sorted(HERE.glob('*.py')):
        files.append({'path':str(path.relative_to(ROOT)),'size_bytes':path.stat().st_size,'sha256':sha(path)})
    for path in [ROOT/'ACTIVE_ANALYSIS.json',ROOT/'README.md',RUN/'FEATURE_SCHEMA_120.csv',RUN/'VALIDATION_FULL.json',RUN/'PRISM_120_Reanalysis_and_Source_Comparison_Report.pdf']:
        files.append({'path':str(path.relative_to(ROOT)),'size_bytes':path.stat().st_size,'sha256':sha(path)})
    manifest={'status':'VALIDATED_ACTIVE','analysis_id':'PRISM_120_20261028','platform':platform.platform(),
        'git_head_before_uncommitted_migration':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        'facekit_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT/'FaceKit',text=True).strip(),
        'pdidb_revision':'d0907a56525dc3b25712d360f9ce3031150eb535','environments':environments,'outputs':files,
        'source_backups':'Separate SHA-256 manifests in archive/historical_125 and private GMDB/PDIDB backups',
        'credentials':'No environment dump or credential file inspected; package metadata only'}
    write('REPRODUCIBILITY_MANIFEST.json',json.dumps(manifest,indent=2))
    claims=[['primary_patients',1845,'ANALYSIS_PLAN.json'],['reference_images',3409,'ANALYSIS_PLAN.json'],
        ['primary_hypotheses',5760,'RESULTS_120/primary_balanced/stats_each_disease_vs_healthy.csv'],
        ['global_BH_discoveries',955,'RESULTS_120/primary_balanced/stats_each_disease_vs_healthy.csv'],
        ['positive_significant_g',557,'REPORT/REPORT_NUMBERS.json'],['negative_significant_g',398,'REPORT/REPORT_NUMBERS.json'],
        ['large_significant_abs_g',482,'REPORT/REPORT_NUMBERS.json'],['pdidb_records',1727,'PDIDB/extraction_summary.json'],
        ['pdidb_accepted',1563,'PDIDB/extraction_summary.json'],['synthetic_primary',0,'UNIFIED_120_MANIFEST.json']]
    with (RUN/'CLAIM_TO_SOURCE.csv').open('w',newline='') as handle:
        writer=csv.writer(handle);writer.writerow(['claim','value','source']);writer.writerows(claims)


if __name__=='__main__':
    assert read_json(RUN/'VALIDATION_FULL.json')['status']=='VALIDATED'
    documentation()
    build_pdf()
    provenance()
    print('Private report, audits, manuscript handoff and reproducibility manifest written')
