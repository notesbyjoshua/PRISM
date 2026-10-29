import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
RUN=ROOT/'revision_results/rebuild_120_20261028/local_only/run_20261028'


def main():
    validation=json.loads((RUN/'VALIDATION_FULL.json').read_text())
    if validation['status']!='VALIDATED':
        raise ValueError('Current full numerical validation is required')
    numbers=json.loads((RUN/'ACTIVE_STUDY_REPORT/REPORT_NUMBERS.json').read_text())
    out=ROOT/'local_review/prism_120_primary_20261028'
    out.mkdir(parents=True,exist_ok=True)
    lines=['# PRISM: statistical investigation using FaceKit','',
        'The sole active study uses 120 pinned FaceKit measurements, 1,845 recorded GMDB patients and 3,409 expanded FairFace reference paths. There are 48 primary syndrome cohorts, 5,760 primary hypotheses and 955 global-BH associations, including 482 with absolute g at least 0.8. Feature extraction is attributed to FaceKit.','',
        '## Validated correction families','',
        '| Analysis | Planned hypotheses | Supported | BH associations |','| --- | ---: | ---: | ---: |']
    for item in numbers['extension_families']:
        lines.append(f"| {item['analysis']} | {item['planned']} | {item['supported']} | {item['significant']} |")
    lines+=['','The adjusted and pairwise families retain prespecified target slots and unsupported cells. Their planned denominators are unchanged. Demographic A/B/C analyses retain separate broad/detailed correction scopes.','',
        '## Robustness and concordance','',
        'The all-eligible reference yields 971 associations, category-by-recorded-sex reference 945 and alternate-photo selection 996. The photo comparison has 918 shared, 37 primary-only and 78 alternate-only associations. One hundred random-photo replicates retain one image per patient. Reference specifications overlap and do not constitute independent replications.','',
        'Only descriptive concordance points are reported. All assessed concordance bootstrap intervals and uncertainty summaries are excluded because the method was biased. The original output is preserved privately for audit. Other primary/within-GMDB bootstrap intervals retain their validated definitions.','',
        '## Source roles and permissions','',
        'PDIDB contains 1,727 separately cataloged synthetic records, of which 1,563 pass QC. It contributes no real-patient inferential observations. No patient-level tables, landmarks or real GMDB photographs are included in manuscript exports; the API returns aggregate histograms. See publication/PI_REVIEW.md for IRB scope, image-reuse permissions, PDIDB lineage/redistribution and existing Git-history remediation questions.','',
        '## Executed verification','',
        f"The independent numerical validator passed {len(validation['components'])} components. Full primary replay, API/privacy checks, package checks and frontend build evidence are in this private folder and FINAL_ACTIVE_STUDY_AUDIT.json. Do not equate computational verification with clinical validation or publication approval."]
    (out/'PRISM_120_STUDY_REPORT.md').write_text('\n'.join(lines)+'\n')
    (out/'MANUSCRIPT_HANDOFF.md').write_text('''# Active manuscript handoff

`manuscript/main.tex` is the sole active paper. It is generated from `revision_results/rebuild_120_20261028/manuscript_120.tex`, with actual current clinical and demographic rows. The four figures and complete current supplementary tables include expanded-reference, demographic, within-GMDB, adjusted, age, correlation and photograph robustness results. Supplementary exports use reference terminology.

Concordance uncertainty is excluded; point estimates are descriptive. The main family and all prespecified extension denominators are unchanged. PRISM is presented as a statistical investigation using FaceKit, with no feature-extractor novelty or diagnostic-accuracy claim. All source roles and PI review items are stated in the manuscript and publication/PI_REVIEW.md. No live Google Doc edit has been made.
''')
    print('Current-only private study report and manuscript handoff written')


if __name__=='__main__':
    main()
