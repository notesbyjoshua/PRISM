# Data Privacy and Repository Audit

Audit date: 2026-10-02

> **PUBLICATION BLOCKER — HUMAN REVIEW REQUIRED.** Committed source tables contain patient- or image-level
> records, identifiers or locators, demographic attributes, diagnoses, genetic/phenotypic metadata, and
> derived facial measurements. Because the repository is already public, treat the present exposure as requiring
> urgent review and do not create a new public release until authorized reviewers confirm what may be redistributed.

This report records repository-relative paths and schema only. It intentionally includes no record values,
patient identifiers, image filenames, image URLs, or machine-local absolute paths.

## Scope and method

The reproducible audit inventoried 138 current or historically committed tabular/database artifacts in
`data/`, `revision_results/`, `publication/`, `manuscript/`, `stats_results/`, and `pds_results_v2/`, plus the
ignored local application database if present. It examined Git tracking/history, file type, approximate row count,
and column names only. The complete file-by-file inventory is `publication/data_privacy_audit.csv`; rerun it with
`python3 publication/scripts/data_privacy_audit.py`.

This is a technical screening, not a legal determination. The classifications are conservative and require
confirmation from people authorized to interpret the study consent, data-use agreement, and institutional policy.

## Blocking findings

- 5 currently tracked artifacts received at least one manual-review flag.
- `data/master_dataset_disease.csv` and `data/master_dataset_healthy.csv` are direct inputs to the frozen analysis
  pipeline, but direct inputs do not need to be publicly redistributed for a manuscript repository to publish
  code, aggregate outputs, checksums, and access instructions.
- The disease cohort appears to derive from the GestaltMatcher Database (GMDB). GMDB describes researcher access
  as proposal-based and subject to Scientific Advisory Board approval, while its study information describes
  pseudonymous sharing under controlled access. That is incompatible with assuming these rows are freely
  redistributable without written authorization.
- The healthy cohort appears to derive from FairFace. Its public repository states a CC BY 4.0 license, but the
  institution must still confirm that this project's redistributed metadata, locators, and derived measurements
  comply with the relevant dataset version, attribution requirements, privacy policy, and journal/repository rules.
- Removing a risky file in a new commit would not remove it from Git history. Any approved remediation plan must
  cover historical objects and existing remote clones or forks.

## High-risk current artifacts

| Path | Approx. rows | Risk-relevant schema (column names only) | Assessment | Recommended status |
|---|---:|---|---|---|
| `data/master_dataset_disease.csv` | 10944 | `disease`, `image_id`, `pose_yaw`, `pose_pitch`, `pose_roll`, `patient_id`, `patient_group`, `internal_syndrome_id`, `internal_syndrome_name`, `gene_names`, `gene_entrez_ids`, `hgvs`, `present_features`, `absent_features`, `ethnicity_category`, `ethnicity_sub_category`, `age_year`, `age_month`, `age_note`, `gender`, `syndrome_id`, `syndrome_name`, `hpo_terms`, `image_dataset_image_id`, `image_dataset_patient_id`, `image_dataset_patient_group`, `image_dataset_internal_syndrome_id`, `image_dataset_internal_syndrome_name`, `image_dataset_gene_names`, `image_dataset_gene_entrez_ids`, `image_dataset_hgvs`, `image_dataset_test_type`, `patient_name`, `image_dataset_person_type`, `image_dataset_distinctiveness`, `image_dataset_present_features`, `image_dataset_absent_features`, `image_dataset_ethnicity_category`, `image_dataset_ethnicity_sub_category`, `image_dataset_age_year`, `image_dataset_age_month`, `image_dataset_age_note`, `image_dataset_gender`, `image_dataset_syndrome_id`, `image_dataset_syndrome_name`, `image_dataset_train`, `image_dataset_OMIM`, `image_dataset_hpo_terms`, `image_dataset_split`, `image_dataset_row_index`, `facekit_image_id`, `facekit_relative_path`, `image_exported`, `image_export_error`; plus 135 other columns (complete schema in CSV) | Patient/image-level disease metadata, identifiers, demographics, diagnosis, genetic/phenotypic metadata, and locator fields. | `REMOVE_FROM_PUBLIC_RELEASE_IF_UNAPPROVED` |
| `data/metadata_with_images.csv` | 10982 | `image_id`, `patient_id`, `patient_group`, `internal_syndrome_id`, `internal_syndrome_name`, `gene_names`, `gene_entrez_ids`, `hgvs`, `present_features`, `absent_features`, `ethnicity_category`, `ethnicity_sub_category`, `age_year`, `age_month`, `age_note`, `gender`, `syndrome_id`, `syndrome_name`, `hpo_terms`, `image_dataset_image_id`, `image_dataset_patient_id`, `image_dataset_patient_group`, `image_dataset_internal_syndrome_id`, `image_dataset_internal_syndrome_name`, `image_dataset_gene_names`, `image_dataset_gene_entrez_ids`, `image_dataset_hgvs`, `image_dataset_test_type`, `patient_name`, `image_dataset_person_type`, `image_dataset_distinctiveness`, `image_dataset_present_features`, `image_dataset_absent_features`, `image_dataset_ethnicity_category`, `image_dataset_ethnicity_sub_category`, `image_dataset_age_year`, `image_dataset_age_month`, `image_dataset_age_note`, `image_dataset_gender`, `image_dataset_syndrome_id`, `image_dataset_syndrome_name`, `image_dataset_train`, `image_dataset_OMIM`, `image_dataset_hpo_terms`, `image_dataset_split`, `image_dataset_row_index`, `facekit_image_id`, `facekit_relative_path`, `image_exported`, `image_export_error`; plus 9 other columns (complete schema in CSV) | Patient/image-level disease metadata with identifier and locator fields. | `REMOVE_FROM_PUBLIC_RELEASE_IF_UNAPPROVED` |
| `data/phenotypes_all_disease.csv` | 10942 | `disease`, `image_id`, `pose_yaw`, `pose_pitch`, `pose_roll`; plus 126 other columns (complete schema in CSV) | Per-image disease labels and facial measurements keyed by an image identifier. | `REMOVE_FROM_PUBLIC_RELEASE_IF_UNAPPROVED` |
| `data/master_dataset_healthy.csv` | 1000 | `file`, `age`, `gender`, `race`, `diagnosis`, `image_id`, `pose_yaw`, `pose_pitch`, `pose_roll`; plus 130 other columns (complete schema in CSV) | Per-image healthy-reference metadata including age category, gender, race, pose, and locator fields. | `REVIEW_LICENSE` |
| `data/phenotypes_all_healthy.csv` | 1000 | `disease`, `image_id`, `pose_yaw`, `pose_pitch`, `pose_roll`; plus 126 other columns (complete schema in CSV) | Per-image healthy-reference measurements keyed by an image identifier. | `REVIEW_LICENSE` |

The high-risk tables should not be copied into a public release merely because they are needed for internal
re-execution. For public verification, publish aggregate feature tables, model/statistical outputs, exact scripts,
environment details, checksums, and a controlled-data access statement. Full reruns can remain conditional on
authorized access to the underlying cohorts.

## Aggregate and publication artifacts

Most files in `stats_results/`, `pds_results_v2/`, `revision_results/`, and `publication/` are classified as
aggregate results, summaries, configuration, or provenance. They are substantially lower risk, but should still
receive a final disclosure-control review for small cells, rare-disease group sizes, free-text columns, and any
unexpected identifier or locator columns before release. The CSV inventory makes that review file-specific.

## Git history review

Git history was searched across all refs by path. The inventory includes the number of commits touching each
artifact and the earliest/latest abbreviated commit hashes. Historical-only data artifacts require review even
when absent from the current checkout.

| Historical-only path | Commits | Latest | Earliest |
|---|---:|---|---|
| `data/master_dataset.csv` | 2 | `bc140f464c33` | `27ba5ce9cdc3` |

If a public repository is approved only after sensitive artifacts are removed, use a coordinated history-rewrite
plan, rotate or retire the affected remote, and obtain explicit confirmation before force-pushing. This audit does
not perform any deletion or history rewrite.

## Data-use and licensing sources

- [GMDB data-access FAQ](https://db.gestaltmatcher.org/faq/topics/data-access-ownership/articles/how-can-i-request-the-data-included-in-gmdb?locale=en)
  describes account, proposal, and Scientific Advisory Board approval requirements.
- [GMDB study information](https://db.gestaltmatcher.org/files/Studyinformation_english.pdf) describes pseudonymous
  sharing with controlled access and identifies privacy risks associated with facial data.
- [GestaltMatcher publication data-availability statement](https://www.nature.com/articles/s41588-021-01010-x)
  describes GMDB as controlled access for approved noncommercial research and educational use.
- [FairFace repository README](https://github.com/joojs/fairface/blob/master/README.md) identifies CC BY 4.0;
  confirm that this is the exact source/version used and preserve required attribution.

## Recommended release boundary

Until authorization is documented, a public release should contain code and reviewed aggregate artifacts only.
Keep patient/image-level source tables and derived per-image tables outside the public repository; describe how
qualified researchers can request the controlled data. Before any publication transition, verify all branches,
tags, release assets, caches, CI artifacts, and remotes—not only the current branch.

## Decisions requiring mentor / institution / data-provider confirmation

1. Does the executed GMDB agreement and participant consent permit public redistribution of derived patient-level
   metadata, GMDB patient/image identifiers, genetic or phenotype fields, or facial measurements in this repository?
2. Are image filenames or file paths controlled data even when the underlying images are not committed?
3. Must all GMDB-derived row-level artifacts be removed from Git history and the existing remote, and who is
   authorized to approve and coordinate that rewrite?
4. May FairFace-derived per-image phenotype rows be redistributed, which FairFace release/license governed this
   cohort, and do the attribution and redistribution method satisfy those terms and institutional policy?
5. Which aggregate outputs may be archived publicly, and what minimum-cell or rare-disease disclosure-control
   threshold does the institution require?
6. Does any applicable data-use agreement require repository access restrictions while remediation is evaluated?
7. Should reproducibility be delivered as code plus checksums and controlled-access instructions rather than bundled
   source data, and what exact data-access language should appear in the manuscript and repository README?
8. Are any other branches, tags, forks, release assets, CI artifacts, or remote caches in scope for remediation?
