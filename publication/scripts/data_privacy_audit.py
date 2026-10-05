#!/usr/bin/env python3
import csv
import gzip
import json
import sqlite3
import subprocess
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CSV_OUT = ROOT / "publication/data_privacy_audit.csv"
MD_OUT = ROOT / "publication/DATA_PRIVACY_AUDIT.md"
DATA_EXTENSIONS = {".csv", ".tsv", ".json", ".jsonl", ".parquet", ".feather", ".xlsx", ".xls", ".db", ".sqlite"}
SCOPES = ("data/", "revision_results/", "publication/", "manuscript/", "stats_results/", "pds_results_v2/")
EXTRA_RELEVANT_SCOPES = ("local_review/bibm2026_evidence/",)
DIRECT_INPUTS = {"data/master_dataset_disease.csv", "data/master_dataset_healthy.csv"}
HIGH_RISK = {
    "data/master_dataset_disease.csv", "data/metadata_with_images.csv", "data/phenotypes_all_disease.csv",
    "data/master_dataset_healthy.csv", "data/phenotypes_all_healthy.csv", "data/master_dataset.csv",
}
FIELDS = [
    "repository_path", "current_exists", "tracked_by_git", "ignored_by_git", "historically_committed",
    "history_commit_count", "latest_commit", "earliest_commit", "file_type", "approximate_rows", "column_names",
    "record_granularity", "sensitive_field_categories", "required_for_publication", "aggregate_replacement_available",
    "recommended_status", "risk_flags", "notes",
]


def git(*args, check=True):
    result = subprocess.run(["git", *args], cwd=ROOT, text=True, capture_output=True, check=check)
    return result.stdout


def is_data_path(path):
    lower = path.lower()
    return lower.endswith(".csv.gz") or Path(lower).suffix in DATA_EXTENSIONS


def under_scope(path):
    return path.startswith(SCOPES)


def bool_text(value):
    return "yes" if value else "no"


def approximate_rows(path):
    opener = gzip.open if path.name.endswith(".gz") else open
    mode = "rt" if path.name.endswith(".gz") else "rb"
    try:
        if mode == "rb":
            with opener(path, mode) as handle:
                return max(sum(chunk.count(b"\n") for chunk in iter(lambda: handle.read(1024 * 1024), b"")) - 1, 0)
        with opener(path, mode, errors="replace") as handle:
            return max(sum(1 for _ in handle) - 1, 0)
    except (OSError, EOFError):
        return ""


def delimited_columns(path, delimiter):
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8-sig", errors="replace", newline="") as handle:
        return next(csv.reader(handle, delimiter=delimiter), [])


def json_columns(path):
    try:
        with open(path, encoding="utf-8") as handle:
            value = json.load(handle)
        if isinstance(value, dict):
            return list(value)
        if isinstance(value, list) and value and isinstance(value[0], dict):
            return list(value[0])
    except (OSError, json.JSONDecodeError):
        pass
    return []


def sqlite_columns(path):
    names = []
    try:
        connection = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
        tables = connection.execute("SELECT name FROM sqlite_master WHERE type='table' ORDER BY name").fetchall()
        for (table,) in tables:
            escaped = table.replace("'", "''")
            for row in connection.execute(f"PRAGMA table_info('{escaped}')"):
                names.append(f"{table}.{row[1]}")
        connection.close()
    except sqlite3.Error:
        pass
    return names


def inspect_file(path):
    lower = path.name.lower()
    if lower.endswith(".csv") or lower.endswith(".csv.gz"):
        return "CSV" if not lower.endswith(".gz") else "CSV.GZ", approximate_rows(path), delimited_columns(path, ",")
    if lower.endswith(".tsv"):
        return "TSV", approximate_rows(path), delimited_columns(path, "\t")
    if lower.endswith((".db", ".sqlite")):
        return "SQLite", "", sqlite_columns(path)
    if lower.endswith((".json", ".jsonl")):
        return lower.rsplit(".", 1)[-1].upper(), "", json_columns(path)
    return path.suffix.lstrip(".").upper(), "", []


def sensitive_categories(columns):
    categories = set()
    for column in columns:
        name = column.lower()
        if any(term in name for term in ("patient_id", "patient_name", "subject_id", "person_id")):
            categories.add("patient identifier")
        if any(term in name for term in ("image_id", "imageid")):
            categories.add("image identifier")
        if any(term in name for term in ("file", "path", "url")):
            categories.add("image/file locator")
        if "age" in name or "birth" in name:
            categories.add("age")
        if any(term in name for term in ("gender", "sex")):
            categories.add("sex/gender")
        if any(term in name for term in ("race", "ethnicity", "ancestry")):
            categories.add("race/ethnicity")
        if any(term in name for term in ("disease", "syndrome", "diagnos", "condition")):
            categories.add("diagnosis")
        if any(term in name for term in ("gene", "hgvs", "variant", "omim")):
            categories.add("genetic/rare-disease metadata")
        if any(term in name for term in ("accession", "database_id", "omim", "entrez")):
            categories.add("accession/database identifier")
        if any(term in name for term in ("hpo", "phenotype", "clinical", "note", "description", "feature")):
            categories.add("clinical/free-text metadata")
        if any(term in name for term in ("institution", "site", "hospital", "country")):
            categories.add("institution/geography")
        if "url" in name:
            categories.add("URL")
        if any(term in name for term in ("date", "timestamp", "created_at", "updated_at")):
            categories.add("date/time")
    return sorted(categories)


def classify(path, columns):
    lower_columns = {column.lower() for column in columns}
    source = path in HIGH_RISK
    has_patient_id = any("patient_id" in column or "patient_name" in column for column in lower_columns)
    has_image_id = any("image_id" in column or "imageid" in column for column in lower_columns)
    locator_names = {"file", "filename", "image_file", "image_path", "relative_path", "facekit_relative_path", "image_url"}
    has_locator = source and any(column in locator_names or column.endswith("_relative_path") for column in lower_columns)
    if source:
        granularity = "patient/image level"
    elif has_patient_id:
        granularity = "patient level"
    elif has_image_id or has_locator:
        granularity = "image/file level"
    else:
        granularity = "aggregate/configuration"

    if path in DIRECT_INPUTS:
        required = "yes — direct analysis input"
    elif path.startswith(("revision_results/final_", "local_review/bibm2026_evidence/")):
        required = "yes — selected evidence input"
    elif path.startswith("publication/"):
        required = "no — publication output or support file"
    elif path.startswith("manuscript/"):
        required = "no — manuscript copy or support file"
    else:
        required = "no — legacy/intermediate analysis output"

    categories = sensitive_categories(columns)
    flags = []
    if granularity in {"patient level", "patient/image level"}:
        flags.append("REVIEW_PATIENT_LEVEL_METADATA")
    if has_locator:
        flags.append("REVIEW_IMAGE_PATH")
    if path.startswith("data/master_dataset_disease") or path.startswith("data/metadata_with_images") or path.startswith("data/phenotypes_all_disease") or path == "data/master_dataset.csv":
        flags.append("REVIEW_CONTROLLED_ACCESS_TERMS")
    if path.startswith(("data/master_dataset_healthy", "data/phenotypes_all_healthy")):
        flags.append("REVIEW_FAIRFACE_LICENSE")

    if "REVIEW_CONTROLLED_ACCESS_TERMS" in flags:
        status = "REMOVE_FROM_PUBLIC_RELEASE_IF_UNAPPROVED"
    elif "REVIEW_FAIRFACE_LICENSE" in flags:
        status = "REVIEW_LICENSE"
    elif granularity != "aggregate/configuration":
        status = "REVIEW_PATIENT_LEVEL_METADATA"
    elif granularity == "aggregate/configuration":
        status = "LIKELY_SAFE_DERIVED"
    else:
        status = "UNKNOWN_NEEDS_HUMAN_REVIEW"
    replacement = "yes for claim verification; no for full rerun" if source else "not needed — already aggregate/configuration"
    notes = "Column names only were inspected; no record values are included in this audit."
    return granularity, categories, required, replacement, status, flags, notes


def history(path):
    commits = [line for line in git("log", "--all", "--format=%H", "--", path).splitlines() if line]
    return commits


def collect_paths(tracked):
    historical = git("log", "--all", "--name-only", "--pretty=format:").splitlines()
    candidates = {path for path in tracked | set(historical) if path and is_data_path(path) and under_scope(path)}
    candidates.update(path for path in tracked if is_data_path(path) and path.startswith(EXTRA_RELEVANT_SCOPES))
    for path in ("backend/data/prism.db",):
        if (ROOT / path).exists():
            candidates.add(path)
    return sorted(candidates)


def make_row(path, tracked):
    full_path = ROOT / path
    exists = full_path.is_file()
    file_type, rows, columns = inspect_file(full_path) if exists else (Path(path).suffix.lstrip(".").upper(), "", [])
    commits = history(path)
    ignored = subprocess.run(["git", "check-ignore", "-q", "--", path], cwd=ROOT).returncode == 0
    granularity, categories, required, replacement, status, flags, notes = classify(path, columns)
    return {
        "repository_path": path,
        "current_exists": bool_text(exists),
        "tracked_by_git": bool_text(path in tracked),
        "ignored_by_git": bool_text(ignored),
        "historically_committed": bool_text(bool(commits)),
        "history_commit_count": len(commits),
        "latest_commit": commits[0][:12] if commits else "",
        "earliest_commit": commits[-1][:12] if commits else "",
        "file_type": file_type,
        "approximate_rows": rows,
        "column_names": "; ".join(columns),
        "record_granularity": granularity,
        "sensitive_field_categories": "; ".join(categories),
        "required_for_publication": required,
        "aggregate_replacement_available": replacement,
        "recommended_status": status,
        "risk_flags": "; ".join(flags),
        "notes": notes,
    }


def write_csv(rows):
    with open(CSV_OUT, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def key_row(rows, path):
    return next((row for row in rows if row["repository_path"] == path), None)


def compact_columns(row):
    columns = row["column_names"].split("; ") if row else []
    terms = (
        "patient", "image_id", "disease", "syndrome", "gene", "hgvs", "hpo", "age", "gender", "race",
        "ethnicity", "diagnos", "present_features", "absent_features", "file", "path", "pose",
    )
    relevant = [column for column in columns if any(term in column.lower() for term in terms)]
    omitted = len(columns) - len(relevant)
    summary = ", ".join(f"`{column}`" for column in relevant)
    return f"{summary}; plus {omitted} other columns (complete schema in CSV)" if omitted else summary


def write_markdown(rows):
    risky = [row for row in rows if row["risk_flags"]]
    tracked_risky = [row for row in risky if row["tracked_by_git"] == "yes"]
    historical_only = [row for row in rows if row["current_exists"] == "no" and row["historically_committed"] == "yes"]
    disease = key_row(rows, "data/master_dataset_disease.csv")
    metadata = key_row(rows, "data/metadata_with_images.csv")
    disease_phenotypes = key_row(rows, "data/phenotypes_all_disease.csv")
    healthy = key_row(rows, "data/master_dataset_healthy.csv")
    healthy_phenotypes = key_row(rows, "data/phenotypes_all_healthy.csv")
    lines = [
        "# Data Privacy and Repository Audit",
        "",
        f"Audit date: {date.today().isoformat()}",
        "",
        "> **PUBLICATION BLOCKER — HUMAN_INPUT_REQUIRED.** Committed source tables contain patient- or image-level",
        "> records, identifiers or locators, demographic attributes, diagnoses, genetic/phenotypic metadata, and",
        "> derived facial measurements. The repository is currently private. Public redistribution remains pending explicit",
        "> mentor, institution, and data-provider approval; private status does not establish that the files may be released.",
        "> This audit performs no deletion, history rewrite, access change, or release action.",
        "",
        "This report records repository-relative paths and schema only. It intentionally includes no record values,",
        "patient identifiers, image filenames, image URLs, or machine-local absolute paths.",
        "",
        "## Scope and method",
        "",
        f"The reproducible audit inventoried {len(rows)} current or historically committed tabular/database artifacts in",
        "`data/`, `revision_results/`, `publication/`, `manuscript/`, `stats_results/`, and `pds_results_v2/`, plus the",
        "ignored local application database if present. It examined Git tracking/history, file type, approximate row count,",
        "and column names only. The complete file-by-file inventory is `publication/data_privacy_audit.csv`; rerun it with",
        "`python3 publication/scripts/data_privacy_audit.py`.",
        "",
        "This is a technical screening, not a legal determination. The classifications are conservative and require",
        "confirmation from people authorized to interpret the study consent, data-use agreement, and institutional policy.",
        "",
        "## Blocking findings",
        "",
        f"- {len(tracked_risky)} currently tracked artifacts received at least one manual-review flag.",
        "- `data/master_dataset_disease.csv` and `data/master_dataset_healthy.csv` are direct inputs to the frozen analysis",
        "  pipeline, but direct inputs do not need to be publicly redistributed for a manuscript repository to publish",
        "  code, aggregate outputs, checksums, and access instructions.",
        "- The disease cohort appears to derive from the GestaltMatcher Database (GMDB). GMDB describes researcher access",
        "  as proposal-based and subject to Scientific Advisory Board approval, while its study information describes",
        "  pseudonymous sharing under controlled access. That is incompatible with assuming these rows are freely",
        "  redistributable without written authorization.",
        "- The healthy cohort appears to derive from FairFace. Its public repository states a CC BY 4.0 license, but the",
        "  institution must still confirm that this project's redistributed metadata, locators, and derived measurements",
        "  comply with the relevant dataset version, attribution requirements, privacy policy, and journal/repository rules.",
        "- Removing a risky file in a new commit would not remove it from Git history. Any approved remediation plan must",
        "  cover historical objects and existing remote clones or forks.",
        "",
        "## High-risk current artifacts",
        "",
        "| Path | Approx. rows | Risk-relevant schema (column names only) | Assessment | Recommended status |",
        "|---|---:|---|---|---|",
    ]
    for row, assessment in (
        (disease, "Patient/image-level disease metadata, identifiers, demographics, diagnosis, genetic/phenotypic metadata, and locator fields."),
        (metadata, "Patient/image-level disease metadata with identifier and locator fields."),
        (disease_phenotypes, "Per-image disease labels and facial measurements keyed by an image identifier."),
        (healthy, "Per-image healthy-reference metadata including age category, gender, race, pose, and locator fields."),
        (healthy_phenotypes, "Per-image healthy-reference measurements keyed by an image identifier."),
    ):
        if row:
            lines.append(
                f"| `{row['repository_path']}` | {row['approximate_rows']} | {compact_columns(row)} | {assessment} | "
                f"`{row['recommended_status']}` |"
            )
    lines += [
        "",
        "The high-risk tables should not be copied into a public release merely because they are needed for internal",
        "re-execution. For public verification, publish aggregate feature tables, model/statistical outputs, exact scripts,",
        "environment details, checksums, and a controlled-data access statement. Full reruns can remain conditional on",
        "authorized access to the underlying cohorts.",
        "",
        "## Aggregate and publication artifacts",
        "",
        "Most files in `stats_results/`, `pds_results_v2/`, `revision_results/`, and `publication/` are classified as",
        "aggregate results, summaries, configuration, or provenance. They are substantially lower risk, but should still",
        "receive a final disclosure-control review for small cells, rare-disease group sizes, free-text columns, and any",
        "unexpected identifier or locator columns before release. The CSV inventory makes that review file-specific.",
        "",
        "## Git history review",
        "",
        "Git history was searched across all refs by path. The inventory includes the number of commits touching each",
        "artifact and the earliest/latest abbreviated commit hashes. Historical-only data artifacts require review even",
        "when absent from the current checkout.",
        "",
        "| Historical-only path | Commits | Latest | Earliest |",
        "|---|---:|---|---|",
    ]
    if historical_only:
        for row in historical_only:
            lines.append(
                f"| `{row['repository_path']}` | {row['history_commit_count']} | `{row['latest_commit']}` | "
                f"`{row['earliest_commit']}` |"
            )
    else:
        lines.append("| None detected in audited extensions/scopes | 0 | — | — |")
    lines += [
        "",
        "If authorized reviewers later require sensitive artifacts to be removed before public release, they must define and",
        "approve a coordinated remediation plan covering current files, Git history, remotes, and any existing clones or",
        "forks. This audit does not recommend or perform deletion, force-pushing, or history rewriting.",
        "",
        "## Data-use and licensing sources",
        "",
        "- [GMDB data-access FAQ](https://db.gestaltmatcher.org/faq/topics/data-access-ownership/articles/how-can-i-request-the-data-included-in-gmdb?locale=en)",
        "  describes account, proposal, and Scientific Advisory Board approval requirements.",
        "- [GMDB study information](https://db.gestaltmatcher.org/files/Studyinformation_english.pdf) describes pseudonymous",
        "  sharing with controlled access and identifies privacy risks associated with facial data.",
        "- [GestaltMatcher publication data-availability statement](https://www.nature.com/articles/s41588-021-01010-x)",
        "  describes GMDB as controlled access for approved noncommercial research and educational use.",
        "- [FairFace repository README](https://github.com/joojs/fairface/blob/master/README.md) identifies CC BY 4.0;",
        "  confirm that this is the exact source/version used and preserve required attribution.",
        "",
        "## Recommended release boundary",
        "",
        "Keep the repository private while mentor, institution, and data-provider review is pending. No source-level,",
        "patient/image-level, per-image derived, aggregate, code, or history content should be represented as approved for",
        "public redistribution until that review defines the release boundary. If a public release is authorized, verify all",
        "branches, tags, release assets, caches, CI artifacts, and remotes—not only the current branch—and provide approved",
        "controlled-data access language for materials that cannot be redistributed.",
        "",
        "## Decisions requiring mentor / institution / data-provider confirmation",
        "",
        "1. Does the executed GMDB agreement and participant consent permit public redistribution of derived patient-level",
        "   metadata, GMDB patient/image identifiers, genetic or phenotype fields, or facial measurements in this repository?",
        "2. Are image filenames or file paths controlled data even when the underlying images are not committed?",
        "3. Must all GMDB-derived row-level artifacts be removed from Git history and the existing remote, and who is",
        "   authorized to approve and coordinate that rewrite?",
        "4. May FairFace-derived per-image phenotype rows be redistributed, which FairFace release/license governed this",
        "   cohort, and do the attribution and redistribution method satisfy those terms and institutional policy?",
        "5. Which aggregate outputs may be archived publicly, and what minimum-cell or rare-disease disclosure-control",
        "   threshold does the institution require?",
        "6. Does any applicable data-use agreement require repository access restrictions while remediation is evaluated?",
        "7. Should reproducibility be delivered as code plus checksums and controlled-access instructions rather than bundled",
        "   source data, and what exact data-access language should appear in the manuscript and repository README?",
        "8. Are any other branches, tags, forks, release assets, CI artifacts, or remote caches in scope for remediation?",
    ]
    MD_OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    tracked = set(git("ls-files").splitlines())
    rows = [make_row(path, tracked) for path in collect_paths(tracked)]
    write_csv(rows)
    write_markdown(rows)
    print(f"Wrote {len(rows)} audit rows to {CSV_OUT.relative_to(ROOT)} and {MD_OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
