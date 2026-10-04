#!/usr/bin/env python3
import hashlib
import os
import platform
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
import numpy
import pandas
import scipy
import statsmodels

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "publication/provenance.yaml"


def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def generation_time():
    epoch = os.environ.get("SOURCE_DATE_EPOCH")
    return datetime.fromtimestamp(int(epoch), timezone.utc) if epoch else datetime.now(timezone.utc)


def yaml_value(value):
    return value if value else "null"


def main():
    sources = [
        "data/master_dataset_disease.csv", "data/master_dataset_healthy.csv",
        "revision_results/final_20260914/exact_feature_names.csv",
        "revision_results/final_20260914/age_3_9/significance_testing/stats_each_disease_vs_healthy.csv.gz",
        "revision_results/final_20260914/age_3_9_id_sensitivity/significance_testing/stats_each_disease_vs_healthy.csv.gz",
        "local_review/bibm2026_evidence/feature_region_family_mapping.csv",
        "local_review/bibm2026_evidence/main_candidate_review.csv",
        "publication/requirements.txt",
    ]
    scripts = {
        "audit": "publication/manuscript_results_audit.py",
        "supplementary_tables": "publication/scripts/generate_supplementary.py",
        "figures": "publication/scripts/generate_figures.py",
        "feature_correlation": "publication/scripts/feature_correlation.py",
        "candidate_evidence": "publication/scripts/generate_candidate_evidence.py",
        "exact_feature_definitions": "publication/scripts/generate_feature_definitions.py",
        "covariate_sensitivity": "publication/covariate_sensitivity/run_covariate_sensitivity.py",
        "provenance": "publication/scripts/generate_provenance.py",
        "reproducibility_verification": "publication/scripts/verify_reproducible_outputs.py",
        "publication_lint": "publication/publication_lint.py",
    }
    source_commit = os.environ.get("PUBLICATION_SOURCE_COMMIT")
    package_commit = os.environ.get("PUBLICATION_PACKAGE_COMMIT")
    lines = [
        "historical_facekit_revision: e474628f3cb17a83edad59ac6e1d9e5438753032",
        "corrected_prism_results_commit: a0952fdf8e4bc57304ed709e091c9e4eca704407",
        "bibm_evidence_commit: a20a1be75512b5d3167acd5edebce0b377b74358",
        f"publication_source_commit: {yaml_value(source_commit)}",
        f"publication_package_commit: {yaml_value(package_commit)}",
        f"generation_utc: \"{generation_time().isoformat()}\"",
        "publication_commit_fields_note: publication_source_commit identifies a committed source state when supplied; publication_package_commit remains null until a release commit or tag is created.",
        "authoritative_results_directory: revision_results/final_20260914",
        "feature_count: 125",
        "main_analysis_age_range: \"3 <= age < 10\"",
        "main_disease_patients: 1850",
        "reference_image_paths: 942",
        "eligible_syndrome_groups: 49",
        "main_tests: 6125",
        "source_sha256:",
    ]
    for source in sources: lines.append(f"  {source}: {digest(ROOT / source)}")
    lines.append("generation_scripts:")
    for name, script in scripts.items():
        path = ROOT / script
        lines.append(f"  {name}:")
        lines.append(f"    path: {script}")
        lines.append(f"    sha256: {digest(path) if path.exists() else 'MISSING'}")
    lines += ["software:", f"  python: \"{platform.python_version()}\"", f"  numpy: \"{numpy.__version__}\"", f"  pandas: \"{pandas.__version__}\"",
              f"  scipy: \"{scipy.__version__}\"", f"  matplotlib: \"{matplotlib.__version__}\"", f"  statsmodels: \"{statsmodels.__version__}\"",
              "known_limitations:",
              "  - Disease images come from GMDB and reference images come from FairFace; source and acquisition are confounded with disease status.",
              "  - Distinct reference paths do not establish distinct reference participants.",
              "  - Recorded age components are used as supplied and do not independently verify age semantics.",
              "  - Small demographic subgroups can yield unstable matched and interaction estimates.",
              "  - Historical extraction provenance identifies source revision but not a per-image execution log or model checksum.",
              "  - Literature classification and clinical interpretation remain unresolved unless explicitly reviewed.",
              "historical_measurement_policy: Reuse corrected 125-feature outputs; do not rerun landmark extraction or substitute newer 120-feature outputs.", ""]
    OUTPUT.write_text("\n".join(lines))
    print(f"Wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__": main()
