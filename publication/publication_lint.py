#!/usr/bin/env python3
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "manuscript/main.tex"
REPORT = ROOT / "publication/PUBLICATION_LINT_REPORT.md"


def issue(rows, severity, check, location, detail): rows.append({"severity": severity, "check": check, "location": location, "detail": detail})


def main():
    rows = []
    tex = MANUSCRIPT.read_text()
    for pattern, name in [(r"\?\?", "Unresolved LaTeX reference"), (r"\bTODO\b", "TODO marker"), (r"\bTBD\b", "TBD marker"),
                          (r"\\cite\{(?:placeholder|citation|ref)\}", "Placeholder citation")]:
        for match in re.finditer(pattern, tex, flags=re.I):
            line = tex.count("\n", 0, match.start()) + 1
            issue(rows, "ERROR", name, f"manuscript/main.tex:{line}", match.group(0))

    figure_blocks = re.findall(r"\\begin\{figure\}.*?\\end\{figure\}", tex, flags=re.S)
    labels = []
    for i, block in enumerate(figure_blocks, 1):
        label = re.search(r"\\label\{([^}]+)\}", block)
        path = re.search(r"\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}", block)
        if not label: issue(rows, "ERROR", "Missing figure label", f"manuscript/main.tex figure block {i}", "Figure has no \\label")
        else: labels.append(label.group(1))
        if not path: issue(rows, "ERROR", "Missing figure path", f"manuscript/main.tex figure block {i}", "Figure has no \\includegraphics")
        elif not (MANUSCRIPT.parent / path.group(1)).exists(): issue(rows, "ERROR", "Broken figure path", f"manuscript/main.tex figure block {i}", path.group(1))
    for label in labels:
        references = len(re.findall(rf"\\ref\{{{re.escape(label)}\}}", tex))
        if references < 1: issue(rows, "WARNING", "Figure never cited", "manuscript/main.tex", label)

    citations = set(re.findall(r"\\cite(?:t|p)?\{([^}]+)\}", tex))
    citations = {key.strip() for group in citations for key in group.split(",") if key.strip()}
    bib_text = (ROOT / "manuscript/references.bib").read_text()
    bib_keys = set(re.findall(r"@\w+\{\s*([^,\s]+)", bib_text))
    for key in sorted(citations - bib_keys): issue(rows, "ERROR", "Citation missing from bibliography", "manuscript/main.tex", key)
    for key in sorted(bib_keys - citations): issue(rows, "WARNING", "Bibliography item never cited", "manuscript/references.bib", key)

    author = re.search(r"\\author\{([^}]*)\}", tex, flags=re.S)
    if not author or "blind" not in author.group(1).lower(): issue(rows, "ERROR", "Blinded author check", "manuscript/main.tex", "Author field is absent or not blinded")
    for match in re.finditer(r"\bnovel phenotype\b", tex, flags=re.I):
        issue(rows, "ERROR", "Unsupported novelty claim", f"manuscript/main.tex:{tex.count(chr(10), 0, match.start()) + 1}", match.group(0))
    for match in re.finditer(r"[^.]*sensitivity[^.]*replication[^.]*\.", tex, flags=re.I):
        sentence = match.group(0)
        if not re.search(r"not (?:an? )?independent replication", sentence, flags=re.I):
            issue(rows, "ERROR", "Sensitivity described as replication", f"manuscript/main.tex:{tex.count(chr(10), 0, match.start()) + 1}", sentence.strip())

    expected = [(r"([\d,]+) disease patients", 1850, "disease patient count"), (r"([\d,]+) distinct reference image paths", 942, "reference path count"),
                (r"([\d,]+)[- ]feature", 125, "feature count"), (r"([\d,]+) syndrome--feature", 6125, "main test count"),
                (r"([\d,]+) globally significant", 811, "global significant count")]
    for pattern, correct, name in expected:
        for match in re.finditer(pattern, tex, flags=re.I):
            value = int(match.group(1).replace(",", ""))
            if value != correct: issue(rows, "ERROR", f"Inconsistent {name}", f"manuscript/main.tex:{tex.count(chr(10), 0, match.start()) + 1}", f"found {value:,}; expected {correct:,}")

    text_files = [p for folder in [ROOT / "publication", ROOT / "manuscript"] for p in folder.rglob("*") if p.is_file() and p.suffix in {".md", ".tex", ".yaml", ".sh"}]
    excluded = {REPORT, ROOT / "publication/REVIEWER_METHODS_AUDIT.md"}
    stale = re.compile(r"44,500|9,002|356 syndromes|579 syndromes")
    newer = re.compile(r"120-feature", re.I)
    for path in text_files:
        if path in excluded: continue
        for number, line in enumerate(path.read_text(errors="replace").splitlines(), 1):
            if stale.search(line): issue(rows, "WARNING", "Legacy manuscript-era count", f"{path.relative_to(ROOT)}:{number}", line.strip())
            if newer.search(line) and not re.search(r"do not|out of scope|not use|substitute|\bno 120-feature", line, flags=re.I):
                issue(rows, "ERROR", "Unqualified newer 120-feature reference", f"{path.relative_to(ROOT)}:{number}", line.strip())

    errors = sum(row["severity"] == "ERROR" for row in rows)
    warnings = sum(row["severity"] == "WARNING" for row in rows)
    lines = ["# Publication lint report", "", f"**Result:** {errors} errors, {warnings} warnings.", ""]
    if rows:
        lines += ["| Severity | Check | Location | Detail |", "|---|---|---|---|"]
        for row in rows: lines.append(f"| {row['severity']} | {row['check']} | `{row['location']}` | {row['detail'].replace('|', '/')} |")
    else:
        lines += ["No lint findings were detected in the current manuscript and publication package.", ""]
    lines += ["", "The institutional and author-confirmation brackets in the manuscript are intentional submission blockers. They are not TODO/TBD markers and require human confirmation before submission.", ""]
    REPORT.write_text("\n".join(lines))
    print(f"Wrote {REPORT.relative_to(ROOT)}: {errors} errors, {warnings} warnings")
    return 1 if errors else 0


if __name__ == "__main__": raise SystemExit(main())
