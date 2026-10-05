#!/usr/bin/env python3
import re
import shutil
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "manuscript/main.tex"
BIBLIOGRAPHY = ROOT / "manuscript/references.bib"
REPORT = ROOT / "publication/PUBLICATION_LINT_REPORT.md"
PLACEHOLDER = re.compile(
    r"\[(?:\\detokenize\{)?(HUMAN_INPUT_REQUIRED|AUTHOR INPUT REQUIRED|AUTHOR CONFIRMATION REQUIRED|AUTHOR AND CLINICAL REVIEW REQUIRED|"
    r"INSTITUTIONAL CONFIRMATION REQUIRED|LITERATURE REVIEW REQUIRED)\}?:\s*([^]]+)\]",
    flags=re.I,
)


def issue(rows, severity, check, location, detail):
    rows.append({"severity": severity, "check": check, "location": location, "detail": detail})


def line_number(text, position):
    return text.count("\n", 0, position) + 1


def section_text(tex, name):
    match = re.search(rf"\\section\{{{re.escape(name)}\}}(.*?)(?=\\section\{{|\\end\{{document\}})", tex, flags=re.S | re.I)
    return match.group(1).strip() if match else ""


def latex_check(rows):
    latexmk, pdflatex = shutil.which("latexmk"), shutil.which("pdflatex")
    if not latexmk and not pdflatex: return "LATEX_COMPILE_NOT_CHECKED"
    with tempfile.TemporaryDirectory(prefix="prism-latex-") as temp:
        workspace = Path(temp) / "manuscript"
        shutil.copytree(MANUSCRIPT.parent, workspace)
        command = ([latexmk, "-pdf", "-interaction=nonstopmode", "-halt-on-error", "main.tex"] if latexmk else
                   [pdflatex, "-interaction=nonstopmode", "-halt-on-error", "main.tex"])
        result = subprocess.run(command, cwd=workspace, text=True, capture_output=True)
        log_path = workspace / "main.log"
        log = log_path.read_text(errors="replace") if log_path.exists() else result.stdout + result.stderr
        if result.returncode:
            issue(rows, "TECHNICAL_ERROR", "LaTeX compilation failed", "manuscript/main.tex", "Compiler returned a nonzero status; inspect the build log locally.")
            return "LATEX_COMPILE_FAILED"
        if re.search(r"LaTeX Warning:.*(?:undefined references|Citation .* undefined)", log, flags=re.I):
            issue(rows, "TECHNICAL_ERROR", "Unresolved LaTeX reference or citation", "manuscript/main.tex", "The compiler reported unresolved references or citations.")
        overfull = [float(value) for value in re.findall(r"Overfull \\[hv]box \(([0-9.]+)pt too wide\)", log)]
        if any(value > 10 for value in overfull):
            issue(rows, "WARNING", "Material overfull box", "manuscript/main.tex", f"Largest detected overflow: {max(overfull):.2f}pt.")
        return "LATEX_COMPILE_PASSED"


def manuscript_checks(rows, tex):
    for pattern, name in [(r"\?\?", "Unresolved LaTeX reference"), (r"\bTODO\b", "TODO marker"), (r"\bTBD\b", "TBD marker"),
                          (r"\\cite\{(?:placeholder|citation|ref)\}", "Placeholder citation")]:
        for match in re.finditer(pattern, tex, flags=re.I):
            issue(rows, "TECHNICAL_ERROR", name, f"manuscript/main.tex:{line_number(tex, match.start())}", match.group(0))
    for match in PLACEHOLDER.finditer(tex):
        issue(rows, "SUBMISSION_BLOCKER", match.group(1).upper(), f"manuscript/main.tex:{line_number(tex, match.start())}", match.group(2).strip())

    figure_blocks = re.findall(r"\\begin\{figure\}.*?\\end\{figure\}", tex, flags=re.S)
    labels = []
    for i, block in enumerate(figure_blocks, 1):
        label = re.search(r"\\label\{([^}]+)\}", block)
        path = re.search(r"\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}", block)
        if not label: issue(rows, "TECHNICAL_ERROR", "Missing figure label", f"manuscript/main.tex figure block {i}", "Figure has no label.")
        else: labels.append(label.group(1))
        if not path: issue(rows, "TECHNICAL_ERROR", "Missing figure path", f"manuscript/main.tex figure block {i}", "Figure has no includegraphics path.")
        elif not (MANUSCRIPT.parent / path.group(1)).exists():
            issue(rows, "TECHNICAL_ERROR", "Broken figure path", f"manuscript/main.tex figure block {i}", path.group(1))
    for label in labels:
        if not re.search(rf"\\ref\{{{re.escape(label)}\}}", tex): issue(rows, "WARNING", "Figure never cited", "manuscript/main.tex", label)

    bib_text = BIBLIOGRAPHY.read_text()
    citations = {key.strip() for group in re.findall(r"\\cite(?:t|p)?\{([^}]+)\}", tex) for key in group.split(",") if key.strip()}
    bib_keys = set(re.findall(r"@\w+\{\s*([^,\s]+)", bib_text))
    for key in sorted(citations - bib_keys): issue(rows, "TECHNICAL_ERROR", "Citation missing from bibliography", "manuscript/main.tex", key)
    for key in sorted(bib_keys - citations): issue(rows, "WARNING", "Bibliography item never cited", "manuscript/references.bib", key)
    if not bib_keys: issue(rows, "SUBMISSION_BLOCKER", "Bibliography is empty", "manuscript/references.bib", "No verified bibliography entries are present.")

    author = re.search(r"\\author\{([^}]*)\}", tex, flags=re.S)
    if not author or "blind" not in author.group(1).lower():
        issue(rows, "TECHNICAL_ERROR", "Blinded author check", "manuscript/main.tex", "Author field is absent or not blinded.")
    for match in re.finditer(r"\bnovel phenotype\b", tex, flags=re.I):
        issue(rows, "TECHNICAL_ERROR", "Unsupported novelty claim", f"manuscript/main.tex:{line_number(tex, match.start())}", match.group(0))
    for match in re.finditer(r"[^.]*sensitivity[^.]*replication[^.]*\.", tex, flags=re.I):
        sentence = match.group(0)
        if not re.search(r"not (?:an? )?independent replication", sentence, flags=re.I):
            issue(rows, "TECHNICAL_ERROR", "Sensitivity described as replication", f"manuscript/main.tex:{line_number(tex, match.start())}", sentence.strip())

    expected = [(r"([\d,]+) disease patients", 1850, "disease patient count"), (r"([\d,]+) distinct reference image paths", 942, "reference path count"),
                (r"([\d,]+)[- ]feature", 125, "feature count"), (r"([\d,]+) syndrome--feature", 6125, "main test count"),
                (r"([\d,]+) globally significant", 811, "global significant count")]
    for pattern, correct, name in expected:
        for match in re.finditer(pattern, tex, flags=re.I):
            value = int(match.group(1).replace(",", ""))
            if value != correct:
                issue(rows, "TECHNICAL_ERROR", f"Inconsistent {name}", f"manuscript/main.tex:{line_number(tex, match.start())}",
                      f"Found {value:,}; expected {correct:,}.")

    for name in ["Introduction", "Discussion", "Conclusion", "Data and Code Availability", "Ethics Statement", "Funding", "Conflicts of Interest", "Author Contributions"]:
        content = section_text(tex, name)
        if not content:
            issue(rows, "SUBMISSION_BLOCKER", f"Unresolved {name}", "manuscript/main.tex", "Section is absent or empty.")
        elif not PLACEHOLDER.search(content) and len(re.findall(r"\b\w+\b", content)) < 20:
            issue(rows, "SUBMISSION_BLOCKER", f"Unresolved {name}", "manuscript/main.tex", "Section is effectively empty.")


def package_checks(rows):
    text_files = [path for folder in [ROOT / "publication", ROOT / "manuscript"] for path in folder.rglob("*")
                  if path.is_file() and path.suffix in {".md", ".tex", ".yaml", ".sh", ".yml"}]
    excluded = {REPORT, ROOT / "publication/REVIEWER_METHODS_AUDIT.md"}
    stale, newer = re.compile(r"44,500|9,002|356 syndromes|579 syndromes"), re.compile(r"120-feature", re.I)
    for path in text_files:
        if path in excluded: continue
        for number, line in enumerate(path.read_text(errors="replace").splitlines(), 1):
            if stale.search(line): issue(rows, "WARNING", "Legacy manuscript-era count", f"{path.relative_to(ROOT)}:{number}", line.strip())
            if newer.search(line) and not re.search(r"do not|out of scope|not use|substitute|newer-only|current 120-feature|\bno 120-feature", line, flags=re.I):
                issue(rows, "TECHNICAL_ERROR", "Unqualified newer 120-feature reference", f"{path.relative_to(ROOT)}:{number}", line.strip())


def main():
    rows = []
    tex = MANUSCRIPT.read_text()
    manuscript_checks(rows, tex)
    package_checks(rows)
    latex_status = latex_check(rows)
    technical_errors = sum(row["severity"] == "TECHNICAL_ERROR" for row in rows)
    warnings = sum(row["severity"] == "WARNING" for row in rows)
    blockers = sum(row["severity"] == "SUBMISSION_BLOCKER" for row in rows)
    ready = technical_errors == 0 and blockers == 0 and latex_status == "LATEX_COMPILE_PASSED"
    lines = ["# Publication lint report", "", f"Technical errors: {technical_errors}", f"Warnings: {warnings}",
             f"Submission blockers: {blockers}", f"LaTeX status: {latex_status}", f"Submission ready: {'YES' if ready else 'NO'}", ""]
    if rows:
        lines += ["| Severity | Check | Location | Detail |", "|---|---|---|---|"]
        for row in rows: lines.append(f"| {row['severity']} | {row['check']} | `{row['location']}` | {row['detail'].replace('|', '/')} |")
    else:
        lines += ["No lint findings were detected."]
    lines += ["", "Submission blockers are unresolved human manuscript tasks, not scientific or software failures. The linter exits",
              "nonzero only for technical errors so CI can pass while reporting `Submission ready: NO`.", ""]
    REPORT.write_text("\n".join(lines), encoding="utf-8")
    print(f"Technical errors: {technical_errors}")
    print(f"Warnings: {warnings}")
    print(f"Submission blockers: {blockers}")
    print(f"LaTeX status: {latex_status}")
    print(f"Submission ready: {'YES' if ready else 'NO'}")
    return 1 if technical_errors else 0


if __name__ == "__main__": raise SystemExit(main())
