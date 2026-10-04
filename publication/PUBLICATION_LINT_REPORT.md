# Publication lint report

Technical errors: 0
Warnings: 0
Submission blockers: 11
LaTeX status: LATEX_COMPILE_NOT_CHECKED
Submission ready: NO

| Severity | Check | Location | Detail |
|---|---|---|---|
| SUBMISSION_BLOCKER | AUTHOR INPUT REQUIRED | `manuscript/main.tex:15` | add the clinical motivation, final interpretation, and journal-length word count. |
| SUBMISSION_BLOCKER | AUTHOR INPUT REQUIRED | `manuscript/main.tex:19` | position the study against the clinical dysmorphology and quantitative phenotyping literature. Do not infer novelty from HPO absence. |
| SUBMISSION_BLOCKER | AUTHOR AND CLINICAL REVIEW REQUIRED | `manuscript/main.tex:59` | interpret recovered established associations, quantitative extensions, and candidate underreported associations without converting significant tests into independent discoveries. |
| SUBMISSION_BLOCKER | AUTHOR INPUT REQUIRED | `manuscript/main.tex:65` | provide a cautious conclusion limited to quantitative association and reproducibility. |
| SUBMISSION_BLOCKER | INSTITUTIONAL CONFIRMATION REQUIRED | `manuscript/main.tex:68` | state which aggregate outputs, code, and any controlled-access data can be shared. |
| SUBMISSION_BLOCKER | INSTITUTIONAL CONFIRMATION REQUIRED | `manuscript/main.tex:71` | insert the approved ethics, consent, exemption, and oversight language. |
| SUBMISSION_BLOCKER | INSTITUTIONAL CONFIRMATION REQUIRED | `manuscript/main.tex:74` | insert funding sources and grant identifiers, or an approved no-funding statement. |
| SUBMISSION_BLOCKER | AUTHOR CONFIRMATION REQUIRED | `manuscript/main.tex:77` | insert the journal-approved disclosure statement. |
| SUBMISSION_BLOCKER | AUTHOR CONFIRMATION REQUIRED | `manuscript/main.tex:80` | insert contribution roles after all authors approve them. |
| SUBMISSION_BLOCKER | LITERATURE REVIEW REQUIRED | `manuscript/main.tex:83` | add verified references to \texttt{references.bib} and cite them in the text. |
| SUBMISSION_BLOCKER | Bibliography is empty | `manuscript/references.bib` | No verified bibliography entries are present. |

Submission blockers are unresolved human manuscript tasks, not scientific or software failures. The linter exits
nonzero only for technical errors so CI can pass while reporting `Submission ready: NO`.
