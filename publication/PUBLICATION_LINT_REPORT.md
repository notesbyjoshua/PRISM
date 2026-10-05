# Publication lint report

Technical errors: 0
Warnings: 0
Submission blockers: 12
LaTeX status: LATEX_COMPILE_NOT_CHECKED
Submission ready: NO

| Severity | Check | Location | Detail |
|---|---|---|---|
| SUBMISSION_BLOCKER | HUMAN_INPUT_REQUIRED | `manuscript/main.tex:15` | add the clinical motivation, author-approved interpretation, and journal-length word count. |
| SUBMISSION_BLOCKER | HUMAN_INPUT_REQUIRED | `manuscript/main.tex:19` | position the study against the clinical dysmorphology and quantitative phenotyping literature using verified citations. Do not infer novelty from HPO absence. |
| SUBMISSION_BLOCKER | HUMAN_INPUT_REQUIRED | `manuscript/main.tex:88` | add a literature-verified and clinician-reviewed interpretation. Do not classify candidate associations, claim novelty, or convert correlated significant measurements into independent discoveries until that review is recorded. |
| SUBMISSION_BLOCKER | HUMAN_INPUT_REQUIRED | `manuscript/main.tex:94` | provide an author-approved conclusion limited to quantitative association and reproducibility. |
| SUBMISSION_BLOCKER | HUMAN_INPUT_REQUIRED | `manuscript/main.tex:97` | insert mentor-, institution-, and data-provider-approved data/code availability language and the authorized public-release boundary. |
| SUBMISSION_BLOCKER | HUMAN_INPUT_REQUIRED | `manuscript/main.tex:100` | insert the institution-approved ethics, consent, exemption, and oversight language. |
| SUBMISSION_BLOCKER | HUMAN_INPUT_REQUIRED | `manuscript/main.tex:103` | insert funding sources and grant identifiers, or an approved no-funding statement. |
| SUBMISSION_BLOCKER | HUMAN_INPUT_REQUIRED | `manuscript/main.tex:106` | insert the author- and journal-approved conflicts-of-interest statement. |
| SUBMISSION_BLOCKER | HUMAN_INPUT_REQUIRED | `manuscript/main.tex:109` | insert contribution roles after all authors approve them. |
| SUBMISSION_BLOCKER | HUMAN_INPUT_REQUIRED | `manuscript/main.tex:112` | add only verified references to \texttt{references.bib} and cite them in the text. |
| SUBMISSION_BLOCKER | HUMAN_INPUT_REQUIRED | `manuscript/main.tex:117` | review and approve the exact historical feature-definition supplement before submission. |
| SUBMISSION_BLOCKER | Bibliography is empty | `manuscript/references.bib` | No verified bibliography entries are present. |

Submission blockers are unresolved human manuscript tasks, not scientific or software failures. The linter exits
nonzero only for technical errors so CI can pass while reporting `Submission ready: NO`.
