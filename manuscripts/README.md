# Companion manuscripts

The papers below share a typography file and contain linked bibliographies.
Definitions, hypotheses, and proofs are included in each paper; explicitly
imported results point to a paper included in this repository or to the
identified public literature.

| Topic | PDF | LaTeX source | Detailed review |
| --- | --- | --- | --- |
| Program overview | [PDF](00-Overview/Program_Overview.pdf) | [Source](00-Overview/Program_Overview.tex) | [Review summary](AUDIT.md) |
| Random-frame RPCholesky | [PDF](01-Random-Frame-RPCholesky/Random_Frame_RPCholesky.pdf) | [Source](01-Random-Frame-RPCholesky/Random_Frame_RPCholesky.tex) | [Audit](01-Random-Frame-RPCholesky/audit.md) |
| Product-probe drift | [PDF](02-Product-Probe-Drift/Product_Probe_Drift_v2.pdf) | [Source](02-Product-Probe-Drift/Product_Probe_Drift_v2.tex) | [Audit](02-Product-Probe-Drift/audit.md) |
| Exact Gaussian range laws | [PDF](03-Exact-Gaussian-Range-Laws/Exact_Gaussian_Range_Laws.pdf) | [Source](03-Exact-Gaussian-Range-Laws/Exact_Gaussian_Range_Laws.tex) | [Audit](03-Exact-Gaussian-Range-Laws/audit.md) |
| SparseStack companion | [PDF](04-SparseStack-Companion/SparseStack_Constants_Nystrom_Adaptive.pdf) | [Source](04-SparseStack-Companion/SparseStack_Constants_Nystrom_Adaptive.tex) | [Audit](04-SparseStack-Companion/audit.md) |
| Two-round SRHT companion | [PDF](05-SRHT-Companion/SRHT_Prescriptions_Nystrom_Adaptive.pdf) | [Source](05-SRHT-Companion/SRHT_Prescriptions_Nystrom_Adaptive.tex) | [Audit](05-SRHT-Companion/audit.md) |
| Approximate and sketched pivots | [PDF](06-Robust-Sketched-RPCholesky/RPCholesky_Approximate_Sketched_Pivots.pdf) | [Source](06-Robust-Sketched-RPCholesky/RPCholesky_Approximate_Sketched_Pivots.tex) | [Audit](06-Robust-Sketched-RPCholesky/audit.md) |

The supporting [Iterated-logarithmic pivot bounds for randomly pivoted
Cholesky](background/Iterated_Logarithmic_RPCholesky.pdf) is included because
the companion arguments explicitly import its finite-count and fractional-moment
statements. It is an unchanged supporting manuscript, rather than an additional
newly audited theorem in this release. The existing graded-operator and hybrid
papers supply the exact operator estimates referenced by the sketch companions.

## Rebuilding the PDFs

Install a TeX distribution with pdfLaTeX and the standard packages used by the
sources, including `amsmath`, `amsthm`, `mathtools`, `lmodern`, `microtype`,
`geometry`, `enumitem`, `booktabs`, `tikz`, `pgfplots`, `hyperref`, `etoolbox`,
and `needspace`. Run from the repository root:

```console
python manuscripts/tools/build.py
```

If pdfLaTeX is not on PATH, use `--pdflatex` followed by its executable path.
The build runs three passes per manuscript and rejects unresolved references.
The shared [style file](manuscript-style.tex) must remain beside the paper folders.

For PDF destination and repository-link checks, install PyMuPDF and run:

```console
python manuscripts/tools/verify_release.py --output release-check.json
python manuscripts/tools/check_urls.py release-check.json link-check.json
```

The second command requires network access. Publisher access restrictions can
produce HTTP 403 despite a valid DOI; accessible preprint alternatives are
provided where available. Neither an HTTP success nor a numerical check certifies
the mathematical novelty or correctness of a theorem.

Each detailed review identifies the numerical checks actually rerun, their
dependencies, and any inherited numerical illustrations that were not rerun.
The proof review is informal and assisted by AI. It is not Lean verification or
external human peer review.

The final [PDF and repository checks](release-check.json) and
[public bibliography checks](link-check.json) accompany this release.
