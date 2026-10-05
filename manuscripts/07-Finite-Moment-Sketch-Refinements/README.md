# Finite-moment refinements for SRHT and Khatri-Rao sketches

Diar Heidary - 5 October 2026.

- [Manuscript (PDF)](Finite_Moment_Sketch_Refinements.pdf)
- [Standalone LaTeX source](Finite_Moment_Sketch_Refinements.tex)
- [Reproducible checks](checks/)

The manuscript proves an ordered-prefix SRHT moment certificate, a
covariance-spectrum boundary refinement for the original Gaussian
Khatri-Rao multiplier, and exact fourth-order calibrations. It also proves
a convex transfer to independent normalized spherical-product probes,
their exact covariance reduction, and a direct spherical fourth-moment
certificate. An appendix gives the conditional spherical linear-pencil
Jacobi comparison.

All new statements have proofs. Imported local estimates and classical
ingredients are explicitly attributed. The numerical comparisons concern
specified finite-moment certificates, not optimal embedding thresholds or
new universal leading constants. The spherical-product rows form a
different distribution from Gaussian-product rows. These new results have
not been formally verified in Lean, and no exhaustive novelty audit is
claimed.

## Reproduce the checks

Python 3.9 or later, using only its standard library:

```sh
python checks/run_checks.py
```

Each program regenerates its adjacent JSON result file. SRHT floating-point
searches only locate candidates; rational enclosures and exact arithmetic
certify the displayed sample counts. The Khatri-Rao comparisons use
rational bounds on the required constants. The spherical example is also
verified independently with an exact covariance-basis calculation.

| Check | Certificate |
| --- | --- |
| `check_srht.py` | Ordered-prefix SRHT bounds |
| `check_srht_terminal.py` | Terminal-grade SRHT comparison and earlier scalar envelopes |
| `check_kr_spectrum.py` | Khatri-Rao covariance-spectrum boundary |
| `check_kr_fourth.py` | Gaussian-product fourth-order bounds |
| `check_spherical.py` | Exact spherical-product fourth moment |

## Build the PDF

The source is standalone and needs no external figures or bibliography file.
Use Tectonic, or run XeLaTeX twice:

```sh
tectonic Finite_Moment_Sketch_Refinements.tex
```

The XeLaTeX path uses CMU Serif when installed, otherwise Latin Modern
Roman, with Noto Sans headings. A pdfLaTeX-compatible Latin Modern fallback
is included. Links use dark green and purple.

The source references pinned versions of the prior public manuscripts.
The checked-in PDF and numerical results correspond to this source version.
