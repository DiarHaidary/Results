# SRHT companion: audit and revision record

Audited 29 September 2026. The audit covers the complete supplied main TeX
source, the repository sources of its imported estimates, its scalar
prescriptions, and the fourth-moment arguments. It is a mathematical review
with numerical checks, not formal verification or independent human peer
review.

## Mathematical findings and corrections

The canonical coefficient-20 proof is valid conditional on the imported
retained row-factor/Bernoulli bound. Its two sampling-density cases, atom
estimate, integer rounding, limited-independence degree count, and confidence
calibration were checked. The rank-only coefficient 11580, the comparison of
the four sufficient coefficients, and the reduced crude-root order were also
checked. No concrete error was found in these scalar proofs.

The selector's heuristic bisection was conflated with an ideal minimum in the
main-results statement. The revision separates those assertions, rechecks
every accepted candidate, supplies a canonical starting point, and explicitly
returns the zero certificate when all rows are retained. It does not claim
that bisection always finds the global minimum.

Global monotonicity of the stated uncapped root is actually false. The audit
found the elementary witness n=4, d=1, q=1, flat factors:
Lambda=25+2 sqrt(66), R(2)=10.2759465..., R(3)=18.3051974.... The revision
proves this comparison exactly and replaces the original monotonicity open
problem with a question about the useful certification region or the stronger
capped allowance. Both witness roots exceed the useful distortion range;
this example does not disprove monotonicity within a feasible region.

The deterministic Nyström comparison, finite-population identity,
Rademacher fourth moments, two-round and one-round energy identities,
covariance map, and trace/diagonal variances were reviewed. Their matrix
orders and finite-population normalizations are consistent. The estimator
widths now respect M<=n; retaining all rows is treated as exact recovery.
The variance discussion now identifies the coefficient 2 explicitly, rather
than suggesting that the displayed bound has the same coefficient as an
importance-sampled matrix-product bound.

The full-rank fixed-range adaptive gate was checked against the ordered
Schur-complement product and spectral-tail denominators. Its factorial
comparison is classical; the revision identifies that attribution. The
fresh-row proof uses one Khintchine envelope and the established size-biased
ratio argument, and does not treat the rows of one sketch as independent
updates.

The sketched-RP section states an imported interface, including an explicit
failure term. It conservatively imposes q>=max(2,k), covering the determinant
path estimate; the independent companion's hybrid path comparison permits
q>=2 and order O(eta^2 k+log(k tr(A)/(epsilon tau))). The text explicitly
recognizes this stronger interface and replaces the obsolete question about
whether a chi-squared comparison is available.
The SRHT row budget still grows quadratically with the pivot budget under
this interface. No new proof of the imported path theorem is claimed here.

The old question asking for a uniform expected tail-relative Nyström bound
has been corrected: the finite-support obstruction already proved in the
SparseStack companion applies to this law when M<n. The revision gives the
spiked-matrix explanation and asks instead for high-probability tail-relative
bounds.

## Exact repository dependencies

All named repository sources are linked in the bibliography. The critical
embedding import is present in existing repository papers:

- `Wick-Hermite Hybrid/srht_hybrid_constants_template.pdf`, Theorem 4.1
  (hybrid subset-effect estimate), and Section 5, including the explicit
  density-sensitive Bernoulli root. This is the exact Lambda_q estimate used
  by the companion.
- `Graded Multiplication Operators for Random-Matrix Moments.pdf`, Theorem
  10.1, Lemmas 3.12, 3.16, 3.17, and Appendix A.2--A.3. These supply the
  positive-selector, kernel, limited-independence and barycenter background.
- `SRHT.pdf` supplies the earlier cubic-allowance Walsh theorem, rather than
  the sharper Lambda_q imported estimate.
- `manuscripts/background/Iterated_Logarithmic_RPCholesky.pdf`, the product
  probe companion, the revised SparseStack companion, and the revised
  robust/sketched-RP companion supply the respective referenced results.

No missing exact source for the central SRHT import was identified. The audit
checked its provenance and subsequent use; it did not formally replay all
underlying kernel proofs or the Lean projects of earlier papers.

## Numerical checks actually run

`checks/audit_srht.py` is an independent Python implementation, using 80
decimal digits for certificate calculations. Its recorded
`checks/audit_results.json` reproduces the main row counts and selected
orders, both boundary certificates for each displayed direct width, the
adaptive-gate table, the nonmonotone-root witness, and binomial atoms for all
n<=100. For explicit PSD fixtures it enumerates every real n=4 sign/subset
law, checks the Gram energy and covariance maps in exact rational arithmetic,
and checks deterministic Nyström comparisons numerically.

The original Julia scripts were also rerun under Julia 1.13.0:
`prescriptions.jl`, `minscan.jl`, `gate.jl`, `scalars.jl`,
`wholesketch.jl`, `oneround_diag.jl`, and `sketched_rp.jl`.
The exhaustive floating-point width scan reproduced minima 58396, 125750,
346315, 1423395 over orders q<=40; their adjacent boundary values were
independently checked in high precision. The full original whole-sketch
enumeration reproduced the displayed n=4 exact-rational identities and n=8
floating-point results, including the rank-three M=3 two-round energy
2512.321428575... and residual 20.012503603.... These finite scans do not
establish a theorem about an unsearched parameter range.

## Novelty and attribution

The revised introduction distinguishes the following claims:

1. The constant-failure optimal row order for fixed real frames under the
   exact two-Walsh law is prior work, including Yuning Yang's Theorem 1; this companion does not
   claim that first result.
2. The retained row-factor estimate, density-sensitive Bernoulli moment,
   exact centered barycenter coupling, and quadratic confidence allowance
   are inherited. The present scalar calibration retains those factors
   together and yields coefficient 20, rather than 72, for this particular
   sufficient condition. It is not an optimality claim or a comparison with
   every prior proof.
3. The whole-sketch calculations provide explicit finite-population
   fourth-moment identities for these shared-row laws. Approximate matrix
   multiplication is the established interpretation of their error scale.
4. The factorial gate is Deshpande--Vempala's classical comparison. Its
   combination with this explicit SRHT budget is the application here.
5. Hutch++, Diag++, and especially Nyström++ precede the estimator
   architecture. The structured-sketch variance calculation, not the
   sketch-and-correct idea or the query order, is the claimed specialization.

Primary sources inspected or located during the check include:

- https://arxiv.org/abs/2602.05394v3
- https://arxiv.org/abs/2609.13946v1
- https://github.com/yuningyang19/OpenProblemsInNLA_TR-01/blob/ed21181197ac839eac95f549404f94e7e3aa6e10/manuscript.pdf
- https://arxiv.org/abs/2508.21189
- https://arxiv.org/abs/1011.1595
- https://arxiv.org/abs/1001.2738v2
- https://eccc.weizmann.ac.il/report/2006/042/download
- https://arxiv.org/abs/2010.09649
- https://arxiv.org/abs/2201.10684
- https://arxiv.org/abs/2109.10659
- https://proceedings.mlr.press/v202/balabanov23a.html
- https://www.cs.princeton.edu/~chazelle/pubs/FJLT-sicomp09.pdf
- https://www.cs.yale.edu/homes/mmahoney/pubs/matrix1_SICOMP.pdf

This was a targeted literature check. Yang's full arXiv PDF was obtained
and pages 1--5 inspected: Theorem 1 (page 2) gives the prescribed width for
every fixed real frame, every Walsh dimension, and 0<epsilon<1, with failure
at most 0.01. Its sampling law matches the two-round law here. Both the
arXiv version and pinned companion-source link are retained. No claim of
exhaustive novelty clearance or independent priority certification is made.

## Exposition and build

The abstract is shortened around the canonical prescription and one principal
whole-sketch consequence. The introduction now explains the quantitative
question, imported assumptions, contribution boundary, and related work.
Long projector, population, covariance, fourth-moment and ratio arguments
are broken into explained stages and separate displays. The paper uses the
shared repository spacing stylesheet. The figures are integrated TikZ and
PGFPlots; there are no external figure-PDF dependencies.

Build from this directory with `pdflatex -interaction=nonstopmode
-halt-on-error SRHT_Prescriptions_Nystrom_Adaptive.tex`, twice after changes
to labels or pagination. Final rendering/page review is recorded below.

The final PDF has 27 pages. All pages were rendered and inspected in four
contact sheets, with larger page views of the abstract, canonical proof,
selector, adaptive gate, deterministic Nyström comparison, fourth-moment
and ratio proofs, numerical tables and bibliography. The final two MiKTeX
passes completed without warnings, undefined citations/references, duplicate
labels, or overfull/underfull boxes. No clipping or collisions were found.
A float barrier keeps the numerical tables together before the open problems
and bibliography. The source's only TeX input is the included repository
stylesheet; all figures are built directly from the manuscript source.

The numerical scripts and their actual outputs remain under `checks/`.
Revision scripts, extracted source texts, compiler intermediates, rendered
QA images and duplicate stand-alone figure PDFs were removed from the
publication directory after review. A read-only cross-review of the revised
SparseStack scalar-certificate and deterministic Nyström arguments found no
concrete mathematical error; a minor clarification about rounding the
reciprocal square-root expression was reported to that manuscript's editor.
