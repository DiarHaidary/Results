# Mathematical audit, exposition revision, and novelty assessment

This release revises six manuscripts and replaces the program overview.
The review read the authoritative complete LaTeX files, rather than relying on
older section drafts or previous audit conclusions. Each paper received a
mathematical and literature review before its substantive rewrite. Additional
cross-reviews examined the random-frame/sketched-pivot connection, the product
drift and budget comparisons, and the two sketch companions' certificates and
energy arguments.

The review is informal and assisted by AI. It is not formal verification, an
external human referee report, or a certificate of priority. The supporting
four-logarithm manuscript is included unchanged and its particular imported
statements were checked for provenance and consistency; its entire proof was
not independently re-audited in this release. The existing repository framework
likewise supplies explicitly identified operator estimates as black boxes.

## What changed mathematically

| Paper | Material corrections and qualifications |
| --- | --- |
| [Random frames](01-Random-Frame-RPCholesky/audit.md) | Removed a false optional asymptotic-decay assertion in the Laplace argument; corrected the supposed necessity of equal Gamma shapes; distinguished the one-step majorant from actual chain monotonicity; credited Colbrook's identical coordinate construction; corrected the description of probability/residual cancellation. |
| [Product drift](02-Product-Probe-Drift/audit.md) | Removed unsupported optimality and “no better” conclusions from upper bounds; made historical budget comparisons self-contained; attributed diagonal correction and the prior probability-sensitive construction; separated sufficient-bound architecture limits from oracle lower bounds; corrected a rounded lower-bound coefficient from 11.85 to the safely downward-rounded 11.84. |
| [Gaussian range laws](03-Exact-Gaussian-Range-Laws/audit.md) | Restricted the complex extension to directly verified Hermitian scalar covariance; corrected the dependent-factor hook example and a missing unit-vector condition; explained the rectangular transpose-query convention; replaced an absent weighted-inverse companion with a complete proof; corrected strict tail-scaled inequalities at zero tail. |
| [SparseStack](04-SparseStack-Companion/audit.md) | Replaced a mismatched private v3 citation by the actual framework theorem and appendix; corrected regression coefficients versus fitted predictions; stated the parameter choice needed for the advertised row rate and the moment-order condition; removed an unsupported full-rank claim; separated unbiased estimation from confidence validation. |
| [Two-round SRHT](05-SRHT-Companion/audit.md) | Separated ideal least-certified width from the heuristic search; seeded the search with the proven canonical prescription; capped estimators at the ambient dimension; corrected inherited-versus-new attribution; replaced a false monotonicity possibility and an already-impossible expected-tail open problem by accurate questions. |
| [Robust/sketched pivots](06-Robust-Sketched-RPCholesky/audit.md) | Added absorbing-state conventions and premature sketch-exhaustion handling; corrected real-sketch/complex-input prefactors and full-transform caps; retained the additive failure term; removed unsupported necessity language; identified the two-sided bracket theorem as already present in the supporting paper, Proposition 8.1. |

The reviewed headline arguments remain in the revised manuscripts with their
conditions and imports made explicit. The corrections do not turn any of the
open unrestricted targets into a solved problem. In particular, original
coordinate RPCholesky at O(r/epsilon), unrestricted adaptive Kronecker trace
complexity, and the uniform fractional warm start remain open in these papers.

## Novelty relative to the sources checked

The strongest novelty assessment concerns a precise theorem, rather than the
whole collection or the reuse of a familiar method. The papers' related-work
sections and detailed reviews give the primary sources and comparisons.

| Paper | Identifiable contribution | What is inherited or already known |
| --- | --- | --- |
| Random-frame RPCholesky | An all-budget comparison between the Haar-averaged sequential RP law and the volume-sampling expected trace benchmark; the ratio-Gamma covariance proof and spectral-chain consequences. No matching all-budget theorem was found in the sources checked. | Volume-sampling identities and optimal tail comparisons; adaptive sampling; the coordinate 2^r sharp example, which is Colbrook's construction. Dirichlet-family monotonicity remains conjectural. |
| Product-probe drift | An order-one, numerator-tilted progress inequality for every PSD residual under product probes, with Nyström and trace-query consequences and a bracket on the best drift exponent. | Product-query access, spherical moment formulas, Hutchinson-style probing, Nyström variance reduction, Holden's diagonal correction and probability-sensitive product probes. The displayed benchmark domination is not domination of all published algorithms. |
| Exact Gaussian range laws | A necessary and sufficient projected-line criterion for a fixed support, structural examples and rank ceilings, and calibrated/frozen-environment reductions to Gaussian algorithm laws. | Gaussian rotation invariance, inverse-Wishart and pseudoinverse means, classical randomized SVD/Nyström/generalized Nyström and block Krylov analysis. The class is restrictive, with an n1-dimensional exact reduction. |
| SparseStack | An exact rational certificate for the displayed sharper sufficient constants; an elementary deterministic Nyström compression inequality and exact signed-hash specialization; restricted fixed-reference frame reuse applications. | CountSketch/OSNAP/SparseStack, optimal asymptotic embedding rates, the graded moment majorant, approximate-matrix-multiplication second moments, stochastic trace/diagonal probing, and the factorial adaptive-versus-volume comparison. This is a refinement, not a new sketch distribution or asymptotic rate. |
| Two-round SRHT | Retaining finite sampling density and the exact scalar factor to reduce the sufficient coefficient from 72 to 20; explicit certified prescriptions and shared-row energy applications. | Fast randomized transforms, the imported row-factor and moment kernels, the exact fixed-size barycenter, Yang's rerandomized SRHT result, and nonadaptive Nyström++ trace estimation. The heuristic selector does not certify a globally minimal width. |
| Approximate/sketched pivots | A functional conversion usable without coordinatewise probability brackets, the one-sided distinction, and a path-local transfer for a single reused sketch with weighted tests and an additive failure allowance. | Coordinate drift, the supporting four-logarithm moment theorem, its existing two-sided probability-bracket theorem, classical projection-DPP/factorial comparisons, and fresh-sketch randomized pivoting. Restart O(r/epsilon) remains conditional on a missing uniform warm start. |

Important primary attribution added or clarified includes:

- Chen–Epperly–Tropp–Webber for original RPCholesky and its analysis;
  Epperly and Divan–Gilles for subsequent pivot analyses.
- Deshpande–Rademacher–Vempala–Wang and Guruswami–Sinop for the volume
  expectation and tail benchmark; Deshpande–Vempala for the factorial
  adaptive-versus-volume likelihood mechanism.
- Colbrook for the prior sharp fixed-coordinate pivoting construction.
- Meyer–Swartworth–Woodruff for the product-query model, Meyer–Avron for
  estimator-specific Kronecker difficulty, and Holden for closely related
  September 2026 probability-sensitive probes and diagonal correction.
- Halko–Martinsson–Tropp, Musco–Musco, Tropp and collaborators, Nakatsukasa,
  and Persson–Cortinovis–Kressner for the established randomized low-rank,
  Krylov, generalized Nyström, and trace-estimation ingredients.
- Kane–Nelson, Nelson–Nguyen, Chenakkod–Dereziński–Dong,
  Camaño–Epperly–Meyer–Tropp, and Mai–Rao for sparse sketch constructions
  and embedding rates; Ailon–Chazelle, Boutsidis–Gittens, and Yang for the
  transform literature.
- Hutchinson, Bekas–Kokiopoulou–Saad, and the XTrace/XDiag work for trace
  and diagonal probing principles.

No blanket “first”, optimality, or independent priority claim is made for the
collection. The precise new inequalities and parameter refinements should be
assessed against the scope of their displayed hypotheses.

## Computational and proof checks

The accompanying detailed audits identify the exact commands, fixtures,
seeds, tolerances, source scripts, and inherited numerical illustrations.
Representative checks actually performed include:

- All 63 SparseStack interval certificates in exact rational arithmetic,
  including integer upward rounding and agreement with the printed table;
  the maximum certified value is 499428/1000000, strictly below 1/2.
- SparseStack ceiling prescriptions and likelihood-price tables; exhaustive
  small signed-hash covariance, drift, trace, and diagonal calculations;
  pointwise full-rank RP set-law comparisons; and log-domain majorant grids.
- Product-probe printed slopes, budgets, crossover values, and exponential
  bases; sequential Nyström identities; independent seed-fixed block/scalar
  simulations with standard errors; and the rounded lower-bound correction.
- Random-frame Laplace integrals for all eight budget-two table entries,
  the nonmonotone-scale counterexample, independent spectral-chain
  simulations, and exact rational coordinate calculations.
- Gaussian-range matrix identities, weighted inverse/order-statistic
  inequalities, the head-only counterexample, and an explicit complex
  Hurwitz countercheck.
- High-precision SRHT prescriptions, integer selector certificates and
  exact small-model covariance/energy fixtures; the original Julia width
  calculation and selected exhaustive whole-sketch experiments were rerun.
- Real and complex sketched-pivot likelihood, residual, reconstruction,
  weighted-state and chi-square identities; parameter tables; and exact
  small gate/union enumerations.

Finite-grid checks and Monte Carlo experiments are numerical evidence, not
proofs of conjectures or uniform inequalities. Large inherited simulation
tables that were not rerun are identified as such in the detailed audits.

## Exposition and repository completeness

The abstracts now foreground the central result and its material hypothesis.
The main theorem statements appear early. Proofs explain the purpose of each
construction and separate long algebraic chains into readable steps and
aligned displays. Shared typography supplies paragraph, theorem, list, and
display spacing; critical diagrams and tables were checked visually and
repaired when needed. The original screenshot's comparison proof is now
spread across explained steps, rather than two compressed paragraphs.

All seven PDFs are rebuilt from the authoritative sources. The complete
rendered pages were reviewed through contact sheets, with full-size inspection
of proof-heavy pages, abstracts, diagrams, and numerical tables. The final
build and link records accompany the release; they are regenerated from the
final versions, rather than from an earlier draft.

Private source-file citations have been removed. The references point to
identified public literature or to exact repository filenames. The existing
framework contains the imported SparseStack and SRHT operator estimates.
The missing supporting RPCholesky paper is included as a PDF. The range-law
paper now contains its own weighted-inverse proof, so the unavailable companion
is not required.

Repository URL checks include exact filename case and URL decoding. Public
bibliography HTTP checks distinguish missing targets from publisher access
restrictions: some valid DOI targets reject automated access with HTTP 403.
Accessible primary preprints are linked alongside those DOIs where available.
The per-paper reviews retain the literature source URLs used for attribution.

See the [manuscript index](README.md) for all PDFs, sources, review reports,
build instructions, and reproducible checking scripts.

## Final release checks

The seven revised PDFs contain 181 pages. The final builds have no unresolved
references, duplicate labels, or overfull boxes. The
[release check](release-check.json) validates internal PDF destinations,
LaTeX imports, relative Markdown links, and all 15 distinct repository URLs
against the exact filenames in the checkout; it reports no errors.

The [bibliography check](link-check.json) records all 67 distinct external
URLs embedded in the revised PDFs. Sixty-three returned HTTP 200. Four DOI
destinations returned publisher HTTP 403 access restrictions; their DOI
metadata was checked, and accessible primary alternatives were added for
the Ailon--Chazelle, Drineas--Kannan--Mahoney, and Kane--Nelson papers.
These access restrictions do not indicate missing repository dependencies.
