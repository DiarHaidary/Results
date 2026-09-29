# Product-Probe Drift: audit and editorial revision

Audit date: 29 September 2026. Scope: the complete source `Product_Probe_Drift_v2.tex`, its proofs and stated numerical comparisons. This was an informal mathematical and bibliographic review, with independently rerun numerical diagnostics; it was not formal verification or external peer review.

## Outcome

The central order-one quotient bound and its stated Nyström and trace-estimation consequences survived review. No blocking mathematical error was found in those results. The source was substantively rewritten, and several unsupported scope or asymptotic claims were corrected. The final compiled PDF has **25 pages**.

## Mathematics checked

- Nyström projection formula, sketch-range monotonicity, and the sequential block-update identity, including singular inputs and pseudoinverses.
- Complex block progress lower bound: the normalized complex vector belongs to the complexification of the real projector range.
- Spherical beta moments, product tensorization by conditional moments and Minkowski, convex moment slopes at order one, and the strict improvement over fourth moments.
- Numerator-weighted change of measure, Jensen and mixed Hölder bound in Lemma 4.1. The zero-denominator convention and kernel condition were checked explicitly. The limit as beta tends to zero is legitimate because the left side is independent of beta.
- Trace cooling: eigenvalue monotonicity, convexity of the scalar tail function, and the decreasing convex potential that gives the rounded stopping budget.
- Compression/Jensen proof for arbitrary nonnegative convex potentials; eigenvalue moment interpolation for Schatten decay; ramp potential for the absolute spectral tail and factor-two regularized preconditioner.
- Diagonal training for product Rademacher queries: isotropy, together with coordinate squares equal to one, suffices for its stated mean-square identity. Conditional unbiasedness removes the variance-of-a-conditional-mean term.
- Real interpolation of complex product responses, nonadaptive query lists, integer rounding, and the real/complex leading-term crossover.
- Algebraic domination of the two explicitly defined benchmark budget formulas, including the uncapped/capped distinction and absorption of rounding allowances.
- The real spike quotient and the complex block Gram eigenvalues in the sharp-exponent lower bound, including the exceptional event and constants 16 and 32.
- Sub-Gaussian fourth-moment bound; fresh SRHT-row conditional independence requirements; CountSketch mixed trace identity and whole-block progress.
- Finite angular laws: symmetry gives exact isotropy, weak convergence plus an equicontinuous unit-direction family gives uniform moment approximation. The result remains existential.

## Corrections and limitations made explicit

1. The previous unidentified local RA-11 report is no longer a dependency. Its two budget formulas are defined in full and treated as source-draft benchmarks; the comparison does not certify the absent report's proofs or claim domination of every published algorithm.
2. The number 1.351600 is a benchmark obtained by substituting a subexponential drift penalty into the current sufficient-bound calculation. It is **not** an algorithmic lower bound or a proved architecture floor.
3. The Chernoff bracket applies to the real one-column and complex two-column laws. The larger retained interpolation span is not covered by that lower-bound proof. Unsupported statements assigning it the same asymptotic exponent were removed.
4. A finite-order experiment does not establish convergence to the conjectured exponent. The conjecture and numerical observations remain separated.
5. The fresh-SRHT result supplies the same dimension-uniform sufficient bound as the Rademacher moment argument. It does not prove equality of optimal fixed-dimension drift constants or that SRHT can do no better.
6. The CountSketch bound counts a whole block as one update. Its identity-matrix example demonstrates a loss of order b; a universal optimality claim is not made.
7. A rounded inequality coefficient in Proposition 7.2 was made conservative: `(52*32)^(1/3) = 11.849984...`, so the displayed chained lower coefficient is now **11.84**, rather than rounding upward to 11.85. An independent read-only review by a second agent also identified this issue.
8. The unseeded inherited Monte Carlo table and unreproducible scan totals were replaced by fresh checks. The retained-span simulation was removed because its numerical setup was not supplied.

## Numerical verification

Run `python verify_numerics.py`. Results are in `verification_results.json`; the script uses paths relative to its own file and requires NumPy and SciPy.

- All seven displayed budget rows recomputed to their reported values.
- Binary real and complex order-one constants: 1.3591409142 and 1.2130613194.
- Binary Chernoff exponents: 0.1091170039 (real), 0.0596601011 (complex).
- Comparison ratios: 18.4417451682 and 44.2301985094.
- All displayed exponential bases independently recomputed.
- Central finite-difference derivative error: at most 1.0614e-10 on the stated grid.
- Exact rounded drift branch crossover was checked for q = 2,...,60, both at epsilon = 0.9^q and epsilon = 0.1; it occurs between q = 30 and q = 31. The first few values of 0.9^q lie outside the catalogue's epsilon < 1/2 range; the numerical budget formula can still be evaluated there.
- Seed 20260929, 20,000 independent complex-spherical samples for each q in {4,6,8,10,12,14}. Scalar and two-column drop means are accompanied by estimated standard errors.
- Sequential Nyström identity on PSD/singular matrices of sizes 4, 8, 12: maximum Frobenius discrepancy 8.91e-15.
- Exhaustive small CountSketch mixed identity and product-Rademacher diagonal-training identity: errors about 7.11e-15.
- Large auxiliary-rank budget optimization searches near the continuous minimizer; it is a floating numerical check, not a certified global integer optimization proof.

## Novelty and fair attribution

The proposed contribution is the **uniform conversion of near-one directional moments to expected residual progress**, its explicit product-law penalty in Nyström and trace-query bounds, and the spike-plus-flat exponent bracket. The literature search did not identify an identical result. This is a bounded search conclusion, not priority certification.

Reused ingredients are credited where they enter: product-probe Hutchinson and scalar-field advantages to Meyer–Avron; oracle distinctions to Meyer–Swartworth–Woodruff; the complex order-one constant, interpolation, diagonal correction and finite angular construction to Holden; residual contraction to Chen–Epperly–Tropp–Webber; convex compression and ramp arguments to Epperly. Nyström itself, Khatri–Rao sketching, and the cited embedding results are not claimed as new.

Primary/public sources inspected or verified:

| Source | What was checked | URL |
|---|---|---|
| Meyer–Swartworth–Woodruff | Oracle, conditioned lower-bound scope | https://arxiv.org/abs/2502.08029v2 |
| Meyer–Avron | Estimator-specific rates; spherical/complex probes | https://arxiv.org/abs/2309.04952v2 |
| Holden, 14 September 2026 | D_n, fractional averaging, Lemma 4.2 interpolation, Sections 6–7 constructions | https://github.com/ajt60gaibb/OpenProblemsInNLA/blob/main/references/holden-further-2026-09-14/RA-11/manuscript.md |
| OpenProblemsInNLA RA-11 | Catalogue model and partial-result status | https://github.com/ajt60gaibb/OpenProblemsInNLA/blob/main/randomized-and-low-rank-approximation/RA-11/README.md |
| Chen–Epperly–Tropp–Webber | RPCholesky reference and contraction context | https://arxiv.org/abs/2207.06503 |
| Epperly | Lemma 4.2 and ramp argument inspected in full text | https://arxiv.org/html/2608.20633v1 |
| Camaño–Epperly–Meyer–Tropp | OSI route, Nyström section, Khatri–Rao discussion | https://arxiv.org/html/2508.21189v1 |
| Saibaba–Verma–Ballard | Latest listed version and low-rank approximation scope | https://arxiv.org/abs/2507.23207v2 |
| Beretta–Musco | Latest version v3, fixed-order OSE scope | https://arxiv.org/abs/2608.28094v3 |
| Divan–Gilles | Randomized/approximate pivoting; different sampling model | https://arxiv.org/html/2609.06287v1 |
| Frangella–Tropp–Udell | Nyström preconditioning attribution | https://arxiv.org/abs/2110.02820v2 |
| Haagerup | Sharp Khintchine reference and publisher record | https://www.impan.pl/en/publishing-house/journals-and-series/studia-mathematica/all/70/3/103223/the-best-constants-in-the-khintchine-inequality |

All proof dependencies needed for the new results are stated or proved in the manuscript, except the explicitly labeled imported sharp Khintchine inequality. External references use linked public records. Links to the included numerical script and JSON use the intended repository's main branch.

## Exposition and PDF quality

- Replaced the long abstract with a concise statement of the progress theorem, consequences, and unresolved scope.
- Added a concrete prior-work/contribution discussion and notation/oracle scope.
- Expanded the ratio, cooling, convex-potential, Schatten, query-count, domination, fractional-overhead, and hard-family proofs into explained steps and aligned displays.
- Main theorem statements now display their core guarantees separately.
- Removed a redundant pipeline diagram, split the dense budget table into two panels, and separated empirical observations from theorems.
- Preserved the original empty author field.
- Built repeatedly with pdfLaTeX. Final pass: no undefined references/citations, no overfull boxes. One harmless bibliography underfull-box warning remains.
- Inspected contact sheets covering all 25 pages and full-size pages for the abstract, Theorem 7.1, Proposition 7.2, the hard-family proof, and the numerical tables. No clipping or overlapping content was observed.
