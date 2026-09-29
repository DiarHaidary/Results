# Mathematical and editorial audit

Reviewed manuscript: `RPCholesky_Approximate_Sketched_Pivots.tex`.
Review date: 29 September 2026. This record reports a mathematical reading,
independent finite numerical checks, source comparisons, and PDF layout
inspection. It is not formal verification or an exhaustive priority search.

## Overall assessment

No fatal gap was found in the functional conversion or in the main
reused-sketch theorem after the corrections below. The headline sketch
guarantee concerns the residual on the original matrix, equivalently an
exact projection onto the selected original columns. The reconstruction
computed from the same sketch has a separate good-event guarantee.

The four-logarithm pivot count relies on the imported initial-moment
estimate from the repository's supporting *Iterated-logarithmic pivot
bounds for randomly pivoted Cholesky*. Its Proposition 7.2, Theorem 6.3,
and the constants on pages 16 and 19 were matched against the imported
formulas. This review did not re-audit that paper's entire entropy proof.
The SparseStack and SRHT constants are imported certificates with exact
repository sources, rather than consequences derived afresh here.

An important attribution correction: the supporting paper already contains
the complete two-sided robust count as Proposition 8.1, on pages 19–21.
The revised abstract, contribution discussion, and two-sided corollary
therefore identify that result as recalled. The distinct developments
here are the functional interface, the complete one-sided trichotomy,
and its application to a single reused sketch.

## Theorem and proof status

| Result | Status and scope |
| --- | --- |
| Order, exhaustion, pivot products, exact path probabilities | The Schur projection proof and determinant telescoping were checked. Zero residuals are absorbing. These identities are established tools, credited to prior pivoting/adaptive-sampling literature. |
| Exact determinant and trace drifts | The determinant-lemma average, the head/tail inequalities, and the trace decrease were checked. The trace drift is prior work. Independent real and complex tests agree with the identities. |
| Scalar completion | Derivatives show the recurrence is increasing and concave. The convex decreasing function `F(u)=1/u-log u` supplies an increment at least `a/r`. The headline statement now specifies `u>epsilon`; no further pivots are needed when the accuracy is already attained. |
| Imported initial moment and Renyi estimate | Conditional inputs. The statements and constants agree with supporting Proposition 7.2 and Theorem 6.3. Their full underlying proofs were not formally verified by this review. |
| Functional robust conversion | All four phases were read and checked: fractional initialization, conditional Jensen cooling, killed determinant contraction above the threshold, and scalar accuracy completion. The statement permits history-dependent kernels and external randomness. |
| Two-sided count | Recalled from supporting Proposition 8.1. The conditional path domination and upper/lower functional comparisons reproduce its guarantee. Its numerical table was independently recomputed. |
| Lower-only bracket | The scalar drift remains valid from any warm start. The bound retains dependence on `log(tr A/tau_r)` from the input. A uniform count without this dependence is still open. |
| Upper-only bracket and finite-budget floor | The first three phases use only upper domination. The tail-only diagonal rule obeys the upper constraint at every state. The matrix is allowed to depend on the fixed budget. The construction now chooses `m>k+r`, ensuring a positive benchmark tail also at `a=1`. All tested edge cases pass. This is not a floor for one fixed finite matrix run until exhaustion. |
| Sketched residual and likelihood identities | Checked by Gram projection and determinant ratios. The selected columns remain independent in both original and sketched spaces until the sketched chain stops. The exact likelihood identity is valid on positive paths. |
| Weighted good-state tests | The three weights are deterministic functions of a fixed reference state. Minkowski and Markov yield the stated fixed-state failure probability without a coordinate union bound. The norm distortion, projection loss, functional comparisons, and per-step chi-square estimate were checked. |
| Main reused-sketch theorem | The hybrid uses the sketched kernel until its first bad state and exact pivoting afterward. Its likelihood second moment is bounded for every fixed sketch. Cauchy–Schwarz is applied jointly to sketch and reference path; the reference states are independent of the sketch. The initial hybrid moment is transferred with half the imported exponent. The terminal error comparison charges only the failure event by `tr A`. The exhaustion convention was made explicit. |
| Same-sketch reconstruction | The QR-factor cancellation and orthogonal decomposition were checked. The factor is `1+eta^2/(2(1-eta)^2)` on the good event, with a corresponding truncated expectation. There is no claim of an unconditional expectation bound for this oblique reconstruction. |
| SparseStack/SRHT row prescriptions | Conditional on the stated moment certificates. SparseStack's fully independent real construction and constants `26.5,179` come from the existing framework's Theorem 10.6 and Appendix B.4–B.5. SRHT's arbitrary-complex certificate comes from the existing hybrid paper's equation (5.6) and Theorem 1.1. The SRHT row count is capped by the padded ambient dimension; a full transform is exact. |
| Pointwise frame gate | The path formula and lower residual spectral tails give the spectrum-dependent factor `Lambda_C`. Both bounds `1<=Lambda_C<=k!` were checked. The fixed frame and fixed target dictionary are essential hypotheses. The flat case is the known projection-DPP/volume-sampling identity; the factorial comparison is due to Deshpande and Vempala. |
| Renyi gate | A valid application of the imported Renyi moment and joint Holder inequality. It is not a new proof of that moment. Numerical crossovers were recomputed. |
| Independent unions and restart corollary | The union projection dominates each individual projection; zero residual columns can be skipped. Independence factors the geometric-mean bound. The unconditional union inequality is valid, but the claimed `O(r/epsilon)` restart method remains conditional on the unproved uniform fractional warm-start conjecture. |

## Corrections made

1. Reduced the abstract to the central reused-sketch guarantee and its
   accuracy/distortion separation. Imported and conjectural parts are
   identified explicitly.
2. Corrected the claim of novelty for the two-sided count: supporting
   Proposition 8.1 already contains it. The corollary now carries that
   attribution in its heading.
3. Added terminal-set and null-transition conventions. An exact reference
   chain may exhaust before the requested budget; its set size is at most
   the step count, rather than always equal to it. The null likelihood
   ratio is one.
4. Corrected the real-sketch/complex-data prefactor. A `2d` real moment
   certificate gives a `2d` complex bound, not literally the definition
   of `MP(d)` with prefactor `d`. The bad-state bound doubles; one extra
   moment order absorbs the resulting square-root factor.
5. Corrected the SRHT prescription when a sufficient threshold exceeds the
   ambient dimension. It uses all rows in that case. The numerical table
   now labels its large values as uncapped sufficient thresholds.
6. Ensured the diagonal floor construction has a positive rank-tail in
   the edge case `a=1`, by requiring `m>k+r`.
7. Replaced statements that a pure chi-square transfer necessarily needs
   `K/epsilon^2` rows with the actual statement: the displayed variance
   estimate certifies the error at that sufficient scale. No lower bound
   on all change-of-measure methods is proved here.
8. Removed the unsupported assertion that the additive failure term is
   universally unavoidable. It is retained because this proof controls
   failed histories only by `tr A`.
9. Clarified that `1+k/(m-k-1)` is the mean for an independent fixed-design
   Gaussian regression sketch, not an exact adaptive same-sketch identity
   or an identity for medians.
10. Distinguished exhaustive finite-history enumeration at finite
    precision from exact arithmetic or a mathematical theorem. Tightened
    the stated fractional-moment numerical range to the rerun values.
11. Replaced the trichotomy figure's crowded lower labels with readable
    conclusions. Its earlier lower-only label also omitted the inverse
    accuracy term; the correct count is in the theorem and caption.
12. Reorganized phase parameters and the weighted good-state proof into
    short explanatory steps and aligned displays. Added shared typography
    and working bibliography links, replacing unnamed private manuscript
    dependencies with repository sources.

## Independent checks actually run

`checks/independent_audit.py` is an independent NumPy implementation and
writes `checks/independent_audit_results.json`. All checks passed:

- 80 real/complex cases for the sketched residual, likelihood, refit,
  determinant drift, and trace drift identities. Maximum relative error
  among the residual/likelihood/refit identities was about `2.52e-15`.
- 80 good-state cases checking all three weighted fourth-moment tests,
  upper and lower functional comparisons, and the chi-square inequality.
- All 24 entries in the robust pivot-count table.
- All eight displayed sketch-parameter rows, including moment orders and
  uncapped row thresholds.
- Five full gate enumerations and five independent-union enumerations.
- All five gate-parameter rows and the Renyi/factorial crossovers at
  ranks 128 and 129 for the two displayed failure probabilities.
- 36 floor-construction edge cases, including `a=1`, zero budget, rank
  one, and higher ranks.

The supplied Julia calculations were also rerun with Julia 1.13.0:

- `identities.jl`: residual and reconstruction identities, and one exact
  likelihood comparison; residual discrepancy about `3.69e-16`.
- `sk_constants.jl`: all robust/sketch count tables, baseline counts
  `71,722,7243,72539,726105`, and the optimized SRHT examples.
- `trichotomy_curves.jl`: every displayed plot coordinate reproduced to
  its printed precision.
- `det01_fractional.jl`: all six scale-separated-family rows reproduced
  by exhaustive enumeration at 4096-bit precision.
- `editor_crossover.jl`: the quoted crossovers and example exponents.
- `sketched_rp.jl`: the 150-trial error/reconstruction table reproduced.
  The largest rerun reconstruction ratio was about 1397.4 in the
  `m=k=40` Gaussian case.
- `sketched_path_quantities2.jl`: the 100-trial path-quantity table
  reproduced to its printed precision.

The latter two reruns are saved separately as `*_review.out`. Empirical
agreement supports the illustrations; it does not certify their good
events, establish statistical significance, or validate an asymptotic
conjecture.

## Prior-work and novelty check

Primary sources examined:

- [Chen, Epperly, Tropp and Webber, RPCholesky](https://arxiv.org/abs/2207.06503):
  prior exact pivoting analysis and trace drift.
- [Divan and Gilles, convergence of randomly pivoted methods](https://arxiv.org/html/2609.06287v1):
  approximate upper-dominated sampling is already in equation (2.7) and
  Theorem 3.1. Section 5.2 analyzes sketchy RPLU using a fresh block per
  iteration. The claimed distinction here is a uniform relative-error
  transfer and one reused sketch for RPCholesky.
- [Epperly, new RPCholesky analysis](https://arxiv.org/html/2608.20633v1):
  prior `O(r/epsilon+r sqrt(log r))` count. The four-log improvement is
  supplied by the supporting repository paper, not newly proved here.
- [Deshpande and Vempala, adaptive sampling](https://eccc.weizmann.ac.il/report/2006/042/):
  the factorial adaptive-sampling/volume-sampling comparison is prior
  work. The source PDF was checked where that comparison appears.
- [Epperly, adaptive randomized pivoting and volume sampling](https://arxiv.org/pdf/2510.02513):
  Theorem 3.3 contains the orthonormal-frame projection-DPP identity.
  Its fixed-frame flat-spectrum case is not a novelty claim here.
- [Camano, Epperly, Meyer and Tropp, structured random matrices](https://arxiv.org/abs/2508.21189):
  credited for the pre-existing SparseStack construction and structured
  sketching context. The numerical constants used here are instead
  attributed to their stated repository certificate.

The manuscript's contribution is stated relative to these inspected
sources. No exhaustive claim of priority is made for the elementary
union inequality, likelihood identities, or general change-of-measure
tools. A spectrum-dependent gate refines an established factorial
comparison; it does not introduce volume sampling.

## Repository dependencies and layout checks

Every non-public supporting proof cited by this manuscript is present in
the repository:

- `manuscripts/background/Iterated_Logarithmic_RPCholesky.pdf`;
- `Graded Multiplication Operators for Random-Matrix Moments.pdf`;
- `Wick-Hermite Hybrid/srht_hybrid_constants_template.pdf`.

The source builds with the shared `manuscripts/manuscript-style.tex`.
The revised PDF bibliography contains real web links rather than local
absolute paths. Links to newly added manuscript/check folders become
public when this revision is uploaded; their targets are checked in the
local repository. PDF annotations were inspected to ensure percent-encoded
repository filenames survived TeX correctly.

The manuscript was compiled repeatedly with `pdflatex` after edits and
all pages were rendered for contact-sheet inspection. Full-size review
covered the abstract, key theorem statements, diagrams, weighted-state
proof, hybrid proof, explicit sketch sizes, and numerical tables. The
final audit distinguishes this layout review from mathematical proof
verification. The final compilation has no unresolved references,
citations, or overfull boxes. The final PDF has 31 pages. The weighted
good-state proof is on pages 16-17, the main theorem and hybrid proof on
pages 18-20, and the explicit sketch prescriptions on pages 20-21.
