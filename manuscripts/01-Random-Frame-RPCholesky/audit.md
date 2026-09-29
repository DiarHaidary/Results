# Random-frame RPCholesky: audit and revision

Scope: the complete current `Random_Frame_RPCholesky.tex`, including every
proof and the extensions through the numerical and open-problem sections.
Old split-section drafts were not treated as authoritative. This is an
informal mathematical review and targeted independent computation, not formal
verification or external peer review.

## Mathematical status

The headline Haar-average comparison has a complete proof in the revised
source. No concrete defect was found in its proof chain after checking the
conditional Haar reduction, Schur mixture identity, Gamma tilt, common-order
covariance argument, and induction over the pivot budget. The main result
averages over the frame and pivots. It does not assert error-distribution
domination, a guarantee in every frame, or a resolution of the original
coordinate RPCholesky pivot-count question.

| Component | Review result |
| --- | --- |
| Conditional Haar law of unused basis vectors | Selection weight depends only on the chosen vector; its normalizer is the residual trace and is independent of the unused basis. The induction is valid. |
| Schur-complement mixture | Checked by the determinant lemma and coefficient comparison; zero eigenvalues are inert. |
| Gamma representation of the energy tilt | Checked homogeneity and independence of the total Gamma mass from its Dirichlet proportions. |
| Lemma G | Checked centering, signed Fubini, separate integrability, both weighted Chebyshev signs, the monotonicity of the reweighting function, and the product of centering shifts. |
| Common monotone order | Checked the deletion identities and Newton inequality; the boundary case with exactly `k` positive eigenvalues has zero numerator. |
| Budget induction | First step is an equality, later steps use the one-step inequality, and rank exhaustion is explicitly handled. |
| Shape endpoints | Checked size-biased Gamma mixture, weak kernel limits, bounded trace, and the deterministic elementary-symmetric recurrence. At infinite shape the profile equals volume sampling's expected profile; the residual laws are not identified. |
| Coordinate contrast | The construction is identical to Colbrook's prior example. Checked leading Cauchy minors, scale separation, prefix limits, determinant/probability cancellation, and the count of `2^r` surviving histories. |
| Higher-degree coordinate extension | Expanded the repeated-tail construction and its determinant cancellation so the `(j+1)^r` claim has an explicit derivation. |
| Tight-frame dilation | Checked the coupled residual recursion and equality of nonzero spectra using `FF*=I`. |
| Rotation-invariant probing | Checked projection onto the unused subspace and independence of radial factors from its uniform direction. |
| Volume-rescaled isotropic sketch identity | Checked Cauchy--Binet expectation under independent isotropic columns, determinant cancellation, and the coefficient count. |
| Population transfer | Checked the common base-space coupling, path likelihoods, almost-sure Gram limits, nonnegative Fatou transfer, and factorial-weighted bound. |
| Inverse-Wishart identity | Checked the Schur-complement degrees of freedom, real/complex Gamma scales, and reciprocal Bartlett product. The underlying Wishart facts are classical imported results. |

The monotonicity in Gamma shape and the proposed flat-prior population limits
remain conjectural or heuristic, as labeled. The numerical evidence is not
used to justify the theorems.

## Corrections made

- Removed the optional assertion that centering by `sum c_i a_i=0` implies
  `m_a(t)=O(t^-2)`. It is false in general; `m_a(t)` can have a `1/t` leading
  term. The preceding integrability argument already suffices.
- Corrected the assertion that common Gamma shapes are essential. The
  covariance proof also adapts to unequal positive shapes with centering
  weights proportional to `alpha_i c_i`. The common-shape statement used for
  spectral chains is unchanged.
- Changed the comparison between adaptive probing and plain Gaussian
  sketching to a comparison of proved upper bounds. The bounds alone do not
  imply pointwise ordering of their actual expected errors.
- Removed a claim of proved asymptotic population sharpness that had only
  been supported by simulations; it is now described as suggested by the
  experiments and heuristic limiting law.
- Corrected the intuition that two coordinate pivots have equal probabilities
  in the sharp example. Two contributions are equal only after normalization
  and cancellation in the history sum; some actual histories are rare.
- Distinguished equality with volume sampling's expected symmetric profile
  from equality of residual laws at infinite Gamma shape.
- Corrected the Hájek reference for the question resolved by Yu to the 1981
  book cited by Yu. Added Dereziński--Warmuth--Hsu's general random-design
  volume-rescaled framework alongside the fixed-design reference.
- Added an elementary proof of the Guruswami--Sinop inequality used here.

## Attribution and scoped novelty assessment

The contribution claimed by this manuscript is the all-budget comparison for
a single Haar-rotated arbitrary PSD input, together with the Gamma-chain
argument and its stated probing/population consequences. No earlier statement
of that comparison was identified in the primary sources inspected. This
search does not establish absolute historical priority.

The following primary sources were inspected or their official publication
records checked:

- Chen, Epperly, Tropp, Webber, *Randomly pivoted Cholesky: Practical
  approximation of a kernel matrix with few entry evaluations*,
  <https://arxiv.org/abs/2207.06503> and
  <https://tropp.caltech.edu/papers/CETW25-Randomly-Pivoted-CPAM.pdf>.
  Their modern analysis and the earlier history of the pivot rule are credited.
- Epperly, *A new analysis of the randomly pivoted Cholesky algorithm*,
  full primary text <https://arxiv.org/html/2608.20633>.
  The coordinate trace-count improvement and its different scope are stated.
- Divan and Gilles, *Convergence rates of randomly pivoted methods for
  low-rank approximation*, full primary text
  <https://arxiv.org/html/2609.06287>. Its spectral-decay and approximate-pivot
  analysis is acknowledged.
- Deshpande, Rademacher, Vempala, Wang, volume-sampling paper,
  <https://theoryofcomputing.org/articles/v002a012/>.
- Guruswami and Sinop, *Optimal Column-Based Low-Rank Matrix Reconstruction*,
  <https://arxiv.org/abs/1104.1732>. The benchmark factor is prior work.
- Yu, *On the inclusion probabilities in some unequal probability sampling
  plans without replacement*, <https://arxiv.org/abs/1005.4107>.
  The diagonal endpoint is credited to the established sampling comparison.
- Dereziński and Warmuth, fixed-design volume sampling,
  <https://arxiv.org/abs/1705.06908>.
- Dereziński, Warmuth, Hsu, *Unbiased estimators for random design regression*,
  <https://jmlr.org/papers/v23/19-571.html>. Volume-rescaled random-design
  sampling and determinant expectation identities are prior work.
- Matthew J. Colbrook, *Sharp worst-case factors for randomly pivoted
  Cholesky and LU*, 11 September 2026. The complete public source was read:
  <https://github.com/ajt60gaibb/OpenProblemsInNLA/blob/main/references/colbrook-random-pivoting-2026-09-11/manuscripts/sharp_random_pivoting.tex>.
  Its exact `L(h)`, `eta=h^(n^2)` construction and `2^r` theorem are reproduced,
  not claimed as new. The bibliography links the accompanying PDF.

The vague authorless iterated-logarithm citation is replaced with the printed
title/author of the supporting manuscript included under
`manuscripts/background/Iterated_Logarithmic_RPCholesky.pdf`. The main paper's
author field remains empty, as in the source. No author was invented.

## Independent checks actually run

`checks/independent_audit.py` and the saved JSON contain the exact seed, values,
sample counts, and quadrature error diagnostics.

- Eight budget-two Laplace integrals reproduce every corresponding table
  entry to its printed rounding. For spectrum `(5,3,2,1)`, the real mean is
  `4.415273304157657`, the complex mean is `4.434463831330548`, and the real
  centered covariance term is `-0.13238865746887965`.
- The nonmonotone-scale example gives `D=+7.428186405745322e-5`, confirming the
  stated sign violation numerically.
- Eight independent 100,000-path spectral-chain experiments cover all four
  main table spectra in real and complex shape. They are consistent with the
  stated means and the theorem within their sampling uncertainty.
- Exact rational coordinate computations reproduce the displayed values
  `3.9956302595719944`, `3.999956443339024`, `7.959093773112895`, and
  `7.999570725226107` for ranks two/three and `h=1/10,1/100`.
- Reciprocal Bartlett products through rank six agree with `1/ell!`.

Original Julia scripts and saved outputs were copied into `checks`, with a
reproduction guide. The complete original multi-million-path Julia runs,
population table, 2,000-frame high-precision experiment, and every exploratory
Monte Carlo count were not rerun. No formal interval certification is claimed
for floating-point quadrature.

## Exposition, build, and visual checks

The abstract now centers the comparison formula and pivot count, states its
averaged scope, credits the prior coordinate example, and supplies the
rank-exhaustion convention. Added motivation and a contribution discussion;
rewrote the central covariance and one-step arguments as explained steps;
split compressed displays; expanded the Wishart calculation; and applied the
shared manuscript typography, A4 pages, 26 mm margins, and T1 Latin Modern fonts. All mathematical statements remain available.

The final PDF has **31 pages**. It was built repeatedly with MiKTeX `pdflatex`
using `-interaction=nonstopmode -halt-on-error`, with sufficient final passes
for stable references. The final log has no undefined citations/references
and no overfull boxes. MiKTeX emits a benign inability to write its external
runtime log; the PDF compilation itself succeeds.

All 31 pages were rendered. Four contact sheets covering every page were
visually inspected, followed by full-page inspection of the abstract, central
covariance proof, one-step proof, sketch comparison, population transfer, and
Wishart calculation. No clipping, overlapping formula, or broken table was
found. The source uses only the shared style file as a build dependency and
contains no unavailable local-manuscript citations. Repository citations use
the final published paths; their availability depends on publishing the
complete repository bundle together.
