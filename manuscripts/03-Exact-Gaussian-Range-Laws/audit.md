# Audit and revision record: exact Gaussian range laws

Reviewed and revised on 30 September 2026. This is a mathematical and
reproducibility review, not formal verification or a claim of exhaustive
literature coverage.

## Scope and outcome

The complete main TeX manuscript was read, including the main results,
projective-space criterion, conditional mechanism, certified families, rank
ceilings, approximation bounds, oracle reduction, generalized Nystrom,
trace estimation, counterexamples, numerical section, and bibliography.
The main real-field results remain supported by their displayed proofs,
after the scope corrections below. No unconditional Gaussian law for
general Khatri-Rao inputs is established. The support must include the
entire spectral tail, and the conditional product class has rank at most
the dimension of the unconditioned factor.

## Corrections and completed dependencies

| Finding | Revision and significance |
| --- | --- |
| The weighted inverse constant was imported from an absent manuscript, which was not the theorem in the repository's existing Khatri-Rao papers. | Replaced that import by Theorem 7.4 and a complete local proof. Convexity of matrix inversion reduces to a squared Gaussian weight; retaining the largest half reduces to a Wishart inverse and a Beta order statistic. The proof explicitly handles both odd and even sample budgets. The constants remain `112/(3*pi)-4/3 = 10.5502357509` and three times this value, `31.6507072526`. |
| The blanket complex-family extension included the real Hurwitz examples. | Restricted the complex statement to Hermitian scalar-covariance certificates. Complex common-Schmidt and single-core common-Schmidt families do transfer. Real symmetric-part Hurwitz identities alone do not. An explicit complex cross covariance of magnitude `1/2` is checked in the independent audit script. |
| The hook construction described its factors as independent. | Explained their shared Gaussian denominator and dependence. It attains the arbitrary-joint-law ceiling, not an independent-factor ceiling. The theorem itself already allowed joint factor laws. |
| A rectangular supported vector was called free through `A^T(A*omega)`. | Explained that a left response is required for the row support. An arbitrary right response is not automatically a legal left product query. The PSD response and left calibration from `Y=A*Omega` remain available without extra queries. |
| The eigenvector lift claimed unit norm without requiring a normalized compressed eigenvector. | Added the condition `norm(y)=1`. |
| Projective notation included the zero support implicitly. | Made `R>=1` explicit; the zero input is handled directly. |
| The eight-bin projective angle illustration stated width `pi/4`. | Corrected the bin width to `pi/8` on `[0,pi)`, matching the actual script. |
| Block Krylov transfer was attributed only to Halko--Martinsson--Tropp. | Added Musco--Musco's primary block Krylov paper. The Gaussian power-method analysis remains credited to HMT. |
| An old trace-exactness floating-point figure varied on rerun. | Updated the reported rerun maximum to `3.2e-12`. Exactness is a theorem in exact arithmetic; the finite-precision numbers are diagnostics. |
| An absent RA-11 report was cited. | Removed it and used the included Product-Probe-Drift companion and the published Kronecker-oracle lower-bound paper. |
| Rounded multiplicative bounds used strict inequality even when the spectral tail vanished. | Changed the weighted-surrogate bound and the main/generalized Nystrom error statements to non-strict inequalities. The constants themselves remain strictly below the rounded constants; the generalized Nystrom error comparison is strict when its tail is positive. Exact-rank inputs now correctly allow zero error on both sides. |

## Mathematical status by component

| Component | Assessment and boundary |
| --- | --- |
| Theorem A, exact criterion | The iid projective lines determine every sketch range. The reverse implication uses the injectivity of the input on its row support. Uniformity on real projective space and almost-sure nonvanishing are both necessary. Equality at one arbitrary sketch size is not asserted to imply equality at size one. |
| Conditional scalar Gaussianity | Standard conditional-law argument gives iid normalized Gaussians independent of their environments. It is sufficient, stronger than projective uniformity, and must hold on the full support. |
| Common-Schmidt, real Hurwitz, and single-core TT families | Their conditional covariances are scalar as claimed. The TT nonzero-scale argument uses its polynomial dependence and isotropy. TT cores regenerated per column and shared/frozen environments are explicitly distinguished. |
| Theorem B, rank ceilings | The smooth image dimension argument supports both `R <= delta+1` and `R <= k+delta` for `k<R`. The latter does not require independent columns. Conditional covariance rank gives `R<=n1`. The attainment example uses dependent factors. |
| Gaussian rangefinder and PSD Nystrom | The pseudoinverse moment, candidate graph subspace, and Ky Fan arguments support the displayed classical coefficients. They are transferred Gaussian bounds, not newly derived general-input embedding bounds. |
| Theorem 7.4, weighted inverse | The added proof covers arbitrary finite positive Schmidt weights summing to one. Positive definiteness and inverse integrability are supplied by the retained Wishart block and the Beta integral. The surrogate bound uses the same nonincreasing row function as the raw regression proof, without a circular dependence in the inverse estimate. |
| Theorem C, product-query reduction | The identity `P_S L L^T v = w*v` proves one-query simulation on the support. The compressed matrix has the same positive spectrum; basis queries reconstruct the PSD input. Query counts assume full response vectors and exact arithmetic. |
| Theorem D, generalized Nystrom | Calibration removes independent left scales exactly. Freezing an environment makes them a common scalar, which cancels. The Gaussian regression identity uses orthogonality of the fitted and residual Gaussian blocks. The `6r+8` call budget is the specified nonadaptive construction and has factor strictly below four. |
| Raw independent left environments | The derivative of the regularized row function is nonpositive. Size-biased monotonicity and the local weighted inverse theorem give the stated conservative coefficient. No automatic two-Gaussian law is claimed. |
| Theorem E, calibrated trace estimator | The Gaussian drift proof, residual telescope, conditional unbiasedness, Gaussian quadratic-form variance, and Chebyshev calculation support the stated bound. The `O(1/epsilon)` principle is classical Gaussian Nystrom++, transferred to legal product queries. The cap by `n1` uses a known dimension, not an unknown rank. |
| Counterexamples | The head-only formula follows from a Cauchy-scale integral and is independently integrated numerically. The left-rescaling example and minimal-norm right-inverse surrogate comparison are valid. |

## Novelty and attribution assessment

The contribution that survives the review is a structured-support
identification: projective uniformity characterizes exact finite-sample
Gaussian range laws; the conditional product certificates expose when this
occurs; a dimension count limits the support; and the same scalar covariance
gives a one-response calibration and exact product-query reduction.
The most substantive algorithmic content is the implementation of Gaussian
procedures by legal product queries and the distinction between calibrated,
frozen, and raw left sketches. The support restriction is essential and
prevents these conclusions from being a solution for unrestricted inputs.

The underlying rotational invariance, Wishart inverse mean, matrix inverse
convexity, Beta order-statistic integral, and Hurwitz identities are classical.
The Gaussian Nystrom coefficient and the trace-estimation rate are not new.
The manuscript now says this explicitly. No searched primary paper was
found stating this combination of exact projective criterion and the
conditional product-query reduction. That is a bounded literature finding,
not a guarantee of priority. The local weighted inverse proof is included to
remove an unpublished dependency; its coefficient is not presented as a
new general embedding result or a sharp weighted-inverse constant.

Primary sources inspected for scope, metadata, and comparison:

- Halko, Martinsson, Tropp: Gaussian rangefinding and power methods:
  <https://arxiv.org/abs/0909.4061>.
- Tropp, Yurtsever, Udell, Cevher: fixed-rank PSD Nystrom, Theorem 4.1:
  <https://arxiv.org/html/1706.05736v1>.
- Tropp, Yurtsever, Udell, Cevher: Gaussian two-sketch reconstruction:
  <https://arxiv.org/abs/1609.00048>.
- Nakatsukasa: generalized Nystrom and stable implementation:
  <https://arxiv.org/abs/2009.11392>.
- Camaño, Epperly, Meyer, Tropp: general-input OSI framework and Remark 2.9:
  <https://arxiv.org/html/2508.21189v1>.
- Bujanović, Grubišić, Kressner, Lam: Khatri-Rao embeddings/eigensolvers:
  <https://arxiv.org/abs/2405.11962>.
- Saibaba, Verma, Ballard: Khatri-Rao approximation on general inputs:
  <https://arxiv.org/html/2507.23207v1>.
- Beretta, Musco: general fixed-subspace Khatri-Rao embeddings, latest
  inspected version 3:
  <https://arxiv.org/abs/2608.28094v3>.
- Meyer, Musco, Musco, Woodruff: Hutch++:
  <https://arxiv.org/abs/2010.09649>.
- Persson, Cortinovis, Kressner: Nystrom++:
  <https://arxiv.org/abs/2109.10659>.
- Musco, Musco: block Krylov approximation:
  <https://arxiv.org/abs/1504.05477>.
- Meyer, Swartworth, Woodruff: unrestricted product-query complexity:
  <https://proceedings.mlr.press/v267/meyer25a.html>.
- Lenzhen, Morier-Genoud, Ovsienko: classical Hurwitz identities:
  <https://arxiv.org/abs/1007.2337>.
- Rakhshan, Rabusseau: tensorized random projection construction:
  <https://proceedings.mlr.press/v108/rakhshan20a.html>.
- Cazeaux, Dupuy, Figueroa Justiniano: tensor-train sketching context:
  <https://arxiv.org/abs/2603.11009>.
- Frangella, Tropp, Udell: Gaussian Nystrom preconditioning:
  <https://arxiv.org/abs/2110.02820>.

## Exposition and repository completeness

The abstract was shortened to about 150 words and centers the exact-law
criterion, the Gaussian Nystrom coefficient, and the thin-class limitation.
A related-work subsection distinguishes the contribution from established
Gaussian facts and general-subspace embedding results. Long proof paragraphs
were reorganized into named steps, with separate displays in the rangefinder,
Nystrom, weighted inverse, oracle reduction, regression, raw reconstruction,
Gaussian drift, and head-only counterexample arguments.

The shared style is loaded through `../manuscript-style.tex`. It is the only
TeX input dependency. Every literature entry now contains a primary-source
URL. The sole companion-paper citation points to
`../02-Product-Probe-Drift/Product_Probe_Drift_v2.pdf` through its canonical
repository URL. That target exists in the release tree. The reproducibility
directory is included beside this paper. Repository web URLs become live on
publication of this prepared tree; the local target and URI spelling were
checked. No absent local report is needed to read or use the theorems.

## Checks actually run

All six supplied Julia scripts were copied into `checks/` and rerun with
Julia 1.13.0, one Julia thread, and one OpenBLAS thread. They all exited
successfully. The `*_rerun_output.txt` files record those runs. Their law
comparisons, dimensional examples, product identities, calibrated/frozen/raw
regression means, exact counterexample values, and trace tables agree with
the reported rounded numbers. The trace exactness rounding figure was
updated as described above.

An independent NumPy/SciPy script, `checks/audit_checks.py`, was also run.
Its output is `checks/audit_checks.json`. It checked:

- 10,000 admissible rank/sample budgets for the even/odd weighted-inverse
  coefficient. The maximum was `10.5502357509`, at `r=1,k=8`.
- 10,000 radial values for the squared-Gaussian CDF inequality.
- Five head-only Cauchy integrals from `y=1` through `y=1e-4`; their
  differences from the closed form were below `2e-15`.
- Scalar conditional covariance, calibration, one-query simulation,
  spectrum preservation, full PSD recovery, and an exact coupled Nystrom
  output. The largest absolute matrix discrepancy was approximately
  `5.1e-14`; the coupled Nystrom relative error was `7.7e-15`.
- The nonzero complex cross covariance in the real Hurwitz example.

The final TeX source was compiled repeatedly with `pdflatex` until references
stabilized. The final log has no undefined-reference, overfull-box, or
underfull-box warning. All pages were rendered for contact-sheet review;
critical theorem/proof/figure pages were also inspected at full-page size.
The final bibliography spacing was adjusted to avoid an isolated final
reference page. No clipped display, table, or figure was observed.
After the zero-tail correction, two additional successful compilations
produced the final 30-page PDF. The affected main-result, weighted-bound,
and generalized Nystrom corollary pages were rendered and inspected again;
the final log still contains none of the warnings listed above.

These finite computations check the stated algebra and reproduce the
illustrations. They do not prove a probability law, establish optimal
constants, classify all admissible supports, or analyze numerical stability.
