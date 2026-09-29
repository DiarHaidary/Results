# SparseStack companion: mathematical and editorial audit

Scope: all statements and proof bodies in `SparseStack_Constants_Nystrom_Adaptive.tex`; the exact certificate table and listing; the parameter table; the original numerical source scripts and recorded outputs; and the imported SparseStack estimates in the actual repository PDFs. This is a targeted mathematical audit and reproducible computational check, not formal verification.

## Main conclusions

The refined constants proof is supported by its stated imported operator majorant. All 63 exact rational certificates reproduce the printed bucket bounds and integer geometric-mean upper bounds. The largest bound is `499428/10^6 < 1/2`. The vacuum-boundary inequality and the `q=1` and ceiling arguments were checked analytically. No error was found in those arguments.

The deterministic Nyström compression bound, signed-hash second-moment identity, finite-support obstruction, block mixed trace and ratio bounds, convex-potential compression argument, bucket-row moment calculation, diagonal/trace variances, and full-rank RP selection guarantee were read in full. The independent checks agree with these arguments. The contribution is incremental: the moment majorant, SparseStack construction, approximate-matrix-multiplication covariance calculation, stochastic probing mechanisms, and factorial adaptive-versus-volume comparison are inherited.

## Corrections made

1. **Import provenance.** The original bibliography named a private SparseStack v3.0 manuscript. The existing repository's `SparseStack.pdf` is an earlier coupled version and does not supply the imported constants. The exact imports are instead present in the repository framework: Theorem 10.6 (page 113), Appendix B.4–B.5 (pages 137–139), equations (224)–(226). The paper now cites that source precisely and distinguishes it from the earlier Lean-checked theorem.
2. **Regression coefficients versus predictions.** The original proof outline wrote `Xhat-Xstar=(Pi Q)^+Pi E`, which omits the selected-column coordinate factor. The revised outline gives the correctly factored prediction identity `F_J Xhat=QQ^T F+QG^{-1}Q^T Pi^T Pi E`; orthogonality then yields the stated excess-error expression.
3. **Advertised row count.** The sketched-RP statement allowed arbitrarily small `eps0 <= eta/2` but advertised a row count proportional to `eta^-2`. The asymptotic specialization now explicitly takes `eps0=eta/2`. It also states the zero-residual completion convention used by the imported theorem.
4. **Block rank claim.** The original remark claimed roughly `n/b` blocks already give a zero residual for the identity. Finite CountSketch matrices can be rank deficient at that budget. The revision replaces the claim with the exact expected occupied-bucket calculation for one block.
5. **Estimation versus validation.** Conditional unbiasedness was described as a validated guarantee. The revised discussion states what is actually proved: a conditionally unbiased estimator and variance control. A confidence certificate requires an additional sampling-error or concentration bound.
6. **Numerical descriptions.** The introductory diagnostic paragraph now distinguishes 20,000–40,000-trial Nyström/diagonal experiments from 150-trial sketched-RP experiments. The matrix `G diag(0.7^j) G^T/30` is labeled by its construction instead of implying its eigenvalues are exactly `0.7^j`.
7. **Certificate robustness.** The published Julia listing now uses exact integer bisection for the final integer bound, removing the floating-point starting guess. The independent Python checker also performs the final comparison entirely with integers and fractions.
8. **Determinant moment order.** The sketched-RP outline used a determinant second moment of degree `2t` with only `q>=2` in the hypotheses. The determinant route requires `q>=t`; the statement now assumes `q>=max{2,k}`. The advertised order `q=O(k+log(1/(eps eta_A)))` is unchanged. Its fixed-state event is aligned with the support paper's weighted fourth-moment argument, with two weight families sufficient for the basic trace recurrence. Their failure bound is `(k+1)(4^-q+2*2^-q) <= (3k+4)2^-q`. The sharper hybrid likelihood argument in the support paper does not require `q>=k`.
9. **Safe radical rounding.** The certificate proof now explicitly rewrites `G_n=(2 sqrt((n-1)/n)+1) sqrt(1/n)` before rounding radicals upward. This matches the implementation and avoids inadvertently rounding a denominator upward.

## Independent check outcomes

Run: `python checks/independent_audit.py`.

| Check | Outcome |
|---|---|
| Exact interval certificates | 63/63 pass; all 63 bucket/T pairs match the LaTeX table |
| Largest certificate integer | 499428 |
| All seven likelihood-price moment orders | Reproduced: 7, 11, 18, 46, 72, 69, 32 |
| Both integer prescriptions in each parameter row | Reproduced |
| Exhaustive small CountSketch model | 64 matrices, n=3, b=2, s=1 |
| Weighted second-moment formula | 6.71, matching the direct enumeration |
| Diagonal MSE and trace variance | 0.21 and 0.42, matching the formulas |
| One-block trace drift | 2.74998476 >= 0.85222222 |
| Convex potential f(x)=x^2 | Mean decrease 5.02254004 >= 1.59966667 |
| Full-rank RP set laws | 48 complex instances; every pointwise bound passes |
| Flat-spectrum equality | Maximum discrepancy 3.34e-16 |
| Majorant maxima, C=81/69/40/35/30 | 0.31861044 / 0.35203741 / 0.44147198 / 0.47654917 / 0.50460940 |

The RP maximum normalized ratio is `1 + 1.6e-14`, a floating-point equality discrepancy in flat cases, below the assertion tolerance. The largest moment orders in the independent majorant check are q=3000 for C=69 and C=35; log-domain propagation avoids underflow. Finite-grid evidence does not prove the C=35 conjecture.

Original Julia Monte Carlo experiments were inspected and their scripts/recorded outputs included. They were not rerun wholesale during the revision. The exact certificate and analytic parameter counts are the only computational ingredients used as proofs.

## Novelty and attribution check

Primary sources consulted during the September 2026 check:

| Claim/mechanism | Prior source and attribution | Contribution claimed here |
|---|---|---|
| SparseStack/OSNAP construction | Kane-Nelson, [JACM 2014](https://doi.org/10.1145/2559902), [accessible preprint](https://arxiv.org/abs/1012.1577); Nelson-Nguyen, [OSNAP](https://arxiv.org/abs/1211.1002); Camaño-Epperly-Meyer-Tropp, [2508.21189](https://arxiv.org/abs/2508.21189) | Same distribution; explicit sufficient constants improved |
| Optimal dimension with near-optimal sparsity | Chenakkod–Dereziński–Dong, [ICALP 2025](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.ICALP.2025.55) | No new asymptotic rate claimed |
| Target Nelson–Nguyen rates | Earlier SparseStack paper, [2609.02978](https://arxiv.org/abs/2609.02978); Mai–Rao, [2609.22548v1](https://arxiv.org/html/2609.22548v1) | A finite all-order constant certificate relative to the repository framework |
| Whole-sketch second moments | Clarkson–Woodruff, [1207.6365](https://arxiv.org/abs/1207.6365); structured-sketch applications in [2508.21189](https://arxiv.org/pdf/2508.21189) | Elementary deterministic Nyström compression argument plus exact signed-hash specialization |
| Factorial adaptive-versus-volume law | Deshpande–Vempala, [ECCC TR06-042, Proposition 7](https://eccc.weizmann.ac.il/report/2006/042/download), page 5 | Retain the already-used spectral denominator as an explicit price and apply it to post-sketch selection |
| Stochastic trace estimation | Hutchinson, [1989](https://doi.org/10.1080/03610918908812806) | Exact SparseStack variance, not a new trace probing principle |
| Stochastic diagonal estimation | Bekas–Kokiopoulou–Saad, [author preprint](https://www-users.cse.umn.edu/~saad/PDF/umsi-2005-082.pdf), published 2007 | Exact SparseStack diagonal variance and independent sampled-diagonal correction |
| Modern variance-reduced probing | Epperly–Tropp–Webber, [XTrace/XDiag](https://arxiv.org/abs/2301.07825) | No claim to the control-variate or diagonal-probing idea |

Mai–Rao's primary HTML explicitly analyzes the same SparseStack distribution by a combinatorial mean bound followed by concentration and coordinate resampling. That is credited as another proof of the target asymptotic rates. Deshpande–Vempala's primary PDF confirms that the factorial comparison and spectral-tail denominators themselves are prior work. No universal claim that the present formulas or elementary compression observation are first in the literature is made; the most clearly identifiable new result is the exact certificate for the displayed constants.

## Repository dependencies and limits

All mathematical black boxes are stated explicitly. Full proofs may be read in papers included in this repository:

- The framework PDF for the operator majorant and baseline direct constants.
- The Product-Probe Drift companion for the general ratio/drift mechanism.
- The Iterated-Logarithmic RPCholesky background paper for the fractional likelihood bound.
- The Robust-Sketched RPCholesky companion for the path-local theorem.

Bibliography links use actual repository target filenames. Publication-time link validation must check these targets after the new files are uploaded; until then, new GitHub `blob/main` targets need not exist online. The fixed-frame gate has a fixed dictionary, fixed V, full-rank selection, fresh RP randomization given the sketch, and a uniform deterministic price bound. The advertised Nyström expectation is rank blind and cannot be turned into a uniform expected spectral-tail guarantee for finite-support sketches with m<n. The C=35 prescription remains a conjecture.

The revised abstract is concise and foregrounds the constant certificate. Proofs now explain their purpose and separate compression, covariance contraction, change of measure, convex-potential reduction, and probability transfer. The final PDF uses the shared manuscript style. Visual QA covered all pages through contact sheets and full-size critical pages: certificate table, parameter envelope, drift proof, and gate diagram. Final compilation and layout details are recorded below.

## Final build and visual check

The final paper is 29 A4 pages with 27 mm margins, Latin Modern fonts and the shared paragraph/display spacing. Two final MiKTeX `pdflatex -interaction=nonstopmode -halt-on-error` passes completed successfully. The final LaTeX log has no overfull/underfull boxes, unresolved citations/references or rerun warnings. The narrow-page interval figure and numerical table were adjusted to fit A4.

Every page was rendered and reviewed in four contact sheets. Full-size checks covered pages 9, 11 and 13 (parameter envelope, interval certificate and deterministic Nyström/covariance proof), page 17 (block moment formula), pages 20-21 (selection guarantee and parameter table), page 23 (revised imported theorem and weighted good event), and page 25 (Monte Carlo table). The final PDF contains 25 external link annotations; text bounds do not cross the safe page boundary. The repository dependency links are validated against local target names. Publisher access restrictions on the Hutchinson and Kane-Nelson DOI pages are not treated as broken scholarly identifiers; Kane-Nelson's accessible arXiv alternative is included.

Temporary inspection dumps, revision scripts, build logs and rendered QA pages were moved outside the publishable repository folder. The retained checks use repository-relative paths. The independent arithmetic/finite-enumeration results are retained separately from the original Julia Monte Carlo records.
