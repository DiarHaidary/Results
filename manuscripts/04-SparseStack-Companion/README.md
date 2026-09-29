# SparseStack companion

[Read the paper](SparseStack_Constants_Nystrom_Adaptive.pdf) · [LaTeX source](SparseStack_Constants_Nystrom_Adaptive.tex) · [Audit](audit.md)

The main result is an exact certificate improving the sufficient universal sparsity/row constants from `(26.5, 179)` to `(25, 154)`, with alternative pairs `(23.5, 163)` and `(22, 172)`. It leaves the SparseStack distribution and the imported operator estimates unchanged. Separate results give a whole-sketch Nyström error bound, CountSketch block drift and exact estimator variances, and an application of the classical adaptive-versus-volume sampling comparison to sketch-dependent frame selection.

The operator dependency is [Theorem 10.6 and Appendix B of the repository framework](../../Graded%20Multiplication%20Operators%20for%20Random-Matrix%20Moments.pdf). The earlier [SparseStack.pdf](../../SparseStack.pdf) uses a coupled proof and different constants. The sharper certificate here is an analytic result with an exact arithmetic check; it is not covered by that earlier paper's Lean verification.

From this directory, build with a current LaTeX distribution:

```text
pdflatex -interaction=nonstopmode -halt-on-error SparseStack_Constants_Nystrom_Adaptive.tex
pdflatex -interaction=nonstopmode -halt-on-error SparseStack_Constants_Nystrom_Adaptive.tex
```

The source includes the shared `../manuscript-style.tex` file. All scholarly dependencies are stated in the paper and linked in its bibliography.

Independent checks (Python and NumPy):

```text
python checks/independent_audit.py
```

The exact 63-case certificate uses only Python standard-library arithmetic. The same script compares every certificate number against the printed LaTeX table; checks the integer parameter prescriptions and likelihood prices; enumerates small sketches and RP pivot paths; and reproduces the numerical majorant maxima. Its results are in [independent_audit_results.json](checks/independent_audit_results.json).

The Julia scripts and `.out` files in `checks/` preserve the original numerical experiments. Run them from `checks/`, because `headroom.jl` and `headroom2.jl` read `diag.jl`. The appendix's revised certificate checker is also supplied as `checks/certificate_listing.jl`; its integer bisection avoids floating point entirely. Recorded Monte Carlo outputs are diagnostics, not mathematical proofs or new simulation runs performed during this revision.
