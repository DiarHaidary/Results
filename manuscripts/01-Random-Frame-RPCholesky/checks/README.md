# Reproducing the checks

`independent_audit.py` is a new, independent targeted check of the current
manuscript. It requires Python, NumPy, SciPy, and SymPy:

```powershell
python checks/independent_audit.py
```

Its seed and complete results are in `independent_audit_results.json`. It runs
eight 100,000-path spectral-chain experiments, reproduces the eight budget-two
Laplace-quadrature entries, checks the nonmonotone-scale covariance example,
computes four coordinate examples with exact rational arithmetic, and checks
the finite Bartlett reciprocal-product coefficients through rank six. A
floating-point quadrature error estimate is a diagnostic, not a rigorous
interval certificate. These checks supplement the proofs; they do not prove
the general theorems.

The `.jl` scripts and `.txt` outputs are the original experiment materials.
The large-run tables were retained from those experiments. The full Julia
runs were not repeated during this revision. The main entry points are:

| Experiment | Julia script | Saved output |
| --- | --- | --- |
| Main real/complex Haar table | `tables_final.jl` | `tables_final.txt` |
| Shape interpolation | `alpha_family.jl` | `alpha_family.txt` |
| Coordinate sharp examples | `det01_coordinate.jl` | `det01_coordinate.txt` |
| Prefix asymptotics | `det01_prefix.jl` | `det01_prefix.txt` |
| Rotated sharp spectra | `det01_rotated.jl` | `det01_rotated_eta1e-8.txt` |
| Population chain | `population.jl` | `population.txt` |
| Tight-frame identities | `dilation.jl` | `dilation.txt` |

For example, from the manuscript directory, run
`julia --threads=auto checks/tables_final.jl`. Monte Carlo output may vary
with the number of threads because the original scripts divide paths among
thread-specific random generators. Entries marked “exact” in the original
floating-point table-generation output mean numerical quadrature; the revised
paper labels them “quadrature.” Exact rational computations are separately
identified in the paper and scripts.
