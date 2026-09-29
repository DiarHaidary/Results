# Reproducing the numerical illustrations

These files accompany the manuscript's numerical section. They are
finite checks and illustrations, not formal proof verification.

The independent checker requires Python 3 and NumPy. Run it from any
directory:

```powershell
python checks/independent_audit.py
```

It writes `independent_audit_results.json` next to the script. The stored
results were obtained with the script's deterministic seed. Its checks
cover real and complex residual/refit identities, fixed-state functional
bounds, the count tables, finite frame/union enumeration, and floor
construction edge cases.

The Julia scripts require only Julia's standard libraries. They were
rerun with Julia 1.13.0. From the manuscript directory, for example:

```powershell
julia checks/sk_constants.jl
julia checks/det01_fractional.jl
julia checks/sketched_rp.jl
julia checks/sketched_path_quantities2.jl
```

| Script | Illustration |
| --- | --- |
| `identities.jl` | One residual, path-likelihood, and refit identity example |
| `sk_constants.jl` | Robust pivot counts and sufficient sketch sizes |
| `trichotomy_curves.jl` | Coordinates used in the trichotomy plots |
| `trichotomy_dp.jl` | Additional diagonal-family dynamic-program checks |
| `det01_fractional.jl` | Exhaustive finite-history enumeration at 4096-bit precision |
| `editor_crossover.jl` | Renyi/factorial gate comparisons |
| `sketched_rp.jl` | 150-trial error and reconstruction illustration |
| `sketched_path_quantities2.jl` | 100-trial path-quantity illustration |

The `.out` files retain the supplied outputs. The two `_review.out`
files contain the audit reruns of the Monte Carlo scripts; the displayed
table entries agree to their printed precision. Numerical agreement does
not establish a uniform conjecture or an adaptive Gaussian expectation
identity.
