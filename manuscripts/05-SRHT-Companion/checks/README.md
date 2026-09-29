# Reproducing the SRHT checks

From `manuscripts/05-SRHT-Companion`:

```text
python checks/audit_srht.py
```

This requires `mpmath` and `numpy`; it writes `checks/audit_results.json`.
It checks scalar certificate boundaries at 80 decimal digits and exact
rational n=4 Gram/covariance identities for explicit test matrices. The
Nyström residual comparisons use numerical pseudoinverses.

The supplied original Julia diagnostics are in `checks/original/`. Run them with
that directory as the working directory, because `gate.jl` and several
other scripts read `prescriptions.jl` relatively:

```text
cd checks/original
julia --startup-file=no prescriptions.jl
julia --startup-file=no minscan.jl
julia --startup-file=no gate.jl
julia --startup-file=no scalars.jl
julia --startup-file=no wholesketch.jl
julia --startup-file=no oneround_diag.jl
julia --startup-file=no sketched_rp.jl
```

The `rerun_*.txt` files record the actual revision-time reruns. The width
scan uses Float64 and a stated 1e-9 log-certificate tolerance; adjacent
thresholds are checked separately in high precision. `wholesketch.jl`
enumerates the entire real Walsh law at n=4 in exact rational arithmetic
for moment identities, and at n=8 in floating point; its Nyström residuals
are numerical. These diagnostics supplement the analytic proofs and do not
formally verify the imported operator estimates.
