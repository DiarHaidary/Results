"""Independent numerical audit diagnostics; these checks are not proofs.

Run with Python, NumPy and SciPy. Writes audit_checks.json beside this file.
"""

from pathlib import Path
import json
import math
import numpy as np
from scipy.integrate import quad
from scipy.special import erf

rng = np.random.default_rng(20260930)
report = {"seed": 20260930, "scope": "numerical identities and constant checks"}

# The added weighted-inverse proof must cover odd as well as even budgets.
constant = 112 / (3 * math.pi) - 4 / 3
coefficients = []
for r in range(1, 101):
    for k in range(4 * (r + 1), 4 * (r + 1) + 100):
        s = (k + 1) // 2
        beta_inverse = k * (k - 1) / ((k - s) * (k - s - 1))
        value = k / (s - r - 1) * (2 * beta_inverse / math.pi - 1 / 3)
        assert value <= constant + 1e-12
        coefficients.append((value, r, k))
report["weighted_order_statistic"] = {
    "grid_size": len(coefficients),
    "maximum": max(coefficients),
    "bound": constant,
}
v = np.geomspace(1e-10, 1e3, 10000)
cdf_squared = erf(np.sqrt(v / 2)) ** 2
radial_rhs = 2 * v / (math.pi * (1 + v / 3))
assert np.all(cdf_squared <= radial_rhs + 1e-14)
report["radial_cdf_inequality"] = {"grid_size": len(v), "passed": True}

# Independently integrate the conditional Cauchy formula in the counterexample.
quad_results = []
for y in [1.0, 0.5, 0.1, 0.01, 0.0001]:
    integral, uncertainty = quad(
        lambda a: 2 / (math.pi * (1 + a * a)) * (a + y) / (1 + a * y),
        0, np.inf, epsabs=1e-11, epsrel=1e-11,
    )
    formula = (2 * y + 2 / math.pi * (1 - y * y) * math.log(1 / y)) / (1 + y * y)
    assert abs(integral - formula) < 1e-9
    quad_results.append({"y": y, "integral": integral, "formula": formula,
                         "reported_quad_error": uncertainty})
report["head_only_quadrature"] = quad_results

# Algebraic support, one-response calibration, one-query recovery, and Nystrom.
n1, n2, rank, ell, width = 18, 4, 6, 2, 4
weights = np.array([0.7, 0.3])
E = np.linalg.qr(rng.normal(size=(n1, n1)))[0]
F = np.linalg.qr(rng.normal(size=(n2, n2)))[0]
V = np.column_stack([
    sum(math.sqrt(weights[b]) * np.kron(E[:, a * ell + b], F[:, b])
        for b in range(ell))
    for a in range(rank)
])
A = (V * np.array([10, 5, 1, 0.5, 0.3, 0.2])) @ V.T
h = rng.normal(size=n2)
L = np.kron(np.eye(n1), h[:, None])
w = sum(weights[b] * (F[:, b] @ h) ** 2 for b in range(ell))
covariance_error = np.linalg.norm(V.T @ L @ L.T @ V - w * np.eye(rank))
supported = A @ np.kron(rng.normal(size=n1), rng.normal(size=n2))
calibration_error = abs(np.linalg.norm(L.T @ supported) ** 2
                        / np.linalg.norm(supported) ** 2 - w)
product_error = np.linalg.norm(A @ supported - A @ L @ (L.T @ supported / w))
B = L.T @ A @ L / w
spectrum_error = np.max(np.abs(np.linalg.eigvalsh(B)[-rank:]
                               - np.linalg.eigvalsh(A)[-rank:]))
recovery = A @ L @ np.linalg.pinv(L.T @ A @ L) @ L.T @ A
recovery_error = np.linalg.norm(recovery - A) / np.linalg.norm(A)
assert max(covariance_error, calibration_error, product_error,
           spectrum_error, recovery_error) < 1e-10

def nystrom(probes):
    Y = A @ probes
    return Y @ np.linalg.pinv(probes.T @ Y) @ Y.T

environments = rng.normal(size=(n2, width))
probes = np.column_stack([np.kron(rng.normal(size=n1), environments[:, j])
                          for j in range(width)])
scales = np.array([sum(weights[b] * (F[:, b] @ environments[:, j]) ** 2
                           for b in range(ell)) for j in range(width)])
G = (V.T @ probes) / np.sqrt(scales)
coupling_error = np.linalg.norm(nystrom(probes) - nystrom(V @ G)) / np.linalg.norm(A)
assert coupling_error < 1e-10
report["matrix_identities"] = dict(
    covariance_error=float(covariance_error), calibration_error=float(calibration_error),
    product_error=float(product_error), spectrum_error=float(spectrum_error),
    recovery_relative_error=float(recovery_error),
    nystrom_coupling_relative_error=float(coupling_error),
)

# The real Hurwitz certificate fails as a complex scalar-covariance certificate.
M1 = np.eye(2) / math.sqrt(2)
M2 = np.array([[0.0, -1.0], [1.0, 0.0]]) / math.sqrt(2)
hc = np.array([1.0, 1.0j]) / math.sqrt(2)
complex_cross = np.vdot(hc, M1.T @ M2 @ hc)
assert abs(complex_cross) > 0.49
report["complex_hurwitz_countercheck"] = {
    "cross_covariance_real": float(complex_cross.real),
    "cross_covariance_imaginary": float(complex_cross.imag),
}

out = Path(__file__).with_name("audit_checks.json")
out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
print(f"Passed independent audit checks; wrote {out.name}")
