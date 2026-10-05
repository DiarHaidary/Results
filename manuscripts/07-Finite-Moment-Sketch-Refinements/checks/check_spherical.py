"""Exact fourth-order spherical-product checks; Python standard library only."""
from fractions import Fraction as F
from pathlib import Path
import json


def mul(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def kron(a, b):
    return [[a[i][j] * b[k][l] for j in range(len(a[0]))
             for l in range(len(b[0]))]
            for i in range(len(a)) for k in range(len(b))]


def trace(a):
    return sum(a[i][i] for i in range(len(a)))


I = [[1, 0], [0, 1]]
X = [[0, 1], [1, 0]]
Z = [[1, 0], [0, -1]]
# The covariance-matched Gaussian coefficients are sqrt(weight)*matrix.
basis = [(kron(a, I), F(1, 2)) for a in [X, Z]]
basis += [(kron(I, a), F(1, 2)) for a in [X, Z]]
basis += [(kron(a, b), F(1, 4)) for a in [X, Z] for b in [X, Z]]
W = [[F(0) for _ in range(4)] for _ in range(4)]
for a, weight in basis:
    aa = mul(a, a)
    for i in range(4):
        for j in range(4):
            W[i][j] += weight * aa[i][j]
assert W == [[F(3 if i == j else 0) for j in range(4)] for i in range(4)]
crossing = sum(wa * wb * trace(mul(mul(mul(a, b), a), b))
               for a, wa in basis for b, wb in basis)
assert crossing == 8
pairing = 2 * trace(mul(W, W)) + crossing
assert pairing == 80
beta = (4 - 1)**4 + 4 - 1
assert beta == 84
eps = F(1, 2)
delta = F(1, 100)
failure = lambda m: F(beta + (m - 1) * pairing, m**3) / eps**4
assert failure(358) <= delta < failure(357)
results = dict(status="passed", covariance_diagonal=3, component_fourth=beta,
               crossing=int(crossing), pairing=int(pairing), sufficient_rows=358,
               failure_upper=float(failure(358)),
               failure_at_previous_integer=float(failure(357)),
               method="Exact rational Pauli covariance basis and integer threshold")
Path(__file__).with_name("spherical_results.json").write_text(
    json.dumps(results, indent=2) + "\n", encoding="utf-8")
print(json.dumps(results, indent=2))
