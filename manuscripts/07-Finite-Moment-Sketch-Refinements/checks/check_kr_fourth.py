from fractions import Fraction as F
from math import isqrt
import json
from pathlib import Path

B = 138612
eps4 = F(1, 16)
delta = F(1, 100)
scale = 10**18
root_floor = isqrt(97 * scale * scale)
root_lo = F(root_floor, scale)
root_hi = F(root_floor + 1, scale)

cases = [
    ('generic_gaussian_envelope', F(2700), F(2700), 2104),
    ('certified_kappa_envelope', 1800 + 30*(9+root_lo), 1800 + 30*(9+root_hi), 1974),
    ('exact_covariance_pairing', F(2252), F(2252), 1928),
]
rows = []
for name, lower, upper, m in cases:
    bound = (B + (m-1)*upper) / m**3 / eps4
    previous_lower = (B + (m-2)*lower) / (m-1)**3 / eps4
    assert bound <= delta
    assert previous_lower > delta
    rows.append(dict(certificate=name, m=m,
                     rational_verification_passed=True,
                     failure_upper=float(bound),
                     failure_lower_at_m_minus_one=float(previous_lower)))
print(json.dumps(rows, indent=2))
Path(__file__).with_name('kr_fourth_results.json').write_text(
    json.dumps(rows, indent=2) + '\n', encoding='utf-8')
