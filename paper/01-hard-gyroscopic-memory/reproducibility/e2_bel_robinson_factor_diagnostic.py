#!/usr/bin/env python3
"""E2-G1 algebraic/affine factor diagnostic.

Diagnostic only; the analytic proof is in
BEST-E2-BEL-ROBINSON-SUPER-ENERGY-BRIDGE-20261004.md.
"""

from fractions import Fraction


def check(a: Fraction, b: Fraction, c: Fraction) -> None:
    assert c > 0

    tr_t2 = 2 * (a * a + b * b)
    psi_abs2 = a * a + b * b

    # Frozen Ferrando-Saez-style convention:
    # B = 1/4 (CC + *C*C)
    b_norm = psi_abs2
    # No-prefactor convention Bhat = 4 B.
    b_raw = 4 * psi_abs2

    assert psi_abs2 * 2 == tr_t2
    assert b_norm * 2 == tr_t2
    assert b_raw == 2 * tr_t2

    # J integrand conversion factors.
    assert tr_t2 == 2 * b_norm
    assert tr_t2 * 2 == b_raw

    # Positive affine rescaling u' = c u:
    # L' = c L, du' = c du, k' = c^-1 k,
    # T' = c^-2 T, B(k'^4) = c^-4 B(k^4).
    lhs_scale = c**3 * c * c**-4
    assert lhs_scale == 1


def main() -> None:
    cases = [
        (Fraction(0), Fraction(0), Fraction(2)),
        (Fraction(1), Fraction(0), Fraction(3)),
        (Fraction(0), Fraction(5, 3), Fraction(7, 2)),
        (Fraction(2, 5), Fraction(-7, 11), Fraction(13, 4)),
    ]
    for case in cases:
        check(*case)
    print("E2-G1 diagnostic PASS:", len(cases), "cases")


if __name__ == "__main__":
    main()
