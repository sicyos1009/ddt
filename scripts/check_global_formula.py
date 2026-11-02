"""Checks for docs/global-formula.md. Requires mpmath.

Run: python scripts/check_global_formula.py
These calculations check identities and normalizations, not the zero list's
completeness or the analytic error estimates in the document.
"""

from fractions import Fraction
from math import factorial

import mpmath as mp


def factor(n):
    result = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            result[p] = result.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        result[n] = result.get(n, 0) + 1
    return result


def prime_coefficient(z, a):
    result = Fraction(1)
    for j in range(a):
        result *= z + j
    return result / factorial(a)


def divisor_coefficient(z, n):
    result = Fraction(1)
    for a in factor(n).values():
        result *= prime_coefficient(z, a)
    return result


def check_fta():
    cases = [
        (Fraction(1, 2), Fraction(1, 2)),
        (Fraction(-1, 2), Fraction(-1, 2)),
        (Fraction(-3, 2), Fraction(1, 3)),
        (Fraction(0), Fraction(2)),
    ]
    count = 0
    for z, w in cases:
        for n in (4, 12, 30, 72, 120, 210):
            # A key d represents the location log(d) before scaling.
            convolution = {1: Fraction(1)}
            for p, a in factor(n).items():
                local = {
                    p**j: prime_coefficient(z, j)
                    * prime_coefficient(w, a - j)
                    for j in range(a + 1)
                }
                next_convolution = {}
                for d, value in convolution.items():
                    for e, weight in local.items():
                        key = d * e
                        next_convolution[key] = (
                            next_convolution.get(key, Fraction(0))
                            + value * weight
                        )
                convolution = next_convolution
            direct = {
                d: divisor_coefficient(z, d) * divisor_coefficient(w, n // d)
                for d in range(1, n + 1)
                if n % d == 0
            }
            assert convolution == direct
            assert sum(direct.values()) == divisor_coefficient(z + w, n)
            count += 1
    print(f"FTA convolution and exact masses: {count} cases passed.")


def ell(s):
    b0 = mp.log(2 * mp.pi) - 1 - mp.euler / 2
    return mp.log(s - 1) + mp.loggamma(1 + s / 2) - b0 * (s - 2)


def psi(rho, s):
    return mp.log(s - rho) - mp.log(2 - rho) + (s - 2) / rho


def check_zero_product():
    # Re(s)=2 puts zeta(s) in the disk |zeta(s)-1| <= zeta(2)-1 < 1.
    # Thus the principal logarithm here is exactly the Euler logarithm.
    s = mp.mpc(2, "2.3")
    alpha = mp.mpf("0.4")
    direct = mp.exp(
        -(mp.log(mp.zeta(s)) + mp.log(mp.zeta(s - 1j * alpha))) / 2
    )
    analytic_factor = mp.exp((ell(s) + ell(s - 1j * alpha)) / 2) / mp.zeta(2)
    log_product = mp.mpc(0)
    first_zero = None
    errors = []
    print("Zero-product check at s=2+2.3i, alpha=0.4 (conjugates included):")
    for k in range(1, 81):
        rho = mp.zetazero(k)
        if first_zero is None:
            first_zero = rho
        for zero in (rho, mp.conj(rho)):
            log_product -= (psi(zero, s) + psi(zero, s - 1j * alpha)) / 2
        if k in (10, 30, 80):
            approximation = analytic_factor * mp.exp(log_product)
            error = abs(approximation / direct - 1)
            errors.append(error)
            print(
                f"  {k:2d} positive zeros, height={mp.nstr(mp.im(rho), 10)}, "
                f"relative error={mp.nstr(error, 12)}"
            )
    # Sanity checks at these particular inputs; no monotonicity theorem is used.
    assert errors[-1] < errors[0]
    assert errors[-1] < mp.mpf("0.05")
    assert abs(ell(2)) < mp.mpf("1e-25")
    assert abs(psi(first_zero, 2)) < mp.mpf("1e-25")
    return first_zero


def check_contour_signs(rho):
    # With alpha=0 there are no square-root branch choices:
    # F_0(s)=1/zeta(s). The rectangle contains s=0 and the first zero pair.
    x, c, a, height = mp.mpf("10.5"), mp.mpf(2), mp.mpf(-1), mp.mpf(20)

    def integrand(s):
        return mp.power(x, s) / (s * mp.zeta(s))

    def segment(start, end):
        return (end - start) * mp.quad(
            lambda u: integrand(start + (end - start) * u),
            [mp.mpf(k) / 8 for k in range(9)],
        ) / (2j * mp.pi)

    right = segment(c - 1j * height, c + 1j * height)
    remainder = (
        segment(c - 1j * height, a - 1j * height)
        + segment(a - 1j * height, a + 1j * height)
        + segment(a + 1j * height, c + 1j * height)
    )
    residue = mp.power(x, rho) / (rho * mp.diff(mp.zeta, rho))
    residue_sum = -2 + 2 * mp.re(residue)
    error = abs(right - remainder - residue_sum)
    print("Contour orientation and s=0 residue error:", mp.nstr(error, 12))
    assert error < mp.mpf("1e-18")


if __name__ == "__main__":
    mp.mp.dps = 30
    check_fta()
    first_zero = check_zero_product()
    check_contour_signs(first_zero)
    print("All checks passed.")
