"""Numerical checks for the Gaussian formula in docs/global-formula.md, section 12.

Run: python scripts/check_smoothed_formula.py
Requires mpmath. These are numerical consistency checks, not interval-certified
quadrature, a certification of the zero list, or a proof of uniform error bounds.
"""

from fractions import Fraction

import mpmath as mp

from check_global_formula import divisor_coefficient


def normal_cdf(t):
    return mp.erfc(-t / mp.sqrt(2)) / 2


def direct_sum(x, h, alpha, cutoff):
    """Compute (12.1), with a bound for the omitted arithmetic tail."""
    weights = [mp.mpf(0)]
    for n in range(1, cutoff + 1):
        value = divisor_coefficient(Fraction(-1, 2), n)
        weights.append(mp.mpf(value.numerator) / value.denominator)
    coefficients = [mp.mpc(0) for _ in range(cutoff + 1)]
    for d in range(1, cutoff + 1):
        twisted = weights[d] * mp.exp(1j * alpha * mp.log(d))
        for m in range(1, cutoff // d + 1):
            coefficients[d * m] += twisted * weights[m]
    total = mp.fsum(
        coefficients[n] * normal_cdf(mp.log(x / n) / h)
        for n in range(1, cutoff + 1)
    )
    # The decreasing Gaussian weight gives sum_{n>N} W(n) <= integral_N^inf W(t)dt.
    # Drop the negative term -N*Phi(-log(N/x)/h) from the exact integral.
    tail_bound = x * mp.exp(h * h / 2) * normal_cdf(
        h - mp.log(cutoff / x) / h
    )
    return total, tail_bound


def check_smoothed_contour():
    x, h, alpha = mp.mpf("10.5"), mp.mpf("0.4"), mp.mpf("0.1")
    c, left, height = mp.mpf(2), mp.mpf(-1), mp.mpf(20)
    rho = mp.zetazero(1)
    max_ratio_argument = mp.mpf(0)

    def continued_F(s, shift):
        nonlocal max_ratio_argument
        zeta_s = mp.zeta(s)
        ratio = zeta_s / mp.zeta(s - 1j * shift)
        argument = abs(mp.arg(ratio))
        max_ratio_argument = max(max_ratio_argument, argument)
        # On the contours used below the ratio stays well away from the
        # principal square root's cut. This realizes continuation from shift=0,
        # where F_0(s)=1/zeta(s), including its negative value at s=0.
        assert argument < mp.pi / 2, "Square-root continuation needs a new path."
        return mp.sqrt(ratio) / zeta_s

    def integrand(s, shift=alpha):
        return mp.exp(s * mp.log(x) + h * h * s * s / 2) * continued_F(
            s, shift
        ) / s

    def segment(start, end):
        return (end - start) * mp.quad(
            lambda u: integrand(start + (end - start) * u),
            [mp.mpf(k) / 8 for k in range(9)],
        ) / (2j * mp.pi)

    def circle(center, shift=alpha, points=128):
        radius = mp.mpf("0.3")
        terms = []
        for k in range(points):
            direction = mp.exp(2j * mp.pi * k / points)
            terms.append(direction * integrand(center + radius * direction, shift))
        return radius * mp.fsum(terms) / points

    direct, arithmetic_tail = direct_sum(x, h, alpha, cutoff=2000)
    right = segment(c - 1j * height, c + 1j * height)
    mellin_tail = (
        mp.power(x, c)
        * mp.exp(h * h * c * c / 2)
        * mp.zeta(c)
        * mp.exp(-h * h * height * height / 2)
        / (mp.pi * h * h * height * height)
    )
    sum_error = abs(direct - right)
    tolerance = mp.mpf("1e-25")
    print("Gaussian divisor sum:", mp.nstr(direct, 14))
    print("Sum vs truncated Mellin error:", mp.nstr(sum_error, 12))
    print("Analytic tail bound:", mp.nstr(arithmetic_tail + mellin_tail, 12))
    assert sum_error <= arithmetic_tail + mellin_tail + tolerance

    # Each circle surrounds both members of its pair. The negative zero's pair
    # is shifted upwards too; its contribution is not the conjugate of the
    # positive zero's contribution at this nonzero alpha.
    centers = [1 + 1j * alpha / 2, rho + 1j * alpha / 2, mp.conj(rho) + 1j * alpha / 2]
    components = [circle(center) for center in centers]
    refined = [circle(center, points=256) for center in centers]
    for coarse, fine in zip(components, refined):
        assert abs(coarse - fine) < tolerance
    boundary = (
        segment(c - 1j * height, left - 1j * height)
        + segment(left - 1j * height, left + 1j * height)
        + segment(left + 1j * height, c + 1j * height)
    )
    origin = continued_F(0, alpha)
    contour_error = abs(right - boundary - origin - mp.fsum(refined))
    print("Nonzero-frequency contour identity error:", mp.nstr(contour_error, 12))
    print("Largest sampled ratio argument:", mp.nstr(max_ratio_argument, 8))
    assert contour_error < tolerance

    weighted_residue = (
        mp.exp(h * h * rho * rho / 2)
        * mp.power(x, rho)
        / (rho * mp.diff(mp.zeta, rho))
    )
    residue_error = abs(circle(rho, shift=mp.mpf(0)) - weighted_residue)
    print("Gaussian-weighted residue error at alpha=0:", mp.nstr(residue_error, 12))
    assert residue_error < tolerance
    assert abs(continued_F(0, mp.mpf(0)) + 2) < tolerance


if __name__ == "__main__":
    mp.mp.dps = 35
    check_smoothed_contour()
    print("All smoothing checks passed.")
