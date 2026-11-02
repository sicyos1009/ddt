"""Checks for section 3 of docs/global-formula.md. Requires mpmath.

Run: python -B scripts/check_desmoothed_formula.py
Tests the lattice smoothing defect and a contour whose upper boundary cuts
through a zero pair. Numerical quadrature is not interval-certified.
"""

import mpmath as mp

from check_smoothed_formula import normal_cdf


def check_lattice_defect():
    print("Half-integer cutoff, h = 1/(X log X):")
    for n in (10, 100, 1000):
        x = mp.mpf(n) + mp.mpf("0.5")
        length = mp.log(x)
        h = 1 / (x * length)
        defect = mp.fsum(
            normal_cdf(-abs(mp.log(x / k)) / h)
            for k in range(1, 2 * n + 2)
        )
        far_tail = x * mp.exp(h * h / 2) * normal_cdf(h - mp.log(2) / h)
        # Gaussian bound for the near-lattice series in (3.2).
        analytic_bound = (
            2 * mp.exp(-length * length / 32)
            / (1 - mp.exp(-length * length / 8))
            + far_tail
        )
        assert defect + far_tail <= analytic_bound
        print(
            f"  X={mp.nstr(x, 8)}, total weight defect <= "
            f"{mp.nstr(defect + far_tail, 12)}"
        )
    assert defect + far_tail < 1 / x


def check_boundary_component():
    x, h, alpha = mp.mpf("10.5"), mp.mpf("0.2"), mp.mpf("0.1")
    rho = mp.zetazero(1)
    zeros = (rho, mp.conj(rho))
    beta, gamma = mp.re(rho), mp.im(rho)
    height = gamma + alpha / 2
    c, left = mp.mpf(2), mp.mpf(-1)
    max_regular_argument = mp.mpf(0)

    def continued_F(s):
        nonlocal max_regular_argument
        shifted = s - 1j * alpha
        zs, zt = mp.zeta(s), mp.zeta(shifted)
        # Remove the pole and the first conjugate pair. This factor has no
        # zero or pole in the region used by this check.
        regular_s = zs * (s - 1) / mp.fprod(s - zero for zero in zeros)
        regular_t = zt * (shifted - 1) / mp.fprod(
            shifted - zero for zero in zeros
        )
        ratio = regular_s / regular_t
        argument = abs(mp.arg(ratio))
        max_regular_argument = max(max_regular_argument, argument)
        assert argument < mp.pi / 2
        value = mp.sqrt(ratio) * mp.sqrt((shifted - 1) / (s - 1)) / zs
        # These rational factors put the cuts on the prescribed vertical
        # segments, unlike the principal square root of the whole zeta ratio.
        for zero in zeros:
            value *= mp.sqrt((s - zero) / (shifted - zero))
        return value

    reference = mp.mpc(2, "2.3")
    euler_value = mp.exp(
        -(mp.log(mp.zeta(reference)) + mp.log(mp.zeta(reference - 1j * alpha)))
        / 2
    )
    tolerance = mp.mpf("1e-22")
    assert abs(continued_F(reference) - euler_value) < tolerance

    def integrand(s):
        return mp.exp(s * mp.log(x) + h * h * s * s / 2) * continued_F(s) / s

    def segment(start, end):
        return (end - start) * mp.quad(
            lambda u: integrand(start + (end - start) * u),
            [mp.mpf(k) / 8 for k in range(9)],
        ) / (2j * mp.pi)

    def circle(center, radius, points=256):
        terms = []
        for k in range(points):
            direction = mp.exp(2j * mp.pi * k / points)
            terms.append(direction * integrand(center + radius * direction))
        return radius * mp.fsum(terms) / points

    cut_at_top = beta + 1j * height
    right = segment(c - 1j * height, c + 1j * height)
    boundary = (
        segment(c - 1j * height, left - 1j * height)
        + segment(left - 1j * height, left + 1j * height)
        # One-sided integrals at the cut crossing are essential.
        + segment(left + 1j * height, cut_at_top)
        + segment(cut_at_top, c + 1j * height)
    )
    pole = circle(1 + 1j * alpha / 2, mp.mpf("0.3"))
    # Both members of the negative pair are inside the outer rectangle.
    # Its circle must fit above the lower boundary.
    negative_pair = circle(mp.conj(rho) + 1j * alpha / 2, mp.mpf("0.08"))

    # The positive pair straddles the upper boundary. Evaluate its open
    # Hankel contribution by a small rectangle with a split top edge.
    radius = mp.mpf("0.3")
    top_left = beta - radius + 1j * height
    top_right = beta + radius + 1j * height
    bottom_left = beta - radius + 1j * (gamma - radius)
    bottom_right = beta + radius + 1j * (gamma - radius)
    open_component = (
        segment(top_left, bottom_left)
        + segment(bottom_left, bottom_right)
        + segment(bottom_right, top_right)
        + segment(top_right, cut_at_top)
        + segment(cut_at_top, top_left)
    )
    without_open_component = right - boundary - continued_F(0) - pole - negative_pair
    error = abs(without_open_component - open_component)
    print("Boundary lies strictly between rho and rho+i*alpha.")
    print("Boundary-component identity error:", mp.nstr(error, 12))
    print("Error if the boundary component is omitted:", mp.nstr(abs(without_open_component), 12))
    print("Largest sampled regular-factor argument:", mp.nstr(max_regular_argument, 8))
    assert error < tolerance
    assert abs(without_open_component) > mp.mpf("1e-5")

    offset = mp.mpf("1e-9")
    f_right = continued_F(cut_at_top + offset)
    f_left = continued_F(cut_at_top - offset)
    jump_error = abs((f_left + f_right) / f_right)
    print("Opposite boundary values across the vertical cut:", mp.nstr(jump_error, 12))
    assert jump_error < mp.mpf("1e-6")


if __name__ == "__main__":
    mp.mp.dps = 30
    check_lattice_defect()
    check_boundary_component()
    print("All desmoothing checks passed.")
