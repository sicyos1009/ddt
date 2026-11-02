"""Checks for section 1 of docs/ng-duality.md. Requires numpy and scipy.

Run: python -B scripts/check_ng_tail.py
Zeros and zeta'(rho) are computed in double precision by the Euler-Maclaurin
formula and Hardy's Z-function. These are consistency checks of the constants
in Lemma 1.3 on low zeros, not a proof.
"""

import math

import numpy as np
from scipy.optimize import brentq
from scipy.special import loggamma

# B_2, B_4, ..., B_20
BERNOULLI = [1 / 6, -1 / 30, 1 / 42, -1 / 30, 5 / 66, -691 / 2730, 7 / 6,
             -3617 / 510, 43867 / 798, -174611 / 330]


def zeta(s):
    """Euler-Maclaurin evaluation of zeta(s) for complex s with moderate |Im s|."""
    s = complex(s)
    n_terms = int(abs(s.imag) / (2 * math.pi)) + 30
    n = np.arange(1, n_terms, dtype=float)
    total = np.sum(np.exp(-s * np.log(n)))
    big_n = float(n_terms)
    total += big_n ** (1 - s) / (s - 1) + 0.5 * big_n ** (-s)
    rising = s
    power = big_n ** (-s - 1)
    for k, b in enumerate(BERNOULLI, start=1):
        total += b / math.factorial(2 * k) * rising * power
        rising *= (s + 2 * k - 1) * (s + 2 * k)
        power /= big_n * big_n
    return total


def hardy_z(t):
    theta = (loggamma(0.25 + 0.5j * t)).imag - 0.5 * t * math.log(math.pi)
    return (np.exp(1j * theta) * zeta(0.5 + 1j * t)).real


def zeros_up_to(height, step=0.02):
    grid = np.arange(10.0, height, step)
    values = np.array([hardy_z(t) for t in grid])
    roots = []
    for i in np.nonzero(np.sign(values[:-1]) != np.sign(values[1:]))[0]:
        roots.append(brentq(hardy_z, grid[i], grid[i + 1], xtol=1e-13))
    return np.array(roots)


def zeta_derivative(s, radius=1e-3, points=16):
    """Cauchy-integral derivative on a small circle."""
    angles = 2 * math.pi * np.arange(points) / points
    nodes = s + radius * np.exp(1j * angles)
    return np.mean([zeta(z) * np.exp(-1j * a) for z, a in zip(nodes, angles)]) / radius


def gallagher_sides(freqs, coeffs, big_y, grid_points=400_000):
    """Return the three sides of (1.4) for S(y)=sum c e^{i lambda y}."""
    y = np.linspace(-big_y, big_y, grid_points)
    s_vals = np.zeros_like(y, dtype=complex)
    for lam, c in zip(freqs, coeffs):
        s_vals += c * np.exp(1j * lam * y)
    lhs = np.trapezoid(np.abs(s_vals) ** 2, y)

    delta = math.pi / big_y
    # A(u)=sum_{u<lambda<=u+delta} c_lambda is constant between breakpoints.
    breaks = np.unique(np.concatenate([freqs, freqs - delta]))
    integral = 0.0
    for left, right in zip(breaks[:-1], breaks[1:]):
        mid = 0.5 * (left + right)
        inside = (freqs > mid) & (freqs <= mid + delta)
        integral += abs(np.sum(coeffs[inside])) ** 2 * (right - left)
    middle = math.pi * big_y ** 2 / 2 * integral

    counts = np.array([np.sum(np.abs(freqs - lam) <= delta) for lam in freqs])
    rhs = math.pi ** 2 * big_y / 2 * np.sum(np.abs(coeffs) ** 2 * counts)
    diagonal = 2 * big_y * np.sum(np.abs(coeffs) ** 2)
    return lhs, middle, rhs, diagonal, counts.max()


def main():
    gammas = zeros_up_to(300.0)
    print(f"zeros found below 300: {len(gammas)} (expected 138)")
    assert len(gammas) == 138
    rhos = 0.5 + 1j * gammas
    derivs = np.array([zeta_derivative(r) for r in rhos])
    residual = max(abs(zeta(r)) for r in rhos[:20])
    print(f"max |zeta(rho)| on first 20 zeros: {residual:.1e}")
    print(f"zeta'(rho_1) = {derivs[0]:.10f}  (reference 0.7832965119+0.1246998297i)")
    assert abs(derivs[0] - (0.7832965119 + 0.1246998297j)) < 1e-8

    a = 1 / (rhos * derivs)
    weights = np.log(gammas + 2) * np.abs(a) ** 2
    print("\npartial sums over 0<gamma<T (the two-sided sum is twice this):")
    for cut in (50, 100, 200, 300):
        mask = gammas < cut
        print(f"  T={cut:3d}: sum|a|^2={np.sum(np.abs(a[mask])**2):.6f}, "
              f"sum log(g+2)|a|^2={np.sum(weights[mask]):.6f}")

    # Lemma 1.3 on the tail P_T - P_{T0}, with both signs of gamma.
    t0 = 30.0
    tail = gammas >= t0
    freqs = np.concatenate([-gammas[tail][::-1], gammas[tail]])
    coeffs = np.concatenate([np.conj(a[tail])[::-1], a[tail]])
    print(f"\nLemma 1.3 on the tail {t0:.0f}<=|gamma|<300:")
    for big_y in (8.0, 20.0, 60.0):
        lhs, middle, rhs, diagonal, max_count = gallagher_sides(freqs, coeffs, big_y)
        ok = lhs <= middle * (1 + 1e-6) and middle <= rhs * (1 + 1e-12)
        print(f"  Y={big_y:4.0f}: int|S|^2={lhs:.4e} <= {middle:.4e} <= {rhs:.4e}"
              f"  [2Y*sum|c|^2={diagonal:.4e}, max N={max_count}]  {'ok' if ok else 'FAIL'}")
        assert ok


if __name__ == "__main__":
    main()
