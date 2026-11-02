"""Checks for section 5 of docs/ng-duality.md. Requires numpy and scipy.

Run: python -B scripts/check_mertens_bridge.py
Zeta values use the Euler-Maclaurin formula in double precision and the
Mertens function is sieved up to 10^7. RH is numerically true at the heights
used, so the RH statements are checked in that range only.
"""

import math

import numpy as np

from check_dual_product import zeta_line
from check_ng_tail import zeros_up_to, zeta_derivative


def local_limit_check():
    """Proposition 17.1 with T = 300: extrapolate sigma -> 1/2."""
    gammas = zeros_up_to(300.0)
    rhos = 0.5 + 1j * gammas
    weights = np.array([1 / abs(r * zeta_derivative(r)) ** 2 for r in rhos])
    target = 2 * weights.sum()
    us, values = [], []
    for sigma in (0.505, 0.51, 0.52):
        width = sigma - 0.5
        t = np.arange(0.0, 300.0, width / 10)
        integrand = 1 / np.abs((sigma + 1j * t) * zeta_line(sigma, t)) ** 2
        values.append((2 * sigma - 1) / math.pi * np.trapezoid(integrand, t))
        us.append(2 * sigma - 1)
    # Quadratic extrapolation in u = 2 sigma - 1 to u = 0.
    coeff = np.polyfit(us, values, 2)
    limit = coeff[-1]
    print("(5.3) local limit, |gamma| < 300:")
    for u, v in zip(us, values):
        print(f"    sigma={0.5 + u / 2:.3f}: {v:.6f}")
    print(f"    extrapolated {limit:.6f} vs sum 1/|rho zeta'(rho)|^2 = {target:.6f}")
    assert abs(limit - target) < 1e-4
    return gammas


def mertens(n_max):
    """M(n) for 0 <= n <= n_max from a Moebius sieve."""
    from check_dual_product import primes_up_to

    mu = np.ones(n_max + 1, dtype=np.int8)
    mu[0] = 0
    for p in primes_up_to(n_max):
        p = int(p)
        mu[p::p] *= -1
        if p * p <= n_max:
            mu[p * p::p * p] = 0
    return np.cumsum(mu.astype(np.int64))


def plancherel_check():
    """(5.4) at sigma = 3/4: int_0^infty phi_sigma^2 dy vs (1/2pi) int |1/(s zeta)|^2 dt."""
    sigma, n_max = 0.75, 10 ** 7
    m = mertens(n_max).astype(float)
    n = np.arange(1, n_max, dtype=float)
    lhs = np.sum(m[1:n_max] ** 2 * (n ** (-2 * sigma) - (n + 1) ** (-2 * sigma))) / (2 * sigma)
    t = np.arange(0.0, 3000.0, 0.01)
    integrand = 1 / np.abs((sigma + 1j * t) * zeta_line(sigma, t)) ** 2
    dual_mean = mean_inverse_square(sigma)
    rhs = np.trapezoid(integrand, t) / math.pi + dual_mean / (math.pi * 3000.0)
    print(f"\n(5.4) Plancherel at sigma={sigma}: int phi^2 dy = {lhs:.5f} (M up to 1e7), "
          f"(1/2pi) int dt/|s zeta|^2 = {rhs:.5f}")
    assert abs(lhs - rhs) < 0.01 * rhs


def mean_inverse_square(sigma):
    """sum mu(n)^2 n^{-2 sigma} = zeta(2 sigma) / zeta(4 sigma)."""
    return float(zeta_line(2 * sigma, [0.0])[0].real / zeta_line(4 * sigma, [0.0])[0].real)


def dual_mean_check():
    """(5.5): mean of |zeta(sigma+it)|^{-2} against central atoms of K^{-1,-1}."""
    sigma = 0.75
    t = np.linspace(10.0, 3010.0, 300001)
    empirical = np.mean(1 / np.abs(zeta_line(sigma, t)) ** 2)
    target = mean_inverse_square(sigma)
    print(f"\n(5.5) mean of |zeta({sigma}+it)|^-2 on [10, 3010]: {empirical:.4f}, "
          f"zeta(2s)/zeta(4s) = {target:.4f}")
    assert abs(empirical - target) < 0.05 * target


def gonek_check(gammas):
    inverse = np.array([1 / abs(zeta_derivative(0.5 + 1j * g)) ** 2 for g in gammas])
    print("\n(5.6) J_{-1}(T) / (3T/pi^3):")
    previous = math.inf
    for big_t in (100, 200, 300):
        ratio = inverse[gammas < big_t].sum() / (3 * big_t / math.pi ** 3)
        print(f"    T={big_t}: {ratio:.3f}")
        assert abs(ratio - 1) < previous
        previous = abs(ratio - 1)


def main():
    gammas = local_limit_check()
    plancherel_check()
    dual_mean_check()
    gonek_check(gammas)


if __name__ == "__main__":
    main()
