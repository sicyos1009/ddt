"""Checks for docs/mertens-psi.md. Requires numpy and scipy.

Run: python -B scripts/check_mertens_psi.py
Zeta values come from the Euler-Maclaurin formula in double precision. The
local kernels are checked at the first zero and at s = 1 only. Zeros for
section 4 are located from sign changes of Hardy's Z up to height 400. These
are consistency checks, not a proof.
"""

import math
import cmath

import numpy as np
from scipy.special import jv, poch, psi as digamma, rgamma

from check_dual_product import (BERNOULLI, divisor_function, primes_up_to,
                                smallest_prime_factor, zeta_line)
from check_ng_tail import zeros_up_to

RHO_1 = 0.5 + 14.134725141734693j


def zeta(s):
    """Euler-Maclaurin zeta for an array of complex s with moderate |Im s|."""
    s = np.atleast_1d(np.asarray(s, dtype=complex))
    big_n = int(np.abs(s.imag).max() / (2 * math.pi)) + 40
    n = np.arange(1, big_n, dtype=float)
    total = np.exp(-np.outer(s, np.log(n))).sum(axis=1)
    total += big_n ** (1 - s) / (s - 1) + 0.5 * big_n ** (-s)
    rising, power = s.copy(), big_n ** (-s - 1)
    for k, b in enumerate(BERNOULLI, start=1):
        total += b / math.factorial(2 * k) * rising * power
        rising = rising * (s + 2 * k - 1) * (s + 2 * k)
        power = power / big_n ** 2
    return total


def beta_transform(z, w, xi, terms=80):
    """B_{z,w}(xi) from the moments (z)_k / Gamma(z+w+k) of main-theorem (4.2)."""
    return sum((1j * xi) ** k / math.factorial(k) * poch(z, k) * rgamma(z + w + k)
               for k in range(terms))


def check_identities(n_max=20000):
    """(1.1)-(1.2): exact convolution identities with Lambda and Lambda_2."""
    spf = smallest_prime_factor(n_max)
    dm = divisor_function(-0.5, n_max, spf).real.tolist()
    dp = divisor_function(0.5, n_max, spf).real.tolist()
    lam, mu = [0.0] * (n_max + 1), [0] * (n_max + 1)
    mu[1] = 1
    for n in range(2, n_max + 1):
        p, m = int(spf[n]), n // int(spf[n])
        mu[n] = 0 if m % p == 0 else -mu[m]
        while m % p == 0:
            m //= p
        if m == 1:
            lam[n] = math.log(p)
    sums = {key: [0.0] * (n_max + 1) for key in ("mu", "eps", "log", "log2", "lamlam")}
    for d in range(1, n_max + 1):
        log_d = math.log(d)
        for m in range(1, n_max // d + 1):
            n = d * m
            sums["mu"][n] += dm[d] * dm[m]
            sums["eps"][n] += dm[d] * dp[m]
            sums["log"][n] += dm[d] * dp[m] * log_d
            sums["log2"][n] += dm[d] * dp[m] * log_d ** 2
            sums["lamlam"][n] += lam[d] * lam[m]
    errors = [0.0] * 4
    for n in range(1, n_max + 1):
        lam2 = lam[n] * math.log(n) + sums["lamlam"][n]
        errors[0] = max(errors[0], abs(sums["mu"][n] - mu[n]))
        errors[1] = max(errors[1], abs(sums["eps"][n] - (n == 1)))
        errors[2] = max(errors[2], abs(sums["log"][n] + lam[n] / 2))
        errors[3] = max(errors[3], abs(sums["log2"][n]
                                       - (lam2 / 4 - 0.75 * lam[n] * math.log(n))))
    print(f"(1.1)-(1.2), n <= {n_max}: max errors "
          + ", ".join(f"{e:.1e}" for e in errors))
    assert max(errors) < 1e-10


def check_coins():
    """(1.4): squarefree divisor sums split into cosine and sine coins."""
    worst = 0.0
    for primes in ([2, 3, 5], [2, 7, 11, 13], [3, 5, 17, 19, 23]):
        n = math.prod(primes)
        for z, w in ((0.5, 0.5), (-0.5, -0.5), (-0.5, 0.5), (0.3, -0.7)):
            for alpha in (0.4, 1.3):
                lhs = 0j
                for mask in range(1 << len(primes)):
                    d = math.prod(p for i, p in enumerate(primes) if mask >> i & 1)
                    k = bin(mask).count("1")
                    lhs += z ** k * w ** (len(primes) - k) * d ** (1j * alpha)
                rhs = n ** (1j * alpha / 2)
                for p in primes:
                    theta = alpha * math.log(p) / 2
                    rhs *= (z + w) * math.cos(theta) + 1j * (z - w) * math.sin(theta)
                worst = max(worst, abs(lhs - rhs))
    print(f"(1.4) coin factorization: max error {worst:.1e}")
    assert worst < 1e-12


def check_kernels():
    """(2.3)-(2.4): Bessel forms, arcsine representation and tails."""
    phi = (np.arange(4000) + 0.5) * (math.pi / 2) / 4000
    u = np.sin(phi) ** 2  # arcsine-distributed when phi is uniform on [0, pi/2]
    worst = 0.0
    for xi in (0.7, 2.5, 6.0):
        e, j0, j1 = cmath.exp(1j * xi / 2), jv(0, xi / 2), jv(1, xi / 2)
        forms = {(0.5, 0.5): e * j0,
                 (-0.5, -0.5): xi / 2 * e * j1,
                 (-0.5, 0.5): -1j * xi / 2 * e * (j0 - 1j * j1),
                 (0.5, -0.5): 1j * xi / 2 * e * (j0 + 1j * j1)}
        for key, value in forms.items():
            worst = max(worst, abs(beta_transform(*key, xi) - value))
        worst = max(worst,
                    abs(1j * xi * np.mean(u * np.exp(1j * xi * u)) - forms[(0.5, -0.5)]),
                    abs(-1j * xi * np.mean((1 - u) * np.exp(1j * xi * u))
                        - forms[(-0.5, 0.5)]),
                    abs(forms[(-0.5, -0.5)]
                        + (forms[(-0.5, 0.5)] + forms[(0.5, -0.5)]) / 2))
    v = np.linspace(0.05, 0.95, 19)
    theta = np.arccos(2 * v - 1)
    tails = [(np.sqrt(v / (1 - v)) / math.pi, 1 / np.tan(theta / 2) / math.pi),
             (-np.sqrt((1 - v) / v) / math.pi, -np.tan(theta / 2) / math.pi),
             ((1 - 2 * v) / (2 * math.pi * np.sqrt(v * (1 - v))),
              -1 / np.tan(theta) / math.pi),
             (1 - 2 / math.pi * np.arcsin(np.sqrt(v)), theta / math.pi)]
    worst = max([worst] + [float(np.abs(a - b).max()) for a, b in tails])
    print(f"(2.3)-(2.4) Bessel forms, arcsine representation, tails: max error {worst:.1e}")
    assert worst < 1e-6


def contour(f, center, radius, points=256):
    angles = 2 * math.pi * np.arange(points) / points
    s = center + radius * np.exp(1j * angles)
    return np.mean(f(s) * radius * np.exp(1j * angles))


def components(big_l, xi):
    """Normalized contour integrals around rho_1 and s = 1 (x = e^L)."""
    alpha = xi / big_l
    radius = (abs(xi) / 2 + 3) / big_l

    def psi_factor(s):
        return np.sqrt(zeta(s) / zeta(s - 1j * alpha))

    def m_factor(s):
        zs = zeta(s)
        return 1 / (zs * np.sqrt(zeta(s - 1j * alpha) / zs))

    def around(point, factor):
        center = point + 1j * alpha / 2
        return contour(lambda s: factor(s) * np.exp(big_l * (s - point)) / s,
                       center, radius)

    return {(0.5, 0.5): around(RHO_1, m_factor) * RHO_1,
            (-0.5, -0.5): around(1, m_factor) * big_l ** 2,
            (-0.5, 0.5): around(1, psi_factor) * big_l,
            (0.5, -0.5): around(RHO_1, psi_factor) * RHO_1 * big_l}


def check_local():
    """(2.1)-(2.2) and the M rows of (2.4): relative error O(1/L)."""
    h = 1e-5
    zeta_prime = (zeta(RHO_1 + h)[0] - zeta(RHO_1 - h)[0]) / (2 * h)
    names = {(0.5, 0.5): "M zero", (-0.5, -0.5): "M pole",
             (-0.5, 0.5): "psi pole", (0.5, -0.5): "psi zero"}
    print("(2.1), (2.4) relative error at rho_1 and s = 1:")
    for xi in (1.0, 4.0):
        previous = None
        for big_l in (25, 50, 100, 200, 400):
            values = components(big_l, xi)
            values[(0.5, 0.5)] *= zeta_prime
            errors = {k: abs(values[k] / beta_transform(*k, xi) - 1) for k in values}
            print(f"    xi={xi}, L={big_l:3d}: "
                  + ", ".join(f"{names[k]} {errors[k]:.2e}" for k in names))
            if previous is not None:
                assert all(0.4 < errors[k] / previous[k] < 0.6 for k in errors)
            previous = errors
    # (2.2): first derivatives at xi = 0 are exact for psi.
    for big_l in (25, 100):
        plus, minus = components(big_l, 1e-4), components(big_l, -1e-4)
        zero_d = (plus[(0.5, -0.5)] - minus[(0.5, -0.5)]) / 2e-4
        pole_d = (plus[(-0.5, 0.5)] - minus[(-0.5, 0.5)]) / 2e-4
        print(f"(2.2) L={big_l}: L rho x^-rho d/dxi I = {zero_d:.8f} (exact i/2), "
              f"L x^-1 d/dxi P = {pole_d:.8f} (exact -i/2)")
        assert abs(zero_d - 0.5j) < 1e-6 and abs(pole_d + 0.5j) < 1e-6


def check_phase(sigma=2.0, alpha_max=3000.0, samples=60000):
    """(3.2): mean of e^{i eta arg zeta} against prod_p 2F1(-eta/2, eta/2; 1; p^-2s)."""
    alphas = np.linspace(0, alpha_max, samples)
    values = np.concatenate([zeta_line(sigma, alphas[i:i + 500])
                             for i in range(0, samples, 500)])
    args = np.angle(values)
    primes = primes_up_to(20000)
    r2 = primes.astype(float) ** (-2 * sigma)
    print(f"(3.2) sigma={sigma}, alpha in [0, {alpha_max:.0f}]:")
    for eta in (2.0, 5.0, 10.0):
        empirical = np.mean(np.exp(1j * eta * args))
        local, term = np.ones_like(r2), np.ones_like(r2)
        for j in range(1, 30):
            term = term * (-eta / 2 + j - 1) * (eta / 2 + j - 1) / j ** 2 * r2
            local = local + term
        product = float(np.prod(local))
        leading = float(np.prod(jv(0, eta * primes.astype(float) ** (-sigma))))
        print(f"    eta={eta:4}: empirical {empirical.real:+.5f}, "
              f"2F1 product {product:+.5f}, J0 only {leading:+.5f}")
        assert abs(empirical - product) < 1e-3
        assert abs(empirical - product) < abs(empirical - leading)


def band_limited_pair(big_p):
    """g(u) = cos^4(pi u / 2X) on |u| < X = log P, and h(t) = int g(u) e^{iut} du."""
    big_x = math.log(big_p)
    terms = ((3 / 8, 0.0), (1 / 2, math.pi / big_x), (1 / 8, 2 * math.pi / big_x))

    def g(u):
        u = np.asarray(u, dtype=float)
        inside = sum(c * np.cos(w * u) for c, w in terms)
        return np.where(np.abs(u) < big_x, inside, 0.0)

    def sinc_part(v):
        v = np.asarray(v, dtype=complex)
        small = np.abs(v) < 1e-12
        return np.where(small, big_x, np.sin(v * big_x) / np.where(small, 1, v))

    def h(t):
        t = np.asarray(t, dtype=complex)
        return sum(c * (sinc_part(t - w) + sinc_part(t + w)) for c, w in terms)

    return g, h


def prime_powers_below(big_p):
    ns, lams = [], []
    for p in primes_up_to(big_p):
        p, q = int(p), int(p)
        while q < big_p:
            ns.append(q)
            lams.append(math.log(p))
            q *= p
    return np.array(ns, dtype=float), np.array(lams)


def theta_prime(t):
    """theta'(t) = Re psi(1/4 + it/2)/2 - log(pi)/2 for the Riemann-Siegel theta."""
    return 0.5 * digamma(0.25 + 0.5j * np.asarray(t)).real - 0.5 * math.log(math.pi)


def check_explicit_formula(gammas, big_p=1000):
    """(4.3) with h_{t0}(t) = h(t - t0) + h(t + t0): zero sum = primes below P."""
    g, h = band_limited_pair(big_p)
    ns, lams = prime_powers_below(big_p)
    t = np.arange(0.0, 3000.0, 0.002)
    print(f"(4.3) P={big_p}, zeros below 400 (tail negligible):")
    for t0 in (0.0, 20.0, 50.0):
        def shifted(s):
            return h(s - t0) + h(s + t0)
        zero_sum = 2 * shifted(gammas).sum().real
        gamma_term = 2 / math.pi * np.trapezoid(shifted(t).real * theta_prime(t), t)
        primes = 2 * np.sum(lams * ns ** -0.5 * 2 * g(np.log(ns)) * np.cos(t0 * np.log(ns)))
        right = (shifted(0.5j) + shifted(-0.5j)).real + gamma_term - primes
        print(f"    t0={t0:4}: zero sum {zero_sum:+.10f}, prime side {right:+.10f}")
        assert abs(zero_sum - right) < 1e-9


def smoothed_zero_density(big_p, ts):
    """Right side of (4.3) for h_{t0}, with the Gamma term replaced by 4 theta'(t0)."""
    g, h = band_limited_pair(big_p)
    ns, lams = prime_powers_below(big_p)
    logs = np.log(ns)
    weights = lams * ns ** -0.5 * g(logs)
    primes = np.concatenate([np.cos(np.outer(ts[i:i + 2000], logs)) @ weights
                             for i in range(0, len(ts), 2000)])
    pole = (h(0.5j - ts) + h(0.5j + ts) + h(-0.5j - ts) + h(-0.5j + ts)).real
    return pole - 4 * primes + 4 * theta_prime(ts), h(0.0).real


def check_zero_detection(gammas):
    """Peaks of the prime-side density against the zeros (section 4 numerics)."""
    print("(4.4) peaks of the smoothed zero density built from primes below P:")
    for big_p, height in ((10 ** 4, 100.0), (10 ** 5, 300.0)):
        ts = np.arange(10.0, height, 0.005)
        density, top = smoothed_zero_density(big_p, ts)
        inner = (density[1:-1] > density[:-2]) & (density[1:-1] > density[2:])
        peaks = ts[1:-1][inner & (density[1:-1] > top / 2)]
        for cut in sorted({min(height, 150.0), height}):
            true = gammas[(gammas > 10.5) & (gammas < cut - 0.5)]
            found = peaks[(peaks > 10.5) & (peaks < cut - 0.5)]
            errors = np.array([np.abs(found - y).min() for y in true])
            missed = true[errors > 0.3]
            spurious = [p for p in found if np.abs(true - p).min() > 0.3]
            print(f"    P=10^{round(math.log10(big_p))}, t<{cut:.0f}: {len(true)} zeros, "
                  f"{len(found)} peaks, max error {errors[errors <= 0.3].max():.4f}, "
                  f"missed {np.round(missed, 3).tolist()}, spurious {np.round(spurious, 2).tolist()}")
            if cut <= 150:
                assert len(missed) == 0 and not spurious and errors.max() < 0.02
            else:
                # The pair 294.965, 295.573 (gap 0.61) is below the resolution 2 pi / log P.
                assert len(missed) == 2 and np.diff(missed)[0] < 2 * math.pi / math.log(big_p) * 1.2


def main():
    check_identities()
    check_coins()
    check_kernels()
    check_local()
    check_phase()
    gammas = zeros_up_to(400.0)
    check_explicit_formula(gammas)
    check_zero_detection(gammas)


if __name__ == "__main__":
    main()
