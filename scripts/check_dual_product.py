"""Checks for section 2 of docs/ng-duality.md. Requires numpy and scipy.

Run: python -B scripts/check_dual_product.py
Zeta values come from the Euler-Maclaurin formula in double precision. These
are consistency checks, not a proof.
"""

import math

import numpy as np

BERNOULLI = [1 / 6, -1 / 30, 1 / 42, -1 / 30, 5 / 66, -691 / 2730, 7 / 6,
             -3617 / 510, 43867 / 798, -174611 / 330]


def zeta_line(sigma, ts):
    """Vectorized Euler-Maclaurin zeta(sigma+it) for an array of t >= 0."""
    ts = np.asarray(ts, dtype=float)
    out = np.empty(ts.shape, dtype=complex)
    for start in range(0, len(ts), 400):
        t = ts[start:start + 400]
        s = sigma + 1j * t
        big_n = int(t.max() / (2 * math.pi)) + 30
        n = np.arange(1, big_n, dtype=float)
        total = np.exp(-np.outer(s, np.log(n))).sum(axis=1)
        total += big_n ** (1 - s) / (s - 1) + 0.5 * big_n ** (-s)
        rising, power = s.copy(), big_n ** (-s - 1)
        for k, b in enumerate(BERNOULLI, start=1):
            total += b / math.factorial(2 * k) * rising * power
            rising = rising * (s + 2 * k - 1) * (s + 2 * k)
            power = power / big_n ** 2
        out[start:start + 400] = total
    return out


def primes_up_to(n):
    sieve = np.ones(n + 1, dtype=bool)
    sieve[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if sieve[p]:
            sieve[p * p::p] = False
    return np.nonzero(sieve)[0]


def smallest_prime_factor(n):
    spf = np.arange(n + 1)
    for p in range(2, int(n ** 0.5) + 1):
        if spf[p] == p:
            block = spf[p * p::p]
            block[block == np.arange(p * p, n + 1, p)] = p
    return spf


def divisor_function(a, n_max, spf):
    """d_a(n) for 1 <= n <= n_max via d_a(p^j) = (a)_j / j!."""
    d = np.zeros(n_max + 1, dtype=complex)
    d[1] = 1
    for n in range(2, n_max + 1):
        p, m, j = spf[n], n, 0
        while m % p == 0:
            m //= p
            j += 1
        local = 1 + 0j
        for i in range(j):
            local *= (a + i) / (i + 1)
        d[n] = d[m] * local
    return d


def local_factor(a, r2, terms=80):
    """2F1(a, a; 1; r2) = sum_j ((a)_j / j!)^2 r2^j, vectorized in r2."""
    r2 = np.asarray(r2, dtype=float)
    total = np.zeros(r2.shape, dtype=complex)
    term = np.ones(r2.shape, dtype=complex)
    for j in range(terms):
        total += term
        term = term * ((a + j) / (j + 1)) ** 2 * r2
    return total


def euler_product(a, sigma, primes):
    return np.prod(local_factor(a, primes ** (-2 * sigma)))


def check_identity(spf):
    """(2.1): symmetric twist of the discrete kernels at sigma=2."""
    sigma, alpha, n_max = 2.0, 0.7, 20000
    print("(2.1) sum_n n^-sigma <K_n, e^{i alpha log n (2u-1)}> vs |zeta|^{2z}:")
    zeta_value = zeta_line(sigma, [alpha])[0]
    for kappa in (-0.5, 0.5, 1.5):
        d = divisor_function(kappa, n_max, spf)
        total = 0j
        for dd in range(1, n_max + 1):
            m = np.arange(1, n_max // dd + 1)
            total += d[dd] * dd ** (-sigma + 1j * alpha) * np.sum(
                d[m] * m ** (-sigma - 1j * alpha))
        target = abs(zeta_value) ** (2 * kappa)
        # Omitted terms dm > N are bounded by sum_{n>N} d_{2|kappa|}(n) n^-2.
        tail = 2 * math.log(n_max) ** max(2 * abs(kappa) - 1, 0) / n_max
        print(f"  z={kappa:+.1f}: {total.real:.10f}{total.imag:+.1e}i vs {target:.10f}"
              f"  (difference {abs(total - target):.1e}, tail bound {tail:.1e})")
        assert abs(total - target) < tail


def check_central_atoms(spf, primes):
    """(2.4): central atoms of K_{d^2} with imaginary exponent = Euler product."""
    d_max = 300000
    print("\n(2.4) sum_d d_a(d)^2 d^(-2 sigma) vs prod_p 2F1(a,a;1;p^(-2 sigma)), a=i*z*eta:")
    for sigma, eta in ((2.0, 3.0), (1.0, 2.0)):
        a = 1j * (-0.5) * eta
        d = divisor_function(a, d_max, spf)
        n = np.arange(1, d_max + 1)
        series = np.sum(d[1:] ** 2 * n ** (-2 * sigma))
        product = euler_product(a, sigma, primes)
        print(f"  sigma={sigma}, eta={eta}: {series:.6f} vs {product:.6f}")
        assert abs(series - product) < 2e-3


def check_mean_square(primes):
    """Lemma 2.2: mean square of zeta*Q_X - 1 on [A, 2A]."""
    sigma, x_cut, big_a = 0.75, 5, 600.0
    t = np.linspace(big_a, 2 * big_a, 120001)
    small = primes[primes <= x_cut]
    q = np.prod([1 - p ** (-(sigma + 1j * t)) for p in small], axis=0)
    g = zeta_line(sigma, t) * q
    mean = np.trapezoid(np.abs(g - 1) ** 2, t) / big_a
    target = math.prod(1 - p ** (-2 * sigma) for p in small) * zeta_2s(2 * sigma) - 1
    print(f"\nLemma 2.2 (sigma={sigma}, X={x_cut}, A={big_a:.0f}): "
          f"mean |zeta Q_X - 1|^2 = {mean:.4f}, limit predicted {target:.4f}")
    assert abs(mean - target) < 0.1 * target


def zeta_2s(x):
    return float(zeta_line(x / 1.0, [0.0])[0].real) if x != 1 else math.inf


def check_empirical(primes):
    """Theorem 1: empirical characteristic function of 2*z*log|zeta(sigma+i alpha)|."""
    sigma, kappa, big_a = 0.75, -0.5, 3000.0
    alpha = np.linspace(0.0, big_a, 300001)
    phi = 2 * kappa * np.log(np.abs(zeta_line(sigma, alpha)))
    print(f"\nTheorem 1 (sigma={sigma}, z={kappa}, alpha in [0,{big_a:.0f}]):")
    print("   eta   empirical            prod_p 2F1          prod_p J0 (leading)")
    from scipy.special import j0
    for eta in (0.5, 1.0, 2.0, 3.0):
        empirical = np.mean(np.exp(1j * eta * phi))
        product = euler_product(1j * kappa * eta, sigma, primes)
        bessel = np.prod(j0(2 * abs(kappa) * eta * primes ** (-sigma)))
        print(f"  {eta:4.1f}  {empirical.real:+.4f}{empirical.imag:+.4f}i   "
              f"{product.real:+.4f}{product.imag:+.4f}i   {bessel:+.4f}")
        assert abs(empirical - product) < 0.03


def main():
    spf = smallest_prime_factor(300000)
    primes = primes_up_to(2_000_000).astype(float)
    check_identity(spf)
    check_central_atoms(spf, primes)
    check_mean_square(primes)
    check_empirical(primes)


if __name__ == "__main__":
    main()
