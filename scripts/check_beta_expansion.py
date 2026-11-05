"""Checks for the J-th order beta expansion (main-theorem.md section 9). Requires numpy and scipy.

Run: python -B scripts/check_beta_expansion.py            # x <= 10^6, about a minute
     python -B scripts/check_beta_expansion.py --x 1e7    # the table of appendix B.2

What is checked:
  * the Taylor coefficients Gt_ab of Gt(s1,s2) = Z(s1)^z Z(s2)^w H(s1,s2), computed by a
    two-variable Cauchy integral, against the closed forms (9.4)-(9.5) of the second order;
  * the arithmetic sums A_x(u^k) (Theorem 1, squarefree model) and A~_x(u^k) (Theorem 2,
    d_z model), computed exactly by a sieve, against the expansion (9.2) in the Li form
    and against its power form (9.3), for J = 1, ..., 6;
  * with --x 1e7, the seven values of the table in appendix B.1 (to 5 digits).
These are numerical consistency checks, not a proof.
"""

import argparse
import math

import numpy as np
from scipy.integrate import quad

EULER_GAMMA = 0.5772156649015329
STIELTJES_1 = -0.07281584548367672
BERNOULLI_2K = [1 / 6, -1 / 30, 1 / 42, -1 / 30, 5 / 66, -691 / 2730, 7 / 6, -3617 / 510]
K_MAX = 4  # moments u^0, ..., u^4


# ---------------------------------------------------------------- zeta and Gt

def z_function(s, n_terms=40):
    """Z(s) = (s-1) zeta(s) by Euler-Maclaurin, for complex arrays s."""
    s = np.asarray(s, dtype=complex)
    flat = s.ravel()
    n = np.arange(1, n_terms, dtype=float)
    head = np.sum(n[:, None] ** (-flat[None, :]), axis=0)
    tail = n_terms ** (1 - flat) / (flat - 1) + 0.5 * n_terms ** (-flat)
    fact, rising = 1.0, flat.copy()
    for k, b in enumerate(BERNOULLI_2K, start=1):
        fact *= 2 if k == 1 else (2 * k - 1) * (2 * k)
        tail = tail + b / fact * rising * n_terms ** (-flat - 2 * k + 1)
        rising = rising * (flat + 2 * k - 1) * (flat + 2 * k)
    return ((flat - 1) * (head + tail)).reshape(s.shape)


def zeta(s):
    s = np.asarray(s, dtype=complex)
    return z_function(s) / (s - 1)


def primes_up_to(n):
    sieve = np.ones(n + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = False
    return np.nonzero(sieve)[0]


def gt_coefficients(z, w, order, model, radius=0.15, points=24, prime_limit=200_000):
    """Gt_ab for a, b <= order by the Cauchy integral on |s1-1| = |s2-1| = radius.

    model "sqfree": H of (3.2), split as in (3.7) so that the prime product converges
    like p^(-3 sigma); model "dz": H = 1 (Theorem 2).
    """
    e = radius * np.exp(2j * np.pi * np.arange(points) / points)
    s1 = (1 + e[:, None] + 0 * e[None, :]).ravel()
    s2 = (1 + 0 * e[:, None] + e[None, :]).ravel()
    log_g = z * np.log(z_function(s1)) + w * np.log(z_function(s2))
    if model == "sqfree":
        e1, e2, e12 = z * (z + 1) / 2, w * (w + 1) / 2, z * w
        log_g += -e1 * np.log(zeta(2 * s1)) - e12 * np.log(zeta(s1 + s2)) - e2 * np.log(zeta(2 * s2))
        logs = np.log(primes_up_to(prime_limit).astype(float))
        for i in range(0, len(logs), 2000):
            lp = logs[i:i + 2000, None]
            x, y = np.exp(-lp * s1[None, :]), np.exp(-lp * s2[None, :])
            log_g += np.sum(np.log(1 + z * x + w * y) + z * np.log(1 - x) + w * np.log(1 - y)
                            - e1 * np.log(1 - x * x) - e12 * np.log(1 - x * y)
                            - e2 * np.log(1 - y * y), axis=0)
    coeff = np.fft.fft2(np.exp(log_g).reshape(points, points)) / points ** 2
    a = np.arange(order + 1)
    return np.real(coeff[:order + 1, :order + 1] / radius ** (a[:, None] + a[None, :]))


def second_order_closed_form(z, w, prime_limit=10 ** 7):
    """(9.4)-(9.5): Gt_ab, a+b <= 2, from kappa, nu, lambda (squarefree model)."""
    r = z + w
    p = primes_up_to(prime_limit).astype(float)
    lp = np.log(p)
    log_c = np.sum(np.log1p(r / p) + r * np.log1p(-1 / p))
    c = math.exp(log_c)
    kappa = EULER_GAMMA + np.sum((r + 1) * lp / ((p - 1) * (p + r)))
    nu = np.sum(lp ** 2 / (p + r) ** 2)
    lam = np.sum(lp ** 2 * ((r + 1) / ((p - 1) * (p + r)) + 1 / (p - 1) ** 2))
    theta = STIELTJES_1 + EULER_GAMMA ** 2 / 2 + lam / 2
    rho2 = kappa ** 2 - nu
    return {(0, 0): c, (1, 0): z * c * kappa, (0, 1): w * c * kappa,
            (2, 0): c * (z * z * rho2 / 2 - z * theta), (1, 1): c * z * w * rho2,
            (0, 2): c * (w * w * rho2 / 2 - w * theta)}


# ---------------------------------------------------------------- expansion (9.2)-(9.3)

def rgamma(x):
    return 0.0 if (x <= 0 and abs(x - round(x)) < 1e-12) else 1 / math.gamma(x)


def rising(a, k):
    out = 1.0
    for i in range(k):
        out *= a + i
    return out


def falling(a, m):
    out = 1.0
    for i in range(m):
        out *= a - i
    return out


def kernel_pairing(a, b, poly):
    """<K_{a,b}, F> for F = sum_k poly[k] u^k, by (4.2)."""
    return sum(c * rising(a, k) * rgamma(a + b + k) for k, c in enumerate(poly) if c)


def li_alpha(alpha, x):
    """Li_alpha(x) = int_2^x (log t)^alpha dt."""
    big_l = math.log(x)
    val, _ = quad(lambda v: math.exp(v - big_l) * v ** alpha, math.log(2), big_l,
                  limit=400, epsabs=0, epsrel=1e-13)
    return x * val


def order_pairings(z, w, gt, poly, n_orders):
    """P_i(F) = sum_{a+b=i} Gt_ab <K_{z-a,w-b}, F>."""
    return [sum(gt[a, i - a] * kernel_pairing(z - a, w - i + a, poly) for a in range(i + 1))
            for i in range(n_orders)]


def predict_li(z, w, pairings, x, order):
    return sum(pairings[i] * li_alpha(z + w - 1 - i, x) for i in range(order))


def power_coefficients(z, w, pairings):
    """<K^(j), F> of (9.3)."""
    r = z + w
    return [sum((-1) ** (j - i) * falling(r - 1 - i, j - i) * pairings[i] for i in range(j + 1))
            for j in range(len(pairings))]


def predict_power(z, w, coeffs, x, order):
    big_l = math.log(x)
    return x * big_l ** (z + w - 1) * sum(coeffs[j] / big_l ** j for j in range(order))


# ---------------------------------------------------------------- arithmetic sums

def arithmetic_sums(z, w, model, x_max, xs, primes, log_n, squarefree):
    """sum_{2<=n<=x} sum_{d|n} weight (log d / log n)^k at the points xs, k <= K_MAX.

    The divisor moments sum_{d|n} weight (log d)^k are built multiplicatively, one prime
    power at a time; weight = z^omega(d) w^omega(n/d) (sqfree) or d_z(d) d_w(n/d) (dz).
    """
    def d_coeff(alpha, a):
        out = 1.0
        for i in range(a):
            out *= (alpha + i) / (i + 1)
        return out

    mom = [np.zeros(x_max + 1) for _ in range(K_MAX + 1)]
    mom[0][:] = 1.0
    for p in primes:
        lp, pa, a = math.log(p), p, 1
        while pa <= x_max and (model == "dz" or a == 1):
            idx = np.arange(pa, x_max + 1, pa)
            if model == "dz":
                idx = idx[idx % (pa * p) != 0]
            old = [m[idx] for m in mom]  # fancy indexing: copies
            for k in range(K_MAX + 1):
                acc = np.zeros(len(idx))
                for j in range(a + 1):
                    c = d_coeff(z, j) * d_coeff(w, a - j)
                    if c:
                        for i in range(k + 1):
                            acc += c * math.comb(k, i) * (j * lp) ** (k - i) * old[i]
                mom[k][idx] = acc
            pa, a = pa * p, a + 1
    mask = np.ones(x_max + 1, dtype=bool) if model == "dz" else squarefree.copy()
    mask[:2] = False
    out = []
    for k in range(K_MAX + 1):
        v = np.zeros(x_max + 1)
        v[mask] = mom[k][mask] / log_n[mask] ** k
        out.append(np.cumsum(v)[xs])
    return np.array(out)


# ---------------------------------------------------------------- main

TESTS = {"u^2": [0, 0, 1, 0, 0], "(u-1/2)^4": [1 / 16, -1 / 2, 3 / 2, -2, 1], "1": [1, 0, 0, 0, 0],
         "u": [0, 1, 0, 0, 0]}
CASES = [("sqfree", 0.5, 0.5, ["u^2", "(u-1/2)^4"]), ("sqfree", 1.0, 1.0, ["u^2"]),
         ("sqfree", 2.0, 1.0, ["1", "u^2"]), ("sqfree", -1.0, 1.0, ["u", "u^2"]),
         ("dz", 0.5, 0.5, ["u^2", "(u-1/2)^4"]), ("dz", 2.0, 1.0, ["u^2"])]
N_ORDERS = 6
# Appendix B.1 at x = 10^7: A / (C x L^(r-1)), observed and second order (2.2).
TABLE_B1 = {("sqfree", 0.5, 0.5, "u^2"): (0.39957, 0.40163),
            ("sqfree", 0.5, 0.5, "(u-1/2)^4"): (0.03173, 0.03343),
            ("sqfree", 1.0, 1.0, "u^2"): (0.43441, 0.43526),
            ("sqfree", 2.0, 1.0, "u^2"): (0.46624, 0.44359),
            ("sqfree", 2.0, 1.0, "1"): (0.87273, 0.84226),
            ("sqfree", -1.0, 1.0, "u"): (-1.07117, -1.06204),
            ("sqfree", -1.0, 1.0, "u^2"): (-0.19026, -0.16535)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--x", type=float, default=1e6)
    x_max = int(parser.parse_args().x)
    xs = np.array(sorted({int(round(10 ** e)) for e in np.arange(4, math.log10(x_max) + 1e-9, 0.5)}))

    print("(9.4)-(9.5) second-order Gt_ab: Cauchy integral vs closed form")
    for z, w in [(0.5, 0.5), (2.0, 1.0), (-1.0, 1.0)]:
        gt = gt_coefficients(z, w, 2, "sqfree")
        closed = second_order_closed_form(z, w)
        err = max(abs(gt[a, b] - v) for (a, b), v in closed.items())
        print(f"  (z,w) = ({z:+.1f},{w:+.1f}): Gt20 = {gt[2, 0]:+.6f}, Gt11 = {gt[1, 1]:+.6f}, "
              f"Gt02 = {gt[0, 2]:+.6f}; max difference {err:.1e}")
        assert err < 1e-5

    primes = primes_up_to(x_max)
    log_n = np.zeros(x_max + 1)
    for p in primes:
        q = p
        while q <= x_max:
            log_n[q::q] += math.log(p)
            q *= p
    squarefree = np.ones(x_max + 1, dtype=bool)
    for p in primes[primes * primes <= x_max]:
        squarefree[p * p::p * p] = False

    print(f"\n(9.2)-(9.3) relative error (observed - expansion) / |C| x L^(r-1) at x = {xs[-1]:.0e}")
    print("  model    (z,w)       F            form   " + "".join(f"{'J=' + str(j):>10}" for j in range(1, N_ORDERS + 1)))
    for model, z, w, names in CASES:
        r = z + w
        gt = gt_coefficients(z, w, N_ORDERS, model)
        sums = arithmetic_sums(z, w, model, x_max, xs, primes, log_n, squarefree)
        x = int(xs[-1])
        scale = abs(gt[0, 0]) * x * math.log(x) ** (r - 1)
        for name in names:
            poly = TESTS[name]
            observed = float(np.dot(poly, sums[:, -1]))
            pairings = order_pairings(z, w, gt, poly, N_ORDERS)
            coeffs = power_coefficients(z, w, pairings)
            rel_li = [(observed - predict_li(z, w, pairings, x, j)) / scale for j in range(1, N_ORDERS + 1)]
            rel_pw = [(observed - predict_power(z, w, coeffs, x, j)) / scale for j in range(1, N_ORDERS + 1)]
            label = f"  {model:<8} ({z:+.1f},{w:+.1f})  {name:<11}"
            print(label + "  power  " + "".join(f"{v:10.1e}" for v in rel_pw))
            print(" " * len(label) + "  Li     " + "".join(f"{v:10.1e}" for v in rel_li))
            key = (model, z, w, name)
            if key in TABLE_B1:
                c_norm = gt[0, 0] * x * math.log(x) ** (r - 1)
                obs_n = observed / c_norm
                second_n = predict_power(z, w, coeffs, x, 2) / c_norm
                print(" " * len(label) + f"  B.1: observed {obs_n:.5f}, second order {second_n:.5f}")
                if x == 10 ** 7:
                    assert abs(obs_n - TABLE_B1[key][0]) < 6e-6
                    assert abs(second_n - TABLE_B1[key][1]) < 6e-6
            # Theorem 1/2 (J = 2) leaves an error of order L^-2; J = 3 must improve on it.
            assert abs(rel_li[2]) < abs(rel_pw[1])
            if (model, z, w, name) == ("sqfree", -1.0, 1.0, "u"):
                # A_x^{-1,1}(u) = -pi(x) exactly, and every order of (9.2) is -(Li_{-1}(x)):
                # what is left is pi(x) - li(x), the oscillation of the zeros (section 9.4).
                pi_x = int(np.searchsorted(primes, x, side="right"))
                assert abs(observed + pi_x) < 1e-6 * pi_x
                print(" " * len(label) + f"  pi(x) = {pi_x}, Li_-1(x) = {li_alpha(-1, x):.1f}")


if __name__ == "__main__":
    main()
