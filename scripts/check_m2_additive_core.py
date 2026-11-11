"""Checks for docs/m2-additive-core.md. Requires numpy and scipy.

Run: python -B scripts/check_m2_additive_core.py

The exact identities of sections 1-3 and 5-8 are checked against the definition
b(n) = sum_{p^a || n} a(a+1)/2 (log p)^2. Numerical agreement is a consistency
check, not a proof.
"""

import math
import random
from itertools import combinations

import numpy as np
from scipy.special import expi

GAMMA = 0.5772156649015329
GAMMA1 = -0.0728158454836767
KAPPA = GAMMA ** 2 + 3 * GAMMA1          # c of m2-divisor-moment (5.2)
N_MAX = 10 ** 6 + 2


def sieve_spf(n):
    spf = np.zeros(n + 1, dtype=np.int64)
    for p in range(2, n + 1):
        if spf[p] == 0:
            spf[p::p][spf[p::p] == 0] = p
    return spf


SPF = sieve_spf(N_MAX)


def factor(n):
    out = {}
    q = 2
    while n > N_MAX:
        if q * q > n:
            out[n] = out.get(n, 0) + 1
            return out
        while n % q == 0:
            out[q] = out.get(q, 0) + 1
            n //= q
        q += 1
    while n > 1:
        q = int(SPF[n])
        a = 0
        while n % q == 0:
            n //= q
            a += 1
        out[q] = a
    return out


def b_of(n):
    return sum(a * (a + 1) / 2 * math.log(q) ** 2 for q, a in factor(n).items())


def w_of(n):
    return sum(a * math.log(q) ** 2 for q, a in factor(n).items())


def r_of(n):
    return sum(a * (a - 1) / 2 * math.log(q) ** 2 for q, a in factor(n).items())


def m2(n):
    return 0.25 + b_of(n) / (4 * math.log(n) ** 2)


def beta(t):
    lt = math.log(t)
    return 0.5 * lt * lt + GAMMA * lt - KAPPA


def p2(n):
    lt = math.log(n)
    return 3 / 8 + GAMMA / (4 * lt) - KAPPA / (4 * lt * lt)


def e2(n):
    return m2(n) - p2(n)


def eps(n):
    return b_of(n) - beta(n)


def h_of(n):
    return w_of(n) - GAMMA * math.log(n)


def saw(t):
    return t - math.floor(t) - 0.5


def check_section1():
    worst = 0.0
    for n in range(2, 5000):
        lt = math.log(n)
        q2p = 3 / 8 * lt * lt + GAMMA / 4 * lt - KAPPA / 4      # Q_2'(n) from (5.3)
        worst = max(worst, abs(p2(n) - q2p / lt ** 2),
                    abs(e2(n) - (b_of(n) - beta(n)) / (4 * lt * lt)))
    print(f"Prop 1.1: p_2 = Q_2'/log^2 and e_2 = (b - beta)/(4 log^2), n < 5000: max error {worst:.1e}")
    assert worst < 1e-13

    def int_p2(x):
        li = expi(math.log(x)) - expi(math.log(2))
        i2 = li - x / math.log(x) + 2 / math.log(2)
        return 3 / 8 * (x - 2) + GAMMA / 4 * li - KAPPA / 4 * i2

    ns = np.arange(2, 10 ** 6 + 1, dtype=float)
    lt = np.log(ns)
    cum = np.cumsum(3 / 8 + GAMMA / (4 * lt) - KAPPA / (4 * lt * lt))
    rows = []
    for x in (1e3 + 0.3, 1e4 + 0.7, 1e5 + 0.1, 1e6 - 0.4):
        k = int(math.floor(x)) - 2
        val = cum[k] - int_p2(x) + saw(x) * p2(x)
        rows.append(f"{x:.1f}: {val:.10f}")
    print("Theorem 1.2: sum p_2 - int p_2 + B1(x) p_2(x) tends to a constant:  " + ";  ".join(rows))


def check_section2():
    worst = 0.0
    for m in range(2, 120):
        for n in range(2, 120):
            if math.gcd(m, n) != 1:
                continue
            lm, ln_ = math.log(m), math.log(n)
            s = lm / (lm + ln_)
            t = 1 - s
            worst = max(worst,
                        abs(eps(m * n) - (eps(m) + eps(n) - lm * ln_ - KAPPA)),
                        abs(e2(m * n) - (s * s * e2(m) + t * t * e2(n) - s * t / 4
                                         - KAPPA / (4 * (lm + ln_) ** 2))),
                        abs((b_of(m * n) - GAMMA * math.log(m * n))
                            - (b_of(m) - GAMMA * lm) - (b_of(n) - GAMMA * ln_)))
    print(f"Prop 2.1 / Cor 2.2: cocycle of eps and e_2, additivity of b - gamma log: max error {worst:.1e}")
    assert worst < 1e-11
    worst = 0.0
    for c in range(2, 300):
        for p in (2, 3, 5, 7, 11, 101, 1009, 10007):
            if c % p == 0:
                continue
            u, v = math.log(c), math.log(p)
            big_l = u + v
            t = u / big_l
            f1 = t * t * e2(c) + (1 - t) * (1 - 3 * t) / 8 - GAMMA * (1 - t) / (4 * big_l)
            f2 = (0.5 * (v - u) ** 2 - u * u + b_of(c) - GAMMA * big_l + KAPPA) / (4 * big_l ** 2)
            worst = max(worst, abs(e2(c * p) - f1), abs(e2(c * p) - f2))
    print(f"Prop 2.3: band formula for e_2(cp), c < 300: max error {worst:.1e}")
    assert worst < 1e-13


def check_section3():
    worst = 0.0
    for n in range(2, 5000):
        lt = math.log(n)
        fac = factor(n)
        v2 = sum(a * (math.log(q) / lt) ** 2 for q, a in fac.items())
        s = r_of(n) / lt ** 2
        worst = max(worst, abs(e2(n) - (0.25 * (v2 - 0.5) + 0.25 * s - GAMMA / (4 * lt)
                                         + KAPPA / (4 * lt * lt))))
    print(f"Prop 3.1: e_2 = (sum V_i^2 - 1/2)/4 + s/4 - gamma/(4L) + kappa/(4L^2), n < 5000: max error {worst:.1e}")
    assert worst < 1e-13
    worst = 0.0
    for n in range(2, 3000):
        tot = 0.0
        for d in range(2, n + 1):
            if n % d == 0:
                fd = factor(d)
                if len(fd) == 1:
                    (q, k), = fd.items()
                    tot += math.log(q) * (k * math.log(q) - GAMMA)
        worst = max(worst, abs(tot - (b_of(n) - GAMMA * math.log(n))))
    print(f"Cor 3.3: b - gamma log = 1 * (Lambda (log - gamma)), n < 3000: max error {worst:.1e}")
    assert worst < 1e-10


def check_section5():
    worst = 0.0
    for d in range(2, 60):
        for m in range(1, 400):
            if math.gcd(d, m * (m + 1)) != 1:
                continue
            core = lambda n: b_of(n) - GAMMA * math.log(n) if n > 1 else 0.0
            worst = max(worst, abs(core(d * m + d) - core(d * m) - (core(m + 1) - core(m))))
            if m >= 2:
                worst = max(worst, abs(eps(d * m + d) - eps(d * m)
                                       - (eps(m + 1) - eps(m) - math.log(d) * math.log(1 + 1 / m))))
    print(f"Prop 5.1 (D1): shifts on dZ reduce to the unit shift: max error {worst:.1e}")
    assert worst < 1e-11


def check_section6():
    worst = 0.0
    for n in range(2, 4000):
        lt = math.log(n)
        a2 = lt * lt * m2(n)
        worst = max(worst, abs(b_of(n) - w_of(n) - r_of(n)),
                    abs(eps(n) - (h_of(n) + r_of(n) - 0.5 * lt * lt + KAPPA)))
        bsq = sum(2 * a * (2 * a + 1) / 2 * math.log(q) ** 2 for q, a in factor(n).items())
        worst = max(worst, abs(bsq - (4 * b_of(n) - w_of(n))))
        eps_sq = bsq - (0.5 * (2 * lt) ** 2 + GAMMA * 2 * lt - KAPPA)
        worst = max(worst, abs(eps_sq - (4 * eps(n) - w_of(n) + 2 * GAMMA * lt - 3 * KAPPA)))
        a2_sq = ((2 * lt) ** 2 + bsq) / 4
        worst = max(worst, abs(a2_sq - (4 * a2 - w_of(n) / 4)))
        hn, hn1 = h_of(n), h_of(n + 1)
        hprod = w_of(n) + w_of(n + 1) - GAMMA * math.log(n * (n + 1))
        worst = max(worst, abs(2 * hn * hn1 - (hprod ** 2 - (2 * hn) ** 2 / 4 - (2 * hn1) ** 2 / 4)))
    print(f"Section 6: b = w + r, eps(n^2), a_2(n^2), polarization identity, n < 4000: max error {worst:.1e}")
    assert worst < 1e-9


def g_weight(q):
    lq = math.log(q)
    return lq * (lq - GAMMA)


def prime_powers(limit):
    primes = np.nonzero(SPF[:limit + 1] == np.arange(limit + 1))[0]
    primes = primes[primes >= 2]
    out = []
    for p in primes:
        p = int(p)
        q = p
        while q <= limit:
            out.append((q, p))
            q *= p
    out.sort()
    return out


def check_section7():
    rnd = random.Random(1)
    bad = 0
    for _ in range(20000):
        d, e = rnd.randint(1, 60), rnd.randint(2, 60)
        if math.gcd(d, e) != 1:
            continue
        x = rnd.uniform(1, 5000)
        cnt = sum(1 for n in range(1, int(x) + 1) if n % d == 0 and (n + 1) % e == 0)
        db = pow(d, -1, e)
        f = x / (d * e) - (saw(x / (d * e) + db / e) - saw(db / e))
        bad += abs(cnt - f) > 1e-9
    print(f"Lemma 7.1: exact count of n <= x with d | n, e | n+1, 20000 random cases: failures {bad}")
    assert bad == 0

    pp = prime_powers(10 ** 6)
    qs = np.array([q for q, _ in pp], dtype=float)
    base = np.array([p for _, p in pp], dtype=float)
    gw = np.log(base) * (np.log(base) - GAMMA) / qs
    cum = np.concatenate([[0.0], np.cumsum(gw)])
    print("Prop 7.2: S(D) = sum_{de <= D, (d,e) = 1} G(d) G(e)/(de);  (S(D) - (log D)^4/24) / (log D)^3")
    for k in (3, 4, 5, 6):
        big_d = 10.0 ** k
        sel = qs <= big_d / 2
        idx = np.searchsorted(qs, big_d / qs[sel], side="right")
        total = float((gw[sel] * cum[idx]).sum())
        same = 0.0
        for (q1, p1) in pp:
            if q1 > big_d / 2:
                break
            q2 = p1
            while q1 * q2 <= big_d:
                same += g_weight(p1) ** 2 / (q1 * q2)
                q2 *= p1
        s_val = total - same
        lam = math.log(big_d)
        print(f"  D = 10^{k}: S = {s_val:.4f}, (log D)^4/24 = {lam ** 4 / 24:.4f}, "
              f"normalized difference {(s_val - lam ** 4 / 24) / lam ** 3:+.4f}")

    x = 10 ** 6
    hv = np.zeros(x + 2)
    for n in range(2, x + 2):
        hv[n] = h_of(n)
    corr = float((hv[1:x + 1] * hv[2:x + 2]).sum())
    sq = float((hv[1:x + 1] ** 2).sum())
    mean = float(hv[1:x + 1].sum()) / x
    lam = math.log(x)
    print(f"Section 7 [numerics], x = 10^6: C(x)/(x log^4 x) = {corr / (x * lam ** 4):.4f}, "
          f"(mean h)^2/log^4 x = {mean ** 2 / lam ** 4:.4f}, sum h^2/(x log^4 x) = {sq / (x * lam ** 4):.4f}")


def check_section8():
    bad = 0
    rnd = random.Random(2)
    for _ in range(5000):
        q = rnd.choice([2, 3, 5, 7, 11])
        j = rnd.randint(1, 3)
        e = q ** j
        t = rnd.uniform(0, 10)
        s = sum(saw(t + a / e) for a in range(1, e) if a % q)
        bad += abs(s - (saw(e * t) - saw(e * t / q))) > 1e-9
        s0 = sum(saw(a / e) for a in range(1, e) if a % q)
        bad += abs(s0) > 1e-9
    print(f"Lemma 8.1: reduced Hermite identity and zero mean of B1(a/e), 5000 cases: failures {bad}")
    assert bad == 0
    x = 10 ** 5 + 0.5
    worst = 0.0
    cum = 0.0
    for n in range(2, int(x) + 1):
        cum += h_of(n)
    pp = prime_powers(int(x))
    rhs = sum(g_weight(p) * (x / q - saw(x / q) - 0.5) for q, p in pp)
    worst = abs(cum - rhs)
    print(f"(8.1): sum_{{n <= x}} h(n) against the sawtooth form, x = 10^5 + 1/2: error {worst:.1e}")
    assert worst < 1e-6 * x


def main():
    check_section1()
    check_section2()
    check_section3()
    check_section5()
    check_section6()
    check_section7()
    check_section8()


if __name__ == "__main__":
    main()
