"""Checks for docs/m2-cp-families.md. Requires numpy and scipy.

Run: python -B scripts/check_cp_families.py
     python -B scripts/check_cp_families.py --figure docs/figures/m2-cp-bands.svg

The identities of sections 1-3 are checked against the definition of M_2 as a
weighted divisor sum. The tables of appendix A are recomputed from scratch.
Numerical agreement is a consistency check, not a proof.
"""

import argparse
import math
from itertools import combinations

import numpy as np
from scipy.integrate import quad

PRIME_LIMIT = 5 * 10 ** 7


def prime_sieve(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    s[4::2] = False
    for i in range(3, int(n ** 0.5) + 1, 2):
        if s[i]:
            s[i * i::2 * i] = False
    return np.nonzero(s)[0]


PRIMES = prime_sieve(PRIME_LIMIT)


def factor(n):
    out, m = {}, n
    for q in PRIMES:
        q = int(q)
        if q * q > m:
            break
        while m % q == 0:
            out[q] = out.get(q, 0) + 1
            m //= q
    if m > 1:
        out[m] = out.get(m, 0) + 1
    return out


def b_of(n):
    """(1.1): b(n) = sum over p^a || n of a(a+1)/2 (log p)^2."""
    return sum(a * (a + 1) / 2 * math.log(q) ** 2 for q, a in factor(n).items())


def m2(n):
    return 0.25 + b_of(n) / (4 * math.log(n) ** 2)


def m2_definition(n):
    """M_2(n) = sum_{d|n} f(d) f(n/d) (log d / log n)^2 with f = d_{1/2}."""
    divs = [(1, 1.0)]
    for q, a in factor(n).items():
        new = []
        for d, wd in divs:
            for k in range(a + 1):
                new.append((d * q ** k, wd))
        divs = new

    def f(m):
        return math.prod(math.comb(2 * a, a) / 4 ** a for a in factor(m).values())

    big_l = math.log(n)
    return sum(f(d) * f(n // d) * (math.log(d) / big_l) ** 2 for d, _ in divs)


def band(c_log, c_b, v):
    """(1.3): F_c(v)."""
    return 0.25 + (v ** 2 + c_b) / (4 * (c_log + v) ** 2)


def primes_for(c, y, lo=2):
    p = PRIMES[(PRIMES >= lo) & (PRIMES <= y)]
    return p[c % p != 0]


def band_mean(c, x):
    """(4.1): E_X(c) and the means m_1, m_2."""
    p = primes_for(c, x // c).astype(float)
    u, b = math.log(c), b_of(c)
    big_l = u + np.log(p)
    m1, m2_ = (1 / big_l).mean(), (1 / big_l ** 2).mean()
    e = (band(u, b, np.log(p))).mean()
    assert abs(e - (0.5 - u / 2 * m1 + (u * u + b) / 4 * m2_)) < 1e-13
    return e, m1, m2_, len(p)


def check_exact():
    worst = 0.0
    for n in range(2, 3001):
        worst = max(worst, abs(m2(n) - m2_definition(n)))
    print(f"(1.1) M_2 = 1/4 + b/(4L^2) against the divisor sum, n <= 3000: max error {worst:.1e}")
    assert worst < 1e-13
    worst_step = worst_rec = 0.0
    for c in range(2, 400):
        u, bc, mc = math.log(c), b_of(c), m2(c)
        for p in PRIMES[:60]:
            p = int(p)
            a = factor(c).get(p, 0)
            worst_step = max(worst_step, abs(b_of(c * p) - bc - (a + 1) * math.log(p) ** 2))
            if a == 0:
                t = u / math.log(c * p)
                rec = (1 - t) ** 2 + t * t * (4 * mc - 1)
                worst_rec = max(worst_rec, abs(4 * m2(c * p) - 1 - rec),
                                abs(m2(c * p) - band(u, bc, math.log(p))))
    print(f"(1.2) b(cp) - b(c) = (a+1)(log p)^2: max error {worst_step:.1e}; "
          f"(1.3)-(1.4): max error {worst_rec:.1e}  (c < 400, first 60 primes)")
    assert worst_step < 1e-10 and worst_rec < 1e-13
    bad = [c for c in range(2, 10 ** 4) if b_of(c) > math.log(c) ** 2 * (1 + 1e-14)]
    eq = [c for c in range(2, 10 ** 4) if abs(b_of(c) - math.log(c) ** 2) < 1e-9]
    is_p = set(int(q) for q in PRIMES[PRIMES < 10 ** 4])
    print(f"Lemma 1.2 for c < 10^4: violations {len(bad)}; equality exactly at primes: {set(eq) == is_p}")
    assert not bad and set(eq) == is_p


def check_shape():
    worst = 0.0
    for c in range(2, 2001):
        u, bc = math.log(c), b_of(c)
        v_star = bc / u
        assert math.exp(v_star) <= c * (1 + 1e-12)
        vs = np.linspace(1e-3, 30, 4001)
        vals = band(u, bc, vs)
        k = int(np.argmin(vals))
        worst = max(worst, abs(vals[k] - (0.25 + bc / (4 * (bc + u * u)))))
        assert np.all(np.diff(vals[vs <= v_star]) < 0) and np.all(np.diff(vals[vs >= v_star]) > 0)
    print(f"Prop 2.1 for c <= 2000: U shape, p* <= c, grid minimum vs closed form within {worst:.1e}")
    c = 194
    u, bc = math.log(c), b_of(c)
    print(f"  c = 194: p* = {math.exp(bc / u):.2f}, minimum {0.25 + bc / (4 * (bc + u * u)):.6f}; "
          + ", ".join(f"M_2(194*{p}) = {m2(194 * p):.4f}" for p in (3, 11, 59, 101)))


def table_families():
    print("Table A.1 / A.2 (X = 20002): G = 1/2 - u m_1/2 + u^2 m_2/4, H = b(c) m_2/4, E_X(c) = G + H")
    for c in (2, 6, 10, 14, 30, 90, 102):
        e, m1, m2_, n = band_mean(c, 20002)
        u, bc = math.log(c), b_of(c)
        g = 0.5 - u / 2 * m1 + u * u / 4 * m2_
        h = bc / 4 * m2_
        print(f"  c = {c:>3}: N = {n:>4}, b = {bc:.6f}, G = {g:.6f}, H = {h:.6f}, E = {e:.6f}")


def check_pairs():
    x = 20002
    cs = list(range(2, 200, 4))
    r_old = []
    r_new = []
    viol = 0
    for c1, c2 in combinations(cs, 2):
        lo, hi = min(c1, c2), max(c1, c2)
        p = PRIMES[(PRIMES <= x // hi)]
        p = p[(lo * hi) % p != 0].astype(float)
        v = np.log(p)
        ul, uh, bl, bh = math.log(lo), math.log(hi), b_of(lo), b_of(hi)
        ll, lh = ul + v, uh + v
        d = band(ul, bl, v) - band(uh, bh, v)
        exact = (v ** 2 + bl) / 4 * (1 / ll ** 2 - 1 / lh ** 2) + (bl - bh) / (4 * lh ** 2)
        assert np.max(np.abs(d - exact)) < 1e-14
        sharp = (v ** 2 + bl) / 4 * (1 / ll ** 2 - 1 / lh ** 2) + abs(bl - bh) / (4 * lh ** 2)
        new = (uh - ul) / (2 * lh) + abs(bl - bh) / (4 * lh ** 2)
        old = (uh - ul) / (2 * ll) + abs(bl - bh) / (4 * lh ** 2)
        viol += int(np.any(np.abs(d) > sharp + 1e-15) or np.any(sharp > new + 1e-15)
                    or np.any(new > old + 1e-15))
        if bl >= bh:
            assert np.max(np.abs(np.abs(d) - sharp)) < 1e-14
        gap = abs(d.mean())
        r_old.append(gap / old.mean())
        r_new.append(gap / new.mean())
        if (lo, hi) == (194, 198):
            print(f"  (194, 198): {len(p)} primes, mean difference {gap:.7f}, "
                  f"bound (3.5) {new.mean():.7f}, averaged (3.6) {old.mean():.7f}")
    print(f"Prop 3.1 on {len(r_new)} pairs (X = 20002): chain |d| <= (3.2) <= (3.3) <= "
          f"(3.6) violated {viol} times")
    print(f"  ratio to averaged bound: (3.5) max {max(r_new):.6f} mean {np.mean(r_new):.6f}; "
          f"averaged (3.6) max {max(r_old):.6f} mean {np.mean(r_old):.6f}")
    assert viol == 0


def check_near_pairs():
    print("Table A.3: pairs with close b(c), odd primes p <= 10^5 dividing neither c")
    for c1, c2 in ((57450, 24634), (89358, 89362)):
        p = PRIMES[(PRIMES >= 3) & (PRIMES <= 10 ** 5)]
        p = p[(c1 % p != 0) & (c2 % p != 0)].astype(float)
        means = []
        for c in (c1, c2):
            means.append(band(math.log(c), b_of(c), np.log(p)).mean())
        print(f"  {c1} vs {c2}: {len(p)} primes, b = {b_of(c1):.9f}, {b_of(c2):.9f}; "
              f"means {means[0]:.9f}, {means[1]:.9f}; difference {means[1] - means[0]:.6e}")


def predictor(c, x):
    """(5.2): odd primes from 3, odd prime factors of c removed (c even)."""
    y = x / c
    u, bc = math.log(c), b_of(c)

    def g(t):
        return band(u, bc, math.log(t))

    qs = [q for q in factor(c) if q % 2 and 3 <= q <= y]
    num = quad(lambda t: g(t) / math.log(t), 3, y, limit=200)[0] - sum(g(q) for q in qs)
    den = quad(lambda t: 1 / math.log(t), 3, y, limit=200)[0] - len(qs)
    return num / den


def check_predictor():
    errs = []
    w_mean = []
    for c in range(2, 200, 4):
        e, m1, m2_, n = band_mean(c, 20002)
        errs.append(predictor(c, 20002) - e)
        w_mean.append(m2_ / 4)
    errs = np.array(errs)
    lam = math.log(20002)
    print(f"Table A.4 (X = 20002, c = 2, 6, ..., 198): MAE {np.abs(errs).mean():.7f}, "
          f"RMSE {math.sqrt((errs ** 2).mean()):.7f}, max {np.abs(errs).max():.7f}")
    print(f"  mean of W_X(c) = m_2/4 over the 50 bands: {np.mean(w_mean):.6f}; "
          f"1/(4 lambda^2) = {1 / (4 * lam ** 2):.6f}, 1/(4 (lambda-1)^2) = {1 / (4 * (lam - 1) ** 2):.6f}")


def check_asymptotics():
    print("Theorem 4.1 / Corollary 4.2: lambda^2 (E - first order), lambda^3 (E - second order)")
    for c in (2, 6, 30):
        u, bc = math.log(c), b_of(c)
        row = []
        for k in (5, 6, 7, 8):
            x = 10 ** k
            if x // c > PRIME_LIMIT:
                continue
            e = band_mean(c, x)[0]
            lam = math.log(x)
            first = 0.5 - u / (2 * lam)
            second = first + (u * u - 2 * u + bc) / (4 * lam ** 2)
            row.append(f"10^{k}: {(e - first) * lam ** 2:+.4f}, {(e - second) * lam ** 3:+.4f}")
        print(f"  c = {c:>2}: " + "; ".join(row))


def check_additivity():
    """Prop 1.3 (1.5) and Prop 3.2 (3.7)-(3.8)."""
    worst = 0.0
    for m in range(2, 150):
        for n in range(2, 150):
            if math.gcd(m, n) != 1:
                continue
            s = math.log(m) / math.log(m * n)
            t = 1 - s
            worst = max(worst, abs(m2(m * n) - (s * s * m2(m) + t * t * m2(n) + s * t / 2)))
    print(f"Prop 1.3 (1.5) on coprime m, n < 150: max error {worst:.1e}")
    assert worst < 1e-13
    x = 20002
    cs = list(range(2, 200, 4))
    worst = 0.0
    viol = 0
    for c1, c2 in combinations(cs, 2):
        lo, hi = min(c1, c2), max(c1, c2)
        p = PRIMES[(PRIMES <= x // hi)]
        p = p[(lo * hi) % p != 0].astype(float)
        v = np.log(p)
        ul, uh, bl, bh = math.log(lo), math.log(hi), b_of(lo), b_of(hi)
        ml, mh = m2(lo), m2(hi)
        tl, th = ul / (ul + v), uh / (uh + v)
        d = band(ul, bl, v) - band(uh, bh, v)
        exact = (th - tl) * (0.5 - (th + tl) * ml) + th * th * (ml - mh)
        worst = max(worst, float(np.max(np.abs(d - exact))))
        bound = (1 - tl) * (uh - ul) / (2 * (uh + v)) + th * th * abs(ml - mh)
        first_old = (uh - ul) / (2 * (uh + v))
        viol += int(np.any(np.abs(d) > bound + 1e-15) or np.any((1 - tl) * first_old > first_old + 1e-15))
    print(f"Prop 3.2 on {len(cs) * (len(cs) - 1) // 2} pairs: (3.7) max error {worst:.1e}, "
          f"(3.8) violated {viol} times")
    assert worst < 1e-13 and viol == 0


def check_shifted_value():
    print("Corollary 4.3: lambda^3 (E_X(c) - F_c(log(Y/e)))")
    for c in (2, 6, 30):
        u, bc = math.log(c), b_of(c)
        row = []
        for k in (5, 6, 7, 8):
            x = 10 ** k
            if x // c > PRIME_LIMIT:
                continue
            e = band_mean(c, x)[0]
            lam = math.log(x)
            f = band(u, bc, lam - u - 1)
            row.append(f"10^{k}: {(e - f) * lam ** 3:+.4f}")
        print(f"  c = {c:>2}: " + "; ".join(row))


def check_small_bands():
    """Prop 7.1 (a): D_c - D~_c = pi(Y) - li(Y) exactly (the excluded prime factors cancel)."""
    from scipy.special import expi
    x = 10 ** 6
    c = np.arange(2, x // 2 + 1)
    y = x / c
    pi_y = np.searchsorted(PRIMES, y, side="right")
    li_y = expi(np.log(y)) - expi(math.log(2))
    total = np.abs(pi_y - li_y).sum()
    print(f"Prop 7.1 (a), x = 10^6: sum_c |D_c - D~_c| / x = {total / x:.4f}")
    assert total / x > 1 / 70


def check_mean():
    x = 10 ** 6
    spf = np.zeros(x + 1, dtype=np.int64)
    for p in PRIMES[PRIMES <= x]:
        p = int(p)
        sel = spf[p::p] == 0
        spf[p::p][sel] = p
    bvals = np.zeros(x + 1)
    sq = np.zeros(x + 1)
    for n in range(2, x + 1):
        m = n
        while m > 1:
            q = spf[m]
            a = 0
            while m % q == 0:
                m //= q
                a += 1
            lq = math.log(q) ** 2
            bvals[n] += a * (a + 1) / 2 * lq
            sq[n] += a * (a - 1) / 2 * lq
    ln = np.log(np.arange(2, x + 1)) ** 2
    mean_m2 = (0.25 + bvals[2:] / (4 * ln)).mean()
    print(f"Prop 6.1 (x = 10^6): mean M_2 = {mean_m2:.6f} (limit 3/8); "
          f"mean square-full correction sum a(a-1)/2 (log p)^2/L^2 = {(sq[2:] / ln).mean():.6f}")


def write_figure(path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager

    names = {f.name for f in font_manager.fontManager.ttflist}
    if "Noto Sans CJK JP" in names:
        plt.rcParams["font.family"] = "Noto Sans CJK JP"
    plt.rcParams.update({"svg.fonttype": "path", "font.size": 10,
                         "axes.edgecolor": "#c3c2b7", "axes.labelcolor": "#52514e",
                         "xtick.color": "#898781", "ytick.color": "#898781"})
    colors = ["#2a78d6", "#d9822b", "#3a9a5b", "#c4473a", "#7d5ba6"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.6, 3.8))
    x = 20002
    for c, col in zip((2, 6, 10, 14, 30), colors):
        p = primes_for(c, x // c).astype(float)
        u, bc = math.log(c), b_of(c)
        ax1.plot(c * p, band(u, bc, np.log(p)), ".", ms=2, color=col, rasterized=True)
        grid = np.geomspace(3 * c, x, 300)
        ax1.plot(grid, band(u, bc, np.log(grid / c)), lw=0.8, color=col, alpha=0.6,
                 label=f"$c={c}$")
    ax1.set_xscale("log")
    ax1.set_xlabel("$n=cp$")
    ax1.set_ylabel("$M_2(cp)$")
    ax1.legend(frameon=False, fontsize=8, loc="lower right")
    c = 194
    u, bc = math.log(c), b_of(c)
    vs = np.linspace(0.3, 9, 400)
    ax2.plot(np.exp(vs), band(u, bc, vs), color="#52514e", lw=1.2, label="$F_{194}(\\log p)$")
    p = primes_for(c, 10 ** 4).astype(float)
    ax2.plot(p, band(u, bc, np.log(p)), ".", ms=3, color="#2a78d6", rasterized=True,
             label="소수 $p\\nmid194$")
    ax2.axvline(math.exp(bc / u), color="#c3c2b7", ls="--", lw=1)
    ax2.annotate("$p^*\\approx58$", (math.exp(bc / u), 0.385), xytext=(-4, 0),
                 textcoords="offset points", ha="right", color="#52514e")
    ax2.axvline(x / c, color="#c4473a", ls=":", lw=1)
    ax2.annotate("$X/c$, $X=20002$", (x / c, 0.405), xytext=(4, 0),
                 textcoords="offset points", color="#c4473a")
    ax2.set_xscale("log")
    ax2.set_xlabel("$p$")
    ax2.set_ylabel("$M_2(194p)$")
    ax2.legend(frameon=False, fontsize=8, loc="upper right")
    for ax in (ax1, ax2):
        ax.grid(axis="y", color="#e1e0d9", lw=1)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, facecolor="white", dpi=200)
    print(f"figure written to {path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--figure", help="write the band figure (SVG)")
    args = parser.parse_args()
    check_exact()
    check_shape()
    table_families()
    check_pairs()
    check_near_pairs()
    check_predictor()
    check_asymptotics()
    check_additivity()
    check_shifted_value()
    check_small_bands()
    check_mean()
    if args.figure:
        write_figure(args.figure)


if __name__ == "__main__":
    main()
