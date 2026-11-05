"""Checks for docs/chowla.md. Requires numpy and scipy.

Run: python -B scripts/check_chowla.py          (N = 2^22, a few minutes)
     python -B scripts/check_chowla.py --full   (N = 2^25, the scale quoted in the document)

The Moebius function and smallest prime factors are sieved up to N. Zeros of
zeta enter only through their known ordinates in the spectral check. Exact
identities are asserted exactly; statistical checks use tolerances in units of
the empirical standard error. These are consistency checks, not a proof.
"""

import itertools
import math
import sys
from fractions import Fraction

import numpy as np
from scipy import sparse, stats
from scipy.special import poch, rgamma

GAMMAS = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
          37.586178, 40.918719, 43.327073, 48.005151, 49.773832]


# ---------------------------------------------------------------- sieves

def sieve(n):
    """mu(0..n) as int8 and the smallest prime factor spf(0..n) as int32 (spf(p) = p)."""
    spf = np.zeros(n + 1, dtype=np.int32)
    for p in range(2, int(n ** 0.5) + 1):
        if spf[p] == 0:
            block = spf[p * p::p]
            block[block == 0] = p
    rest = np.nonzero(spf == 0)[0]
    spf[rest] = rest.astype(np.int32)
    mu = np.ones(n + 1, dtype=np.int8)
    mu[0] = 0
    for p in np.nonzero(spf[:n + 1] == np.arange(n + 1))[0]:
        if p < 2:
            continue
        mu[p::p] *= -1
        if p * p <= n:
            mu[p * p::p * p] = 0
    return mu, spf


def prime_list(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.nonzero(s)[0]


def densities(h_max):
    """Mirsky density S2(h) and primitive density delta_flat(h) for 0 <= h <= h_max."""
    ps = prime_list(max(h_max, 10 ** 6))
    base = float(np.prod(1 - 2.0 / ps.astype(float) ** 2))
    s2 = np.full(h_max + 1, base)
    flat = np.full(h_max + 1, base)
    for p in ps[ps <= h_max]:
        flat[::p] *= (1 - 1 / p) / (1 - 2 / p ** 2)
        if p * p <= h_max:
            s2[::p * p] *= (1 - 1 / p ** 2) / (1 - 2 / p ** 2)
    return s2, flat, base


def prime_factors(n):
    out, p = [], 2
    while p * p <= n:
        if n % p == 0:
            out.append(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        out.append(n)
    return out


# ---------------------------------------------------------------- section 1

def check_exact_identities(mu):
    """(1.2) dilation identity and (1.3) gcd decomposition, exactly."""
    m = mu.astype(np.int64)
    x = 200_000
    for p in (2, 3, 5):
        for h in (1, 2, 6, 7, 12, 30):
            n = np.arange(p, x + 1, p)
            left = int((m[n] * m[n + p * h]).sum())
            k = np.arange(1, x // p + 1)
            ok = (k % p != 0) & ((k + h) % p != 0)
            right = int((m[k] * m[k + h] * ok).sum())
            assert left == right, (p, h)
    for h in (1, 2, 6, 12, 30, 60, 210):
        lhs = int((m[1:x + 1] * m[1 + h:x + 1 + h]).sum())
        rhs = 0
        for g in range(1, h + 1):
            if h % g or any(g % (q * q) == 0 for q in range(2, int(g ** 0.5) + 1)):
                continue
            k = h // g
            mm = np.arange(1, x // g + 1)
            ok = (np.gcd(mm, k) == 1) & (np.gcd(g, mm * (mm + k)) == 1)
            rhs += int((m[mm] * m[mm + k] * ok).sum())
        assert lhs == rhs, h
    print("(1.2) dilation identity, (1.3) gcd decomposition: exact for all tested h")


def check_square_coincidences(n_max=5000, h_max=60):
    """Proposition 1.4: n(n+h) = n'(n'+h') has at most tau(|h^2-h'^2|)/2 solutions."""
    n = np.arange(1, n_max + 1, dtype=np.int64)
    values = {h: n * (n + h) for h in range(1, h_max + 1)}
    worst, total = 0.0, 0
    for h in range(1, h_max + 1):
        for g in range(h + 1, h_max + 1):
            count = len(np.intersect1d(values[h], values[g]))
            diff = g * g - h * h
            tau = sum(1 for d in range(1, int(diff ** 0.5) + 1) if diff % d == 0) * 2
            assert count <= tau // 2, (h, g)
            worst = max(worst, count / (tau / 2))
            total += count
    print(f"(1.4) n(n+h)=n'(n'+h'), h<h'<={h_max}, n<={n_max}: {total} solutions in all, "
          f"max count/(tau/2) = {worst:.2f}")


def autocorrelation(a, h_max, x_cut):
    """C_h = sum_{n <= x_cut} a(n) a(n+h) for 0 <= h <= h_max (a indexed from 0)."""
    size = 1 << int(np.ceil(np.log2(2 * len(a))))
    short = np.zeros(len(a))
    short[:x_cut + 1] = a[:x_cut + 1]
    out = np.fft.irfft(np.conj(np.fft.rfft(short, size)) * np.fft.rfft(a, size), size)
    return np.rint(out[:h_max + 1])


def check_variance_law(mu, h_max):
    """Mirsky variance law and the square-class covariance lattice (section 1)."""
    big_n = len(mu) - 1
    x = big_n - h_max
    c = autocorrelation(mu.astype(np.float64), h_max, x)
    s2, _, _ = densities(h_max)
    h = np.arange(h_max + 1)
    z = c / np.sqrt(s2 * x)
    z[0] = np.nan
    sel = h >= 1
    print(f"(1.5) N={big_n}, h<={h_max}: Var C_h/sqrt(S2 x) = {np.nanvar(z):.3f}, "
          f"kurtosis {stats.kurtosis(z[sel]):+.3f}")
    classes = {"h odd": h % 2 == 1, "4|h": h % 4 == 0, "9|h": h % 9 == 0,
               "3||h": (h % 3 == 0) & (h % 9 != 0)}
    for name, mask in classes.items():
        mask = mask & sel
        raw = np.mean((c[mask] / np.sqrt(x)) ** 2)
        ratio = raw / s2[mask].mean()
        print(f"    {name:5s}: E[C^2/x] = {raw:.4f}, E[S2] = {s2[mask].mean():.4f}, ratio {ratio:.3f}")
        assert abs(ratio - 1) < 4 * math.sqrt(2 / mask.sum()) + 0.02

    def lam(p, k):
        if k % (p * p) == 0:
            return (1 - 1 / p) / (1 - 1 / p ** 2)
        if k % p == 0:
            return (1 - 1 / p) / (1 - 2 / p ** 2)
        return (1 - 2 / p) / (1 - 2 / p ** 2)

    print("(1.6) covariance lattice E[z_ak z_bk]: empirical / predicted")
    for a, b in ((1, 2), (1, 3), (2, 3), (1, 5), (1, 4), (2, 9)):
        ks = np.arange(1, h_max // max(a, b) + 1)
        emp = float(np.mean(z[a * ks] * z[b * ks]))
        squarefree = all(a % (q * q) and b % (q * q) for q in (2, 3, 5, 7))
        if squarefree and math.gcd(a, b) == 1:
            pf = prime_factors(a * b)
            pred = float(np.mean([s2[k] * np.prod([lam(p, k) for p in pf]) / max(a, b)
                                  / math.sqrt(s2[a * k] * s2[b * k]) for k in ks]))
        else:
            pred = 0.0
        se = 1 / math.sqrt(len(ks))
        print(f"    (a,b)=({a},{b}): {emp:+.4f} / {pred:+.4f}   (SE {se:.4f})")
        assert abs(emp - pred) < 4 * se + 0.015
    return c


def primitive_correlations(mu, h_max):
    """C_flat_h = sum_{n <= N-h_max, (n,h)=1} mu(n)mu(n+h) via Moebius-inverted FFTs."""
    big_n = len(mu) - 1
    x = big_n - h_max
    muf = mu.astype(np.float64)
    out = np.zeros(h_max + 1)
    for d in range(1, h_max + 1):
        if mu[d] == 0:
            continue
        m_all, m_short, k_max = big_n // d, x // d, h_max // d
        b = muf[d * np.arange(1, m_all + 1)]
        if m_all > 4096:
            size = 1 << int(np.ceil(np.log2(2 * m_all)))
            ac = np.fft.irfft(np.conj(np.fft.rfft(b[:m_short], size)) * np.fft.rfft(b, size),
                              size)[:k_max + 1]
        else:
            ac = np.array([np.dot(b[:m_short], b[k:m_short + k]) for k in range(k_max + 1)])
        out[d * np.arange(1, k_max + 1)] += int(mu[d]) * np.rint(ac[1:k_max + 1])
    return out, x


def check_primitive(mu, h_max):
    """Section 3: primitive parts have variance delta_flat x and no cross-shift covariance."""
    m = mu.astype(np.int64)
    cp, x = primitive_correlations(mu, h_max)
    for h in (1, 6, 30, 210):
        n = np.arange(1, x + 1)
        direct = int((m[1:x + 1] * m[1 + h:x + 1 + h] * (np.gcd(n, h) == 1)).sum())
        assert direct == int(cp[h])
    _, flat, _ = densities(h_max)
    h = np.arange(h_max + 1)
    z = cp / np.sqrt(flat * x)
    z[0] = np.nan
    zz = z[1:]
    print(f"(3.1) primitive parts, h<={h_max}: Var {np.var(zz):.4f}, kurtosis {stats.kurtosis(zz):+.3f}")
    for name, mask in {"h odd": h % 2 == 1, "4|h": h % 4 == 0, "3|h": h % 3 == 0,
                       "30|h": h % 30 == 0}.items():
        mask = mask & (h >= 1)
        v = float(np.mean(z[mask] ** 2))
        print(f"    {name:5s}: E[z^2] = {v:.4f}")
        assert abs(v - 1) < 4 * math.sqrt(2 / mask.sum()) + 0.03
    for a in (2, 3, 4, 5):
        ks = np.arange(1, h_max // a + 1)
        v = z[ks] * z[a * ks]
        se = v.std() / math.sqrt(len(ks))
        print(f"    E[z_k z_{a}k] = {v.mean():+.4f} (SE {se:.4f})")
        assert abs(v.mean()) < 4 * se + 0.005


def trajectories(mu, h_max, primitive, grid=600, x_min=1e3):
    """C_h(x) on a log grid of x, for 1 <= h <= h_max (optionally gcd(n,h)=1)."""
    big_n = len(mu) - 1
    xs = np.unique(np.rint(np.exp(np.linspace(math.log(x_min), math.log(big_n - h_max - 1),
                                              grid))).astype(np.int64))
    starts = np.concatenate([[0], xs[:-1]])
    out = np.zeros((h_max, len(xs)))
    for h in range(1, h_max + 1):
        prod = mu[1:big_n + 1 - h] * mu[1 + h:big_n + 1]
        if primitive:
            for p in prime_factors(h):
                prod[p - 1::p] = 0
        seg = np.add.reduceat(prod[:xs[-1]], starts, dtype=np.int64)
        out[h - 1] = np.cumsum(seg)
    return xs.astype(float), out


def spectrum(f, y, ts):
    w = 0.5 - 0.5 * np.cos(2 * math.pi * (y - y[0]) / (y[-1] - y[0]))
    dy = np.gradient(y)
    g = (f - (f * w * dy).sum(-1, keepdims=True) / (w * dy).sum()) * w * dy
    return np.abs(g @ np.exp(-1j * np.outer(y, ts))) ** 2


def check_no_zero_lines(mu, h_max):
    """Section 2: zero frequencies appear in M(x)/sqrt(x) but not in C_h(x)."""
    xs, c = trajectories(mu, h_max, primitive=False)
    s2, _, _ = densities(h_max)
    y = np.log(xs)
    m_x = np.cumsum(mu.astype(np.int64))[xs.astype(np.int64)]
    ts = np.linspace(0.5, 52, 1031)
    pc = spectrum(c / np.sqrt(s2[1:h_max + 1, None] * xs[None, :]), y, ts).mean(0)
    pm = spectrum(m_x / np.sqrt(xs), y, ts)

    def ratio(p, g):
        sel = np.abs(ts - g) < 3
        for gg in GAMMAS:
            sel &= np.abs(ts - gg) > 1
        return p[np.argmin(np.abs(ts - g))] / np.median(p[sel])

    rm = [ratio(pm, g) for g in GAMMAS[:6]]
    rc = [ratio(pc, g) for g in GAMMAS[:6]]
    print(f"(2.3) line/background at gamma_1..6: M {np.round(rm, 1).tolist()}")
    print(f"      Chowla mean over h<={h_max}: {np.round(rc, 2).tolist()}")
    assert min(rm) > 10 and max(rc) < 1.6


def check_ou(mu, h_max, seed=1):
    """Section 3: primitive trajectories have independent increments and Brownian maxima."""
    xs, c = trajectories(mu, h_max, primitive=True)
    _, flat, _ = densities(h_max)
    h = np.arange(1, h_max + 1)
    y = np.log(xs)
    z = c / np.sqrt(flat[h][:, None] * xs[None, :])
    x0 = xs[0] * 10
    edges = sorted({int(np.argmin(np.abs(xs - x0 * 2 ** j))) for j in range(40)
                    if x0 * 2 ** j <= xs[-1]})
    xe = xs[edges]
    inc = np.diff(c[:, edges], axis=1) / np.sqrt(flat[h][:, None] * np.diff(xe)[None, :])
    adj = float(np.mean(inc[:, :-1] * inc[:, 1:]))
    se = 1 / math.sqrt(inc[:, :-1].size)
    m_x = np.cumsum(mu.astype(np.int64))[xs.astype(np.int64)]
    m_inc = np.diff(m_x[edges]) / np.sqrt(np.diff(xe))
    m_adj = float(np.corrcoef(m_inc[:-1], m_inc[1:])[0, 1])
    print(f"(3.2) adjacent dyadic increments: primitive corr {adj:+.4f} (SE {se:.4f}), "
          f"M(x) corr {m_adj:+.3f}")
    assert abs(adj) < 4 * se
    j0 = int(np.searchsorted(y, y[0] + 2.3))
    step = float(np.mean(np.diff(y)))
    rows = []
    for dy in (0.5, 1.0, 2.0, 3.0):
        k = int(round(dy / step))
        a, b = z[:, j0:-k], z[:, j0 + k:]
        r = np.mean(a * b) / math.sqrt(np.mean(a * a) * np.mean(b * b))
        per_h = np.mean(a * b, axis=1) / math.sqrt(np.mean(a * a) * np.mean(b * b))
        rows.append((dy, r, per_h.std() / math.sqrt(len(per_h))))
    print("      log-x autocorrelation vs exp(-dy/2): "
          + ", ".join(f"{dy}: {r:.3f}/{math.exp(-dy / 2):.3f} (SE {se:.3f})" for dy, r, se in rows))
    for dy, r, se in rows:
        assert abs(r - math.exp(-dy / 2)) < 4 * se + 0.02
    rng = np.random.default_rng(seed)
    dx = np.diff(np.concatenate([[0.0], xs]))
    sims = np.concatenate([np.abs(np.cumsum(rng.standard_normal((1000, len(xs))) * np.sqrt(dx),
                                            axis=1)).max(1) for _ in range(20)]) / math.sqrt(xs[-1])
    emp = np.abs(c).max(1) / np.sqrt(flat[h] * xs[-1])
    pval = stats.ks_2samp(emp, sims).pvalue
    print(f"      sup|C_flat|/sqrt(delta x): mean {emp.mean():.3f} (Brownian {sims.mean():.3f}), "
          f"KS p = {pval:.2f}")
    assert pval > 0.001


# ---------------------------------------------------------------- section 2

def check_pole_functional():
    """(2.2): <K_{-1/2,-1/2}, u^k> = -(1/4) <K''_{3/2,3/2}, u^k> for all k."""
    worst = 0.0
    for k in range(0, 25):
        left = poch(-0.5, k) * rgamma(k - 1)
        right = -0.25 * k * (k - 1) * (poch(1.5, k - 2) if k >= 2 else 0.0) * rgamma(k + 1)
        worst = max(worst, abs(left - right))
    print(f"(2.2) K_(-1/2,-1/2) = -K''_(3/2,3/2)/4 on moments k<25: max error {worst:.1e}")
    assert worst < 1e-12


# ---------------------------------------------------------------- section 4

def check_short_interval_lemma(seed=5):
    """Theorem 4.1: |T log(n1 n2/n3 n4)/(2 pi)| <= 1/2 iff m1+m2 = m3+m4, else >= 1."""
    rng = np.random.default_rng(seed)
    for x, y, full in ((4000, 31, True), (10 ** 7, 1581, False), (10 ** 12, 500_000, False)):
        if full:
            m = np.arange(1, y + 1)
            g = np.stack(np.meshgrid(m, m, m, m, indexing="ij"), -1).reshape(-1, 4)
        else:
            g = rng.integers(1, y + 1, (2_000_000, 4))
            g[:1_000_000, 3] = g[:1_000_000, 0] + g[:1_000_000, 1] - g[:1_000_000, 2]
            g = g[(g[:, 3] >= 1) & (g[:, 3] <= y)]
        n = (x + g).astype(np.float64)
        same = g[:, 0] + g[:, 1] == g[:, 2] + g[:, 3]
        lam = np.abs(np.log1p((n[:, 0] * n[:, 1] - n[:, 2] * n[:, 3]) / (n[:, 2] * n[:, 3])))
        for t in (8 * math.pi * x / 3 * 1.001, 3 * math.pi * x, 4 * math.pi * x):
            v = t * lam / (2 * math.pi)
            assert v[same].max() <= 0.5 and v[~same].min() >= 1.0
        print(f"(4.1) X={x}, Y={y}: {len(g)} quadruples, dichotomy holds for T in (8piX/3, 4piX]")


def intervals(big_n, x0, y):
    return np.arange(x0, big_n - y - 1, y)


def check_m4(mu, x0, seed=3):
    """Theorem 4.1 numerics, (4.4) Wick formula with the mu^2 spectrum, and (4.5)."""
    big_n = len(mu) - 1
    y = int(math.sqrt(x0) / 2)
    starts = intervals(big_n, x0, y)
    muf = mu.astype(np.float64)
    named = [0.0, 0.25, 0.5, 1 / 9, 1 / 3, 1 / 36, 1 / 25, 0.2, 1 / 6, 0.125, 0.318309886]
    thetas = np.concatenate([named, np.linspace(0.001, 0.999, 499)])
    k = np.arange(-(y - 1), y)
    e_k = np.exp(-2j * np.pi * np.outer(thetas, k))
    e_m = np.exp(2j * np.pi * np.outer(thetas, np.arange(1, y + 1)))
    size = 1 << int(np.ceil(np.log2(2 * y)))
    r_sum = np.zeros(len(thetas))
    p_sum = np.zeros(len(thetas))
    d2 = np.zeros(len(thetas))
    cs = []
    for s in starts:
        a = muf[s + 1:s + y + 1]
        f = np.fft.rfft(a, size)
        ac = np.fft.irfft(f * np.conj(f), size)
        c = np.concatenate([ac[size - (y - 1):], ac[:y]])
        q = c[y - 1]
        c2 = c ** 2
        c2[y - 1] = 0
        r = (e_k @ c2).real / q ** 2
        p = (np.abs(e_m @ (a * a)) ** 2 - q) / q ** 2
        r_sum += r
        p_sum += p
        d2 += (r - p) ** 2
        if len(cs) < 6000:
            cs.append(c)
    n = len(starts)
    r_mean, p_mean = r_sum / n, p_sum / n
    se = np.sqrt(np.maximum(d2 / n - (r_mean - p_mean) ** 2, 1e-30) / n)
    print(f"(4.4) Y={y}, {n} intervals: R(theta) = sum_(k!=0) C_k^2 e(k theta)/Q^2 vs "
          f"(|mu2^(theta)|^2 - Q)/Q^2")
    labels = ["0", "1/4", "1/2", "1/9", "1/3", "1/36", "1/25", "1/5", "1/6", "1/8", "0.318"]
    for i, lab in enumerate(labels):
        print(f"    theta={lab:5s}: {r_mean[i]:+.5f} / {p_mean[i]:+.5f}")
        assert abs(r_mean[i] - p_mean[i]) < 4 * se[i] + 2e-4
    zg = (r_mean[len(named):] - p_mean[len(named):]) / se[len(named):]
    print(f"    dense grid of 499 theta: rms z = {np.sqrt(np.mean(zg ** 2)):.2f}, "
          f"max|z| = {np.abs(zg).max():.2f}")
    assert np.sqrt(np.mean(zg ** 2)) < 1.5
    cross = {}
    for gap in (1, 2, 5):
        v = []
        for i in range(len(cs) - gap):
            ci, cj = cs[i], cs[i + gap]
            prod = ci * cj
            prod[y - 1] = 0
            v.append(prod.sum() / (ci[y - 1] * cj[y - 1]))
        v = np.array(v)
        cross[gap] = (v.mean(), v.std() / math.sqrt(len(v)))
    print("(4.5) disjoint intervals, sum_(k!=0) C_k(I)C_k(J)/(Q_I Q_J): "
          + ", ".join(f"gap {g}: {m:+.4f} (SE {s:.4f})" for g, (m, s) in cross.items()))
    for m, s in cross.values():
        assert abs(m) < 4 * s
    # Dirichlet polynomial side: fourth moment and an echo at a height offset.
    rng = np.random.default_rng(seed)
    rows = []
    for x in rng.integers(max(x0, big_n // 2), big_n - y - 1, 12):
        nn = np.arange(x + 1, x + y + 1)
        a = muf[nn]
        keep = a != 0
        nn, a = nn[keep], a[keep]
        q = a @ a
        logn = np.log(nn.astype(float))
        t = rng.uniform(3 * math.pi * x, 6 * math.pi * x, 2000)
        p0 = np.exp(-1j * np.outer(t, logn)) @ a
        p1 = np.exp(-1j * np.outer(t + 2 * math.pi * x / 4, logn)) @ a
        rows.append((np.mean(np.abs(p0) ** 2) / q, np.mean(np.abs(p0) ** 4) / q ** 2,
                     np.mean(np.abs(p0) ** 2 * np.abs(p1) ** 2) / q ** 2 - 1))
    rows = np.array(rows)
    print(f"(4.1) P_I side: <|P|^2>/Q = {rows[:, 0].mean():.3f}, <|P|^4>/Q^2 = {rows[:, 1].mean():.3f}"
          f" (Wick 2), echo at offset 2 pi X/4: {rows[:, 2].mean():+.3f} (predicted ~0.11)")
    assert abs(rows[:, 0].mean() - 1) < 0.05 and abs(rows[:, 1].mean() - 2) < 0.25


# ---------------------------------------------------------------- section 4.4

D_POLY = {2: (2, -1, 0), 3: (6, -9, 4, 0), 4: (24, -72, 82, -33, 0)}


def diag_count(q, k):
    """D_k(Q) = (k!)^2 [z^k] (sum_{m<=k} z^m/(m!)^2)^Q, computed exactly."""
    base = [Fraction(1, math.factorial(m) ** 2) for m in range(k + 1)]
    out = [Fraction(1)] + [Fraction(0)] * k
    for _ in range(q):
        out = [sum(out[i] * base[j - i] for i in range(j + 1)) for j in range(k + 1)]
    return int(out[k] * math.factorial(k) ** 2)


def diag_poly(q, k):
    return np.polyval(np.array(D_POLY[k], dtype=np.float64), np.asarray(q, dtype=np.float64))


def check_moment_ladder(mu):
    """Theorem 4.5: rung r separates equal power sums p_1..p_r; top rung is a rearrangement."""
    cases = ((2, 1, 4000, 31), (3, 1, 4000, 29), (4, 1, 20000, 45),
             (3, 2, 10 ** 6, 50), (4, 2, 10 ** 7, 60), (4, 3, 10 ** 9, 25))
    for k, r, x, y in cases:
        ms = np.array(list(itertools.combinations_with_replacement(range(1, y + 1), k)), dtype=np.int64)
        base = x ** k
        v = [math.prod(int(x + m) for m in row) - base for row in ms]
        keys = np.stack([(ms ** j).sum(1) for j in range(1, r + 1)], 1)
        eta = k * r * y ** (r + 1) / ((r + 1) * (x - y))
        kept, spread = {}, 0.0
        for key, val in zip(map(tuple, keys.tolist()), v):
            lo, hi = kept.get(key, (val, val))
            kept[key] = (min(lo, val), max(hi, val))
        for lo, hi in kept.values():
            spread = max(spread, math.log1p((hi - lo) / (lo + base)))
        order = sorted(range(len(v)), key=v.__getitem__)
        ks = keys[order]
        vs = [v[i] for i in order]
        cut = np.nonzero(np.any(ks[1:] != ks[:-1], axis=1))[0]
        gap = min(math.log1p((vs[i + 1] - vs[i]) / (vs[i] + base)) for i in cut)
        print(f"(4.7) k={k}, r={r}, X={x}, Y={y}: eta={eta:.3f}; equal p_1..p_r: |lambda| X^r <= "
              f"{spread * x ** r:.2e} (bound {eta / r:.2e}); otherwise >= {gap * x ** r:.3f} "
              f"(bound {(1 - eta) / r:.3f})")
        assert eta < 1 / 3 and spread * x ** r <= eta / r and gap * x ** r >= (1 - eta) / r
    x, y, k = 10 ** 6, 500, 3
    vals = [int(n) for n in range(x + 1, x + y + 1) if mu[n] != 0]
    seen, bad = {}, 0
    for c in itertools.combinations_with_replacement(vals, k):
        p = math.prod(c)
        bad += p in seen
        seen[p] = c
    print(f"(4.8) X={x}, Y={y} (Y^3 > X): {len(vals)} squarefree n, {bad} products n1n2n3 = n4n5n6 "
          f"that are not rearrangements")
    assert bad == 0


def brute_sign_mean(a_set, k, r):
    """E_eps V_{k,r} for independent signs on a_set: count (n,n') with equal p_1..p_r, even multiplicities."""
    idx = {v: i for i, v in enumerate(a_set)}
    groups = {}
    for t in itertools.product(a_set, repeat=k):
        par = [0] * len(a_set)
        for v in t:
            par[idx[v]] ^= 1
        key = (tuple(sum(v ** j for v in t) for j in range(1, r + 1)), tuple(par))
        groups[key] = groups.get(key, 0) + 1
    return sum(c * c for c in groups.values())


def additive_counts(a_set):
    """E*(A) (distinct a+b=c+d, ordered) and AP3(A) = #{(a,b,c): b != c, b+c = 2a}."""
    a = np.zeros(max(a_set) + 1, dtype=np.int64)
    a[list(a_set)] = 1
    r2 = np.convolve(a, a)
    q = len(a_set)
    ap3 = int(sum(r2[2 * v] for v in a_set)) - q
    return int((r2 ** 2).sum()) - (2 * q * q - q) - 2 * ap3, ap3


def shifted_model3(a2, thetas):
    """Proposition 4.6 (2): perm G - sum |G_jl|^2 + 4Q with G_jl = sum a2_m e(m(theta_j - theta_l))."""
    m = np.arange(a2.shape[-1])
    g = [[a2 @ np.exp(2j * np.pi * m * (tj - tl)) for tl in thetas] for tj in thetas]
    perm = sum(g[0][s[0]] * g[1][s[1]] * g[2][s[2]] for s in itertools.permutations(range(3)))
    return perm.real - sum(np.abs(g[j][l]) ** 2 for j in range(3) for l in range(3)) + 4 * a2.sum(-1)


def check_sign_model(seed=2):
    """Proposition 4.6: exact independent-sign means, by enumeration over small sets."""
    for q, k in ((0, 2), (1, 3), (5, 4), (9, 3), (12, 4)):
        assert diag_count(q, k) == round(diag_poly(q, k))
    rng = np.random.default_rng(seed)
    for _ in range(5):
        a_set = sorted(rng.choice(np.arange(1, 41), int(rng.integers(5, 9)), replace=False).tolist())
        q = len(a_set)
        e_star, ap3 = additive_counts(a_set)
        for k, r in ((2, 1), (3, 1), (3, 2), (4, 1), (4, 2)):
            want = diag_count(q, k) + (9 * e_star + 6 * ap3 if (k, r) == (4, 1) else 0)
            assert brute_sign_mean(a_set, k, r) == want, (a_set, k, r)
    worst = 0.0
    for _ in range(4):
        y, size = 24, 128
        a2 = np.zeros(y)
        a2[rng.choice(y, 9, replace=False)] = 1
        supp = np.nonzero(a2)[0]
        th = [0.0] + rng.random(2).tolist()
        mean = 0.0
        for signs in itertools.product((-1.0, 1.0), repeat=len(supp)):
            a = np.zeros(y)
            a[supp] = signs
            vals = np.ones(size)
            for t in th:
                vals *= np.abs(np.fft.fft(a * np.exp(2j * np.pi * np.arange(y) * t), size)) ** 2
            mean += vals.mean()
        mean /= 2 ** len(supp)
        worst = max(worst, abs(mean - shifted_model3(a2, th)) / mean)
    print(f"(4.9) independent signs: E V_(k,r) = D_k(Q) (+9E*+6AP3 at k=4, r=1) on 5 random sets; "
          f"shifted k=3 formula exact to {worst:.0e}")
    assert worst < 1e-12


def sign_models(mu, spf, seed):
    """mu, independent signs on squarefree n, and a Rademacher random multiplicative function, on [0, N]."""
    big_n = len(mu) - 1
    rng = np.random.default_rng(seed)
    sq = mu != 0
    indep = rng.choice(np.array([-1, 1], dtype=np.int8), size=big_n + 1) * sq
    eps_p = rng.choice(np.array([-1, 1], dtype=np.int8), size=big_n + 1)
    mult = np.ones(big_n + 1, dtype=np.int8)
    m = np.arange(big_n + 1, dtype=np.int64)
    m[0] = 1
    while True:
        act = np.nonzero(m > 1)[0]
        if len(act) == 0:
            break
        p = spf[m[act]].astype(np.int64)
        mult[act] *= eps_p[p]
        m[act] //= p
    return {"mu": mu, "independent": indep.astype(np.int8), "multiplicative": (mult * sq).astype(np.int8)}


TRIPLES = {"0,0,0": (0.0, 0.0, 0.0), "0,0,1/4": (0.0, 0.0, 0.25), "0,1/4,1/2": (0.0, 0.25, 0.5),
           "0,1/9,2/9": (0.0, 1 / 9, 2 / 9), "0,1/4,1/9": (0.0, 0.25, 1 / 9),
           "generic": (0.0, 0.318309886, 0.123456789)}


def check_m2k(mu, spf, x0, seed=7):
    """Conjecture 4.7: 2k-th moments (k=2,3,4) and shifted sixth moments against the sign model."""
    big_n = len(mu) - 1
    y = int(math.sqrt(x0) / 2)
    starts = intervals(big_n, x0, y)
    seqs = sign_models(mu, spf, seed)
    size = 1 << int(np.ceil(np.log2(4 * y)))
    m = np.arange(y)
    thetas = sorted({t for tr in TRIPLES.values() for t in tr})
    phase = {t: np.exp(2j * np.pi * m * t) for t in thetas}
    sq = mu != 0
    names = list(seqs)
    rows = {(nm, key): [] for nm in names for key in ["k2", "k3", "k4"] + list(TRIPLES)}
    model = {key: [] for key in ["k2", "k3", "k4"] + list(TRIPLES)}
    qs = []
    for b0 in range(0, len(starts), 400):
        idx = starts[b0:b0 + 400, None] + 1 + m[None, :]
        a2 = sq[idx].astype(np.float64)
        q = a2.sum(1)
        qs.extend(q.tolist())
        f2 = np.fft.rfft(a2, size, axis=1)
        r2 = np.rint(np.fft.irfft(f2 * f2, size, axis=1))
        ap3 = np.array([r2[i, 2 * m[a2[i] > 0]].sum() for i in range(len(q))]) - q
        e_star = (r2 ** 2).sum(1) - (2 * q * q - q) - 2 * ap3
        model["k2"].extend(diag_poly(q, 2).tolist())
        model["k3"].extend(diag_poly(q, 3).tolist())
        model["k4"].extend((diag_poly(q, 4) + 9 * e_star + 6 * ap3).tolist())
        for key, tr in TRIPLES.items():
            model[key].extend(shifted_model3(a2, tr).tolist())
        for nm in names:
            a = seqs[nm][idx].astype(np.float64)
            s2 = {t: np.abs(np.fft.fft(a * phase[t], size, axis=1)) ** 2 for t in thetas}
            for k in (2, 3, 4):
                rows[(nm, f"k{k}")].extend((s2[0.0] ** k).mean(1).tolist())
            for key, tr in TRIPLES.items():
                rows[(nm, key)].extend((s2[tr[0]] * s2[tr[1]] * s2[tr[2]]).mean(1).tolist())
    qs = np.array(qs)
    print(f"(4.10) Y={y}, {len(starts)} intervals: int|S_I|^(2k) / exact sign-model mean "
          "(mean, SE, sd over intervals)")
    for k in (2, 3, 4):
        mk = np.array(model[f"k{k}"])
        line = f"    k={k}: "
        for nm in names:
            v = np.array(rows[(nm, f"k{k}")]) / mk
            se = v.std() / math.sqrt(len(v))
            line += f"{nm} {v.mean():.5f} (SE {se:.5f}, sd {v.std():.3f})  "
            assert abs(v.mean() - 1) < 4 * se + 1e-3
        print(line)
    print("(4.10) shifted sixth moments: model/Q^3 and (int prod|S(a+theta_j)|^2 - model)/Q^3")
    for key in TRIPLES:
        mk = np.array(model[key])
        line = f"    theta=({key:9s}): model/Q^3 {np.mean(mk / qs ** 3):.4f};"
        for nm in names:
            d = (np.array(rows[(nm, key)]) - mk) / qs ** 3
            se = d.std() / math.sqrt(len(d))
            line += f"  {nm} {d.mean():+.5f} (SE {se:.5f})"
            assert abs(d.mean()) < 4 * se + 1e-3
        print(line)


def rung1_ratio(seq, starts, y, k):
    """int |S_I|^(2k) / D_k(Q_I) on (s, s+y] for every start s (exact integer coefficients)."""
    size = 1 << int(np.ceil(np.log2(k * y)))
    out = []
    for b0 in range(0, len(starts), 20000):
        idx = starts[b0:b0 + 20000, None] + 1 + np.arange(y)[None, :]
        a = seq[idx].astype(np.float64)
        c = np.rint(np.fft.irfft(np.fft.rfft(a, size, axis=1) ** k, size, axis=1))
        out.append((c * c).sum(1) / diag_poly((a != 0).sum(1), k))
    return np.concatenate(out)


def check_ladder_k3(mu, spf, x0, y, n_max, scan, seed=7):
    """Theorem 4.5 at k=3 on (X, X+Y]: rung 1 int|S|^6, rung 2 (equal p1, p2), rung 3 = D_3(Q)."""
    big_n = len(mu) - 1
    assert 6 * y ** 3 < x0 - y and y ** 3 <= x0
    seqs = sign_models(mu, spf, seed)
    seqs = {"mu": seqs["mu"], "mu^2": (mu != 0).astype(np.int8), "independent": seqs["independent"]}
    starts = intervals(big_n, x0, y)
    sub_starts = starts[::max(1, len(starts) // n_max)]
    ms = np.array(list(itertools.combinations_with_replacement(range(y), 3)), dtype=np.int64)
    perms = np.array([6 // math.prod(math.factorial(c) for c in np.unique(r, return_counts=True)[1])
                      for r in ms], dtype=np.float64)
    _, gid, gsize = np.unique(ms.sum(1) * (3 * y * y + 1) + (ms ** 2).sum(1),
                              return_inverse=True, return_counts=True)
    sub = np.nonzero(gsize[gid] > 1)[0]
    _, gsub = np.unique(gid[sub], return_inverse=True)
    incidence = sparse.csr_matrix((np.ones(len(sub)), (gsub, np.arange(len(sub)))))
    print(f"(4.11) k=3 ladder, Y={y}: rung 1 on all {len(starts)} intervals, rung 2 on {len(sub_starts)}; "
          f"{len(ms)} multisets, {len(sub)} in (p1,p2) classes with more than one member")
    for nm, seq in seqs.items():
        r1 = rung1_ratio(seq, starts, y, 3)
        rung2, qs = [], []
        for b0 in range(0, len(sub_starts), 200):
            idx = sub_starts[b0:b0 + 200, None] + 1 + np.arange(y)[None, :]
            a = seq[idx].astype(np.float64)
            qs.extend((a != 0).sum(1).tolist())
            x = perms[sub] * a[:, ms[sub, 0]] * a[:, ms[sub, 1]] * a[:, ms[sub, 2]]
            rung2.extend((((incidence @ x.T) ** 2).sum(0) - (x * x).sum(1)).tolist())
        r2 = np.array(rung2) / diag_poly(np.array(qs), 3)
        se1, se2 = r1.std() / math.sqrt(len(r1)), r2.std() / math.sqrt(len(r2))
        print(f"    {nm:11s}: rung 1 V31/D3 = {r1.mean():.4f} (SE {se1:.4f}); rung 2 (V32-D3)/D3 = "
              f"{r2.mean():+.1e} (SE {se2:.1e}); rung 3 = 1")
        if nm == "mu^2":
            assert r2.mean() > 0.1
        else:
            assert abs(r1.mean() - 1) < 4 * se1 + 1e-3 and abs(r2.mean()) < 4 * se2 + 1e-5
    line = "    rung 1 on all intervals, other Y (mu / independent):"
    for yy in scan:
        st = intervals(big_n, x0, yy)
        v = {nm: rung1_ratio(seqs[nm], st, yy, 3) for nm in ("mu", "independent")}
        line += f"  Y={yy}: " + " / ".join(f"{x.mean():.4f}" for x in v.values()) \
            + f" (SE {v['mu'].std() / math.sqrt(len(st)):.4f})"
        for x in v.values():
            assert abs(x.mean() - 1) < 4 * x.std() / math.sqrt(len(st)) + 1e-3
    print(line)


# ---------------------------------------------------------------- section 5

def d_half(a):
    """d_{-1/2}(p^a) = (-1/2)_a / a!."""
    return poch(-0.5, a) / math.factorial(a)


def check_kernel_coins():
    """(5.1)-(5.2): K_n = mu(n) nu_n for squarefree n, cosine product, sin^2 coin at p^2."""
    worst = 0.0
    for primes in ((2, 3, 5), (2, 7, 11, 13), (3, 5, 17, 19, 23), (101, 103)):
        n = math.prod(primes)
        r = len(primes)
        log_n = math.log(n)
        for xi in (0.7, 3.0, 6.0):
            lhs = 0j
            for mask in range(1 << r):
                d = math.prod(p for i, p in enumerate(primes) if mask >> i & 1)
                lhs += d_half(1) ** bin(mask).count("1") * d_half(1) ** (r - bin(mask).count("1")) \
                    * np.exp(1j * xi * math.log(d) / log_n)
            mu_n = (-1) ** r
            rhs = mu_n * np.exp(1j * xi / 2) * np.prod([math.cos(xi * math.log(p) / (2 * log_n))
                                                        for p in primes])
            worst = max(worst, abs(lhs - rhs))
    for p, q in ((2, 3), (5, 7), (3, 101)):
        n = p * p * q
        log_n = math.log(n)
        for xi in (1.0, 4.0):
            lhs = 0j
            for i in range(3):
                for j in range(2):
                    w = d_half(i) * d_half(2 - i) * d_half(j) * d_half(1 - j)
                    lhs += w * np.exp(1j * xi * (i * math.log(p) + j * math.log(q)) / log_n)
            phi = xi * math.log(p) / log_n
            coin_p = 0.5 * np.exp(1j * phi) * math.sin(phi / 2) ** 2
            coin_q = -np.exp(1j * xi * math.log(q) / (2 * log_n)) * math.cos(xi * math.log(q) / (2 * log_n))
            worst = max(worst, abs(lhs - coin_p * coin_q))
    print(f"(5.1)-(5.2) kernel = mu nu_n, cosine product, sin^2 coin at p^2: max error {worst:.1e}")
    assert worst < 1e-12


def cosine_shape(lo, hi, spf, xi):
    """prod_{p|n} cos(xi log p / (2 log n)) for lo <= n < hi."""
    n = np.arange(lo, hi, dtype=np.int64)
    log_n = np.log(n.astype(np.float64))
    c = np.ones(len(n))
    m = n.copy()
    while True:
        act = m > 1
        if not act.any():
            break
        p = spf[m[act]].astype(np.int64)
        c[act] *= np.cos(xi * np.log(p.astype(np.float64)) / (2 * log_n[act]))
        m[act] //= p
    return c


def dirichlet_limit(xi, samples=4_000_000, seed=0):
    """E cos(xi S (2B-1)), S ~ Beta(1/2,1/2), B ~ Beta(1/4,1/4): D1-D3 for Dir(1/4,1/2,1/4)."""
    rng = np.random.default_rng(seed)
    s = rng.beta(0.5, 0.5, samples)
    b = rng.beta(0.25, 0.25, samples)
    return float(np.mean(np.cos(xi * s * (2 * b - 1))))


def check_dirichlet_limit(mu, spf):
    """(5.3): mean of prod cos^2 over squarefree n tends to the Dir(1/4,1/2,1/4) value."""
    big_n = len(mu) - 1
    for xi in (3.0, 6.0):
        limit = dirichlet_limit(xi)
        ls, vals = [], []
        k = 3
        while 2 * 10 ** k <= big_n:
            lo, hi = 10 ** k, 2 * 10 ** k
            c = cosine_shape(lo, hi, spf, xi)
            sq = mu[lo:hi] != 0
            vals.append(float(np.mean(c[sq] ** 2)))
            ls.append(math.log(1.5 * lo))
            k += 1
        ls, vals = np.array(ls), np.array(vals)
        design = np.stack([np.ones_like(ls), 1 / ls, 1 / ls ** 2], 1)
        coef = np.linalg.lstsq(design, vals, rcond=None)[0]
        print(f"(5.3) xi={xi}: " + ", ".join(f"L={l:.1f}: {v:.4f}" for l, v in zip(ls, vals))
              + f"; fit a+b/L+c/L^2 gives a = {coef[0]:.4f}, Dirichlet limit {limit:.4f}")
        assert abs(coef[0] - limit) < 0.03
        assert np.all(np.diff(np.abs(vals - limit)) < 0)


def wick_excess(seq, starts, y):
    size = 1 << int(np.ceil(np.log2(2 * y)))
    out = []
    for s in starts:
        a = seq[s:s + y]
        f = np.fft.rfft(a, size)
        c = np.fft.irfft(f * np.conj(f), size)[:y]
        g = np.fft.rfft(a * a, size)
        w = np.fft.irfft(g * np.conj(g), size)[:y]
        d = (c ** 2 - w) / c[0] ** 2
        d[0] = 0
        out.append(2 * d[1:].sum())
    out = np.array(out)
    return out.mean(), out.std() / math.sqrt(len(out))


def check_pole_return(mu, spf, x0):
    """Section 5.3: at xi != 0 the Wick excess is the residue-class (major arc) model."""
    big_n = len(mu) - 1
    y = int(math.sqrt(x0) / 2)
    nn = np.arange(x0, big_n + 1)
    m = mu[x0:].astype(np.float64)
    starts = np.arange(0, len(m) - y - 1, y)
    isprime = spf[x0:] == nn
    for xi in (6.0, 3.0):
        a = m * cosine_shape(x0, big_n + 1, spf, xi)
        raw, raw_se = wick_excess(a, starts, y)
        r = nn % 30
        blocks = np.minimum((np.log(nn / x0) / np.log(big_n / x0) * 4).astype(int), 3)
        key = r * 4 + blocks
        means = np.bincount(key, weights=a, minlength=120) / np.maximum(
            np.bincount(key, minlength=120), 1)
        after, after_se = wick_excess(a - means[key], starts, y)
        cls1 = r == 1
        print(f"(5.4) xi={xi}: mean {a.mean():+.5f}; class 1 mod 30 mean {a[cls1].mean():+.4f} "
              f"(primes alone {np.mean(a[cls1] * isprime[cls1]):+.4f}); Wick excess {raw:+.4f} "
              f"(SE {raw_se:.4f}) -> {after:+.4f} (SE {after_se:.4f}) after the mod 30 model")
        if xi == 6.0:
            cm = np.bincount(r, weights=a, minlength=30) / np.bincount(r, minlength=30)
            coprime = [c for c in range(30) if math.gcd(c, 30) == 1]
            print(f"      class means mod 30: coprime classes {cm[coprime].min():+.3f}..{cm[coprime].max():+.3f}, "
                  f"others {cm[[c for c in range(30) if c not in coprime]].min():+.3f}.."
                  f"{cm[[c for c in range(30) if c not in coprime]].max():+.3f}")
            two = {k: float(np.dot(a[:-k], a[k:]) / (len(a) - k)) for k in (1, 29, 30, 31, 60, 210)}
            print("      E[a_n a_(n+k)]: " + ", ".join(f"k={k}: {v:+.1e}" for k, v in two.items()))
        if xi == 6.0:
            # The excess grows roughly like Y; at N = 2^22 (Y = 500) it is about 3 SE.
            assert raw > 2.5 * raw_se
        assert abs(after) < 3 * after_se + 2e-4


def check_kernel_m4(mu, spf, x0):
    """(5.3): kernel-resolved Wick formula for pairs (xi, eta) at theta = 0, 1/4, 1/9."""
    big_n = len(mu) - 1
    y = int(math.sqrt(x0) / 2)
    m = mu[x0:].astype(np.float64)
    seqs = {0.0: m, 3.0: m * cosine_shape(x0, big_n + 1, spf, 3.0),
            6.0: m * cosine_shape(x0, big_n + 1, spf, 6.0)}
    starts = np.arange(0, len(m) - y - 1, y)
    thetas = [0.0, 0.25, 1 / 9]
    k = np.arange(-(y - 1), y)
    e_k = np.exp(-2j * np.pi * np.outer(thetas, k))
    e_m = np.exp(2j * np.pi * np.outer(thetas, np.arange(y)))
    size = 1 << int(np.ceil(np.log2(2 * y)))

    def auto(a):
        f = np.fft.rfft(a, size)
        ac = np.fft.irfft(f * np.conj(f), size)
        return np.concatenate([ac[size - (y - 1):], ac[:y]])

    pairs = [(0.0, 3.0), (3.0, 3.0), (6.0, 6.0), (3.0, 6.0)]
    acc = {p: np.zeros((3, len(thetas))) for p in pairs}
    for s in starts:
        blocks = {xi: seq[s:s + y] for xi, seq in seqs.items()}
        cs = {xi: auto(b) for xi, b in blocks.items()}
        for (x1, x2) in pairs:
            a, b = blocks[x1], blocks[x2]
            qa, qb = a @ a, b @ b
            prod = cs[x1] * cs[x2]
            prod[y - 1] = 0
            lhs = (e_k @ prod).real / (qa * qb)
            st = a * b
            pred = (np.abs(e_m @ st) ** 2 - np.sum(st ** 2)) / (qa * qb)
            acc[(x1, x2)] += np.stack([lhs, pred, (lhs - pred) ** 2])
    n = len(starts)
    print(f"(5.3) kernel Wick formula, Y={y}, {n} intervals: lhs / formula at theta = 0, 1/4, 1/9")
    for p in pairs:
        lhs, pred, sq = acc[p] / n
        se = np.sqrt(np.maximum(sq - (lhs - pred) ** 2, 1e-30) / n)
        print(f"    (xi,eta)=({p[0]:.0f},{p[1]:.0f}): "
              + "   ".join(f"{lhs[i]:+.4f}/{pred[i]:+.4f}" for i in range(3))
              + "   max|diff|/SE " + f"{np.max(np.abs(lhs - pred) / se):.1f}")
        if p == (0.0, 3.0):
            assert np.all(np.abs(lhs - pred) < 4 * se + 2e-4)


def check_sign_models(mu, spf, x0, seed=11):
    """Section 5.3: the xi = 6 excess belongs to mu, not to the sign models with the same shapes."""
    big_n = len(mu) - 1
    y = int(math.sqrt(x0) / 2)
    rng = np.random.default_rng(seed)
    nn = np.arange(x0, big_n + 1, dtype=np.int64)
    sq = mu[x0:] != 0
    shape = cosine_shape(x0, big_n + 1, spf, 6.0)
    eps_p = rng.choice([-1.0, 1.0], size=big_n + 1)
    mult = np.ones(len(nn))
    m = nn.copy()
    while True:
        act = m > 1
        if not act.any():
            break
        p = spf[m[act]].astype(np.int64)
        mult[act] *= eps_p[p]
        m[act] //= p
    models = {"mu": mu[x0:].astype(np.float64),
              "independent signs": rng.choice([-1.0, 1.0], size=len(nn)) * sq,
              "random multiplicative": mult * sq}
    starts = np.arange(0, len(nn) - y - 1, y)
    out = {}
    for name, sign in models.items():
        out[name] = wick_excess(sign * shape, starts, y)
    print("(5.4) xi=6 Wick excess with the same cosine shapes: "
          + ", ".join(f"{k} {v[0]:+.4f} (SE {v[1]:.4f})" for k, v in out.items()))
    assert out["mu"][0] > 2.5 * out["mu"][1]
    for name in ("independent signs", "random multiplicative"):
        assert abs(out[name][0]) < 4 * out[name][1] + 2e-4


# ---------------------------------------------------------------- main

def main():
    full = "--full" in sys.argv
    big_n = 2 ** 25 if full else 2 ** 22
    mu, spf = sieve(big_n)
    print(f"N = {big_n}, M(N) = {int(mu.astype(np.int64).sum())}")
    check_exact_identities(mu)
    check_square_coincidences()
    check_variance_law(mu, 200_000 if full else 20_000)
    check_pole_functional()
    check_no_zero_lines(mu, 1000 if full else 300)
    check_primitive(mu, 100_000 if full else 10_000)
    check_ou(mu, 1000 if full else 300)
    check_short_interval_lemma()
    check_m4(mu, 10_000_000 if full else 1_000_000)
    check_moment_ladder(mu)
    check_sign_model()
    check_m2k(mu, spf, 10_000_000 if full else 1_000_000)
    check_ladder_k3(mu, spf, 10_000_000 if full else 1_000_000, 80 if full else 50, 20_000,
                    (20, 200) if full else (20, 120))
    check_kernel_coins()
    check_dirichlet_limit(mu, spf)
    x_kernel = 10_000_000 if full else 1_000_000
    check_kernel_m4(mu, spf, x_kernel)
    check_pole_return(mu, spf, x_kernel)
    check_sign_models(mu, spf, x_kernel)


if __name__ == "__main__":
    main()
