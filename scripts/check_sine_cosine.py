"""Checks and figures for docs/sine-cosine-rh.md. Requires numpy and scipy;
matplotlib only with --figures.

Run: python -B scripts/check_sine_cosine.py [--high] [--figures]

Zeta values off the critical line come from the Euler-Maclaurin formula with
about |Im s|/2 terms. On the critical line above height 200 the Riemann-Siegel
formula with the corrections C0, C1, C2 is used. Zeros are located from sign
changes of Hardy's Z. Logarithmic singularities at zeros are integrated with
exact cell averages. --high adds heights 10^3-10^5 (about three minutes) and
the window statistics up to 10^8. These are consistency checks, not a proof.
"""

import argparse
import math
import os

import numpy as np
from scipy.integrate import quad, simpson
from scipy.special import loggamma

from check_ng_tail import BERNOULLI, zeros_up_to

SIGMA_1 = 3.0
FIG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "docs", "figures")


# --------------------------------------------------------------------------
# zeta and Hardy's Z
# --------------------------------------------------------------------------

def zeta(s, chunk=2000):
    """Euler-Maclaurin zeta with about |Im s|/2 terms, so that |s|/(2 pi N) <= 1/pi
    and ten Bernoulli corrections suffice also at heights of several thousand."""
    s = np.atleast_1d(np.asarray(s, dtype=complex))
    out = np.empty(s.shape, dtype=complex)
    for i in range(0, len(s), chunk):
        ss = s[i:i + chunk]
        big_n = int(np.abs(ss.imag).max() / 2) + 30
        log_n = np.log(np.arange(1, big_n, dtype=float))
        total = np.zeros(len(ss), dtype=complex)
        for j in range(0, len(log_n), 400):
            total += np.exp(-np.outer(ss, log_n[j:j + 400])).sum(axis=1)
        total += big_n ** (1 - ss) / (ss - 1) + 0.5 * big_n ** (-ss)
        rising, power = ss.copy(), float(big_n) ** (-ss - 1)
        for k, b in enumerate(BERNOULLI, start=1):
            total += b / math.factorial(2 * k) * rising * power
            rising = rising * (ss + 2 * k - 1) * (ss + 2 * k)
            power = power / big_n ** 2
        out[i:i + chunk] = total
    return out


def theta(t):
    """Riemann-Siegel theta."""
    t = np.asarray(t, dtype=float)
    asym = t / 2 * np.log(t / (2 * math.pi) + (t <= 50)) - t / 2 - math.pi / 8 \
        + 1 / (48 * t + (t == 0)) + 7 / (5760 * (t + (t == 0)) ** 3)
    exact = np.imag(loggamma(0.25 + 0.5j * t)) - 0.5 * t * math.log(math.pi)
    return np.where(t > 50, asym, exact)


# C0 of the Riemann-Siegel remainder as a series in z = 2p - 1 (Gabcke);
# C1 and C2 follow from its derivatives. check_rs() compares C0 with the
# closed form cos(2 pi (p^2 - p - 1/16)) / cos(2 pi p).
_C0 = [0.38268343236508977, 0.43724046807752045, 0.13237657548034352,
       -0.01360502604767419, -0.01356762197010358, -0.00162372532314446,
       0.00029705353733380, 0.00007943300879521, 0.00000046556124614,
       -0.00000143272516310, -0.00000010354847112, 0.00000001235792708,
       0.00000000178810839, -0.00000000003391414, -0.00000000001632663]
_P0 = np.polynomial.Polynomial(np.ravel(np.column_stack([_C0, np.zeros(len(_C0))])))
_P2, _P3, _P6 = _P0.deriv(2), _P0.deriv(3), _P0.deriv(6)


def hardy_z_rs(t, chunk=200000):
    t = np.asarray(t, dtype=float)
    out = np.empty_like(t)
    for i in range(0, len(t), chunk):
        tt = t[i:i + chunk]
        a = np.sqrt(tt / (2 * math.pi))
        big_n = np.floor(a).astype(int)
        z = 2 * (a - big_n) - 1
        th = theta(tt)
        acc = np.zeros_like(tt)
        for n in range(1, big_n.max() + 1):
            m = big_n >= n
            acc[m] += np.cos(th[m] - tt[m] * math.log(n)) / math.sqrt(n)
        w = np.sqrt(2 * math.pi / tt)
        c0 = _P0(z)
        c1 = -_P3(z) / (12 * math.pi ** 2)
        c2 = _P6(z) / (288 * math.pi ** 4) + _P2(z) / (16 * math.pi ** 2)
        out[i:i + chunk] = 2 * acc + (-1.0) ** (big_n - 1) * np.sqrt(w) * (c0 + c1 * w + c2 * w * w)
    return out


def hardy_z(t):
    t = np.asarray(t, dtype=float)
    out = np.empty_like(t)
    low = t < 200
    out[low] = np.real(np.exp(1j * theta(t[low])) * zeta(0.5 + 1j * t[low]))
    out[~low] = hardy_z_rs(t[~low])
    return out


def check_rs():
    p = np.array([0.03, 0.17, 0.41, 0.58, 0.66, 0.91])
    closed = np.cos(2 * math.pi * (p * p - p - 1 / 16)) / np.cos(2 * math.pi * p)
    t = np.array([210.3, 500.7, 1234.5, 3000.2])
    em = np.real(np.exp(1j * theta(t)) * zeta(0.5 + 1j * t))
    e0, e1 = np.abs(_P0(2 * p - 1) - closed).max(), np.abs(hardy_z_rs(t) - em).max()
    print(f"Riemann-Siegel: C0 series {e0:.1e}, Z(t) vs Euler-Maclaurin at t<=3000 {e1:.1e}")
    assert e0 < 1e-12 and e1 < 1e-5


# --------------------------------------------------------------------------
# helpers on the critical line
# --------------------------------------------------------------------------

def cell_log_line(t0, t1, h):
    """Midpoints u of cells of width h on [t0, t1] and U = log|f(1/2+iu)|,
    f = (s-1) zeta(s) / s^2, with exact cell averages of the logarithmic
    singularity near each zero; also the zeros g (sign changes of Hardy Z,
    linear interpolation)."""
    u = np.arange(t0 + h / 2, t1, h)
    zu = hardy_z(u)
    k = np.nonzero(np.signbit(zu[1:]) != np.signbit(zu[:-1]))[0]
    g = u[k] + h * zu[k] / (zu[k] - zu[k + 1])
    U = np.log(np.abs(zu)) - 0.5 * np.log(0.25 + u * u)

    def F(x):
        return x * np.log(np.abs(x) + 1e-300) - x

    for d in range(-3, 5):
        c = np.clip(k + d, 0, len(u) - 1)
        a = u[c] - h / 2 - g
        U[c] += (F(a + h) - F(a)) / h - np.log(np.abs(u[c] - g))
    return u, U, g


def arg_f_line(t, n_count):
    """arg f(1/2+it) = arg(s-1) - 2 arg s + pi (N(t)-1) - theta(t)  (Riemann-von Mangoldt)."""
    s = 0.5 + 1j * np.asarray(t)
    return np.angle(s - 1) - 2 * np.angle(s) + math.pi * (n_count - 1) - theta(t)


def conjugate_transform(t_targets, u, U, h, big_t, near=150.0, block=0.1):
    """V(t) = (2t/pi) int_0^infty [U(u)-U(t)]/(u^2-t^2) du with U ~ -log u beyond big_t."""
    nb = int(big_t / block)
    idx = np.minimum((u // block).astype(int), nb - 1)
    mid = (np.arange(nb) + 0.5) * block
    m0 = np.bincount(idx, weights=U, minlength=nb) * h
    m1 = np.bincount(idx, weights=U * (u - mid[idx]), minlength=nb) * h
    cnt = np.bincount(idx, minlength=nb) * h
    out = []
    for t in t_targets:
        ut = np.interp(t, u, U)
        sel = np.abs(u - t) <= near
        un = u[sel]
        integrand = (U[sel] - ut) / (un * un - t * t)
        integrand[np.argmin(np.abs(un - t))] = 0.0
        val = h * integrand.sum()
        far = np.abs(mid - t) > near + block
        mf = mid[far]
        val += np.sum((m0[far] - ut * cnt[far]) / (mf * mf - t * t)
                      - 2 * mf / (mf * mf - t * t) ** 2 * m1[far])
        # beyond big_t: U ~ -log u; substitute u = big_t / x, log x carried by the quadrature weight
        kern = lambda x: big_t / (big_t ** 2 - t * t * x * x)
        tail = quad(kern, 0, 1, weight="alg-loga", wvar=(0, 0))[0] \
            - (math.log(big_t) + ut) * math.atanh(t / big_t) / t
        out.append(2 * t / math.pi * (val + tail))
    return np.array(out)


def log_f_interior(s0):
    """log f(s0) continued horizontally from the right; real on the real axis."""
    sig = np.linspace(40.0, s0.real, 40001)
    s = sig + 1j * s0.imag
    v = (s - 1) * zeta(s) / s ** 2
    return math.log(abs(v[-1])) + 1j * np.unwrap(np.angle(v))[-1]


# --------------------------------------------------------------------------
# section 1: logarithmic coins
# --------------------------------------------------------------------------

def check_log_coin(n_max=200000):
    """(1.1)-(1.2): the logarithm of D(s-bar, s) splits into cosine and sine coins."""
    n = np.arange(2, n_max + 1)
    lam = np.zeros(n_max + 1)
    sieve = np.ones(n_max + 1, dtype=bool)
    for p in range(2, n_max + 1):
        if sieve[p]:
            sieve[p * p::p] = False
            q = p
            while q <= n_max:
                lam[q] = 1.0 / math.log(q) * math.log(p)
                q *= p
    c = lam[2:]
    worst = 0.0
    for s in (2 + 3j, 2.5 - 7j, 2 + 20j):
        sig, al = s.real, s.imag
        cos_coin = np.sum(c * n ** -sig * np.cos(al * np.log(n)))
        sin_coin = np.sum(c * n ** -sig * np.sin(al * np.log(n)))
        zs = zeta(s)[0]
        for z, w, exact in ((-0.5, -0.5, -math.log(abs(zs))), (-0.5, 0.5, 1j * np.angle(zs)),
                            (0.3, 0.9, 1.2 * math.log(abs(zs)) + 0.6j * np.angle(zs))):
            coin = (z + w) * cos_coin + 1j * (z - w) * sin_coin
            worst = max(worst, abs(coin - exact))
    print(f"(1.1)-(1.2) logarithmic coins at sigma>=2: max error {worst:.1e}")
    assert worst < 1e-4


# --------------------------------------------------------------------------
# section 2: criterion 1 (Littlewood, M-psi form)
# --------------------------------------------------------------------------

def littlewood_parts(sigma0, t_values, h=2e-3):
    """Cosine part int_0^T [C(sigma0+ia) - C(sigma1+ia)] da and sine part
    int_{sigma0}^{sigma1} S(sigma+iT) dsigma, C=log|Z|, S=arg Z, Z=(s-1)zeta(s)."""
    t_max = max(t_values)
    a = np.arange(h / 2, t_max, h)

    def log_abs_z(s):
        return np.log(np.abs((s - 1) * zeta(s)))

    if sigma0 == 0.5:      # zeros on the edge: exact cell averages of log|t-gamma|
        u, U, _ = cell_log_line(0.0, t_max, h)
        a = u
        left = U + np.log(0.25 + u * u)          # log|(s-1)zeta(s)| on the line
    else:
        left = log_abs_z(sigma0 + 1j * a)
    edge = left - log_abs_z(SIGMA_1 + 1j * a)
    cum = np.concatenate([[0.0], np.cumsum(edge) * h])
    grid = np.concatenate([[0.0], a + h / 2])
    sig = sigma0 + np.concatenate([[0.0], np.logspace(-7, math.log10(SIGMA_1 - sigma0), 3000)])
    cos_part, sin_part = [], []
    for t in t_values:
        v = (sig + 1j * t - 1) * zeta(sig + 1j * t)
        arg = np.unwrap(np.angle(v[::-1]))[::-1]
        sin_part.append(simpson(arg, x=sig))
        cos_part.append(np.interp(t, grid, cum))
    return np.array(cos_part), np.array(sin_part)


def check_littlewood(gammas):
    """(2.1)-(2.2) with sigma1 = 3."""
    worst = 0.0
    for sigma0 in (0.6, 0.5, 0.3, 0.1):
        ts = [20.0, 35.0, 50.0]
        c, s = littlewood_parts(sigma0, ts)
        for t, ci, si in zip(ts, c, s):
            pred = 2 * math.pi * max(0.5 - sigma0, 0.0) * np.sum(gammas <= t)
            worst = max(worst, abs(ci + si - pred))
            print(f"    sigma0={sigma0:.1f} T={t:2.0f}: cosine {ci:+9.5f}, sine {si:+9.5f}, "
                  f"sum {ci + si:9.5f}, 2 pi sum(beta-sigma0)+ {pred:9.5f}")
    print(f"(2.1)-(2.2) Littlewood in cosine/sine form: max error {worst:.1e}")
    assert worst < 5e-3


# --------------------------------------------------------------------------
# section 3: criterion 2 (conjugate transform)
# --------------------------------------------------------------------------

def check_schwarz():
    """(3.1): log f(s0) equals the Schwarz integral of the cosine data (no zero off the line)."""
    big_t, h = 400.0, 0.0025
    u, U, _ = cell_log_line(0.0, big_t, h)
    worst = 0.0
    for s0 in (1.0 + 0j, 2 + 3j, 1.2 + 10j, 0.7 + 20j, 0.55 + 30j):
        a, al = s0.real - 0.5, s0.imag
        kernel = 1 / (s0 - 0.5 - 1j * u) + 1 / (s0 - 0.5 + 1j * u)
        val = h * np.sum(U * kernel) / math.pi
        val += -(2 * a + 2j * al) * (math.log(big_t) + 1) / (math.pi * big_t)
        exact = 0j if s0 == 1 else log_f_interior(s0)
        worst = max(worst, abs(val - exact))
        print(f"    s0={s0}: Schwarz integral {val.real:+.4f}{val.imag:+.4f}i, "
              f"log f(s0) {exact.real:+.4f}{exact.imag:+.4f}i")
    print(f"(3.1) Schwarz representation without Blaschke factor: max error {worst:.1e}")
    assert worst < 3e-3


def critical_line_data(t_max=60.0, step=0.05, h=0.002, big_t=600.0):
    u, U, g = cell_log_line(0.0, big_t, h)
    t = np.arange(0.5, t_max, step)
    t = t[np.min(np.abs(t[:, None] - g[None, g < t_max + 1]), axis=1) > 0.03]
    n_t = np.searchsorted(g, t)
    v = conjugate_transform(t, u, U, h, big_t)
    sine = arg_f_line(t, n_t)
    return t, v, sine, np.interp(t, u, U), g[g < t_max]


def check_critical_line():
    """(3.3): on the critical line the sine side is the conjugate transform of the cosine side."""
    t, v, sine, _, _ = critical_line_data()
    err = np.abs(v - sine)
    print(f"(3.3) critical line t<=60: median |H[cos]-sin| {np.median(err):.1e}, max {err.max():.1e}")
    assert np.median(err) < 2e-3 and err.max() < 0.05


def check_defect_energy():
    """(3.6): int (pi sgn u - 2 arctan(u/delta))^2 du = 8 pi log 2 * delta."""
    for delta in (0.1, 1.0):
        val = 2 * quad(lambda x: (math.pi - 2 * math.atan(x / delta)) ** 2, 0, np.inf, limit=400)[0]
        print(f"(3.6) delta={delta}: energy {val:.10f}, 8 pi log2 delta {8 * math.pi * math.log(2) * delta:.10f}")
        assert abs(val - 8 * math.pi * math.log(2) * delta) < 1e-7


# --------------------------------------------------------------------------
# section 4: heights and windows (--high)
# --------------------------------------------------------------------------

def check_heights():
    """Section 4: (3.3) at heights 10^3, 10^4, 10^5."""
    for t0, h in ((1e3, 0.003), (1e4, 0.003), (1e5, 0.006)):
        big_t = 1.6 * t0 + 500
        u, U, g = cell_log_line(0.0, big_t, h)
        n0 = int(np.searchsorted(g, t0))
        t = np.arange(t0, t0 + 10.0, 0.02)
        t = t[np.min(np.abs(t[:, None] - g[None, (g > t0 - 1) & (g < t0 + 11)]), axis=1) > 0.05]
        v = conjugate_transform(t, u, U, h, big_t)
        err = np.abs(v - arg_f_line(t, np.searchsorted(g, t)))
        print(f"    t0=1e{round(math.log10(t0))}: N(t0)={n0}, median {np.median(err):.4f}, "
              f"95% {np.quantile(err, 0.95):.4f}, max {err.max():.4f}")
        assert np.quantile(err, 0.95) < 0.05


def check_windows(heights):
    """Section 4: equal variance and zero correlation of the two sides in windows."""
    for t0 in heights:
        span = 600.0 if t0 >= 1e8 else 1000.0
        t = np.arange(t0, t0 + span, 0.002)
        z = hardy_z_rs(t)
        x = np.log(np.abs(z))
        base = np.cumsum(np.concatenate([[0], np.signbit(z[1:]) != np.signbit(z[:-1])])) \
            - theta(t) / math.pi - 1
        y = math.pi * (base + np.round(-base.mean()))   # Turing-style offset: mean S ~ 0
        ratio, corr = x.var() / y.var(), np.corrcoef(x, y)[0, 1]
        print(f"    T=1e{round(math.log10(t0))}: Var(cos) {x.var():.3f}, Var(sin) {y.var():.3f}, "
              f"ratio {ratio:.4f}, corr {corr:+.4f}, (1/2)loglogT {0.5 * math.log(math.log(t0)):.3f}")
        assert abs(ratio - 1) < 5e-3 and abs(corr) < 1e-2


def check_second_order():
    """(4.1): the contour integral of (log zeta)^2 over [1/2,3]x[T1,T2] vanishes."""
    def row(t):
        x = np.concatenate([[0], np.logspace(-8, math.log10(SIGMA_1 - 0.5), 8000)])
        sig = 0.5 + x
        v = zeta(sig + 1j * t)
        return sig, np.log(np.abs(v)) + 1j * np.unwrap(np.angle(v[::-1]))[::-1]

    def left(t1, t2, n1, h):
        e = np.arange(t1, t2 + h / 2, h)
        t = 0.5 * (e[1:] + e[:-1])
        ze, zt = hardy_z(e), hardy_z(t)
        flip = np.signbit(ze[1:]) != np.signbit(ze[:-1])
        k = np.nonzero(flip)[0]
        g = e[k] + h * ze[k] / (ze[k] - ze[k + 1])
        gp = np.log(np.abs((ze[k + 1] - ze[k]) / h))
        cnt = np.cumsum(flip)
        vals = (np.log(np.abs(zt)) + 1j * (math.pi * (n1 + cnt - 1) - theta(t))) ** 2 * h

        def i0(length, c):
            out = np.zeros_like(c)
            m = length > 0
            lg = np.log(length[m]) + c[m]
            out[m] = length[m] * (lg * lg - 2 * lg + 2)
            return out

        for d in range(-3, 4):
            c = np.clip(k + d, 0, len(t) - 1)
            a, b = e[c] - g, e[c + 1] - g
            vr = math.pi * (n1 + cnt[c] - 1) - theta(t[c])
            cr, cl = gp + 1j * vr, gp + 1j * (vr - math.pi * (d == 0))
            vals[c] = np.where(b <= 0, i0(-a, cl) - i0(-b, cl),
                               np.where(a >= 0, i0(b, cr) - i0(a, cr), i0(-a, cl) + i0(b, cr)))
        return vals.sum()

    for t1, t2 in ((1000.0, 1500.0), (5000.0, 5600.0)):
        sb, lb = row(t1)
        st, lt = row(t2)
        n1 = int(round(theta(np.array([t1]))[0] / math.pi + 1 + lb[0].imag / math.pi))
        gl = 2 * left(t1, t2, n1, 0.0005) - left(t1, t2, n1, 0.001)     # error is O(h)
        tr = np.linspace(t1, t2, 400001)
        gr = simpson(np.log(zeta(SIGMA_1 + 1j * tr)) ** 2, x=tr)
        total = simpson(lb ** 2, x=sb) + 1j * gr - simpson(lt ** 2, x=st) - 1j * gl
        print(f"    [{t1:.0f},{t2:.0f}]: int(cos^2-sin^2) {gl.real:+.3f}, 2 int cos*sin {gl.imag:+.3f}, "
              f"contour {abs(total):.1e}")
        assert abs(total) < 0.05


# --------------------------------------------------------------------------
# figures (--figures)
# --------------------------------------------------------------------------

INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#ffffff"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"


def _style(plt):
    plt.rcParams.update({
        "font.family": ["DejaVu Sans", "Noto Sans CJK KR"], "font.size": 10,
        "axes.edgecolor": GRID, "axes.labelcolor": INK2, "axes.linewidth": 0.8,
        "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True, "grid.color": GRID,
        "grid.linewidth": 0.6, "axes.facecolor": SURFACE, "figure.facecolor": SURFACE,
        "legend.frameon": False, "lines.linewidth": 1.6, "svg.hashsalt": "sine-cosine",
        "axes.spines.top": False, "axes.spines.right": False,
    })


def _save(fig, name):
    path = os.path.join(FIG_DIR, name)
    fig.savefig(path, format="svg", bbox_inches="tight", metadata={"Date": None})
    print(f"    wrote {os.path.relpath(path)}")


def make_figures(gammas):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _style(plt)

    # Figure 1: criterion 1. Labels are plain Unicode (no mathtext) so that the
    # Hangul fallback font is used for the Korean words.
    ts = np.arange(0.25, 60.0, 0.25)
    c5, s5 = littlewood_parts(0.5, ts)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.4, 5.6), sharex=True)
    ax1.plot(ts, c5, color=BLUE, label="코사인 쪽 (M):  ∫₀ᵀ [C(½+iα) − C(3+iα)] dα")
    ax1.plot(ts, -s5, color=ORANGE, ls=(0, (4, 3)), label="사인 쪽 (ψ)의 부호 반전:  −∫ S(σ+iT) dσ")
    ax1.set_title("σ₀ = ½ : 두 쪽이 모든 높이에서 겹친다", loc="left", color=INK, fontsize=10)
    ax1.set_ylim(-5.6, 0.9)
    ax1.legend(loc="upper right", fontsize=8.5)
    for sigma0, color, ls in ((0.5, BLUE, "-"), (0.3, ORANGE, (0, (4, 3))), (0.1, AQUA, (0, (1, 2)))):
        c, s = (c5, s5) if sigma0 == 0.5 else littlewood_parts(sigma0, ts)
        ax2.plot(ts, c + s, color=color, ls=ls, label=f"σ₀ = {sigma0}")
        n = np.array([np.sum(gammas <= t) for t in ts])
        ax2.plot(ts, 2 * math.pi * (0.5 - sigma0) * n, color=INK2, lw=0.8, alpha=0.6)
    ax2.set_title("코사인 쪽 + 사인 쪽 = 2π Σ (β − σ₀)⁺   (회색: 영점 개수로 계산한 값)", loc="left",
                  color=INK, fontsize=10)
    ax2.set_xlabel("T")
    ax2.legend(loc="upper left", fontsize=8.5)
    for ax in (ax1, ax2):
        for g in gammas[gammas < 60]:
            ax.axvline(g, color=GRID, lw=0.6, zorder=0)
    fig.tight_layout()
    _save(fig, "sine-cosine-littlewood.svg")
    plt.close(fig)

    # Figure 2: criterion 2 on the critical line.
    t, v, sine, U, zeros = critical_line_data()
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.4, 5.2), sharex=True,
                                   gridspec_kw={"height_ratios": [1, 1.25]})
    ax1.plot(t, np.maximum(U, -6), color=BLUE)
    ax1.set_ylim(-6.2, 1)
    ax1.set_title("코사인 쪽  log|f(½+it)|   (영점에서 −∞, −6에서 자름)", loc="left", color=INK, fontsize=10)
    ax2.plot(t, sine, color=ORANGE, label="사인 쪽  arg f(½+it)")
    ax2.plot(t[::6], v[::6], ls="none", marker="o", ms=3.2, color=BLUE, mec=SURFACE, mew=0.8,
             label="코사인 쪽의 켤레 변환 (3.3)")
    ax2.set_title("사인 쪽: 영점마다 π 계단, 그 사이 기울기 −θ′(t)", loc="left", color=INK, fontsize=10)
    ax2.set_xlabel("t")
    ax2.legend(loc="lower left", fontsize=8.5)
    for ax in (ax1, ax2):
        for g in zeros:
            ax.axvline(g, color=GRID, lw=0.6, zorder=0)
    fig.tight_layout()
    _save(fig, "sine-cosine-critical-line.svg")
    plt.close(fig)

    # Figure 3: defect of an off-line pair.
    x = np.linspace(-3, 3, 2401)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.4, 4.8), sharex=True)
    ax1.plot(x, math.pi * np.sign(x), color=INK, lw=1.2, label="사인 쪽 계단  π·sgn(u)")
    for delta, color, ls in ((0.1, BLUE, "-"), (0.3, ORANGE, (0, (4, 3))), (1.0, AQUA, (0, (1, 2)))):
        ax1.plot(x, 2 * np.arctan(x / delta), color=color, ls=ls,
                 label=f"코사인 쪽 켤레 변환  2·arctan(u/δ), δ = {delta}")
        ax2.plot(x, math.pi * np.sign(x) - 2 * np.arctan(x / delta), color=color, ls=ls,
                 label=f"δ = {delta}")
    ax1.set_title("선에서 δ만큼 떨어진 영점과 그 거울상   (u = t − γ)", loc="left", color=INK, fontsize=10)
    ax1.legend(loc="upper left", fontsize=8)
    ax2.set_title("결손  π·sgn(u) − 2·arctan(u/δ),   에너지 8π·log 2·δ", loc="left", color=INK, fontsize=10)
    ax2.set_xlabel("u")
    ax2.legend(loc="upper right", fontsize=8.5)
    fig.tight_layout()
    _save(fig, "sine-cosine-defect.svg")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--high", action="store_true", help="heights 10^3-10^5 and windows up to 10^8")
    parser.add_argument("--figures", action="store_true", help="write docs/figures/sine-cosine-*.svg")
    args = parser.parse_args()
    gammas = zeros_up_to(60.0)
    check_rs()
    check_log_coin()
    check_littlewood(gammas)
    check_schwarz()
    check_critical_line()
    check_defect_energy()
    print("Section 4 windows:")
    check_windows([1e3, 1e4, 1e5] + ([1e6, 1e7, 1e8] if args.high else []))
    if args.high:
        print("Section 4 heights:")
        check_heights()
        print("(4.1) second-order contour:")
        check_second_order()
    if args.figures:
        print("Figures:")
        make_figures(gammas)


if __name__ == "__main__":
    main()
