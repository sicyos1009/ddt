"""Checks for docs/divisor-cumsum.md. Requires numpy and scipy.

Run: python -B scripts/check_divisor_cumsum.py
     python -B scripts/check_divisor_cumsum.py --figure docs/figures/divisor-cumsum.svg

The per-divisor shares f_x(d) are sums of numbers +-2^-k, so for x <= 10^6
every partial sum is exact in double precision and the identities of section 1
are checked with ==. The comparisons of section 2 are consistency checks, not a
proof; the cumulative prediction (2.4) is not covered by a theorem.
"""

import argparse
import math

import numpy as np
from scipy.special import betainc, gamma, hyp2f1, poch, rgamma

EULER_GAMMA = 0.5772156649015329
G1 = -EULER_GAMMA / 2          # main-theorem (6.1) at z = w = -1/2
G2 = -EULER_GAMMA / 2 - 1


def mobius_weights(x):
    """mu(n) and w(n) = mu(n)/tau(n) = mu(n) 2^-omega(n) for n <= x."""
    mu = np.ones(x + 1, dtype=np.int64)
    omega = np.zeros(x + 1, dtype=np.int64)
    is_prime = np.ones(x + 1, dtype=bool)
    is_prime[:2] = False
    for p in range(2, x + 1):
        if is_prime[p]:
            is_prime[2 * p::p] = False
            mu[p::p] *= -1
            omega[p::p] += 1
            mu[p * p::p * p] = 0
    mu[0] = 0
    w = np.zeros(x + 1)
    sf = mu != 0
    w[sf] = mu[sf] * 2.0 ** (-omega[sf])
    return mu, w


def shares(w, x):
    """(1.1): f_x(d) = sum of w(n) over multiples n <= x of d."""
    f = np.zeros(x + 1)
    for d in range(1, x + 1):
        f[d] = w[d::d].sum()
    return f


def kernel_cdf(z, w, v):
    """(2.3): <1_[0,v] K_{z,w}>, the continued incomplete beta integral."""
    return v ** z * hyp2f1(z, 1 - w, z + 1, v) / (gamma(z + 1) * gamma(w))


def prediction(v, big_l, order=1):
    """(2.4): M_x(v) / (x/L^2) to first order in 1/L, closed form."""
    s = v * (1 - v)
    lead = (v - 0.5) / (math.pi * np.sqrt(s))
    if order == 0:
        return lead
    poly = EULER_GAMMA / 2 + (EULER_GAMMA + 1) * (3 * v - 12 * v ** 2 + 8 * v ** 3)
    return lead - poly / (4 * math.pi * big_l * s ** 1.5)


def check_exact(cases):
    """(1.2)-(1.4): total, first moment and factorization."""
    for x, (mu, w, f) in cases.items():
        m_x = int(mu[1:].sum())
        assert f[1:].sum() == m_x
        n = np.arange(1, x + 1)
        lhs = (f[1:] * np.log(n)).sum()
        rhs = 0.5 * (mu[1:] * np.log(n)).sum()
        assert abs(lhs - rhs) < 1e-6 * max(1.0, abs(rhs))
        print(f"(1.2) x = {x:>7}: sum f_x(d) = {f[1:].sum():.0f} = M(x); "
              f"(1.3) error {abs(lhs - rhs):.1e}")
    x = 10 ** 4
    mu, w, f = cases[x]
    worst = 0.0
    for d in range(1, x + 1):
        if mu[d] == 0:
            worst = max(worst, abs(f[d]))
            continue
        m = np.arange(1, x // d + 1)
        inner = w[m][np.gcd(m, d) == 1].sum()
        worst = max(worst, abs(f[d] - w[d] * inner))
    print(f"(1.4) factorization, x = {x}: max error {worst:.1e}")
    assert worst == 0.0


def check_first_share(cases):
    """(1.5): f_x(1) against main-theorem (2.2) with (z, w) = (-1/2, 0), F = 1."""
    mu, w, f = cases[max(cases)]
    p = np.nonzero(w == -0.5)[0].astype(float)          # omega(n) = 1, squarefree
    r = -0.5
    c_r = math.exp(np.sum(np.log1p(r / p) + r * np.log1p(-1 / p)))
    kappa = EULER_GAMMA + np.sum((r + 1) * np.log(p) / ((p - 1) * (p + r)))
    c0, c1 = c_r / gamma(r), c_r * (r * kappa - 1) / gamma(r - 1)
    print(f"(1.5) C(-1/2)/Gamma(-1/2) = {c0:.4f}; f_x(1) L^(3/2)/x vs c0 + c1/L:")
    for x, (mu, w, f) in cases.items():
        big_l = math.log(x)
        obs = f[1] * big_l ** 1.5 / x
        print(f"  x = {x:>7}: {obs:+.4f} vs {c0 + c1 / big_l:+.4f}")


def check_kernels():
    """(2.3): continuation of the incomplete beta and the closed forms."""
    v = np.linspace(0.02, 0.98, 49)
    errors = []
    for z, w in ((0.5, 0.5), (2.0, 1.0), (1.5, 0.7)):
        errors.append(np.abs(kernel_cdf(z, w, v) - betainc(z, w, v) * rgamma(z + w)).max())
    errors.append(np.abs(kernel_cdf(-0.5, -0.5, v) - prediction(v, 1.0, order=0)).max())
    tail = -(1 - 2 * v) / (2 * math.pi * np.sqrt(v * (1 - v)))      # -T of main (7.5)
    errors.append(np.abs(kernel_cdf(-0.5, -0.5, v) - tail).max())
    errors.append(np.abs(kernel_cdf(-0.5, 0.5, v)
                         - (1 / math.pi) * np.sqrt((1 - v) / v)).max())
    for big_l in (5.0, 13.8):
        first = (kernel_cdf(-0.5, -0.5, v)
                 + (G1 * kernel_cdf(-1.5, -0.5, v) + G2 * kernel_cdf(-0.5, -1.5, v)) / big_l)
        errors.append(np.abs(first - prediction(v, big_l)).max())
    print("(2.3)-(2.4) kernel CDFs and closed forms: max error "
          f"{max(errors):.1e}")
    assert max(errors) < 1e-10


def moment_prediction(k, big_l):
    """(2.2): <K + L^-1 (G1 K_{-3/2,-1/2} + G2 K_{-1/2,-3/2}), u^k> via main (4.2)."""
    m = lambda z, w: poch(z, k) * rgamma(z + w + k)
    return m(-0.5, -0.5) + (G1 * m(-1.5, -0.5) + G2 * m(-0.5, -1.5)) / big_l


def check_moments(cases):
    """(2.2): smooth test functions u^2, u^3 against the first-order law."""
    print("(2.2) sum f_x(d) (log d / L)^k / (x/L^2):")
    for k in (2, 3):
        for x, (mu, w, f) in cases.items():
            big_l = math.log(x)
            u = np.log(np.arange(1, x + 1)) / big_l
            obs = (f[1:] * u ** k).sum() / (x / big_l ** 2)
            pred = moment_prediction(k, big_l)
            print(f"  k = {k}, x = {x:>7}: observed {obs:+.4f}, predicted {pred:+.4f}, "
                  f"(obs - pred) L^2 = {(obs - pred) * big_l ** 2:+.2f}")


def check_cumulative(cases):
    """(2.4): cumulative shares against the closed-form prediction."""
    print("(2.4) M_x(v) / (x/L^2) at d <= x^v:")
    for x, (mu, w, f) in cases.items():
        big_l = math.log(x)
        c = np.cumsum(f[1:])
        scale = x / big_l ** 2
        row = [f"f_x(1) {f[1] / scale:+.3f}"]
        for v in (0.2, 0.5, 0.8):
            obs = c[int(x ** v + 1e-9) - 1] / scale
            row.append(f"v={v}: {obs:+.3f} ({prediction(v, big_l, 0):+.3f}, "
                       f"{prediction(v, big_l):+.3f})")
        row.append(f"M(x) {int(mu[1:].sum())}")
        print(f"  x = {x:>7}: " + "; ".join(row))


def write_figure(path, cases):
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
    x = max(cases)
    mu, w, f = cases[x]
    big_l, scale = math.log(x), x / math.log(x) ** 2
    c = np.cumsum(f[1:]) / scale
    idx = np.unique(np.concatenate([np.arange(1, 3001),
                                    np.geomspace(3000, x, 3000).astype(int)])) - 1
    v = np.log(idx + 1) / big_l
    grid = np.linspace(0.12, 0.88, 300)

    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    ax.axhline(0, color="#c3c2b7", lw=1)
    ax.grid(axis="y", color="#e1e0d9", lw=1)
    ax.set_axisbelow(True)
    ax.step(v, c[idx], where="post", color="#2a78d6", lw=2, label=f"$\\mathcal{{M}}_x(v)$, $x=10^{{{round(math.log10(x))}}}$")
    ax.plot(grid, prediction(grid, big_l), color="#52514e", lw=1.5, ls="--", label="예측 (2.4), $O(L^{-1})$까지")
    ax.plot(grid, prediction(grid, big_l, 0), color="#898781", lw=1.2, ls=":", label="주항 $P_{-1/2,-1/2}(v)$")
    for vv, yy, text, ha, dy in ((0, c[0], "$f_x(1)$", "left", -11), (1, c[-1], "$M(x)$", "right", 0)):
        ax.plot([vv], [yy], "o", ms=7, color="#2a78d6", mec="white", mew=2, zorder=5)
        ax.annotate(text, (vv, yy), xytext=(6 if ha == "left" else -8, dy),
                    textcoords="offset points", ha=ha, va="center", color="#52514e")
    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(c.min() - 0.3, max(c.max(), prediction(grid, big_l).max()) + 0.15)
    ax.set_xlabel("$v=\\log d/\\log x$")
    ax.set_ylabel("$\\mathcal{M}_x(v)\\,/\\,(x/L^2)$")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.legend(frameon=False, loc="upper left", fontsize=9)
    fig.tight_layout()
    fig.savefig(path, facecolor="white")
    print(f"figure written to {path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--figure", help="write the cumulative-share figure (SVG)")
    args = parser.parse_args()
    cases = {}
    for x in (10 ** 4, 10 ** 5, 10 ** 6):
        mu, w = mobius_weights(x)
        cases[x] = (mu, w, shares(w, x))
    check_exact(cases)
    check_first_share(cases)
    check_kernels()
    check_moments(cases)
    check_cumulative(cases)
    if args.figure:
        write_figure(args.figure, cases)


if __name__ == "__main__":
    main()
