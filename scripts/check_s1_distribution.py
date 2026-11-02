"""Checks for section 3 of docs/ng-duality.md. Requires numpy and scipy.

Run: python -B scripts/check_s1_distribution.py
Zeta values use the Euler-Maclaurin formula in double precision; zeros below
height 300 are located with Hardy's Z-function. These are consistency checks
at low height, not a proof.
"""

import math

import numpy as np
from scipy.special import j0, loggamma

from check_dual_product import primes_up_to, zeta_line
from check_ng_tail import zeros_up_to

# Gauss-Legendre nodes for the sigma-integral of log|zeta(sigma+it)|.
_GL_X, _GL_W = np.polynomial.legendre.leggauss(20)


def _gauss(a, b):
    return 0.5 * (b - a) * _GL_X + 0.5 * (a + b), 0.5 * (b - a) * _GL_W


def prime_powers(limit):
    out = []
    for p in primes_up_to(limit):
        q, k = int(p), 1
        while q <= limit:
            out.append((q, math.log(p), k))
            q *= int(p)
            k += 1
    return out


def integrated_log_zeta(ts):
    """I(t) = int_{1/2}^infty log|zeta(sigma+it)| d sigma for an array of t >= 10."""
    ts = np.asarray(ts, dtype=float)
    total = np.zeros_like(ts)
    # sigma = 1/2 + w^2 on [1/2, 3/4]: removes the logarithmic endpoint singularity.
    w, wt = _gauss(0.0, 0.5)
    for wi, wti in zip(w, wt):
        total += 2 * wi * wti * np.log(np.abs(zeta_line(0.5 + wi * wi, ts)))
    for a, b in ((0.75, 1.25), (1.25, 2.0)):
        nodes, weights = _gauss(a, b)
        for s, wgt in zip(nodes, weights):
            total += wgt * np.log(np.abs(zeta_line(s, ts)))
    # sigma >= 2: int_2^infty Re log zeta = Re sum Lambda(n) n^{-2-it} / log^2 n.
    for n, logp, k in prime_powers(20000):
        log_n = k * logp
        total += logp / log_n ** 2 * n ** -2.0 * np.cos(ts * log_n)
    return total


def integrated_log_zeta_real_axis():
    """I(0) = int_{1/2}^infty log|zeta(sigma)| d sigma (log singularity at sigma=1)."""
    from scipy.integrate import quad
    from check_ng_tail import zeta

    def integrand(s):
        return math.log(abs(zeta(complex(s, 0.0)).real))

    pieces = [(0.5, 0.9), (0.9, 1.0), (1.0, 1.1), (1.1, 2.0), (2.0, 40.0)]
    return sum(quad(integrand, a, b, limit=200)[0] for a, b in pieces)


def riemann_siegel_theta(t):
    return loggamma(0.25 + 0.5j * t).imag - 0.5 * t * math.log(math.pi)


def check_littlewood_identity():
    """pi*S_1(t) - I(t) is constant (= -I(0)) where RH has been verified."""
    gammas = zeros_up_to(300.0)
    probe = np.array([20.0, 50.0, 90.5, 140.0, 199.9, 250.0, 299.0])
    # S(u) = N(u) - 1 - theta(u)/pi, so int_0^t N = sum_{gamma<=t} (t - gamma) exactly
    # and only the smooth theta-integral needs quadrature.
    u = np.linspace(0.0, 300.0, 600001)
    theta_integral = np.concatenate([[0.0], np.cumsum(
        0.5 * (riemann_siegel_theta(u[1:]) + riemann_siegel_theta(u[:-1])) * np.diff(u))])
    s1_probe = np.array([np.sum(np.clip(t - gammas, 0, None)) for t in probe]) \
        - probe - np.interp(probe, u, theta_integral) / math.pi
    difference = math.pi * s1_probe - integrated_log_zeta(probe)
    print("(3.2) pi*int_0^t S - I(t) at heights", probe.tolist())
    print("   ", np.round(difference, 5).tolist())
    spread = difference.max() - difference.min()
    i0 = integrated_log_zeta_real_axis()
    print(f"    spread {spread:.1e}; common value {difference.mean():.5f}; "
          f"-I(0) computed directly {-i0:.5f}")
    assert spread < 2e-3 and abs(difference.mean() + i0) < 2e-3


def prime_factor_data(primes, theta_points=4096, k_max=60):
    """Per-prime Y_p(theta) = sum_k cos(k theta) / (k^2 p^{k/2} log p) on a theta grid."""
    theta = 2 * math.pi * np.arange(theta_points) / theta_points
    ks = np.arange(1, k_max + 1)
    cos_k = np.cos(np.outer(ks, theta))
    data = []
    for p in primes:
        coeff = p ** (-ks / 2.0) / (ks ** 2 * math.log(p))
        data.append(coeff @ cos_k)
    return data


def model_statistics(etas):
    """Characteristic function, variance and third cumulant of sum_p Y_p."""
    small = primes_up_to(10000).astype(float)
    large = primes_up_to(2_000_000)[len(small):].astype(float)
    ys = prime_factor_data(small)
    char = np.ones(len(etas), dtype=complex)
    var, k3 = 0.0, 0.0
    for y in ys:
        char *= np.array([np.mean(np.exp(1j * eta * y)) for eta in etas])
        var += np.mean(y ** 2)
        k3 += np.mean(y ** 3)
    a1 = 1 / (np.sqrt(large) * np.log(large))
    for i, eta in enumerate(etas):
        char[i] *= np.prod(j0(eta * a1))
        # primes beyond 2e6: Gaussian tail with variance sum 1/(2 p log^2 p).
        char[i] *= math.exp(-eta ** 2 / (8 * math.log(2e6) ** 2))
    var += np.sum(a1 ** 2) / 2 + 1 / (4 * math.log(2e6) ** 2)
    return char, var, k3


def check_distribution():
    etas = [0.5, 1.0, 2.0, 3.0]
    ts = np.linspace(10.0, 3010.0, 60001)
    values = integrated_log_zeta(ts)
    centered = values - values.mean()
    emp_var = np.mean(centered ** 2)
    emp_k3 = np.mean(centered ** 3)
    char, var, k3 = model_statistics(etas)
    print(f"\nTheorem 2: I(t) for t in [10, 3010] against the prime model sum_p Y_p")
    print(f"  mean {values.mean():+.4f} (model 0), variance {emp_var:.4f} (model {var:.4f}), "
          f"third cumulant {emp_k3:+.4f} (model {k3:+.4f})")
    print("   eta   empirical            model")
    for eta, model in zip(etas, char):
        empirical = np.mean(np.exp(1j * eta * values))
        print(f"  {eta:4.1f}  {empirical.real:+.4f}{empirical.imag:+.4f}i   "
              f"{model.real:+.4f}{model.imag:+.4f}i")
        assert abs(empirical - model) < 0.02
    assert abs(emp_var - var) < 0.05 and emp_k3 > 0 and k3 > 0


def check_l1_decay():
    """(3.6): (1/T) int |I - I_X| decreases with X; compare with 1/log X."""
    ts = np.linspace(1500.0, 3000.0, 30001)
    values = integrated_log_zeta(ts)
    primes = primes_up_to(20000)
    print("\n(3.6) mean |I - I_X| on t in [1500, 3000]:")
    previous = math.inf
    for x_cut in (2, 13, 101, 1009, 10007):
        model = np.zeros_like(ts)
        for p in primes[primes <= x_cut]:
            logp = math.log(p)
            for k in range(1, 60):
                amplitude = p ** (-k / 2) / (k * k * logp)
                if amplitude < 1e-12:
                    break
                model += amplitude * np.cos(k * ts * logp)
        error = np.mean(np.abs(values - model))
        print(f"  X={x_cut:5d}: {error:.4f}   (1/log X = {1 / math.log(x_cut):.4f})")
        assert error < previous and error < 1 / math.log(x_cut)
        previous = error


def main():
    check_littlewood_identity()
    check_l1_decay()
    check_distribution()


if __name__ == "__main__":
    main()
