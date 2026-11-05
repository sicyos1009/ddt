"""M과 ψ의 영점 성분: 선분 표현 (2.5)–(2.6)과 1/L 전개 (2.7)의 검산.

numpy와 scipy만 사용한다. ζ는 Euler–Maclaurin 합으로 계산한다.

    python -B scripts/check_mertens_segment.py

확인 항목
  1. (2.5): 첫 영점 ρ₁에서 Euler 가지를 따라 계산한 윤곽 적분 I^M_ρ, I^ψ_ρ와
     선분 위의 Gauss–Chebyshev 적분이 α = ξ/L ≤ 5 까지 일치한다.
  2. (2.7): I^M_ρ/c_ρ에서 B_{1/2,1/2}와 1/L 보정을 뺀 나머지가 O(L^-2)이다.
     또한 Re κ_ρ = -ϑ_RS'(γ) 를 확인한다.
"""
import math

import numpy as np
from scipy.optimize import brentq
from scipy.special import bernoulli, digamma, j0, j1, loggamma

# ---------------------------------------------------------------- zeta
_B = bernoulli(2 * 14)
_FACT = np.array([float(math.factorial(k)) for k in range(2 * 14 + 1)])


def zeta(s, N=40, K=12):
    """Euler–Maclaurin, |Im s| <= 40 정도에서 기계 정밀도."""
    s = np.asarray(s, dtype=complex)
    flat = s.ravel()
    n = np.arange(1, N, dtype=float)
    out = np.exp(-np.outer(np.log(n), flat)).sum(axis=0)
    out = out + N ** (1 - flat) / (flat - 1) + 0.5 * N ** (-flat)
    poch = flat.copy()
    for k in range(1, K + 1):
        out = out + _B[2 * k] / _FACT[2 * k] * poch * N ** (-flat - 2 * k + 1)
        poch = poch * (flat + 2 * k - 1) * (flat + 2 * k)
    return out.reshape(s.shape)


def theta(t):
    t = np.asarray(t, dtype=float)
    return np.imag(loggamma(0.25 + 0.5j * t)) - 0.5 * t * np.log(np.pi)


def hardy_z(t):
    return np.real(np.exp(1j * theta(t)) * zeta(0.5 + 1j * np.asarray(t, dtype=float)))


assert abs(zeta(2.0) - np.pi ** 2 / 6) < 1e-13
assert abs(zeta(0.5) + 1.4603545088095868) < 1e-12

GAMMA = brentq(hardy_z, 14.0, 14.3, xtol=1e-15)
RHO = 0.5 + 1j * GAMMA
GAP = 21.022039638771555 - GAMMA  # 다음 영점까지의 거리; |α| < GAP 이어야 함


def taylor_at(center, kmax=6, r=0.05, M=96):
    th = 2 * np.pi * np.arange(M) / M
    vals = zeta(center + r * np.exp(1j * th))
    return np.array([np.mean(vals * np.exp(-1j * k * th)) / r ** k for k in range(kmax + 1)])


TC = taylor_at(RHO)          # TC[k] = ζ^{(k)}(ρ)/k!
ZP = TC[1]                   # ζ'(ρ)
KAPPA = TC[2] / TC[1]        # ζ''(ρ)/(2ζ'(ρ))


# ------------------------------------------- f(s) = ζ(s)/(s-ρ) on the critical line
def f_line(t):
    t = np.asarray(t, dtype=float)
    h = 1j * (t - GAMMA)
    out = np.empty(t.shape, dtype=complex)
    near = np.abs(h) < 0.02
    out[~near] = zeta(0.5 + 1j * t[~near]) / h[~near]
    out[near] = sum(TC[k] * h[near] ** (k - 1) for k in range(1, len(TC)))
    return out


# f(1/2+it) = -i e^{-iθ(t)} Z(t)/(t-γ) 이고 Z(t)/(t-γ)는 (0, γ₂)에서 부호가 일정하다.
# 따라서 arg f + θ 는 상수이고, 아래 ell은 log f 의 연속 가지이며 exp(ell(γ)) = ζ'(ρ).
PHI0 = np.angle(ZP * np.exp(1j * theta(GAMMA)))


def ell(t):
    return np.log(np.abs(f_line(t))) + 1j * (PHI0 - theta(t))


assert abs(np.exp(ell(np.array([GAMMA])))[0] - ZP) < 1e-12


# ------------------------------------------------------------------ (2.6)
def R_M(alpha, u):
    s_u = RHO + 1j * alpha * u
    return RHO / s_u * ZP * np.exp(-0.5 * (ell(GAMMA + alpha * u) + ell(GAMMA - alpha * (1 - u))))


def R_psi(alpha, u):
    s_u = RHO + 1j * alpha * u
    return RHO / s_u * np.exp(0.5 * (ell(GAMMA + alpha * u) - ell(GAMMA - alpha * (1 - u))))


# Gauss–Chebyshev: ∫_0^1 a(u) φ(u) du = mean φ(u_j)
NQ = 600
UQ = 0.5 * (1 - np.cos((2 * np.arange(1, NQ + 1) - 1) * np.pi / (2 * NQ)))


def segment_M(xi, L):
    """I^M_ρ / c_ρ,  c_ρ = x^ρ/(ρ ζ'(ρ))."""
    a = xi / L
    return np.mean(np.exp(1j * xi * UQ) * R_M(a, UQ))


def segment_psi(xi, L):
    """I^ψ_ρ / (x^ρ/ρ)."""
    a = xi / L
    return 1j * xi / L * np.mean(UQ * np.exp(1j * xi * UQ) * R_psi(a, UQ))


# ------------------------------------------- direct contour, Euler branch
def euler_log_zeta(s_end, n=4000):
    s = np.linspace(3.0, s_end.real, n) + 1j * s_end.imag
    z = zeta(s)
    return np.log(np.abs(z[-1])) + 1j * np.unwrap(np.angle(z))[-1]


def direct(xi, L, eA, eB, n=24000, pad=0.3):
    """(1/2πi)∮ ζ(s-iα)^eA ζ(s)^eB x^{s-ρ}/s ds 를 선분 [ρ, ρ+iα] 둘레 원에서 계산.
    가지는 원의 오른쪽 끝에서 Euler 가지로 정하고 원을 따라 연속으로 잇는다."""
    a = xi / L
    r = abs(a) / 2 + pad
    th = 2 * np.pi * np.arange(n) / n
    s = RHO + 0.5j * a + r * np.exp(1j * th)
    logs = []
    for shift in (a, 0.0):
        z = zeta(s - 1j * shift)
        ph = np.unwrap(np.angle(z))
        ph += euler_log_zeta(s[0] - 1j * shift).imag - ph[0]
        logs.append(np.log(np.abs(z)) + 1j * ph)
    g = np.exp(eA * logs[0] + eB * logs[1] + L * (s - RHO)) / s
    ds = 1j * (s - RHO - 0.5j * a)
    return np.mean(g * ds) / 1j  # = (1/2πi) ∮ g ds  (trapezoid, dθ = 2π/n)


# ------------------------------------------------------------- kernels
def B_half_half(xi):
    return np.exp(0.5j * xi) * j0(xi / 2)


def B_half_mhalf(xi):   # B_{1/2,-1/2}
    return 0.5j * xi * np.exp(0.5j * xi) * (j0(xi / 2) + 1j * j1(xi / 2))


def B_mhalf_mhalf(xi):  # B_{-1/2,-1/2}
    return 0.5 * xi * np.exp(0.5j * xi) * j1(xi / 2)


def main():
    print(f"γ₁ = {GAMMA:.15f}")
    print(f"ζ'(ρ₁) = {ZP:.12f}")
    print(f"κ_ρ₁ = ζ''/(2ζ') = {KAPPA:.12f}")
    print()

    print("1. (2.5) 선분 표현 vs Euler 가지 윤곽 적분")
    print(f"{'L':>5} {'ξ':>6} {'α':>6}   {'|I^M 상대차|':>12} {'|I^ψ 상대차|':>12}   {'I^M/c_ρ':>26}")
    worst = 0.0
    for L, xi in [(25, 1.0), (25, 4.0), (25, 10.0), (10, 15.0), (5, 12.5), (5, 25.0), (5, -25.0)]:
        sM, dM = segment_M(xi, L), direct(xi, L, -0.5, -0.5) * RHO * ZP
        sP, dP = segment_psi(xi, L), direct(xi, L, -0.5, 0.5) * RHO
        eM, eP = abs(sM - dM) / abs(dM), abs(sP - dP) / abs(dP)
        worst = max(worst, eM, eP)
        print(f"{L:5d} {xi:6.1f} {xi / L:6.2f}   {eM:12.1e} {eP:12.1e}   {sM.real:+.9f}{sM.imag:+.9f}i")
    assert worst < 1e-8, worst
    assert abs(segment_M(0.0, 25) - 1) < 1e-14
    print(f"   최대 상대차 {worst:.1e}")
    print()

    print("2. (2.7) 1/L 전개의 나머지 (×L² 가 일정하면 O(L^-2))")
    print(f"{'ξ':>5} {'L':>5}   {'|나머지|':>10} {'×L²':>8}   {'|B 만 뺀 차|×L':>14}")
    for xi in (1.0, 4.0):
        prev = None
        for L in (25, 50, 100, 200, 400):
            v = segment_M(xi, L)
            first = v - B_half_half(xi)
            rem = first - (KAPPA * B_mhalf_mhalf(xi) - B_half_mhalf(xi) / RHO) / L
            print(f"{xi:5.1f} {L:5d}   {abs(rem):10.3e} {abs(rem) * L * L:8.4f}   {abs(first) * L:14.6f}")
            if prev is not None:
                assert 3.6 < prev / abs(rem) < 4.4
            prev = abs(rem)
    theta_prime = 0.5 * np.real(digamma(0.25 + 0.5j * GAMMA)) - 0.5 * np.log(np.pi)
    d = abs(KAPPA.real + theta_prime)
    print(f"   Re κ + ϑ_RS'(γ₁) = {d:.1e}")
    assert d < 1e-12
    print()
    print("모든 확인 통과")


if __name__ == "__main__":
    main()
