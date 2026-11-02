# 영점의 국소 기여와 모멘트·복원

제곱 인수 없는 가중 약수합에서 영점의 중복도가 어느 모멘트부터 나타나는지 계산하고, 보조 통계로 영점 높이와 제타함수의 미분값을 복원하는 방법을 살핀다. **[국소 계산]**은 고정된 영점의 기여를 뜻하며, 무한 영점합의 수렴을 포함하지 않는다. RH는 리만 가설이다.

- [§1. 반사된 베타핵](#section-1)
- [§2. 아크사인과 이차 모멘트](#section-2)
- [§3. 모멘트 차수와 영점 중복도](#section-3)
- [§4. 영점 높이와 미분값의 복원](#section-4)
- [부록 A. 영점 기여의 수치 기록](#appendix-a)

<a id="notation"></a>

**기호와 전제.** $r=z+w$, $L=\log x$로 두고, $\mu,\omega,\Lambda$는 뫼비우스 함수, 서로 다른 소인수의 개수, 폰 망골트 함수이다. $\rho=\beta_\rho+i\gamma$는 비자명 영점이며 중복도는 $m$이다. 다루는 지수합과 산술 인자는

$$
S_x(\xi)=\sum_{\substack{dm\le x\\\mu^2(dm)=1}}
z^{\omega(d)}w^{\omega(m)}e^{i\xi\log d/L},\qquad
H(s,t)=\prod_p(1+zp^{-s}+wp^{-t})(1-p^{-s})^z(1-p^{-t})^w
$$

이다. $\mathcal K_{z,w}(u)=g_z(u)g_w(1-u)$, $g_z(v)=v_+^{z-1}/\Gamma(z)$를 전체해석적 분포로 읽고, $\mathcal B_{z,w}(\xi)=\langle\mathcal K_{z,w},e^{i\xi u}\rangle$로 쓴다. 아래의 “S 수준” 모멘트는 로그 좌표의 분모를 $\log x$로 고정한 모멘트이다. 기본 핵 항등식과 이동 변환은 [주정리 §4–5](main-theorem.md#section-4)를 사용한다.

이 문서의 국소 기여와 별개로, 모든 정수의 $d_z$ 이산핵에 대한 정확한 영점 표현은 [전역 공식](global-formula.md)에 있다. 보조 통계 $C_2$의 확률적 극한과 평균제곱 증명은 [극한분포](limiting-distribution.md)에서 다룬다.

<a id="section-1"></a>

## 1. 반사된 핵

**[국소 계산]** 중복도 $m$인 비자명 영점 $\rho=\beta_\rho+i\gamma$에서 $H$가 $(\rho,\rho)$ 근방에 정칙이라고 하자. 이 정칙성은 $\beta_\rho\ge\tfrac12$이면 [주정리 (3.7)](main-theorem.md#equation-3-7)로 확보된다. $b_\rho=\zeta^{(m)}(\rho)/m!$로 두고 가지를 고정하면 $\zeta(s_1)^z\zeta(s_2)^w\approx b_\rho^r(s_1-\rho)^{mz}(s_2-\rho)^{mw}$이다. $s-\rho=q/L$로 옮긴 국소 기여는 고정된 $\rho,\xi$에서

<a id="equation-1-1"></a>

$$
x^\rho L^{-mr-1}\,\frac{b_\rho^rH(\rho,\rho)}{\rho}\,\mathcal B_{-mz,-mw}(\xi)
+O_{\rho,z,w,\xi}\!\left(x^{\beta_\rho}L^{-m\Re r-2}\right)
\tag{1.1}
$$

이다. 반사된 핵은 $\mathcal K_{-mz,-mw}$이며, $m=1$이면 $\mathcal K_{-z,-w}$이다. 이 국소식에는 RH가 필요 없다. RH 아래에서는 $x$의 거듭제곱이 $x^{1/2}$이고, 로그 인자를 포함한 기준 규모는 $x^{1/2}L^{-m\Re r-1}$이다.

<a id="section-2"></a>

## 2. 아크사인: 2차 모멘트부터 영점이 보인다

$z=w=\tfrac12$이고 $\Re\rho\ge\tfrac12$인 단순 영점에서는 $H(\rho,\rho)=1/\zeta(2\rho)$이고 핵은 $\mathcal K_{-1/2,-1/2}=h_2'$이다. [주정리 §7](main-theorem.md#section-7)의 표에서 $\langle h_2',1\rangle=\langle h_2',u\rangle=0$, $\langle h_2',u^2\rangle=-\tfrac14$이다. 따라서 이 영점의 기여는 2차 모멘트부터 나타난다. 실제로 제곱 인수 없는 수의 개수의 급수 $\zeta(s)/\zeta(2s)$는 $\rho$에서 특이점이 없다.

이 핵은 아크사인의 이차 보정 [주정리 (2.5)](main-theorem.md#equation-2-5)에도 나타난다. 대칭 핵에서
$z(\mathcal K_{z-1,z}+\mathcal K_{z,z-1})=\frac{z}{z-1}\mathcal K_{z-1,z-1}$ ($z\ne1$)이므로,
보정 핵과 반사된 핵의 지수가 일치하는 경우는 $z-1=-z$, 즉 $z=\tfrac12$이다.

**정확한 확인.** 제곱 인수 없는 $n$의 균등한 무작위 약수를 $D$라 하면 $\log D=\sum_{p\mid n}\varepsilon_p\log p$이고, $\varepsilon_p$는 독립인 공정한 동전이다. 따라서

<a id="equation-2-1"></a>

$$
2^{-\omega(n)}\sum_{d\mid n}(\log d)^2=\frac14\Big(\log^2n+\sum_{p\mid n}\log^2p\Big).
\tag{2.1}
$$

$T_2(x)=\sum_{n\le x}\mu^2(n)\sum_{p\mid n}\log^2p$의 디리클레 급수는 $\dfrac{\zeta(s)}{\zeta(2s)}\sum_p\dfrac{\log^2p\,p^{-s}}{1+p^{-s}}$이다.

- **영점에서:** $\sum_p\log^2p\,p^{-s}=(\zeta'/\zeta)'(s)-\sum_{k\ge2}k\sum_p\log^2p\,p^{-ks}\approx-(s-\rho)^{-2}$이므로, 단순극과 유수 $-\zeta'(\rho)/\zeta(2\rho)$가 생긴다. 이를 $4L^2$로 나눈 S 수준 $u^2$ 통계의 영점 항은 (1.1)에 $\langle h_2',u^2\rangle=-\tfrac14$를 넣은 것과 정확히 같다.
- **$s=\tfrac12$에서:** $-3\sum_p\log^2p\,p^{-2s}\approx-\tfrac34(s-\tfrac12)^{-2}$이고 $\zeta(s)/\zeta(2s)\approx2\zeta(\tfrac12)(s-\tfrac12)$이다. 그래서 진동 없는 편향 $-3\zeta(\tfrac12)\sqrt x=4.381\sqrt x$가 생긴다. S 수준 $u^2$로는 $-\tfrac34\zeta(\tfrac12)x^{1/2}L^{-2}=4\zeta(\tfrac12)\langle\mathcal T,u^2\rangle x^{1/2}L^{-2}$로, [고차 전개 (3.2)](higher-order.md#equation-3-2)의 세 점 핵과 일치한다.
- **더 낮은 층:** $s=\tfrac13$에서는 $x^{1/3}$에 일차 로그 다항식이 곱해진 항이 나온다. $\zeta(2s)$의 영점 $s=\rho/2$도 특이점 후보이다. **RH와 해당 $\rho$의 단순성 아래** 이 점의 삼중극은 진동수 $\gamma/2$, 규모 $x^{1/4}(\log x)^2$의 항을 만든다.

**보조 통계.** $c_n=\sum_{d\mid n}\Lambda(d)\log d$로 두면 제곱 인수 없는 $n$에서 $c_n=\sum_{p\mid n}\log^2p$, 즉 $\log D$ 분산의 4배이다. 디리클레 급수 $Q(s)=\zeta(s)(\zeta'/\zeta)'(s)$에는 $s=\tfrac12$ 항도 $\zeta(2s)$ 층도 없다. 절단된 명시적 공식을

<a id="equation-2-2"></a>

$$
\sum_{n\le x}c_n=xP_2(\log x)-2\,\mathrm{Re}\sum_{0<\gamma\le T}\frac{\zeta'(\rho)}{\rho}\,x^\rho+R_2(x,T)
\tag{2.2}
$$

로 쓴다. $xP_2(\log x)=\operatorname{Res}_{s=1}Q(s)x^s/s$이고, 합은 서로 다른 비자명 영점에 걸친다. 중복도 $m\ge2$이면 $Q(s)=O((s-\rho)^{m-2})$로 정칙이고 $\zeta'(\rho)=0$이므로 그 항은 사라진다. $R_2$의 RH 아래 평균제곱 추정은 [극한분포 §3](limiting-distribution.md#section-3)에서 증명한다. **[수치]** [부록 A](#appendix-a).

<a id="section-3"></a>

## 3. 모멘트의 차수와 영점의 중복도

공정한 동전의 누적량은 $\kappa_2=\tfrac14$, $\kappa_3=0$, $\kappa_4=-\tfrac18$이고, $\kappa_k(\log D)=\kappa_k(\varepsilon)\sum_{p\mid n}(\log p)^k$이다.

**명제 1.** 정수 $k\ge1$과 $\Re\rho\ge\tfrac12$인 중복도 $m$의 비자명 영점에 대해, $\sum_n\mu^2(n)\sum_{p\mid n}(\log p)^k\,n^{-s}$는 $k>m$이면 $k-m$차 극을 가지고, $k\le m$이면 정칙이다.

**증명.** 해당 급수는 $\zeta(s)/\zeta(2s)$와 $\sum_p(\log p)^kp^{-s}/(1+p^{-s})$의 곱이다. 두 번째 인수에서 $(-d/ds)^{k-1}(-\zeta'/\zeta)(s)\sim-m(k-1)!(s-\rho)^{-k}$를 빼면 $\rho$ 근방에 정칙이다. 소수 제곱 항도 $\zeta(2\rho)\ne0$에 의해 이 근방으로 연장된다. 첫 번째 인수는 정확히 $m$차 영점이므로 결론이 따른다. $\square$

진폭이 $\zeta'(\rho)$이므로 2차 모멘트는 단순 영점만 본다. 특히 $C_2$의 오차에서 얻는 영점 배제도 단순 영점에 한정된다. 아크사인에서는 홀수 누적량이 $0$이라 이중 영점은 4차부터 보인다.

<a id="section-4"></a>

## 4. 영점 높이와 미분값의 복원

$M(x)=\sum_{n\le x}\mu(n)$, $C_2(x)=\sum_{n\le x}c_n$, $\psi(x)=\sum_{n\le x}\Lambda(n)$, $y=\log x$로 둔다.
각 잔차를 $e^{-y/2}$로 정규화하면, 단순 영점 $\rho=\beta_\rho+i\gamma$의 기여는 $2\Re(a_\rho e^{(\rho-1/2)y})$이다. 세 디리클레 급수의 국소 유수를 비교하면 RH 없이

| 잔차 | $a_\rho$ |
|:--|:--|
| $M(x)$ | $1/(\rho\zeta'(\rho))$ |
| $C_2(x)-xP_2(\log x)$ | $-\zeta'(\rho)/\rho$ |
| $\psi(x)-x$ | $-1/\rho$ |

이다. $\psi$는 중복도 $m_\rho$인 영점도 계수 $-m_\rho/\rho$로 포착한다.
따라서 영점은 $\psi$에서 찾고, $C_2$의 복소 계수로 미분값을 추정하는 방안을 택한다.
두 통계는 뫼비우스 반전으로 연결되며, 단순 영점에서

<a id="equation-4-1"></a>

$$
(\mu*c)(n)=\Lambda(n)\log n,\qquad
\frac{a_\rho^{(2)}}{a_\rho^{(\psi)}}=\zeta'(\rho).
\tag{4.1}
$$

공통 컨볼루션 필터의 응답이 복소 주파수 $\nu_\rho=\gamma-i(\beta_\rho-\tfrac12)$에서 정의되고 영이 아니면 계수비가 보존된다. 조화 역산은 $\nu_\rho$를 맞추어 $\beta_\rho$도 추정할 수 있다. 임계선 밖 후보는 별도의 영점 검증으로 확인한다.

높은 영점의 후속 실험에서는 목표 높이로 주파수를 이동하고 해당 대역의 진폭을 보존하는 필터를 사용한다. 절단오차와 경계효과를 평가하며 FFT와 조화 역산(필터 대각화; Main 등)을 비교한다. 성능은 영점 위치를 입력하지 않고 얻은 추정 오차와 누락 수로 측정한다.

<a id="appendix-a"></a>

## 부록 A. 영점 기여의 수치 기록

아래 표와 실험 설명은 초기 원고의 보고값이다. 계산 코드와 원자료가 없어 재현하지 않았으며, [전역 공식의 검산](global-formula.md#checks)과는 검증 범위가 다르다.

**보고된 방법.** $T_2$와 $c_n$은 각각 $5\cdot10^7$, $3\cdot10^7$까지 계산했다. 영점과 $\zeta^\prime(\rho)$, $\zeta(2\rho)$는 Euler–Maclaurin 공식으로 계산했으며, 높이 $1000$까지의 영점 649개를 사용했다.

**알려진 영점으로 잔차 재현.** (2.2)을 적분한 $\int_0^x C_2(t)\,dt$의 잔차를 $x^{3/2}$로 나눈 것과, 알려진 영점 $N$개의 합 $-2\mathrm{Re}\sum\zeta'(\rho)x^{\rho+1}/(\rho(\rho+1))$를 $x\in[10^6,3\cdot10^7]$에서 비교했다.

| 입력 영점 수 | 1 | 10 | 79 | 649 |
|:--|--:|--:|--:|--:|
| 상관 | 0.594 | 0.928 | 0.991 | 0.995 |

잔차의 표준편차는 0.0087, 차이의 표준편차는 0.0009이다. 평활화하지 않은 합에서는 상관이 0.126, 0.383, 0.654, 0.837이다.

**스펙트럼.** 잔차를 $\log x$에 대해 푸리에 변환한 앞의 세 봉우리는 다음과 같다.

| 봉우리 | 14.20 | 21.08 | 25.04 |
|:--|--:|--:|--:|
| $\gamma$ | 14.13 | 21.02 | 25.01 |

**아크사인 설정의 $T_2$.** $x(AL^2+BL+C)+c\sqrt x$로 피팅한 편향은 $c=4.335$로, 예측 $-3\zeta(\tfrac12)=4.381$과 맞는다. 선도 계수 $A=0.30394$도 $1/(2\zeta(2))=0.30396$과 맞는다. 진동은 $\zeta(2s)$ 층과 $s=\tfrac13$ 항이 섞여, 영점 649개로 상관 0.70까지만 나온다.

<a id="references"></a>

## 참고 문헌

- NIST DLMF, [§5.9](https://dlmf.nist.gov/5.9): Hankel 적분. [§13.6](https://dlmf.nist.gov/13.6): Bessel 함수. [§25.16](https://dlmf.nist.gov/25.16): $\psi$의 명시적 공식.
- N. Ng, [The distribution of the summatory function of the Möbius function](https://arxiv.org/abs/math/0310381), Proc. London Math. Soc. (2004): $M(x)$의 극한분포와 $J_{-1}(T)$ 가설.
- J. Main, V. A. Mandelshtam, G. Wunner, H. S. Taylor, [Harmonic inversion as a general method for periodic orbit quantization](https://arxiv.org/abs/chao-dyn/9709009) (1998): 필터 대각화를 이용한 영점 복원.
