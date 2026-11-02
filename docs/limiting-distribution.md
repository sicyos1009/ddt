# RH 아래 보조 통계의 극한분포

보조 통계 $C_2$의 정규화한 오차가 RH 아래 극한분포를 갖는다는 명제와, 그 증명에 필요한 평균제곱 추정을 한 문서에 모았다. RH는 리만 가설, LI는 명제에서 지정한 영점 높이들의 유리수 위 선형독립이다. 영점의 단순성은 가정하지 않는다.

- [§1. 통계와 명시적 공식](#section-1)
- [§2. 극한분포 명제](#section-2)
- [§3. 절단오차의 평균제곱 추정](#section-3)

<a id="section-1"></a>

## 1. 통계와 명시적 공식

$\Lambda$를 폰 망골트 함수라 하고

$$
c_n=\sum_{d\mid n}\Lambda(d)\log d,\qquad
C_2(x)=\sum_{n\le x}c_n,\qquad
Q(s)=\sum_{n\ge1}\frac{c_n}{n^s}
=\zeta(s)(\zeta'/\zeta)'(s)\quad(\Re s>1)
$$

로 둔다. 주항은 $xP_2(\log x)=\operatorname{Res}_{s=1}Q(s)x^s/s$로 정한다. $P_2$는 이차 다항식이다. $y=\log x$, $\rho=\beta_\rho+i\gamma$를 쓰고 절단오차 $R_2$를

$$
C_2(x)=xP_2(\log x)
-2\Re\sum_{0<\gamma\le T}\frac{\zeta'(\rho)}{\rho}x^\rho
+R_2(x,T)
$$

로 정의한다. 합은 서로 다른 비자명 영점에 걸친다. 중복도 $m\ge2$인 영점에서는 $Q(s)=O((s-\rho)^{m-2})$가 정칙이고 $\zeta'(\rho)=0$이므로 해당 항이 사라진다. 이 공식의 산술적 배경은 [영점 기여 §2](zero-contributions.md#section-2)에 있으며, 필요한 절단오차는 이 문서 §3에서 증명한다.

$\mathcal E(y,T)=e^{-y/2}R_2(e^y,T)$로 쓰고, $N(T)$는 양의 높이 $T$ 이하의 영점을 중복도와 함께 세는 함수로 둔다. $J_0$는 제1종 Bessel 함수이며 $\widehat\nu(\xi)=\int e^{i\xi t}\,d\nu(t)$이다.

<a id="section-2"></a>

## 2. RH 아래 극한분포

**명제 1.** RH 아래
$$
E_2(y)=e^{-y/2}\bigl(C_2(e^y)-e^yP_2(y)\bigr)
$$
는 $B^2$-거의 주기함수이며 극한분포 $\nu$를 가진다. 추가로 단순 영점의 서로 다른 양의 높이가 $\mathbb Q$ 위에서 선형독립(LI)이면

<a id="equation-2-1"></a>

$$
\widehat\nu(\xi)=\prod_{\substack{\gamma>0\\\rho\text{ simple}}}
J_0\!\left(\frac{2|\zeta'(\rho)|\,\xi}{|\rho|}\right).
\tag{2.1}
$$

**증명.** §3은 §1의 명시적 공식의 정규화한 나머지 $\mathcal E(y,T)=e^{-y/2}R_2(e^y,T)$에 대해
$$
\lim_{Y\to\infty}\frac1Y\int_{\log2}^{Y}|\mathcal E(y,e^Y)|^2\,dy=0
$$
을 준다. 볼록성 상계와 Cauchy 추정으로 $\zeta'(\tfrac12+i\gamma)\ll_\varepsilon\gamma^{1/4+\varepsilon}$이므로 $r_\gamma=-2\zeta'(\rho)/\rho\ll_\varepsilon\gamma^{-3/4+\varepsilon}$이다. 영점 개수 $N(T+1)-N(T)\ll\log T$와 함께 Akbary–Ng–Shahabi의 따름정리 1.3(a)를 적용하면 $B^2$-거의 주기성과 극한분포의 존재가 따른다. LI 아래 곱 공식은 같은 논문의 정리 1.9이다. $\square$

중복 영점은 계수가 $0$이므로 단순영점 가정이 필요 없다. [영점 기여 §4](zero-contributions.md#section-4)의 잔차 계수는 단순 영점에서 $a_\rho^{(M)}=-(a_\rho^{(\psi)})^2/a_\rho^{(2)}$이지만, 그 제곱합의 제어에는 여전히 역미분값의 추정이 필요하다. 이 결과로 Ng의 $J_{-1}$ 가정을 제거할 수는 없다.

<a id="section-3"></a>

## 3. 절단오차의 평균제곱 추정

$Q(s)=\zeta(s)(\zeta'/\zeta)'(s)$, $a=\tfrac14$로 둔다. 소인수분해에서
$$
0\le c_n=\frac12\sum_{p^h\parallel n}h(h+1)(\log p)^2\le(\log n)^2.
$$
따라서 $2\le x\le T$, $c=1+1/\log x$에서 유효 Perron 공식의 오차는
<a id="equation-3-1"></a>

$$
O\!\left(\frac{x\log^3(2x)}T+\log^2(2x)\right)
\tag{3.1}
$$
이다. $n\asymp x$에서는 $\min(1,x/(T|n-x|))$를 합하고, 나머지는 절대수렴으로 처리한다. 정수 끝점의 반가중치 차이도 이 오차에 포함된다.

RH 아래 $\Re s=a$는 모든 비자명 영점과 거리가 적어도 $1/4$이다. 로그미분의 부분분수 전개와 $N(t+1)-N(t)\ll\log(t+2)$에서 $(\zeta'/\zeta)'(a+it)\ll\log(|t|+2)$를 얻는다. 볼록성 상계를 곱하면, 고정된 $0<\varepsilon<1/8$에 대해
<a id="equation-3-2"></a>

$$
\Phi(t):=\frac{Q(a+it)}{a+it}
\ll_\varepsilon(1+|t|)^{-5/8+\varepsilon},
\qquad \Phi\in L^2(\mathbb R).
\tag{3.2}
$$

$X=e^Y$라 두고 모든 영점 높이에서 $\gg1/\log X$ 떨어진 $T\in[X,X+1]$를 고른다. 높이 $\pm T$에서는 $a\le\sigma\le c$에 걸쳐 $(\zeta'/\zeta)'(\sigma\pm iT)\ll\log^3T$이다. Perron 경로를 $\Re s=a$로 옮길 때 두 수평선의 적분은 볼록성 상계로
$$
\ll_\varepsilon T^\varepsilon\left(\frac xT+\frac{x^{1/4}}{T^{5/8}}\right).
$$
지나는 극은 $s=1$과 단순 영점뿐이며 그 유수가 §1의 명시적 공식의 주항과 영점 합을 준다. 왼쪽 수직선 적분은 $x^{1/4}V_T(\log x)$이고
<a id="equation-3-3"></a>

$$
V_T(y)=\frac1{2\pi}\int_{-T}^{T}\Phi(t)e^{ity}\,dt,
\qquad
\int_{\mathbb R}|V_T(y)|^2\,dy
=\frac1{2\pi}\int_{-T}^{T}|\Phi(t)|^2\,dt\ll1
\tag{3.3}
$$
이다. 마지막 식은 Plancherel 항등식이다. (3.1)과 수평선 추정을 합치면 $\log2\le y\le Y$에서
$$
\mathcal E(y,T)=e^{-y/4}V_T(y)
+O_\varepsilon\!\left((1+y)^2e^{-y/2}+(1+Y)^3e^{-(1/2-\varepsilon)Y}\right).
$$
따라서
<a id="equation-3-4"></a>

$$
\frac1Y\int_{\log2}^{Y}|\mathcal E(y,T)|^2\,dy
\ll_\varepsilon\frac1Y+(1+Y)^6e^{-(1-2\varepsilon)Y}\longrightarrow0.
\tag{3.4}
$$
$T$를 $X=e^Y$로 바꿀 때 추가되는 영점은 $O(\log X)$개이고, 각 정규화한 계수는 $O_\varepsilon(X^{-3/4+\varepsilon})$이다. 이 차이도 균등하게 $0$으로 가므로 §2의 평균제곱 조건이 성립한다.

<a id="references"></a>

## 참고 문헌

- N. Ng, [The distribution of the summatory function of the Möbius function](https://arxiv.org/abs/math/0310381), Proc. London Math. Soc. (2004): $M(x)$의 극한분포와 $J_{-1}(T)$ 가설.
- A. Akbary, N. Ng, M. Shahabi, [Limiting distributions of the classical error terms of prime number theory](https://arxiv.org/pdf/1306.1657), Q. J. Math.: 따름정리 1.3(a)의 극한분포 기준과 정리 1.9의 LI 아래 곱 공식.
