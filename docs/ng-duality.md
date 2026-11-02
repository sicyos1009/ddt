# Ng 공식과 소수·약수 쪽 쌍대

Mertens 함수의 극한분포에 관한 Ng의 결과를 두 방향에서 다룬다. §1은 Ng 쪽의 가정을 약화하고, §2–3은 같은 구조를 소수와 약수 쪽에서 가정 없이 세운다. §4는 두 쪽을 비교하고, §5는 둘을 Mertens 함수의 푸리에 쌍으로 다시 잇는다. 절·정리·수식 번호는 이 문서 안에서 독립적이다.

기호는 [주정리 §8](main-theorem.md#section-8)의 일반화된 약수계수 $d_z$와 이산핵 $\mathscr K_n^{z,w}$, 그리고 [전역 공식 §1.1](global-formula.md#section-1-1)의 소수 거듭제곱 측도 $\lambda_{p,a}^{z,w}$를 따른다. **RH**는 리만 가설, **LI**는 양의 영점 높이들이 $\mathbb Q$ 위에서 선형독립이라는 가설이다.

| 결과 | 대상 | 가정 |
|:--|:--|:--|
| [명제 1.1](#section-1) | $e^{-y/2}M(e^y)$의 극한분포 | RH와 $\sum_\gamma\log(\lvert\gamma\rvert+2)/\lvert\rho\zeta'(\rho)\rvert^2<\infty$, 곱 공식에는 LI |
| [정리 1](#section-2-2) | $\sigma>\tfrac12$에서 $\log\lvert\zeta(\sigma+i\alpha)\rvert$의 극한분포 | 없음 |
| [정리 2](#section-3-2) | $S_1(t)$의 극한분포 | 없음 |
| [명제 5.1](#section-5-2) | 국소 극한 공식: $\lvert1/(s\zeta(s))\rvert^2$의 영점 봉우리 넓이 | 없음 |
| [명제 5.2](#section-5-3) | Ng의 분산과 Mertens 평균제곱 | RH, 등식에는 (1.1) |

- [§1. Ng 쪽: 로그 가중 조건 아래의 극한분포](#section-1)
- [§2. 소수·약수 쪽 쌍대: $\sigma>\tfrac12$의 아크사인 곱](#section-2)
- [§3. 임계선의 쌍대: 누적량 $S_1$](#section-3)
- [§4. Ng 공식과 쌍대의 대응](#section-4)
- [§5. Mertens 함수와의 연결](#section-5)
- [검산](#checks)
- [참고 문헌](#references)

<a id="section-1"></a>

## 1. Ng 쪽: 로그 가중 조건 아래의 극한분포

[극한분포 §2](limiting-distribution.md#section-2)의 방법을 $M(x)$ 자체에 적용하여, Ng의 $J_{-1}(T)\ll T$를 로그 가중 제곱합의 수렴으로 약화한다. 이 절의 논증은 이산핵을 쓰지 않는 고전적인 해석이다. 단순 영점 $\rho=\tfrac12+i\gamma$에 대해

$$
a_\gamma=\frac1{\rho\zeta'(\rho)},\qquad
P_T(y)=\sum_{|\gamma|<T}a_\gamma e^{i\gamma y}
$$

로 두고, 합은 양·음의 높이를 모두 포함한다. $a_{-\gamma}=\overline{a_\gamma}$이므로 $P_T$는 실수값이다. RH 아래 서로 다른 영점은 서로 다른 높이를 가진다.

**명제 1.1 (로그 가중 조건 아래 Ng 분포).** RH와

<a id="equation-1-1"></a>

$$
\sum_\gamma\frac{\log(|\gamma|+2)}{|\rho\zeta'(\rho)|^2}<\infty
\tag{1.1}
$$

를 가정한다. $\zeta'(\rho)=0$인 항은 $+\infty$로 읽으므로, (1.1)은 모든 영점이 단순하다는 조건을 포함한다. 그러면 $\varphi(y)=e^{-y/2}M(e^y)$는 $B^2$-근사주기 함수이고 극한분포 $\nu$를 가진다. $\nu$는 $P_{T_0}$의 극한분포 $\nu_{T_0}$들의 $T_0\to\infty$ 약극한이다. LI를 추가하면

<a id="equation-1-2"></a>

$$
\widehat\nu(\xi)=\prod_{\gamma>0}J_0\!\left(\frac{2\xi}{|\rho\zeta'(\rho)|}\right)
\tag{1.2}
$$

이고, 이 곱은 절대수렴한다.

**가정의 위치.** 어떤 $\theta<2$에 대해 $J_{-1}(T)\ll T^\theta$이면 (1.1)이 성립한다. $|\rho|\ge\gamma$와 부분합산으로

$$
\sum_{0<\gamma\le T}\frac{\log\gamma}{|\rho\zeta'(\rho)|^2}
\le\int_{14}^{T}\frac{\log t}{t^2}\,dJ_{-1}(t)
\ll T^{\theta-2}\log T+\int_{14}^{T}t^{\theta-3}\log t\,dt
$$

가 유계이기 때문이다. 따라서 명제 1.1은 Ng의 $\theta=1$, Akbary–Ng–Shahabi의 $\theta<3-\sqrt3$, Ng(2025)의 $\theta<2$ 조건을 포함한다. 반대로 $J_{-1}(T)\asymp T^2(\log T)^{-3}$이면 (1.1)은 성립하지만 거듭제곱 조건은 성립하지 않는다.

한편 (1.2)의 곱이 특성함수가 되려면 $\sum_\gamma|a_\gamma|^2<\infty$가 필요하다. 발산하면 $\xi\ne0$에서 곱이 $0$으로 가기 때문이다. 같은 조건은 $\varphi$가 계수 $a_\gamma$를 갖는 $B^2$-근사주기 함수일 때 베셀 부등식에서도 나온다. 따라서 필요조건과의 차이는 $\log(|\gamma|+2)$ 가중치뿐이다. 이 가중치의 출처와 개선은 아래 덧붙임에 있다.

**보조정리 1.2 (명시적 공식의 절단).** RH와 모든 영점의 단순성을 가정하고 $0<\varepsilon\le\tfrac14$를 고정한다. 충분히 큰 정수 $\nu$마다 $T_\nu\in[\nu,\nu+1]$이 존재하여, 모든 $x\ge2$에서

<a id="equation-1-3"></a>

$$
M(x)=\sum_{|\gamma|<T_\nu}\frac{x^\rho}{\rho\zeta'(\rho)}
+O_\varepsilon\!\left(\frac{x\log x}{T_\nu^{1-\varepsilon}}+1\right).
\tag{1.3}
$$

**증명.** RH 아래에서는 충분히 큰 각 $\nu$에 대해 $T_\nu\in[\nu,\nu+1]$을 골라, $\tfrac12\le\sigma\le2$에서 $|\zeta(\sigma\pm iT_\nu)|^{-1}\ll_\varepsilon T_\nu^\varepsilon$이 되게 할 수 있다(Titchmarsh 정리 14.16). $-1\le\sigma<\tfrac12$에서는 $\zeta(s)=\chi(s)\zeta(1-s)$, $|\chi(\sigma+it)|\asymp|t|^{1/2-\sigma}$와 $\zeta(\bar s)=\overline{\zeta(s)}$에 의해 $|\zeta(\sigma\pm iT_\nu)|^{-1}\ll_\varepsilon T_\nu^{\sigma-1/2+\varepsilon}\le T_\nu^\varepsilon$이다.

$c=1+1/\log x$로 두면 $|\mu(n)|\le1$이므로 유효 Perron 공식은

$$
M(x)=\frac1{2\pi i}\int_{c-iT_\nu}^{c+iT_\nu}\frac{x^s}{s\zeta(s)}\,ds
+O\!\left(\frac{x\log x}{T_\nu}+1\right)
$$

를 준다. 정수 $x$에서 끝점의 반가중치 차이도 $O(1)$에 들어가며, 계산은 [극한분포 (3.1)](limiting-distribution.md#equation-3-1)과 같다. 경로를 직사각형 $-1\le\Re s\le c$, $|\Im s|\le T_\nu$의 나머지 세 변으로 옮긴다. 안쪽의 극은 $s=0$과 비자명 영점뿐이다. $1/\zeta$는 $s=1$에서 정칙이고, 자명한 영점은 $\Re s\le-2$에 있다. $s=0$의 유수는 $1/\zeta(0)=-2$이고, 단순 영점의 유수는 $x^\rho/(\rho\zeta'(\rho))$이다. 두 수평변의 기여는

$$
\ll T_\nu^{\varepsilon-1}\int_{-1}^{c}x^\sigma\,d\sigma\ll\frac{x\,T_\nu^{\varepsilon-1}}{\log x}
$$

이다. 왼쪽 변에서는 $|\zeta(-1+it)|^{-1}\ll(1+|t|)^{-3/2}$이므로 기여가 $O(x^{-1})$이다. 모두 모으면 (1.3)을 얻는다. $\square$

(1.3)은 Titchmarsh §14.27과 Ng(2004)가 사용한 절단 공식과 같은 형태이며, $J_{-1}$을 쓰지 않는다.

**보조정리 1.3 (Gallagher형 평균제곱).** 실수 $\lambda$의 유한집합 위에서 $S(y)=\sum_\lambda c_\lambda e^{i\lambda y}$라 하고, $Y>0$, $\delta=\pi/Y$로 둔다. 그러면

<a id="equation-1-4"></a>

$$
\int_{-Y}^{Y}|S(y)|^2\,dy
\le\frac{\pi Y^2}{2}\int_{\mathbb R}\Big|\sum_{u<\lambda\le u+\delta}c_\lambda\Big|^2du
\le\frac{\pi^2Y}{2}\sum_\lambda|c_\lambda|^2N_\lambda
\tag{1.4}
$$

이다. 여기서 $N_\lambda=\#\{\lambda':|\lambda'-\lambda|\le\delta\}$이다.

**증명.** $A(u)=\sum_{u<\lambda\le u+\delta}c_\lambda$로 두면

$$
\int_{\mathbb R}A(u)e^{iuy}\,du
=\sum_\lambda c_\lambda\int_{\lambda-\delta}^{\lambda}e^{iuy}\,du
=S(y)\,\frac{1-e^{-i\delta y}}{iy}.
$$

$|y|\le Y=\pi/\delta$이면 $|\sin(\delta y/2)|\ge\delta|y|/\pi$이므로, 마지막 인수의 절댓값은 $2\delta/\pi$ 이상이다. Plancherel 항등식 $\int|\widehat A|^2=2\pi\int|A|^2$에서

$$
\int_{-Y}^{Y}|S|^2
\le\Big(\frac{\pi}{2\delta}\Big)^2\,2\pi\int_{\mathbb R}|A|^2
=\frac{\pi Y^2}{2}\int_{\mathbb R}|A|^2
$$

이다. Cauchy–Schwarz로 $|A(u)|^2\le\#\{\lambda\in(u,u+\delta]\}\sum_{u<\lambda\le u+\delta}|c_\lambda|^2$이다. 각 $\lambda$는 $u\in[\lambda-\delta,\lambda)$에서만 나타나고, 이때 $(u,u+\delta]\subset(\lambda-\delta,\lambda+\delta]$이므로 개수는 $N_\lambda$ 이하이다. 따라서 $\int|A|^2\le\delta\sum_\lambda|c_\lambda|^2N_\lambda$이다. $\square$

이는 Gallagher(1970) 보조정리 1의 상수를 명시한 형태이다.

**명제 1.1의 증명.**

*(i) 꼬리의 평균제곱.* $Y\ge\pi$이면 $\delta\le1$이고, Riemann–von Mangoldt 공식에서 $N_\gamma\le\#\{\gamma':|\gamma'-\gamma|\le1\}\le C_0\log(|\gamma|+2)$이다. 보조정리 1.3을 $S=P_T-P_{T_0}$에 적용하면, $2\le T_0<T$와 $Y\ge\pi$에서

<a id="equation-1-5"></a>

$$
\frac1Y\int_{-Y}^{Y}|P_T(y)-P_{T_0}(y)|^2\,dy\le C_1\tau(T_0),
\qquad
\tau(T_0)=\sum_{|\gamma|\ge T_0}|a_\gamma|^2\log(|\gamma|+2)
\tag{1.5}
$$

이다. $C_1=\pi^2C_0/2$는 절대상수이고, 상계는 $T$와 $Y$에 무관하다. (1.1)에 의해 $\tau(T_0)\to0$이다.

*(ii) 절단의 평균제곱.* 보조정리 1.2를 $\varepsilon=\tfrac14$로 쓴다. 큰 $Y$에 대해 $\nu=\lceil e^Y\rceil$로 두면 $T_\nu\ge e^Y$이다. $\log2\le y\le Y$에서 $x=e^y$로 (1.3)을 정규화하면

$$
\varphi(y)=P_{T_\nu}(y)+R_Y(y),\qquad
|R_Y(y)|\ll y\,e^{y/2}e^{-3Y/4}+e^{-y/2}
$$

이므로 $Y^{-1}\int_{\log2}^{Y}|R_Y|^2\,dy\ll Y^2e^{-Y/2}+Y^{-1}\to0$이다.

*(iii) $B^2$-근사주기성.* (i)과 (ii)에서 각 $T_0$에 대해

$$
\limsup_{Y\to\infty}\frac1Y\int_{\log2}^{Y}|\varphi(y)-P_{T_0}(y)|^2\,dy\le2C_1\tau(T_0)
$$

이다. 오른쪽은 $T_0\to\infty$에서 $0$으로 가고 $P_{T_0}$는 실수값 삼각다항식이므로, $\varphi$는 $B^2$-근사주기 함수이다. $[0,\log2)$에서 $\varphi$는 유계이므로 구간의 시작점은 결과에 영향을 주지 않는다.

*(iv) 극한분포.* $P_{T_0}$의 주파수는 유한개이므로, Kronecker–Weyl 정리에 의해 $P_{T_0}$는 극한분포 $\nu_{T_0}$를 가진다. 서로 다른 주파수의 직교성에서 $\int t^2\,d\nu_{T_0}=\sum_{|\gamma|<T_0}|a_\gamma|^2\le\sum_\gamma|a_\gamma|^2$이므로 $\{\nu_{T_0}\}$는 조밀(tight)하다. $|f|\le1$이고 Lipschitz 상수가 $1$ 이하인 $f$에 대해 $I_Y(g)=Y^{-1}\int_0^Yf(g(y))\,dy$로 두면, Cauchy–Schwarz와 (iii)에서

$$
\limsup_{Y\to\infty}\Big|I_Y(\varphi)-\int f\,d\nu_{T_0}\Big|\le\eta(T_0):=\bigl(2C_1\tau(T_0)\bigr)^{1/2}
$$

이다. 따라서 $I_Y(\varphi)$의 상극한과 하극한의 차이는 모든 $T_0$에서 $2\eta(T_0)$ 이하이므로 극한 $\ell(f)$가 존재하고, $|\ell(f)-\int f\,d\nu_{T_0}|\le\eta(T_0)$이다. 조밀성에 의해 어떤 부분열이 확률측도 $\nu$로 약수렴하며, 위 부등식으로 모든 유계 Lipschitz 함수 $f$에서 $\int f\,d\nu=\ell(f)$이다. 유계 Lipschitz 함수는 약수렴을 결정하므로 $\nu$는 $\varphi$의 극한분포이고 유일하며, 전체 족에서 $\nu_{T_0}\to\nu$이다. 이는 Akbary–Ng–Shahabi 정리 2.9의 일반적인 사실을 이 경우에 맞춰 쓴 것이다.

*(v) LI 아래 곱 공식.* LI 아래 $0<\gamma<T_0$의 높이들은 $\mathbb Q$ 위에서 선형독립이므로 $(\gamma y\bmod2\pi)_{0<\gamma<T_0}$는 전체 토러스에서 균등분포한다. $P_{T_0}(y)=\sum_{0<\gamma<T_0}2|a_\gamma|\cos(\gamma y+\arg a_\gamma)$이므로, $\nu_{T_0}$는 독립 균등 위상 $\theta_\gamma$에 대한 $\sum_{0<\gamma<T_0}2|a_\gamma|\cos\theta_\gamma$의 분포이다. $\mathbb Ee^{i\xi r\cos\theta}=J_0(r\xi)$이므로

$$
\widehat{\nu_{T_0}}(\xi)=\prod_{0<\gamma<T_0}J_0(2|a_\gamma|\xi)
$$

이다. $\nu_{T_0}\to\nu$에서 특성함수가 점별수렴하여 (1.2)를 얻는다. $|1-J_0(t)|\le t^2/4$와 $\sum_\gamma|a_\gamma|^2<\infty$에서 곱은 절대수렴한다. $\square$

**덧붙임 (로그 가중치의 출처와 개선).** (1.5)의 $\log(|\gamma|+2)$는 폭 $2\delta=2\pi/Y$인 창의 영점 수 $N_\gamma$를 단위 창의 영점 수로 어림한 데서만 나온다. RH 아래 $S(t)\ll\log t/\log\log t$(Titchmarsh 정리 14.13)이고, Riemann–Siegel 세타함수 $\vartheta$에 대해 $N(t)=\pi^{-1}\vartheta(t)+1+S(t)$이므로

$$
N_\gamma\le\frac{\delta}{\pi}\log\frac{|\gamma|}{2\pi}
+O\!\left(\frac{\log|\gamma|}{\log\log|\gamma|}\right)
$$

이다. (ii)에서 쓰는 영점은 $|\gamma|<T_\nu\le e^Y+2$를 만족하므로 첫 항은 $O(1)$이다. 따라서 (1.1)의 가중치는 $\log(|\gamma|+3)/\log\log(|\gamma|+3)$로 바꿀 수 있다.

같은 이유로 $\delta\le\pi/\log(|\gamma|/2)$이다. 그러므로 각 영점 주변의 길이 $\asymp1/\log|\gamma|$인 창, 즉 평균 간격 규모의 창에 들어가는 영점 수가 균등하게 유계라면 $N_\gamma=O(1)$이고, 필요조건 $\sum_\gamma|a_\gamma|^2<\infty$만으로 명제 1.1이 성립한다. 필요조건과의 차이는 짧은 구간에 영점이 몰리는 정도로 정확히 측정된다.

**[수치]** [검산 스크립트](../scripts/check_ng_tail.py)는 높이 $300$ 아래의 영점 $138$개와 $\zeta'(\rho)$를 배정밀도로 계산한다. 꼬리 $30\le|\gamma|<300$에서 (1.4)의 세 변을 $Y=8,20,60$에 대해 비교하며, 부등식이 모두 성립한다. $Y=60$에서 왼쪽 변은 $0.8819$, 대각 근사 $2Y\sum|c_\lambda|^2$는 $0.8820$, 가운데 변은 $2.176$이다. (1.1)의 $\gamma>0$ 부분합은 $T=50,100,200,300$에서 $0.0381,0.0421,0.0445,0.0453$이다.

    python -B scripts/check_ng_tail.py

<a id="section-2"></a>

## 2. 소수·약수 쪽 쌍대: $\sigma>\tfrac12$의 아크사인 곱

Ng의 정리와 §1에서는 $e^{-y/2}M(e^y)$가 $y=\log x$의 근사주기 함수이고, 그 주파수는 영점의 높이이다. 극한분포의 특성함수 $\prod_\gamma J_0(2\eta/|\rho\zeta'(\rho)|)$를 얻으려면 영점의 위치(RH), 진폭의 제곱합($J_{-1}$ 류), 높이들의 독립성(LI)을 가정해야 한다.

이 절에서는 같은 구조를 소수와 약수 쪽에서 세운다. 변수는 약수의 대칭 주파수 $\alpha$이고, 주파수는 $\log p$이다. 독립성은 FTA에서, 진폭의 제곱합은 $\sum_pp^{-2\sigma}<\infty$에서 나오므로 **아무 가정도 쓰지 않는다**. 결과는 Bohr–Jessen 정리의 $\log|\zeta|$ 부분과 같으며, 여기서는 그것을 이산핵의 언어로 다시 쓰고 증명한다.

<a id="section-2-1"></a>

### 2.1. 대칭 비틀기와 소수별 인수

실수 $z$와 $\sigma>1$, $\alpha\in\mathbb R$에 대해 [주정리 §8의 이산핵](main-theorem.md#section-8)을 중심 좌표 $2u-1$로 비틀어

$$
\mathfrak F_\sigma^z(\alpha)
:=1+\sum_{n\ge2}n^{-\sigma}\left\langle\mathscr K_n^{z,z},e^{i\alpha\log n\,(2u-1)}\right\rangle
=\sum_{d,m\ge1}\frac{d_z(d)d_z(m)}{(dm)^\sigma}\Big(\frac dm\Big)^{i\alpha}
$$

로 둔다. $u=\log d/\log n$이면 $e^{i\alpha\log n(2u-1)}=(d/m)^{i\alpha}$이다. [주정리 명제 8.2](main-theorem.md#section-8-2)에서 $s=\sigma-i\alpha$, $t=\sigma+i\alpha$로 두고 $\zeta(\bar s)=\overline{\zeta(s)}$를 쓰면

<a id="equation-2-1"></a>

$$
\boxed{
\mathfrak F_\sigma^z(\alpha)=\zeta(\sigma-i\alpha)^z\zeta(\sigma+i\alpha)^z=|\zeta(\sigma+i\alpha)|^{2z}
}
\tag{2.1}
$$

이다. 두 거듭제곱은 서로 켤레인 Euler 곱 로그로 정해지므로 곱은 양수이다. 이는 반사 대칭 $R_*\mathscr K_n^{z,z}=\mathscr K_n^{z,z}$([주정리 명제 8.1](main-theorem.md#section-8-1))의 생성함수 판이다. $\tfrac12<\sigma\le1$에서는 (2.1)의 오른쪽으로 $\mathfrak F_\sigma^z$를 정의하고

$$
\mathscr L_\sigma^z(\alpha):=\log\mathfrak F_\sigma^z(\alpha)=2z\log|\zeta(\sigma+i\alpha)|
$$

로 둔다. 이 값은 직선 $\Re s=\sigma$ 위의 영점, 즉 가산개의 $\alpha$에서만 $-\infty$이다.

**소수별 인수.** [전역 공식 (1.1)](global-formula.md#equation-1-1)의 소수 거듭제곱 측도 $\lambda_{p,a}^{z,z}$는 원자 $j\log p$에 $d_z(p^j)d_z(p^{a-j})$를 놓는다. 대칭 비틀기에서 이 원자의 위상은 $e^{i\theta(2j-a)}$, $\theta=\alpha\log p$이다. 따라서 $r=p^{-\sigma}$에 대해

<a id="equation-2-2"></a>

$$
\sum_{a\ge0}r^a\sum_{j=0}^ad_z(p^j)d_z(p^{a-j})e^{i\theta(2j-a)}
=(1-re^{i\theta})^{-z}(1-re^{-i\theta})^{-z}
=|1-re^{i\theta}|^{-2z}
\tag{2.2}
$$

이다. 특히 $z=\tfrac12$의 일차 원자 $\lambda_{p,1}=\tfrac12(\delta_0+\delta_{\log p})$, 즉 DDT의 동전 하나는 $\tfrac12(e^{-i\theta}+e^{i\theta})=\cos\theta$가 된다. 위상 $\theta$가 균등하면 이것이 아크사인 분포이다.

<a id="section-2-2"></a>

### 2.2. 정리

**정리 1 (Ng 공식의 쌍대).** $\sigma>\tfrac12$이고 $z$가 실수라 하자. $\alpha$를 $[0,A]$에서 균등하게 고를 때, $A\to\infty$이면 $\mathscr L_\sigma^z(\alpha)$는 극한분포 $\nu_{\sigma,z}$를 가진다. $\nu_{\sigma,z}$는 독립 균등 위상 $\theta_p$에 대한

<a id="equation-2-3"></a>

$$
\sum_pY_p,\qquad Y_p=-2z\log|1-p^{-\sigma}e^{i\theta_p}|
\tag{2.3}
$$

의 분포이며, 이 급수는 거의 확실히 수렴한다. 특성함수는

<a id="equation-2-4"></a>

$$
\boxed{
\widehat{\nu_{\sigma,z}}(\eta)
=\prod_p{}_2F_1\!\left(iz\eta,iz\eta;1;p^{-2\sigma}\right)
=\sum_{d\ge1}\frac{d_{iz\eta}(d)^2}{d^{2\sigma}}
}
\tag{2.4}
$$

이고, 곱과 급수는 절대수렴한다. 각 인수는

$$
{}_2F_1\!\left(iz\eta,iz\eta;1;p^{-2\sigma}\right)=J_0\!\left(\frac{2z\eta}{p^\sigma}\right)+O_{z,\eta}\!\left(p^{-4\sigma}\right)
$$

이므로 소수마다 진폭 $2|z|p^{-\sigma}$의 아크사인 인수가 붙고, 고차항은 소수 거듭제곱의 보정이다. RH, LI, $J_{-1}$ 류의 가정이나 영점 밀도 추정은 쓰지 않는다.

**(2.4) 오른쪽의 해석.** $a=iz\eta$이면 $\log\zeta(\sigma-i\alpha)=\overline{\log\zeta(\sigma+i\alpha)}$에서 복소 지수에서도 $\zeta(\sigma-i\alpha)^a\zeta(\sigma+i\alpha)^a=|\zeta(\sigma+i\alpha)|^{2a}$이므로 $e^{i\eta\mathscr L_\sigma^z}=\mathfrak F_\sigma^{a}$이다. $\sigma>1$에서 $\mathfrak F_\sigma^{a}(\alpha)$의 $\alpha$-평균은 상수항이고, 상수항은 $d=m$, 즉 $n=d^2$의 중심 원자 $u=\tfrac12$에서만 나온다. 따라서

$$
\widehat{\nu_{\sigma,z}}(\eta)
=\sum_{d\ge1}d^{-2\sigma}\,\mathscr K_{d^2}^{a,a}\big(\{\tfrac12\}\big),
\qquad a=iz\eta
$$

이다($d=1$ 항은 $1$). 극한분포의 특성함수가 **허수 지수 이산핵의 중심 원자 합**으로 주어진다. 정리 1은 이 등식을 $\sigma>\tfrac12$까지 연장한다.

**보조정리 2.1 (FTA와 위상의 독립).** $\{\log p\}$는 $\mathbb Q$ 위에서 선형독립이다. 따라서 각 $X$에서 $(\alpha\log p\bmod2\pi)_{p\le X}$는 $\alpha\in[0,A]$, $A\to\infty$일 때 토러스 $\mathbb T^{\pi(X)}$에서 균등분포한다. 그러므로

$$
\mathscr L_{\sigma,X}^z(\alpha):=-2z\sum_{p\le X}\log|1-p^{-\sigma-i\alpha}|
$$

는 $\sum_{p\le X}Y_p$의 분포를 극한분포로 가진다. 또한 $\mathbb Ee^{i\eta Y_p}={}_2F_1(iz\eta,iz\eta;1;p^{-2\sigma})$이다.

**증명.** 정수 $c_p$가 모두 $0$은 아닌데 $\sum c_p\log p=0$이면 $\prod_{c_p>0}p^{c_p}=\prod_{c_p<0}p^{-c_p}$가 서로 다른 두 소인수분해가 되어 FTA에 모순이다. 균등분포는 Kronecker–Weyl 정리이다. 인수는 $a=iz\eta$로 $(1-re^{\pm i\theta})^{-a}=\sum_jd_a(p^j)r^je^{\pm ij\theta}$를 곱하고, $\theta$에 대한 직교성으로 $\mathbb E|1-re^{i\theta}|^{-2a}=\sum_jd_a(p^j)^2r^{2j}$를 얻는다. 이것이 ${}_2F_1(a,a;1;r^2)$이다. $\square$

**보조정리 2.2 (평균제곱 근사).** $\sigma>\tfrac12$, $X\ge2$, $Q_X(s)=\prod_{p\le X}(1-p^{-s})$라 하면

<a id="equation-2-5"></a>

$$
\lim_{A\to\infty}\frac1A\int_A^{2A}\left|\zeta(\sigma+it)Q_X(\sigma+it)-1\right|^2dt
=\varepsilon_X:=\sum_{\substack{m>1\\P^-(m)>X}}m^{-2\sigma}
=\zeta(2\sigma)\prod_{p\le X}(1-p^{-2\sigma})-1
\tag{2.5}
$$

이고, $X\to\infty$이면 $\varepsilon_X\to0$이다. 여기서 $P^-(m)$은 $m$의 가장 작은 소인수이다.

**증명.** Titchmarsh 정리 4.11을 $x=2A$로 쓰면 $A\le t\le2A$에서 $\zeta(s)=\sum_{n\le2A}n^{-s}+O(A^{-\sigma})$이다. $P_X=\prod_{p\le X}p$에 대해 $Q_X(s)=\sum_{k\mid P_X}\mu(k)k^{-s}$이므로 곱 $\sum_{n\le2A}n^{-s}Q_X(s)=\sum_mb_mm^{-s}$의 계수는 다음과 같다. $m\le2A$이면 $b_m=\sum_{k\mid(m,P_X)}\mu(k)=\mathbf 1_{P^-(m)>X}$이고, $2A<m\le2AP_X$이면 $|b_m|\le2^{\pi(X)}$이다.

Montgomery–Vaughan 평균값 정리 $\int_A^{2A}|\sum_mc_mm^{-it}|^2dt=\sum_m|c_m|^2(A+O(m))$를 적용한다. $1<m\le2A$ 부분은 $A\sum_{1<m\le2A,\,P^-(m)>X}m^{-2\sigma}+O(\sum_{m\le2A}m^{1-2\sigma})=A\varepsilon_X+o(A)$를 준다. $m>2A$ 부분은 $\ll_XA^{2-2\sigma}=o(A)$이고, $O(A^{-\sigma}|Q_X|)$ 항과 교차항은 Cauchy–Schwarz로 $o(A)$이다. 마지막으로 $\sum_pp^{-2\sigma}<\infty$이므로 $\varepsilon_X=\prod_{p>X}(1-p^{-2\sigma})^{-1}-1\to0$이다. $\square$

이 보조정리는 영점에 관한 정보를 쓰지 않는다. 직선 $\Re s=\sigma$ 근처에 영점이 있으면 그곳에서 $|\zeta Q_X-1|$이 $1$에 가까워질 뿐이며, 그런 $t$의 밀도가 (2.5)로 제어된다.

**정리 1의 증명.** $g_X=\zeta Q_X$를 $\sigma+i\alpha$에서 평가하면 $\mathscr L_\sigma^z-\mathscr L_{\sigma,X}^z=2z\log|g_X|$이다. $0<\delta\le\tfrac12$이고 $|g_X-1|\le\delta$이면 $|\log|g_X||\le2\delta$이다.

*(i) 나쁜 집합.* $[1,A]$를 이진 구간 $[B,2B]$로 나누고 각 구간에 Chebyshev 부등식과 (2.5)를 쓰면

$$
\limsup_{A\to\infty}\frac1A\,\mathrm{meas}\{\alpha\in[0,A]:|g_X-1|>\delta\}\le\frac{\varepsilon_X}{\delta^2}
$$

이다. $\delta=\varepsilon_X^{1/3}$로 잡으면, 위쪽 밀도가 $\varepsilon_X^{1/3}$ 이하인 집합 밖에서 $|\mathscr L_\sigma^z-\mathscr L_{\sigma,X}^z|\le4|z|\varepsilon_X^{1/3}$이다.

*(ii) 극한분포.* $|f|\le1$이고 Lipschitz 상수가 $1$ 이하인 $f$에 대해, (i)과 보조정리 2.1에서

$$
\limsup_{A\to\infty}\left|\frac1A\int_0^Af\big(\mathscr L_\sigma^z(\alpha)\big)\,d\alpha-\mathbb Ef\Big(\sum_{p\le X}Y_p\Big)\right|
\le(4|z|+2)\,\varepsilon_X^{1/3}
$$

이다. Jensen 공식에서 $\mathbb EY_p=0$이고, $\mathbb E(\log|1-re^{i\theta}|)^2=\tfrac12\sum_kr^{2k}/k^2$이므로 $\sum_p\operatorname{Var}Y_p<\infty$이다. Kolmogorov의 두 급수 정리로 $\sum_pY_p$는 거의 확실히 수렴하고, $\mathbb Ef(\sum_{p\le X}Y_p)\to\mathbb Ef(\sum_pY_p)$이다. $X\to\infty$로 보내면 모든 유계 Lipschitz 함수 $f$에서 시간평균이 $\mathbb Ef(\sum_pY_p)$로 수렴하므로, $\nu_{\sigma,z}$는 (2.3)의 분포이다.

*(iii) 특성함수.* 독립성에서 특성함수는 소수별 인수의 곱이다. $|(a)_j/j!|\le(|a|)_j/j!$이므로 ${}_2F_1(a,a;1;r^2)-1=a^2r^2+O_{|a|}(r^4)$이고, $\sum_pp^{-2\sigma}<\infty$에서 곱이 절대수렴한다. 같은 상계 $|d_a(d)|\le d_{|a|}(d)$에서 $\sum_dd_a(d)^2d^{-2\sigma}$도 절대수렴하며, 곱셈성으로 곱과 같다. $a=iz\eta$이면 $a^2r^2=-z^2\eta^2r^2$이고 $J_0(2z\eta r)=1-z^2\eta^2r^2+O(r^4)$이므로 인수의 선도항이 아크사인이다. $\square$

실수 지수 $k$에서는 같은 인수 ${}_2F_1(k,k;1;p^{-2\sigma})$가 무작위 Euler 곱의 모멘트 $\mathbb E|\zeta(\sigma,\mathbb X)|^{2k}$를 주고, 선도항은 $J_0$ 대신 수정 Bessel 함수 $I_0$가 된다. 정리 1의 특성함수는 그 허수 지수 판이다.

<a id="section-3"></a>

## 3. 임계선의 쌍대: 누적량 $S_1$

§2의 곱은 $\sigma=\tfrac12$에서 퇴화한다(§4.3). 진폭 $p^{-1/2}$의 제곱합 $\sum_p1/p$가 발산하기 때문이며, 이는 영점과 무관한 소수의 성질이다. 반면 Ng 쪽이 $x^{-1/2}$ 정규화에서 퇴화하지 않는 것은 $M(x)$가 누적합이라 Mellin 변환에 $1/s$가 붙고, 그것이 진폭 $1/(\rho\zeta'(\rho))$에 주파수의 역수 $1/\rho$를 주기 때문이다. 쌍대 쪽에서 같은 역할을 하는 것은 $\sigma$ 방향의 적분이다.

$$
\int_{1/2}^\infty p^{-\sigma}\,d\sigma=\frac1{\sqrt p\,\log p},
\qquad
\sum_p\frac1{p\log^2p}<\infty.
$$

이 절은 이 누적량이 $S(t)=\pi^{-1}\arg\zeta(\tfrac12+it)$의 적분 $S_1$과 같고, 임계선에서 퇴화하지 않는 소수별 아크사인 곱을 극한분포로 가짐을 보인다. **RH를 포함하여 아무 가정도 쓰지 않는다.**

<a id="section-3-1"></a>

### 3.1. 누적량과 Littlewood 공식

$t\in\mathbb R$에 대해

<a id="equation-3-1"></a>

$$
\mathscr I(t)=\int_{1/2}^\infty\log|\zeta(\sigma+it)|\,d\sigma
\tag{3.1}
$$

로 둔다. $\sigma\ge2$에서 $\log|\zeta(\sigma+it)|=O(2^{-\sigma})$이고 영점과 극점에서의 로그 특이점은 적분 가능하므로, 적분은 모든 $t$에서 절대수렴한다. $\mathscr I(0)=\int_{1/2}^\infty\log|\zeta(\sigma)|\,d\sigma$는 상수이다.

$S(t)$는 $2$에서 $2+it$, 다시 $\tfrac12+it$로 가는 선분을 따라 연속적으로 정한 편각이며, 영점 높이에서는 좌우 극한의 평균이다. 문헌의 관례대로

$$
S_1(t)=\int_0^tS(u)\,du+\frac{\mathscr I(0)}{\pi}
$$

로 두면, 모든 $t$에서

<a id="equation-3-2"></a>

$$
\boxed{\pi S_1(t)=\mathscr I(t)}
\tag{3.2}
$$

이다. 직사각형 $\tfrac12\le\sigma\le A$, $0\le\Im s\le t$의 경계에서 $\log\zeta$를 적분하면 된다. $s=1$의 극은 위로 비켜 가고, $\log\zeta$는 오른쪽 변에서 수평으로 이어 붙인다. Littlewood 보조정리에서 경계 적분은 $-2\pi i\int\nu(\sigma)\,d\sigma$로 순허수이므로, 실수부를 취하면 영점의 위치와 무관하게 (3.2)가 나온다(Titchmarsh §9.9). 이 식은 RH가 필요 없다. 따라서 $S_1(t)=O(\log t)$도 무조건 성립한다.

$\pi S_1$의 Dirichlet 급수형 전개는 $\Re\log\zeta(s)=\Re\sum_n\Lambda(n)n^{-s}/\log n$을 $\sigma$로 적분한 것이다.

$$
\pi S_1(t)\ \approx\ \Re\sum_{n\ge2}\frac{\Lambda(n)}{n^{1/2+it}\log^2n}
=\sum_p\sum_{k\ge1}\frac{\cos(kt\log p)}{k^2p^{k/2}\log p}.
$$

이 급수는 절대수렴하지 않으므로, 아래에서는 이 형태를 평균의 의미에서 정당화한다.

<a id="section-3-2"></a>

### 3.2. 정리

**정리 2 (임계선의 쌍대).** $t$를 $[0,T]$에서 균등하게 고르면, $T\to\infty$일 때 $\pi S_1(t)=\mathscr I(t)$는 극한분포 $\mu$를 가진다. $\mu$는 독립 균등 위상 $\theta_p$에 대한

<a id="equation-3-3"></a>

$$
\sum_pY_p,\qquad
Y_p=-\int_{1/2}^\infty\log|1-p^{-\sigma}e^{i\theta_p}|\,d\sigma
=\sum_{k\ge1}\frac{\cos(k\theta_p)}{k^2p^{k/2}\log p}
\tag{3.3}
$$

의 분포이며, 이 급수는 거의 확실히 수렴한다. 특성함수는

<a id="equation-3-4"></a>

$$
\boxed{
\widehat\mu(\eta)=\prod_p\mathbb Ee^{i\eta Y_p},
\qquad
\mathbb Ee^{i\eta Y_p}=J_0\!\left(\frac{\eta}{\sqrt p\,\log p}\right)+O_\eta(p^{-2})
}
\tag{3.4}
$$

이고 곱은 절대수렴한다. $\mu$의 평균은 $0$, 분산은 $\sum_p\sum_k(2k^4p^k\log^2p)^{-1}$이며, 3차 누율은 양수이다. RH, LI, $J_{-1}$ 류의 가정이나 영점 밀도 추정은 쓰지 않는다.

$Y_p$는 §2의 소수별 변수 $-\log|1-p^{-\sigma}e^{i\theta_p}|$ ($z=\tfrac12$)를 $\sigma$로 적분한 것이다. 각 소수의 선도항 $\cos\theta_p/(\sqrt p\log p)$는 아크사인 분포를 따른다. §2의 중심 원자 공식 (2.4)는 $\sigma$ 적분이 곱셈 구조를 깨므로 여기서는 성립하지 않는다.

**비대칭.** $a_k=(k^2p^{k/2}\log p)^{-1}>0$이고 $\mathbb E[\cos j\theta\cos k\theta\cos l\theta]=\tfrac14(\mathbf 1_{j+k=l}+\mathbf 1_{k+l=j}+\mathbf 1_{l+j=k})\ge0$이므로 $\mathbb EY_p^3\ge\tfrac34a_1^2a_2>0$이다. 따라서 $\mu$는 오른쪽으로 치우친다. 원인은 소수 제곱 항 $k=2$이며, Chebyshev 편향과 같은 출처이다. Ng 쪽의 분포는 LI 아래 대칭이므로 이 현상은 쌍대 쪽에만 있다.

<a id="section-3-3"></a>

### 3.3. 증명

$X\ge2$, $P_X=\prod_{p\le X}p$, $Q_X(s)=\prod_{p\le X}(1-p^{-s})$, $g_X=\zeta Q_X$로 둔다. $Q_X$는 $\Re s>0$에서 영이 아니므로, $\Re s\ge\tfrac12$에서 $g_X$의 영점은 $\zeta$의 영점과 같다. 또

$$
\mathscr I_X(t):=\sum_{p\le X}Y_p(t\log p)=-\int_{1/2}^\infty\log|Q_X(\sigma+it)|\,d\sigma,
\qquad
\mathscr I(t)-\mathscr I_X(t)=\int_{1/2}^\infty\log|g_X(\sigma+it)|\,d\sigma
$$

이다. $Y_p$는 $\theta$에 대해 짝함수이므로 위상의 부호는 상관없다.

**보조정리 3.1 (임계선 근처의 평균제곱).** $\sigma=\tfrac12+\delta$, $0<\delta\le\tfrac32$, $T\ge2$이면

$$
\frac1T\int_T^{2T}|g_X(\sigma+it)-1|^2\,dt\le3\varepsilon_X(\sigma)+C_X(\log T)\,T^{-2\min(\delta,1/2)},
$$

$$
\varepsilon_X(\sigma)=\prod_{p>X}(1-p^{-2\sigma})^{-1}-1
\ll\Big(1+\frac1{\delta\log X}\Big)e^{-2\delta\log X}+\frac1X
$$

이다. 따라서

<a id="equation-3-5"></a>

$$
\int_{1/2}^{2}\Big(\frac1T\int_T^{2T}|g_X(\sigma+it)-1|^2\,dt\Big)^{1/2}d\sigma
\ll\frac1{\log X}+\frac{C_X'}{\sqrt{\log T}}
\tag{3.5}
$$

이다.

**증명.** 보조정리 2.2의 증명을 $\sigma$에 대한 의존성을 남긴 채 반복한다. Titchmarsh 정리 4.11에서 $T\le t\le2T$, $\tfrac12\le\sigma\le2$이면 $\zeta(s)=\sum_{n\le2T}n^{-s}+O(T^{-\sigma})$이다. 곱 $\sum_{n\le2T}n^{-s}Q_X(s)=\sum_mb_mm^{-s}$의 계수는 $m\le2T$에서 $\mathbf 1_{P^-(m)>X}$, $2T<m\le2TP_X$에서 절댓값 $2^{\pi(X)}$ 이하이다. Montgomery–Vaughan 평균값 정리로 첫 부분의 평균제곱은 $\varepsilon_X(\sigma)+O(T^{-1}\sum_{m\le2T}m^{1-2\sigma})$ 이하이고, 둘째 부분은 $\ll_X\sum_{2T<m\le2TP_X}m^{-2\sigma}(1+m/T)$이다. 둘 다 $\ll_X(\log T)T^{-2\min(\delta,1/2)}$에 들어간다.

$\varepsilon_X$의 상계에는 Mertens 정리 $\sum_{p\le u}1/p=\log\log u+M+O(1/\log u)$를 부분합산에 쓴다. 그러면 $y=2\delta\log X$에 대해 $\sum_{p>X}p^{-1-2\delta}\le E_1(y)+O(X^{-2\delta}/\log X)$이고, $E_1(y)\le e^{-y}\log(1+1/y)$이다. 소수 거듭제곱 $k\ge2$의 기여는 $O(1/X)$이다. 따라서 $\varepsilon_X=e^V-1$에서 $V\le e^{-y}\log(1+1/y)+O(e^{-y}/\log X+1/X)$이고, $y<1$이면 $\varepsilon_X\ll1/y$, $y\ge1$이면 $\varepsilon_X\ll e^{-y}+1/X$이다.

(3.5)는 Cauchy–Schwarz $\sqrt{a+b}\le\sqrt a+\sqrt b$와 $\int_0^\infty\sqrt{(1+1/y)e^{-y}}\,dy<\infty$, 그리고 $\int_0^\infty T^{-\min(\delta,1/2)}d\delta\ll1/\log T+T^{-1/2}$에서 나온다. $\varepsilon_X$의 적분이 작은 이유는 $\varepsilon_X$가 큰 영역이 폭 $1/\log X$의 띠뿐이고, 그 안에서도 $\varepsilon_X^{1/2}\ll(\delta\log X)^{-1/2}$가 적분 가능하기 때문이다. $\square$

**보조정리 3.2 (Littlewood 보조정리의 양성).** $\tfrac12\le\sigma_0\le2$이고 $T$, $2T$가 영점의 높이가 아니면

$$
\int_T^{2T}\big|\log|g_X(\sigma_0+it)|\big|\,dt
\le2\int_T^{2T}|g_X(\sigma_0+it)-1|\,dt+C\Big(\frac TX+\log T+\pi(X)\Big)
$$

이다. $C$는 절대상수이다.

**증명.** 직사각형 $[\sigma_0,2]\times[T,2T]$에 Littlewood 보조정리를 적용하면

$$
\int_T^{2T}\log|g_X(\sigma_0+it)|\,dt
=\int_T^{2T}\log|g_X(2+it)|\,dt+\int_{\sigma_0}^{2}\big(\arg g_X(\sigma+iT)-\arg g_X(\sigma+2iT)\big)d\sigma
+2\pi\int_{\sigma_0}^{2}\nu(\sigma)\,d\sigma
$$

이다. 여기서 $\nu(\sigma)\ge0$은 직사각형 안에서 실수부가 $\sigma$보다 큰 영점의 개수이다. 마지막 항이 음이 아니므로 오른쪽의 앞 두 항을 $B$라 하면 $\int\log|g_X|\ge B$이다. 따라서

$$
\int\big|\log|g_X|\big|=2\int\log^+|g_X|-\int\log|g_X|\le2\int\log^+|g_X|+|B|
$$

이다. $\log^+|g|\le\log(1+|g-1|)\le|g-1|$이다. $\sigma=2$에서는 $|g_X-1|\le\sum_{m>X}m^{-2}<1/X$이므로 첫 항은 $O(T/X)$이다. 편각에는 $|\arg\zeta(\sigma+iT)|\ll\log T$ ($\sigma\ge\tfrac12$, Titchmarsh §9.4)를 쓴다. 또 $\Re(1-p^{-s})>0$이므로 $|\arg Q_X(\sigma+iT)|\le\tfrac\pi2\pi(X)$이다. $\square$

영점은 $\log|\zeta|$를 국소적으로 $-\infty$로 보내지만, 이 보조정리에서는 평균을 **키우는** 쪽으로만 들어간다. Jensen 공식과 같은 원리이다. 그래서 임계선 밖의 영점에 대해 아무것도 알 필요가 없다.

**정리 2의 증명.** $T$와 $2T$를 길이 $1$ 이내로 옮겨 영점의 높이를 피해도 평균은 $O(\log T/T)$만 변한다. $\mathscr I=\pi S_1=O(\log t)$이기 때문이다. Fubini 정리와 보조정리 3.2를 $\sigma_0\in[\tfrac12,2]$에서 적분하고, $\sigma\ge2$에서는 $|\log|g_X||\le2|g_X-1|\le4X^{1-\sigma}$를 쓴다. 그러면

$$
\frac1T\int_T^{2T}|\mathscr I(t)-\mathscr I_X(t)|\,dt
\le2\int_{1/2}^{2}\Big(\frac1T\int_T^{2T}|g_X-1|^2\Big)^{1/2}d\sigma
+O\Big(\frac1X+\frac{\log T+X}T\Big)
$$

이다. (3.5)에서

<a id="equation-3-6"></a>

$$
\limsup_{T\to\infty}\frac1T\int_T^{2T}|\mathscr I(t)-\mathscr I_X(t)|\,dt\ll\frac1{\log X}
\tag{3.6}
$$

이다. $\mathscr I_X(t)$는 $(t\log p\bmod2\pi)_{p\le X}$의 연속함수이므로, 보조정리 2.1(FTA와 Kronecker–Weyl)에 의해 $\sum_{p\le X}Y_p$의 분포를 극한분포로 가진다.

이후는 정리 1의 증명 (ii)와 같다. 유계 Lipschitz 함수 $f$에 대해 시간평균의 차이는 (3.6)으로 제어되고, $[0,T]$는 이진 구간으로 나눈다. $\mathbb EY_p=0$이고 $\operatorname{Var}Y_p=\sum_ka_k^2/2\ll1/(p\log^2p)$이므로 Kolmogorov 정리로 $\sum_pY_p$가 거의 확실히 수렴하며, $\mu$는 그 분포이다.

특성함수는 독립성에서 소수별 인수의 곱이다. $Y_p=a_1\cos\theta_p+a_2\cos2\theta_p+O(p^{-3/2})$이고 $\mathbb E[\cos2\theta\,e^{iz\cos\theta}]=-J_2(z)=O(z^2)$이므로 $\mathbb Ee^{i\eta Y_p}=J_0(\eta a_1)+O_\eta(a_1^2a_2+a_2^2)=J_0(\eta a_1)+O_\eta(p^{-2})$이다. $|1-\mathbb Ee^{i\eta Y_p}|\le\eta^2\mathbb EY_p^2/2$에서 곱은 절대수렴한다. 평균, 분산, 3차 누율은 독립 합의 누율이 더해진다는 사실과 위의 비대칭 계산에서 나온다. $\square$

**RH의 역할.** RH를 가정하면 Selberg의 근사(Titchmarsh §14.20–14.21)로 $\mathscr I$를 길이 $T^\varepsilon$의 소수 다항식에 직접 붙일 수 있어 수렴 속도를 개선할 수 있다. 그러나 극한분포의 존재와 곱 공식에는 RH가 필요 없다.

**[수치]** [검산 스크립트](../scripts/check_s1_distribution.py)는 다음을 확인한다.

- 높이 $300$ 아래의 영점 $138$개로 $S(u)=N(u)-1-\vartheta(u)/\pi$를 정확히 적분하여, $\pi\int_0^tS-\mathscr I(t)$가 $t=20,\dots,299$의 일곱 점에서 $3.6\times10^{-6}$ 안에서 상수임을 확인한다. 그 값 $-2.56779$는 (3.1)을 $t=0$에서 직접 적분한 $-\mathscr I(0)$과 일치한다.
- $t\in[1500,3000]$에서 (3.6)의 평균 오차 $\frac1T\int|\mathscr I-\mathscr I_X|$를 잰다. $X=2,13,101,1009,10007$에서 각각 $0.405$, $0.139$, $0.077$, $0.043$, $0.027$로 줄어들며, 모두 $1/\log X$보다 작다.
- $t\in[10,3010]$의 $6\times10^4$개 점에서 $\mathscr I(t)$를 계산하여 소수 모형 (3.3)과 비교한다. 평균은 $0.0002$, 분산은 $0.780$(모형 $0.783$), 3차 누율은 $+0.168$(모형 $+0.179$)이다. $\eta=0.5,1,2,3$에서 경험적 특성함수와 (3.4)의 차이는 $0.005$ 이하이다.

    python -B scripts/check_s1_distribution.py

<a id="section-4"></a>

## 4. Ng 공식과 쌍대의 대응

<a id="section-4-1"></a>

### 4.1. 구조의 대응

| 항목 | Ng (§1) | 정리 1 ($\sigma>\tfrac12$, §2) | 정리 2 ($\sigma=\tfrac12$, §3) |
|:--|:--|:--|:--|
| 변수 | $y=\log x$ | 대칭 약수 주파수 $\alpha$ | 높이 $t$ |
| 대상 | $e^{-y/2}M(e^y)$ | $2z\log\lvert\zeta(\sigma+i\alpha)\rvert$ | $\pi S_1(t)=\mathscr I(t)$ |
| 생성함수의 곱 | $1/(s\zeta(s))$의 Hadamard 곱 | $\zeta(\sigma-i\alpha)^z\zeta(\sigma+i\alpha)^z$의 Euler 곱 | $\log\zeta$의 Euler 곱을 $\sigma$로 적분 |
| 급수 | $\sum_\gamma a_\gamma e^{i\gamma y}$ | $-2z\sum_p\log\lvert1-p^{-\sigma-i\alpha}\rvert$ | $\Re\sum_n\Lambda(n)n^{-1/2-it}/\log^2n$ |
| 주파수 | 영점의 높이 $\gamma$ | $\log p$ | $\log p$ |
| 진폭 | $a_\gamma=1/(\rho\zeta'(\rho))$ | $2zp^{-\sigma}$ | $1/(\sqrt p\log p)$ |
| 위상의 독립성 | LI (추측) | FTA | FTA |
| 진폭의 제곱합 | (1.1) 류 (추측) | $\sum_pp^{-2\sigma}<\infty$ | $\sum_p(p\log^2p)^{-1}<\infty$ |
| 영점의 위치 | RH | 필요 없음 | 필요 없음 |
| 국소 구조 | 단순영점 필요 | 소수 거듭제곱은 ${}_2F_1$의 고차항 | 소수 제곱 항이 양의 왜도를 만듦 |
| 절단 | 높이 $T_\nu$와 Gallagher 부등식 | $\zeta Q_X$의 평균제곱 | $\zeta Q_X$의 평균제곱과 Littlewood 보조정리의 양성 |
| 소수별·영점별 변수 | $2\lvert a_\gamma\rvert\cos\theta_\gamma$ | $-2z\log\lvert1-p^{-\sigma}e^{i\theta_p}\rvert$ | $\sum_k\cos(k\theta_p)/(k^2p^{k/2}\log p)$ |
| 인수 | $J_0(2\eta\lvert a_\gamma\rvert)$ (정확) | ${}_2F_1(iz\eta,iz\eta;1;p^{-2\sigma})$ | $J_0(\eta/(\sqrt p\log p))+O(p^{-2})$ |
| 곱의 닫힌 꼴 | 알려지지 않음 | 중심 원자 합 $\sum_dd_{iz\eta}(d)^2d^{-2\sigma}$ | 알려지지 않음 |
| 아크사인의 출처 | 위상 $\gamma y$의 균등분포 | 약수 동전의 위상 $\alpha\log p$ | 위상 $t\log p$ |
| 약수 좌표의 핵 | 영점별 $e^{i\xi/2}J_0(\xi/2)$ ([전역 공식 (1.18)](global-formula.md#equation-1-18)) | 소수별 $\cos\theta_p$ (2.2) | (2.2)를 $\sigma$로 적분한 것 |
| 분포의 대칭 | 대칭 (LI 아래) | $z$와 같은 부호의 왜도 | 양의 왜도 |

정리 1과 정리 2의 변수 $\alpha$와 $t$는 같은 역할을 한다. (2.1)에서 대칭 비틀기 $(d/m)^{i\alpha}$는 높이 $\alpha$에서의 $|\zeta|$를 주기 때문이다. $z=-\tfrac12$는 $\mu=d_{-1/2}*d_{-1/2}$의 핵, $z=\tfrac12$는 균등 약수의 아크사인 핵에 해당한다.

<a id="section-4-2"></a>

### 4.2. Ng의 결과와 쌍대 결과

| Ng (RH와 $J_{-1}$ 류 가정) | 쌍대 ($S_1$) | 쌍대 쪽 근거 |
|:--|:--|:--|
| $M(x)\ll x^{1/2}(\log x)^{3/2}$ | $S_1(t)\ll\log t$, RH 아래 $\ll\log t/(\log\log t)^2$ | Titchmarsh §9.9, Littlewood |
| 약한 Mertens 가설 $\int_2^X(M(x)/x)^2dx\ll\log X$ | $\int_0^TS_1(t)^2dt\ll T$ | Littlewood, 점근식은 Selberg(무조건) |
| $\int_1^xM(u)u^{-3/2}du=o(\log x)$ | $\int_0^TS_1(t)\,dt=o(T)$ | 정리 2와 2차 모멘트 상계 |
| 극한분포의 존재 (명제 1.1) | 극한분포의 존재 | 정리 2 (무조건) |
| LI 아래 $\prod_\gamma J_0$ | LI 없이 소수별 곱 (3.4) | 정리 2 |
| 밀도의 수치 계산에 영점과 $\zeta'(\rho)$ 필요 | 소수만으로 계산 | (3.3) |
| 대칭인 극한분포 | 양의 왜도 | 정리 2 |
| 꼬리의 크기와 $\Omega$ 추측 | $S_1$의 $\Omega_\pm$ 결과 | Tsang |

왼쪽 열의 $M(x)$에 대한 결과는 Ng(2004)가 RH와 $J_{-1}(T)\ll T$에서 얻은 것이고, 극한분포는 §1에서 조건 (1.1)로 약화했다. 오른쪽 열은 Littlewood의 점별 상계 개선을 빼면 모두 무조건이다. 평균에 관한 행은 정리 2의 분포 수렴에 Selberg의 2차 모멘트 상계(균등 적분 가능성)를 더한 것이다.

<a id="section-4-3"></a>

### 4.3. 임계선에서의 퇴화와 누적

Ng의 정규화 $x^{-1/2}$는 쌍대 쪽의 $\sigma=\tfrac12$에 대응한다. $\sigma=\tfrac12$에서는 $\sum_pp^{-1}=\infty$이므로, $\eta\ne0$에서 (2.4)의 곱이 $0$으로 가고 퇴화하지 않은 극한분포가 없다. 이 경우에는 Selberg의 중심극한정리에 따라 $\log|\zeta(\tfrac12+i\alpha)|/\sqrt{\tfrac12\log\log A}$가 표준정규분포로 수렴하며, 이는 RH를 가정해도 바뀌지 않는다. 퇴화의 원인은 영점이 아니라 소수 진폭의 제곱합이 발산한다는 사실이다.

누적하면 퇴화가 사라진다. 두 쪽을 미분 단계와 누적 단계로 나누면 다음과 같다.

| 단계 | Ng 쪽 | 쌍대 쪽 |
|:--|:--|:--|
| 미분 | 측도 $\sum\mu(n)\delta_n$ (함수로서의 극한분포 없음) | $S(t)$, $\log\lvert\zeta(\tfrac12+it)\rvert$ (가우스, Selberg) |
| 누적 | $M(x)$: 영점별 아크사인 곱 (조건부) | $S_1(t)$: 소수별 아크사인 곱 (무조건) |

누적은 두 쪽 모두에서 진폭을 주파수로 나눈다(§3 도입부). 그 결과 Ng 쪽에는 진폭의 제곱합 (1.1)과 RH, LI가 가정으로 남지만, 쌍대 쪽에서는 각각 $\sum_p(p\log^2p)^{-1}<\infty$, Littlewood 보조정리의 양성, FTA가 그 자리를 채운다.

<a id="section-5"></a>

## 5. Mertens 함수와의 연결

§2–3의 쌍대 정리는 $\zeta$에 관한 진술이다. 이 절에서는 그것을 Mertens 함수 $M(x)$에 직접 연결한다. Ng 쪽과 쌍대 쪽은 같은 함수 $M$을 시간 $y=\log x$와 주파수 $t$에서 본 한 쌍이며, 두 쪽은 $\sigma\downarrow\tfrac12$라는 하나의 극한에서 만난다.

<a id="section-5-1"></a>

### 5.1. 푸리에 쌍

$y\ge0$에서 $\varphi_\sigma(y)=e^{-\sigma y}M(e^y)$로 두면, $\varphi_{1/2}$는 §1의 $\varphi$, 즉 Ng의 함수이다. $s=\sigma+it$에 대해

<a id="equation-5-1"></a>

$$
\widehat{\varphi_\sigma}(t):=\int_0^\infty\varphi_\sigma(y)e^{-ity}\,dy=\int_1^\infty M(x)x^{-s-1}\,dx=\frac1{s\,\zeta(s)}
\tag{5.1}
$$

이다. $\sigma>1$에서는 Abel 합산 $s\int_1^\infty M(x)x^{-s-1}dx=\sum_n\mu(n)n^{-s}$으로 성립한다. RH 아래에서는 $M(x)\ll x^{1/2+\varepsilon}$이므로 적분이 $\sigma>\tfrac12$에서 절대수렴하고, 해석적 연장으로 같은 식이 성립한다. RH 없이 $\sigma\le1$에서 $\widehat{\varphi_\sigma}$라고 쓰면 $1/(s\zeta(s))$의 연장을 뜻한다.

이로써 쌍대 정리들은 $M$에 대한 진술이 된다.

<a id="equation-5-2"></a>

$$
\mathscr L_\sigma^{-1/2}(t)=-\log|\zeta(\sigma+it)|=\log\big|(\sigma+it)\,\widehat{\varphi_\sigma}(t)\big|,
\qquad
\pi S_1(t)=-\int_{1/2}^\infty\log\big|(\sigma+it)\,\widehat{\varphi_\sigma}(t)\big|\,d\sigma.
\tag{5.2}
$$

Ng는 $\varphi_{1/2}$의 값을 시간 $y$에서 보고, 정리 1과 정리 2는 $\widehat{\varphi_\sigma}$의 로그 크기를 주파수 $t$에서 본다. 이 절의 식 (5.1)–(5.2)는 $1/\zeta=s\int Mx^{-s-1}dx$를 다시 쓴 것으로, 새로운 내용은 다음 절들에 있다.

**중간선의 양성.** [전역 공식 (2.1)](global-formula.md#equation-2-1)의 계수 $b_\alpha(n)$으로 비튼 Mertens 함수 $M_\alpha(x)=\sum_{n\le x}b_\alpha(n)$을 두자. 이는 전역 공식의 $\mathfrak M_x(\xi)$를 $\alpha=\xi/\log x$에서 본 것이며, $M_0=M$이다. 그 생성함수 $F_\alpha(s)=\zeta(s-i\alpha)^{-1/2}\zeta(s)^{-1/2}$는 두 분기점 $1$, $1+i\alpha$의 중간선 $\Im s=\alpha/2$에서

$$
s\int_1^\infty M_\alpha(x)x^{-s-1}dx\Big|_{s=\sigma+i\alpha/2}=F_\alpha(\sigma+\tfrac{i\alpha}2)=\big|\zeta(\sigma+\tfrac{i\alpha}2)\big|^{-1}>0
\qquad(\sigma>1)
$$

이다. 이는 (2.1)을 $z=-\tfrac12$에서 높이 $\alpha/2$로 옮긴 것이다. 따라서 정리 1은 비튼 Mertens 함수족의 중간선 Mellin 변환에 대한 분포 정리로도 읽힌다.

<a id="section-5-2"></a>

### 5.2. 국소 극한 공식

**명제 5.1 (국소 극한).** $T>0$이 영점의 높이가 아니라 하자. 그러면

<a id="equation-5-3"></a>

$$
\boxed{
\lim_{\sigma\downarrow\frac12}\frac{2\sigma-1}{2\pi}\int_{-T}^{T}\frac{dt}{|(\sigma+it)\zeta(\sigma+it)|^2}
=\sum_{\substack{\Re\rho=1/2\\|\gamma|<T}}\frac1{|\rho\zeta'(\rho)|^2}
}
\tag{5.3}
$$

이다. 오른쪽 합은 임계선 위의 영점에 걸치며, 그중 중복 영점이 있으면 $+\infty$로 읽는다. 임계선 밖의 영점은 기여하지 않는다. RH를 쓰지 않는다.

**증명.** $|\gamma|<T+1$인 영점은 유한개이다. 임계선 밖 영점의 실수부 중 최솟값을 $\tfrac12+\eta_0$라 하고 $\tfrac12<\sigma<\tfrac12+\eta_0/2$로 둔다. 임계선 위의 각 영점 $\rho=\tfrac12+i\gamma$ 주위에 반지름 $r$의 서로소인 구간을 잡으면, 그 밖에서 피적분함수는 $\sigma$에 균등하게 유계이다. 따라서 그 부분의 기여는 $O(2\sigma-1)\to0$이다.

중복도 $m$인 영점 근처에서는 $1/\zeta(s)=c_\rho(s-\rho)^{-m}(1+O(|s-\rho|))$, $c_\rho=m!/\zeta^{(m)}(\rho)$이다. $w=\sigma-\tfrac12$로 두면

$$
\frac{2\sigma-1}{2\pi}\int_{|t-\gamma|<r}\frac{dt}{|s\zeta(s)|^2}
=\frac{w}\pi\cdot\frac{|c_\rho|^2}{|\rho|^2}\int_{-r}^{r}\frac{(1+O(r))\,du}{(w^2+u^2)^m}
$$

이다. $m=1$이면 $\int_{-r}^r du/(w^2+u^2)=(2/w)\arctan(r/w)\to\pi/w$이므로 기여는 $(1+O(r))/|\rho\zeta'(\rho)|^2$로 간다. $m\ge2$이면 적분이 $\asymp w^{1-2m}$이라 기여가 $\infty$로 간다. 마지막으로 $r\to0$으로 보낸다. $\square$

각 영점은 직선 $\Re s=\sigma$ 위에서 $|1/\zeta|^2$에 폭 $\sigma-\tfrac12$인 Lorentz 봉우리를 만들고, 그 넓이가 정확히 Ng의 진폭 제곱 $|a_\gamma|^2$이다.

<a id="section-5-3"></a>

### 5.3. Plancherel과 Ng의 분산

RH 아래 $\varphi_\sigma\in L^2(0,\infty)$ ($\sigma>\tfrac12$)이므로 (5.1)과 Plancherel 정리에서

<a id="equation-5-4"></a>

$$
\int_0^\infty\varphi_\sigma(y)^2\,dy=\int_1^\infty\frac{M(x)^2}{x^{2\sigma+1}}\,dx
=\frac1{2\pi}\int_{\mathbb R}\frac{dt}{|(\sigma+it)\zeta(\sigma+it)|^2}
\tag{5.4}
$$

이다. 왼쪽은 Mertens 함수의 평균제곱이고, 오른쪽은 쌍대 쪽 양 $e^{-2\log|\zeta|}$에 $|s|^{-2}$의 가중치를 준 적분이다. $\varphi_\sigma^2=e^{-(2\sigma-1)y}\varphi_{1/2}^2$이므로 왼쪽에 $2\sigma-1$을 곱한 것은 $\varphi_{1/2}^2$의 Abel 평균이다.

**명제 5.2 (Ng의 분산).**

1. RH 아래
$$
\sum_\gamma\frac1{|\rho\zeta'(\rho)|^2}\le\liminf_{\sigma\downarrow\frac12}\,(2\sigma-1)\int_0^\infty\varphi_\sigma(y)^2\,dy
$$
이다. 오른쪽이 유한하면 모든 영점은 단순하다. 특히 약한 Mertens 가설 $\int_0^Y\varphi_{1/2}^2\,dy\ll Y$는 단순영점과 $\sum_\gamma|a_\gamma|^2<\infty$를 함의한다.
2. RH와 (1.1) 아래 극한이 존재하고
$$
\sum_\gamma\frac1{|\rho\zeta'(\rho)|^2}
=\lim_{Y\to\infty}\frac1Y\int_0^Y\varphi_{1/2}^2\,dy
=\lim_{\sigma\downarrow\frac12}\frac{2\sigma-1}{2\pi}\int_{\mathbb R}\frac{dt}{|(\sigma+it)\zeta(\sigma+it)|^2}
$$
이다.

**증명.** (1) RH 아래 모든 영점이 임계선 위에 있다. 각 $T$에서 (5.3)의 왼쪽은 (5.4)의 오른쪽에 $2\sigma-1$을 곱한 것보다 작으므로, $\sigma\downarrow\tfrac12$ 다음 $T\to\infty$를 취하면 부등식이 나온다. 약한 Mertens 가설이 $\int_0^Y\varphi_{1/2}^2\le CY$이면, 부분적분으로 $(2\sigma-1)\int_0^\infty e^{-(2\sigma-1)y}\varphi_{1/2}^2\,dy=(2\sigma-1)^2\int_0^\infty e^{-(2\sigma-1)Y}\big(\int_0^Y\varphi_{1/2}^2\big)dY\le C$이다.

(2) 명제 1.1의 증명에서 $\varphi_{1/2}$는 Bohr–Fourier 계수가 $a_\gamma$인 $B^2$-근사주기 함수이다. $B^2$의 Parseval 등식으로 평균 $\frac1Y\int_0^Y\varphi_{1/2}^2$은 $\sum|a_\gamma|^2$로 수렴한다. Cesàro 평균이 수렴하면 Abel 평균도 같은 값으로 수렴하고, (5.4)로 마지막 등식을 얻는다. $\square$

(1)의 둘째 주장은 고전적인 사실(Titchmarsh §14.29)과 같으며, 여기서는 (5.3)에서 바로 나온다. 명제 5.1은 $T$를 고정하고 $\sigma\downarrow\tfrac12$로 보낸 것이고, 명제 5.2는 $T=\infty$에서 같은 극한을 취한 것이다. 일반적으로는 (1)의 부등식만 성립하며, **Ng의 가정 (1.1)이 하는 일은 두 극한의 순서를 바꿔도 되게 하는 것**이다. 즉 $J_{-1}$ 문제는 $|1/(s\zeta(s))|^2$의 질량이 $\sigma\downarrow\tfrac12$에서 영점 봉우리 밖으로 새어 나가지 않느냐는 문제이다.

<a id="section-5-4"></a>

### 5.4. 쌍대 평균과 Gonek 상수 (휴리스틱)

RH 아래 고정된 $\sigma>\tfrac12$에서 $|\zeta(\sigma+it)|^{-2}$의 $t$-평균은

<a id="equation-5-5"></a>

$$
\lim_{T\to\infty}\frac1T\int_T^{2T}\frac{dt}{|\zeta(\sigma+it)|^2}=\sum_n\frac{\mu(n)^2}{n^{2\sigma}}=\frac{\zeta(2\sigma)}{\zeta(4\sigma)}
\tag{5.5}
$$

이다(Titchmarsh 14장의 평균값 정리). $d_{-1}=\mu$이므로 오른쪽은 이산핵 $\mathscr K_{d^2}^{-1,-1}$의 중심 원자 합, 즉 (2.4)의 오른쪽을 $a=-1$에서 형식적으로 읽은 것이다. $\sigma\downarrow\tfrac12$이면 (5.5)는 $\sim1/(\zeta(2)(2\sigma-1))$로 커진다.

이 질량을 명제 5.1처럼 영점 봉우리가 모두 나른다고 가정하자. 높이 $t$에서 영점의 밀도는 $(2\pi)^{-1}\log(t/2\pi)$이고 각 봉우리의 넓이는 $\pi/(|\zeta'(\rho)|^2(\sigma-\tfrac12))$이다. 그러면 높이 $t$ 근처의 평균 $|\zeta'(\rho)|^{-2}$는 $6/(\pi^2\log t)$이 되고, 이를 적분하면

<a id="equation-5-6"></a>

$$
J_{-1}(T)=\sum_{0<\gamma\le T}\frac1{|\zeta'(\rho)|^2}\sim\frac{3}{\pi^3}\,T
\tag{5.6}
$$

를 얻는다. 이는 Gonek의 추측과 같은 상수이다. 이 계산은 $\sigma\downarrow\tfrac12$와 $T\to\infty$의 순서를 바꾼 것으로, §5.3에서 Ng의 가정이 하는 일과 정확히 같다. 따라서 증명이 아니라, 쌍대 쪽의 가정 없는 평균 (5.5)가 Gonek 상수를 가리킨다는 관찰이다.

**[수치]** [검산 스크립트](../scripts/check_mertens_bridge.py)는 다음을 확인한다.

- (5.3), $T=300$: $\sigma=0.505,0.51,0.52$의 값 $0.04083$, $0.05282$, $0.07564$를 $2\sigma-1$의 이차식으로 외삽하면 $0.028463$이다. 오른쪽 합 $0.028441$과 $2\times10^{-5}$ 안에서 일치한다.
- (5.4), $\sigma=\tfrac34$: $10^7$까지 체로 구한 $M(x)$에서 왼쪽은 $0.71866$, $|t|\le3000$의 수치 적분과 꼬리 보정에서 오른쪽은 $0.71868$이다.
- (5.5), $\sigma=\tfrac34$: $t\in[10,3010]$의 평균은 $2.1774$, $\zeta(\tfrac32)/\zeta(3)=2.1733$이다.
- (5.6): $J_{-1}(T)/(3T/\pi^3)$는 $T=100,200,300$에서 $1.115$, $1.057$, $1.022$이다.

이 범위에서는 RH가 수치로 확인되어 있으므로, RH를 쓰는 식도 그 범위에서 검산한 것이다.

    python -B scripts/check_mertens_bridge.py

<a id="checks"></a>

## 검산

네 스크립트는 numpy와 scipy를 사용하며, 제타함수와 영점을 배정밀도로 직접 계산한다.

- [Ng 쪽 검산](../scripts/check_ng_tail.py): 보조정리 1.3의 평균제곱 부등식과 (1.1)의 부분합 (§1).
- [쌍대 곱 검산](../scripts/check_dual_product.py): 대칭 비틀기 항등식 (2.1), 중심 원자 합 (2.4), 평균제곱 (2.5), $-\log|\zeta(\tfrac34+i\alpha)|$의 경험적 특성함수 (§2). $\alpha\in[0,3000]$에서 경험적 특성함수와 (2.4)의 차이는 $0.01$ 안팎이며, 선도항 $\prod_pJ_0$만 쓰면 차이가 더 커서 소수 거듭제곱 보정이 보인다.
- [Mertens 연결 검산](../scripts/check_mertens_bridge.py): 국소 극한 (5.3), Plancherel (5.4), 쌍대 평균 (5.5), Gonek 비율 (5.6) (§5).
- [누적량 검산](../scripts/check_s1_distribution.py): Littlewood 공식 (3.2), 근사 (3.6)의 오차 감소, $\mathscr I(t)$의 분포·왜도 (§3).

이 수치 확인은 증명을 대체하지 않으며, 낮은 높이의 영점과 유한한 구간에서의 계산이다.

    python -B scripts/check_ng_tail.py
    python -B scripts/check_dual_product.py
    python -B scripts/check_s1_distribution.py
    python -B scripts/check_mertens_bridge.py

<a id="references"></a>

## 참고 문헌

- N. Ng, [The distribution of the summatory function of the Möbius function](https://arxiv.org/abs/math/0310381), Proc. London Math. Soc. 89 (2004): RH와 $J_{-1}(T)\ll T$ 아래의 극한분포, LI 아래의 곱 공식, 약한 Mertens 가설.
- A. Akbary, N. Ng, M. Shahabi, [Limiting distributions of the classical error terms of prime number theory](https://arxiv.org/abs/1306.1657), Q. J. Math. 65 (2014): 정리 2.9의 $B^2$-근사주기 함수의 극한분포, $J_{-1}(T)\ll T^\theta$ ($\theta<3-\sqrt3$) 조건.
- N. Ng, [Prime Number Error Terms](https://arxiv.org/abs/2505.11295), arXiv:2505.11295 (2025): Akbary–Ng–Shahabi 따름정리 1.3(b)의 조건을 $\theta<2$로 개선.
- E. C. Titchmarsh, *The Theory of the Riemann Zeta-Function*, 2nd ed. (1986): 정리 4.11($\zeta$의 근사식), §9.4($\arg\zeta$의 상계), §9.9(Littlewood 보조정리와 $S_1$), 정리 14.13($S(t)$의 RH 상계), 정리 14.16(높이 $T_\nu$의 선택), §14.20–14.21(Selberg의 근사), §14.27($M(x)$의 명시적 공식).
- P. X. Gallagher, A large sieve density estimate near $\sigma=1$, Invent. Math. 11 (1970), 329–339: 보조정리 1(지수합의 평균제곱).
- H. L. Montgomery, R. C. Vaughan, Hilbert's inequality, J. London Math. Soc. (2) 8 (1974), 73–82: Dirichlet 다항식의 평균값 정리.
- H. Bohr, B. Jessen, Über die Werteverteilung der Riemannschen Zetafunktion, Acta Math. 54 (1930), 1–35; 58 (1932), 1–55: $\sigma>\tfrac12$에서 $\log\zeta(\sigma+it)$의 극한분포.
- Y. Lamzouri, [On the distribution of extreme values of zeta and L-functions in the strip 1/2<σ<1](https://arxiv.org/abs/1005.4640), IMRN (2011): 무작위 Euler 곱 모형 $\zeta(\sigma,\mathbb X)$.
- A. Selberg, Contributions to the theory of the Riemann zeta-function, Arch. Math. Naturvid. 48 (1946): $S(t)$, $S_1(t)$의 모멘트와 $\log\zeta$의 소수합 근사.
- K.-M. Tsang, Some $\Omega$-theorems for the Riemann zeta-function, Acta Arith. 46 (1986): $S_1(t)$의 $\Omega_\pm$ 결과.
- S. M. Gonek, On negative moments of the Riemann zeta-function, Mathematika 36 (1989), 71–88: $J_{-1}(T)$와 상수 $3/\pi^3$의 추측.
- M. Radziwiłł, K. Soundararajan, Selberg's central limit theorem for $\log|\zeta(1/2+it)|$, Enseign. Math. 63 (2017), 1–19: §4.3의 임계선 극한.
