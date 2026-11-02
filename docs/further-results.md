# 후속 결과: 뫼비우스 경우, 고차 전개와 제타함수 영점

**2026-11-01 (개정판)**

「일반화된 베타 핵과 가중 약수합」의 후속 결과와 부록 B.2–B.4 및 부록 C를 모았다. 2026-11-03에 주정리 문서에 이산핵을 다루는 §7을 추가하면서, 기존 §7–9를 이 문서의 §8–10으로 옮기고 정리·수식 번호와 참조를 함께 조정했다.

기본 기호 $A_x^{z,w}(F)$, $r=z+w$, $L=\log x$, $C(r)$와 핵 $\mathcal K_{z,w}$는 [주정리 문서의 §1](main-theorem.md#section-1)을 따른다. 오일러 곱과 핵의 성질은 같은 문서의 [§2](main-theorem.md#section-2), [§3](main-theorem.md#section-3), [§6](main-theorem.md#section-6)에 있다.

이 문서의 $A_x$와 $H$는 기존 제곱 인수 없는 모형에 대한 것이다. 모든 정수의 $d_z$ 가중치를 사용하는 $\widetilde A_x$와 $\widetilde H=1$은 [주정리 문서의 §7](main-theorem.md#section-7)에서 정의한다.

**[개요]**는 증명의 뼈대만 제시한 부분이며, **RH**는 리만 가설을 뜻한다. 각 결과에 붙은 가정을 그대로 유지했다.

- [§8. 뫼비우스 경우](#section-8): 명제 8.1과 관련 예제
- [§9. 베타 급수와 수렴 반경](#section-9): 명제 9.1, 9.2와 국소 기여에 관한 논의
- [§10. 영점 층과 복원](#section-10): 명제 10.1, 10.2와 복원 방법
- [부록 B.2–B.4. 수치 확인](#appendix-b)
- [부록 C. RH 아래 평균제곱 추정](#appendix-c)

부록 B.1의 수치 확인은 [주정리 문서](main-theorem.md#appendix-b-1)에 있다.

<a id="section-8"></a>

## 8. 뫼비우스 경우

<a id="section-8-1"></a>

### 8.1. 정수 지수

**명제 8.1.** $z,w\in\{0,-1,-2,\dots\}$, $F\in C^\infty([0,1])$이면 모든 $k$에 대해 $A_x^{z,w}(F)=O_{F,z,w,k}(xL^{-k})$이다. $F$가 다항식이면 $O(xe^{-c\sqrt L})$이다.

**증명.** 지수가 정수이므로 $\mathfrak D_{\xi,L}(s)=\zeta(s-i\alpha)^z\zeta(s)^wH$에는 분기점이 없고, 영점 없는 영역 $1-\delta\le\Re s\le c$, $|\Im s|\le T$에서 해석적이다. [(A.5)](main-theorem.md#equation-a-5)에서 $|\zeta|^{\pm1}\ll(\log\tau)^{O(1)}$이므로 [(A.6)](main-theorem.md#equation-a-6)도 그대로 성립한다. 따라서 Perron 경로를 우회 없이 $\Re s=1-\delta$로 옮길 수 있고, $|\xi|\le L^{1/4}$, $j=0,1$에서 $\partial_\xi^jS_x(\xi)\ll xe^{-c\sqrt L}$이다. [(5.2)](main-theorem.md#equation-5-2)에서 $|\xi|\le L^{1/8}$에 대해 $\Theta_x(\xi)\ll xe^{-c\sqrt L}(1+|\xi|)+|\xi|\sqrt xL^{k_0}$이다. [§5.3](main-theorem.md#section-5-3)의 고주파 꼬리는 $xL^{k_0-1-(N-1)/8}$이고 $N$은 임의로 크게 잡을 수 있다. $F=u^k$이면 $\sum_{n\le t}\sum_{d\mid n}(\text{가중치})(\log d)^k$의 디리클레 급수 $(-\partial_{s_1})^kD(s_1,s)|_{s_1=s}$도 같은 영역에서 해석적이므로, $(\log t)^{-k}$에 대한 부분합산으로 $xe^{-c\sqrt L}$을 얻는다. $\square$

<a id="section-8-2"></a>

### 8.2. $(-1,-1)$: 뫼비우스의 곱

제곱 인수 없는 $n=dm$에서 $d,m$은 서로소이므로 쌍 가중치는 $(-1)^{\omega(d)}(-1)^{\omega(m)}=\mu(d)\mu(m)=\mu(n)$이고

$$
A_x^{-1,-1}(F)=\sum_{2\le n\le x}\mu(n)\sum_{d\mid n}F\!\left(\frac{\log d}{\log n}\right).
$$

$F=1$이면 $\sum_{2\le n\le x}\mu(n)\tau(n)$이다. 주항은 서로 독립인 두 이유로 사라진다.

- **산술:** $C(-2)=0$이다. $p=2$의 국소 인수 $1+r/p=1-2/2$가 $0$이기 때문이다.
- **해석:** $\mathcal K_{-1,-1}=0$이다. 지수가 정수이기 때문이다.

$D(s,s)=\prod_p(1-2p^{-s})=\zeta(s)^{-2}H(s,s)$는 $s=1$에서 3중 영점을 가진다($\zeta^{-2}$에서 2, $p=2$에서 1).

서로소 조건을 풀면([§2.1](main-theorem.md#section-2-1)) 쌍 가중치는 $\mu(d)\mu(m)$이고, 디리클레 급수는 정확히 $1/(\zeta(s)\zeta(t))$, 즉 $H^\natural\equiv1$이다. 상수는 $C(-1)^2=1$이 되지만 핵은 여전히 $0$이다. 따라서 서로소 조건을 풀어도 명제 8.1의 소멸 추정은 유지된다.

<a id="section-8-3"></a>

### 8.3. $(-1,1)$: 한쪽만 뫼비우스

쌍 가중치는 $\mu(d)$이다. $C(0)=1$이고, $g_1(1-u)=\mathbf 1_{u<1}$이 $u=0$ 근처에서 $1$이므로 $\mathcal K_{-1,1}=\delta_0'$, $\mathcal K_{-2,1}=\delta_0''$, $\mathcal K_{-1,0}=0$이다. [(1.4)](main-theorem.md#equation-1-4)에서

<a id="equation-8-1"></a>

$$
A_x^{-1,1}(F)=-\frac xL\Big[F'(0)+\frac{F'(0)+\kappa(0)F''(0)}{L}\Big]+O\!\left(\frac x{L^3}\right),
\qquad \kappa(0)=\gamma+\sum_p\frac{\log p}{p(p-1)}=1.3325822\ldots
\tag{8.1}
$$

뫼비우스 반전은 이 차수까지 $F$의 원점에서의 1차, 2차 테일러 계수만 남긴다.

- **$F=u$:** $\sum_{d\mid n}\mu(d)\log d=-\Lambda(n)$이므로 좌변은 정확히 $-\pi(x)$이다. (8.1)은 $\pi(x)=x/L+x/L^2+O(x/L^3)$, 즉 $\mathrm{li}(x)$의 앞 두 항이다.
- **$F=u^2$:** 제곱 인수 없는 $n$에 대해 $\sum_{d\mid n}\mu(d)(\log d/\log n)^2$는 $\omega(n)\le2$일 때만 영이 아니고, $n=p$에서 $-1$, $n=pq$에서 $2ab$이다($a=\log p/\log n$, $b=\log q/\log n$). 따라서
<a id="equation-8-2"></a>

$$
\sum_{\substack{p<q\\pq\le x}}\frac{2\log p\log q}{\log^2(pq)}=\pi(x)-\frac{2\kappa(0)\,x}{L^2}+O\!\left(\frac x{L^3}\right).
\tag{8.2}
$$

**서로소 해제 변형.** 쌍 가중치 $\mu(d)\mu^2(m)$이면 [(2.6)](main-theorem.md#equation-2-6)에서 상수는 $C(-1)C(1)=6/\pi^2$이고, $H_{-1}\equiv1$이므로 $\eta(-1)=0$, $\kappa(-1)=\gamma$이다.

<a id="equation-8-3"></a>

$$
A_x^{\natural,-1,1}(F)=-\frac{6}{\pi^2}\frac xL\Big[F'(0)+\frac{F'(0)+\gamma F''(0)}{L}\Big]+O\!\left(\frac x{L^3}\right).
\tag{8.3}
$$

$6/\pi^2$은 $\mu*\mu^2$가 제곱수에만 놓인다는 사실($\sum(\mu*\mu^2)(n)n^{-s}=1/\zeta(2s)$)에서 온다. **[수치]** [B.1](main-theorem.md#appendix-b-1), B.2.

<a id="section-8-4"></a>

### 8.4. $(-\tfrac12,-\tfrac12)$

제곱 인수 없는 $n=dm$에서 쌍 가중치는 $\mu(n)/\tau(n)$이다. 따라서 $F=1,u$에 대한 합은 정확히 $M(x)-1$, $(M(x)-1)/2$이며, [§6](main-theorem.md#section-6) 표의 영인 주계수와 일치한다. 여기서 $M(x)=\sum_{n\le x}\mu(n)$이다.

<a id="section-9"></a>

## 9. 베타 급수의 수렴 반경과 $s=\tfrac12$

<a id="section-9-1"></a>

### 9.1. 베타 급수

[부록 A.3](main-theorem.md#appendix-a-3)의 국소 인수 $G$를 전개하여

<a id="equation-9-1"></a>

$$
G(s_1,s_2)=\sum_{a,b\ge0}G_{ab}(s_1-1)^a(s_2-1)^b
\tag{9.1}
$$

로 쓴다. $G_{00}=C$, $G_{10}=G_1$, $G_{01}=G_2$이다.

**명제 9.1.** 각 $J\ge1$에 대해, 같은 범위의 $\xi$에서

<a id="equation-9-2"></a>

$$
S_x(\xi)=xL^{r-1}\sum_{a+b<J}G_{ab}L^{-a-b}\mathcal B_{z-a,w-b}(\xi)+O\!\left(xL^{\beta-1-J}(1+|\xi|)^{B_J}\right).
\tag{9.2}
$$

$J=2$가 [(5.1)](main-theorem.md#equation-5-1)이고, 증명은 [부록 A.3](main-theorem.md#appendix-a-3)의 테일러 전개를 $J$차까지 늘린 것이다. 즉 베타핵 전개는 $G$의 $s=1$에서의 테일러 전개이다.

$j$차 핵들의 $k$차 모멘트를 함께 모으면 $\sum_{a+b=j}G_{ab}\langle\mathcal K_{z-a,w-b},u^k\rangle=m_{j,k}/\Gamma(r-j+k)$이고

<a id="equation-9-3"></a>

$$
m_{j,k}=\sum_{a+b=j}G_{ab}\,(z-a)_k.
\tag{9.3}
$$

$m_{j,0}$은 $G(1+t,1+t)$의 테일러 계수이다. $m_{j,1}$, $m_{j,2}$에는 대각선 위의 $\partial_1G$, $\partial_1^2G$가 들어간다.

<a id="section-9-2"></a>

### 9.2. 수렴 반경과 핵들 사이의 상쇄

[(2.7)](main-theorem.md#equation-2-7)에서 가장 가까운 실축 쪽 특이점 후보는 $s_1=\tfrac12$, $s_1+s_2=1$, $s_2=\tfrac12$이다. 소거되지 않는 후보는 해당 방향의 테일러 수렴 반경을 $\tfrac12$로 제한한다. 더 왼쪽에는 고차 소수 거듭제곱의 층과 $1/s_2$의 극이 있다.

대각선에서는 [(2.8)](main-theorem.md#equation-2-8)에 의해 세 특이점이 $\zeta(2s)^{-r(r+1)/2}$로 합쳐진다.

- **$r(r+1)/2$가 0 이상의 정수인 경우:** 합친 결과가 $s=\tfrac12$에서 해석적이므로 $m_{j,0}$의 반경은 다음 특이점까지 늘어난다. $r=0,\pm1$에서는 $s=0$의 극까지 거리 $1$이다. 즉 같은 차수의 핵들이 총질량에서 서로 지워진다.
- **비상수 모멘트:** $\partial_1^2G$ 등은 대각선에서도 $s=\tfrac12$에서 특이하다. 그래서 상쇄가 불완전하다.

**[개요]** $r\notin\mathbb Z$에서 $1/\Gamma(r-j)$의 계승 성장은, 반경 $\tfrac12$인 모멘트 급수의 최소항이 $j\approx L/2$에 있고 크기가 대략 $x^{1/2}$임을 시사한다. (9.2)는 고정된 $J$의 추정이므로, 이를 최적 절단 오차로 주장하려면 상수와 나머지의 $J$에 대한 균등 추정을 별도로 증명해야 한다.

<a id="section-9-3"></a>

### 9.3. 아크사인의 세 점 핵

$z=w=\tfrac12$이면 $e_1=e_2=\tfrac38$, $zw=\tfrac14$이고, 지수 합이 $1$이다. 대각선에서 $H(s,s)=1/\zeta(2s)$는 $s=\tfrac12$에서 오히려 영점이므로 $m_{j,0}$에는 $s=\tfrac12$이 보이지 않는다. $m_{j,1}$도 $d\leftrightarrow m$ 대칭 때문에 보이지 않는다. 반면 $m_{j,2}$의 생성함수에는 대각선에서 다음 단순극이 생긴다.

$$
\partial_1^2\big[(2s_1-1)^{3/8}(s_1+s_2-1)^{1/4}\big](2s_2-1)^{3/8}\Big|_{s_1=s_2=s}=-\frac34(2s-1)^{-1}.
$$

따라서

<a id="equation-9-4"></a>

$$
m_{j,2}\sim c\,(-2)^j,\qquad c=-\frac3{16}\cdot2Z(\tfrac12)=-\frac{3}{16}\bigl(-\zeta(\tfrac12)\bigr)=-0.27382\ldots
\tag{9.4}
$$

**[수치]** B.3.

비틀림이 있으면 $s-\tfrac12=q/L$에서 $H\approx\frac2L(q-i\xi)^{3/8}(q-\tfrac{i\xi}2)^{1/4}q^{3/8}$이다($\widetilde H(\tfrac12,\tfrac12)=1$). $s=1$에서 뻗은 두 절단을 위아래로 감싸는 Hankel 고리 위에서는 $\zeta(s_1)^{1/2}\zeta(s_2)^{1/2}\approx\zeta(\tfrac12)$이다. 세 점 핵을

<a id="equation-9-5"></a>

$$
\langle\mathcal T,e^{i\xi u}\rangle=\frac1{2\pi i}\int_{\mathcal U_Y}e^q(q-i\xi)^{3/8}\Big(q-\frac{i\xi}2\Big)^{1/4}q^{3/8}\,dq
\tag{9.5}
$$

로 정의하면 $\mathcal T$는 $[0,1]$에 지지된 분포이고

$$
\langle\mathcal T,1\rangle=0,\qquad\langle\mathcal T,u\rangle=0,\qquad\langle\mathcal T,u^2\rangle=-\frac3{16}.
$$

질량은 $q$의 총 차수가 $1$이라서, 1차 모멘트는 상수항의 Hankel 적분이라서 $0$이다. 2차 모멘트는 $(1-t)^{3/8}(1-t/2)^{1/4}$의 $t^2$ 계수 $-3/32$에서 온다.

**[개요]** 고정된 $\xi$에서 $s=\tfrac12$ 근방의 국소 기여는 $4\zeta(\tfrac12)\,x^{1/2}L^{-2}\langle\mathcal T,e^{i\xi u}\rangle+O_\xi(x^{1/2}L^{-3})$이다. $F=u^2$(S 수준)이면 $-\tfrac34\zeta(\tfrac12)x^{1/2}L^{-2}$이고, 이는 §10.2의 정확한 계산과 일치한다. 다음 층은 소수 세제곱 항의 $s=\tfrac13$이다.

$r\in\mathbb Z$이면 $1$과 $1+i\alpha$의 절단을 두 점 사이의 선분으로 묶을 수 있다. 두 점을 함께 도는 모노드로미가 $e^{-2\pi ir}=1$이므로 $s=\tfrac12$ 근방의 가지가 단일값이고, 그 극한은 $\zeta(\tfrac12)^r$이다. 이는 $r=-1$에도 적용된다. 일반 정수 $r$의 세 점 핵은 (9.5)의 지수들을 $e_1,zw,e_2$로 바꾸어 얻으며, 총차수는 $r(r+1)/2$이다. 비정수 $r$에서는 절단의 도약을 따로 처리해야 한다.

<a id="section-9-4"></a>

### 9.4. 임계띠와의 관계

실축 근처에서는 무조건적으로 $\Re s=\tfrac12$까지 연장된다. $\zeta$는 $(0,1)$ 위에서, 그리고 $|t|<14$에서 영이 아니기 때문이다. 세 점 핵을 더하면 $\tfrac13$까지 간다. 베타 급수는 실축 근처의 특이점만 보므로, 전체 오차를 막는 $|\Im s|\ge14.13$ 쪽 경로는 핵을 어떻게 재배열해도 건드릴 수 없다.

**명제 9.2.** 모든 $\varepsilon>0$에 대해 $\sum_{n\le x}(\mu*\mu)(n)\ll x^{\theta+\varepsilon}$이면 $\sigma>\theta$에서 $\zeta(s)\ne0$이다. 역으로 RH이면 $\sum_{n\le x}(\mu*\mu)(n)\ll x^{1/2+\varepsilon}$이다.

**증명.** 부분합산으로 $1/\zeta(s)^2=\sum(\mu*\mu)(n)n^{-s}$가 $\sigma>\theta$에서 수렴하고 해석적이다. 그 영역의 영점은 극이 되므로 있을 수 없다. 역은 RH 아래 $\sigma\ge\tfrac12+\varepsilon$에서 $1/\zeta(s)\ll|t|^\varepsilon$인 것과 Perron 공식으로 따른다. $\square$

따라서 뫼비우스 경우의 거듭제곱 오차 추정은 영점 없는 영역과 직결된다.

<a id="section-10"></a>

## 10. 영점 층과 복원

<a id="section-10-1"></a>

### 10.1. 반사된 핵

**[개요]** 중복도 $m$인 비자명 영점 $\rho=\beta_\rho+i\gamma$에서 $H$가 $(\rho,\rho)$ 근방에 정칙이라고 하자. 이 정칙성은 $\beta_\rho\ge\tfrac12$이면 [(2.7)](main-theorem.md#equation-2-7)로 확보된다. $b_\rho=\zeta^{(m)}(\rho)/m!$로 두고 가지를 고정하면 $\zeta(s_1)^z\zeta(s_2)^w\approx b_\rho^r(s_1-\rho)^{mz}(s_2-\rho)^{mw}$이다. $s-\rho=q/L$로 옮긴 국소 기여는 고정된 $\rho,\xi$에서

<a id="equation-10-1"></a>

$$
x^\rho L^{-mr-1}\,\frac{b_\rho^rH(\rho,\rho)}{\rho}\,\mathcal B_{-mz,-mw}(\xi)
+O_{\rho,z,w,\xi}\!\left(x^{\beta_\rho}L^{-m\Re r-2}\right)
\tag{10.1}
$$

이다. 반사된 핵은 $\mathcal K_{-mz,-mw}$이며, $m=1$이면 $\mathcal K_{-z,-w}$이다. 이 국소식에는 RH가 필요 없다. RH 아래에서는 $x$의 거듭제곱이 $x^{1/2}$이고, 로그 인자를 포함한 기준 규모는 $x^{1/2}L^{-m\Re r-1}$이다.

<a id="section-10-2"></a>

### 10.2. 아크사인: 2차 모멘트부터 영점이 보인다

$z=w=\tfrac12$이고 $\Re\rho\ge\tfrac12$인 단순 영점에서는 $H(\rho,\rho)=1/\zeta(2\rho)$이고 핵은 $\mathcal K_{-1/2,-1/2}=h_2'$이다. [§6](main-theorem.md#section-6)의 표에서 $\langle h_2',1\rangle=\langle h_2',u\rangle=0$, $\langle h_2',u^2\rangle=-\tfrac14$이다. 따라서 이 영점의 기여는 2차 모멘트부터 나타난다. 실제로 제곱 인수 없는 수의 개수의 급수 $\zeta(s)/\zeta(2s)$는 $\rho$에서 특이점이 없다.

이 핵은 아크사인의 이차 보정 [(1.5)](main-theorem.md#equation-1-5)에도 나타난다. 대칭 핵에서
$z(\mathcal K_{z-1,z}+\mathcal K_{z,z-1})=\frac{z}{z-1}\mathcal K_{z-1,z-1}$ ($z\ne1$)이므로,
보정 핵과 반사된 핵의 지수가 일치하는 경우는 $z-1=-z$, 즉 $z=\tfrac12$이다.

**정확한 확인.** $D$를 $n$의 균등한 무작위 약수라 하면 $\log D=\sum_{p\mid n}\varepsilon_p\log p$이고, $\varepsilon_p$는 독립인 공정한 동전이다. 따라서

<a id="equation-10-2"></a>

$$
2^{-\omega(n)}\sum_{d\mid n}(\log d)^2=\frac14\Big(\log^2n+\sum_{p\mid n}\log^2p\Big).
\tag{10.2}
$$

$T_2(x)=\sum_{n\le x}\mu^2(n)\sum_{p\mid n}\log^2p$의 디리클레 급수는 $\dfrac{\zeta(s)}{\zeta(2s)}\sum_p\dfrac{\log^2p\,p^{-s}}{1+p^{-s}}$이다.

- **영점에서:** $\sum_p\log^2p\,p^{-s}=(\zeta'/\zeta)'(s)-\sum_{k\ge2}k\sum_p\log^2p\,p^{-ks}\approx-(s-\rho)^{-2}$이므로, 단순극과 유수 $-\zeta'(\rho)/\zeta(2\rho)$가 생긴다. 이를 $4L^2$로 나눈 S 수준 $u^2$ 통계의 영점 항은 (10.1)에 $\langle h_2',u^2\rangle=-\tfrac14$를 넣은 것과 정확히 같다.
- **$s=\tfrac12$에서:** $-3\sum_p\log^2p\,p^{-2s}\approx-\tfrac34(s-\tfrac12)^{-2}$이고 $\zeta(s)/\zeta(2s)\approx2\zeta(\tfrac12)(s-\tfrac12)$이다. 그래서 진동 없는 편향 $-3\zeta(\tfrac12)\sqrt x=4.381\sqrt x$가 생긴다. S 수준 $u^2$로는 $-\tfrac34\zeta(\tfrac12)x^{1/2}L^{-2}=4\zeta(\tfrac12)\langle\mathcal T,u^2\rangle x^{1/2}L^{-2}$로, (9.5)의 세 점 핵과 일치한다.
- **더 낮은 층:** $s=\tfrac13$의 $x^{1/3}(\log)$ 항, 그리고 $\zeta(2s)$ 영점에서 오는 진동수 $\gamma/2$의 $x^{1/4}(\log x)^2$ 항이 있다.

**보조 통계.** $c_n=\sum_{d\mid n}\Lambda(d)\log d$로 두면 제곱 인수 없는 $n$에서 $c_n=\sum_{p\mid n}\log^2p$, 즉 $\log D$ 분산의 4배이다. 디리클레 급수 $Q(s)=\zeta(s)(\zeta'/\zeta)'(s)$에는 $s=\tfrac12$ 항도 $\zeta(2s)$ 층도 없다. 절단된 명시적 공식을

<a id="equation-10-3"></a>

$$
\sum_{n\le x}c_n=xP_2(\log x)-2\,\mathrm{Re}\sum_{0<\gamma\le T}\frac{\zeta'(\rho)}{\rho}\,x^\rho+R_2(x,T)
\tag{10.3}
$$

로 쓴다. $xP_2(\log x)=\operatorname{Res}_{s=1}Q(s)x^s/s$이고, 합은 서로 다른 비자명 영점에 걸친다. 중복도 $m\ge2$이면 $Q(s)=O((s-\rho)^{m-2})$로 정칙이고 $\zeta'(\rho)=0$이므로 그 항은 사라진다. $R_2$의 RH 아래 평균제곱 추정은 부록 C에서 증명한다. **[수치]** B.4.

<a id="section-10-3"></a>

### 10.3. 모멘트의 차수와 영점의 중복도

공정한 동전의 누적량은 $\kappa_2=\tfrac14$, $\kappa_3=0$, $\kappa_4=-\tfrac18$이고, $\kappa_k(\log D)=\kappa_k(\varepsilon)\sum_{p\mid n}(\log p)^k$이다.

**명제 10.1.** $\Re\rho\ge\tfrac12$인 중복도 $m$의 비자명 영점에서, $\sum_n\mu^2(n)\sum_{p\mid n}(\log p)^k\,n^{-s}$는 $k>m$이면 $k-m$차 극을 가지고, $k\le m$이면 정칙이다.

**증명.** 해당 급수는 $\zeta(s)/\zeta(2s)$와 $\sum_p(\log p)^kp^{-s}/(1+p^{-s})$의 곱이다. 두 번째 인수에서 $(-d/ds)^{k-1}(-\zeta'/\zeta)(s)\sim-m(k-1)!(s-\rho)^{-k}$를 빼면 $\rho$ 근방에 정칙이다. 소수 제곱 항도 $\zeta(2\rho)\ne0$에 의해 이 근방으로 연장된다. 첫 번째 인수는 정확히 $m$차 영점이므로 결론이 따른다. $\square$

진폭이 $\zeta'(\rho)$이므로 2차 모멘트는 단순 영점만 본다. 특히 $C_2$의 오차에서 얻는 영점 배제도 단순 영점에 한정된다. 아크사인에서는 홀수 누적량이 $0$이라 이중 영점은 4차부터 보인다.

<a id="section-10-4"></a>

### 10.4. 영점 높이와 미분값의 복원

$C_2(x)=\sum_{n\le x}c_n$, $\psi(x)=\sum_{n\le x}\Lambda(n)$, $y=\log x$로 둔다.
각 잔차를 $e^{-y/2}$로 정규화하면, 단순 영점 $\rho=\beta_\rho+i\gamma$의 기여는 $2\Re(a_\rho e^{(\rho-1/2)y})$이다. 세 디리클레 급수의 국소 유수를 비교하면 RH 없이

| 잔차 | $a_\rho$ |
|:--|:--|
| $M(x)$ | $1/(\rho\zeta'(\rho))$ |
| $C_2(x)-xP_2(\log x)$ | $-\zeta'(\rho)/\rho$ |
| $\psi(x)-x$ | $-1/\rho$ |

이다. $\psi$는 중복도 $m_\rho$인 영점도 계수 $-m_\rho/\rho$로 포착한다.
따라서 영점은 $\psi$에서 찾고, $C_2$의 복소 계수로 미분값을 추정하는 방안을 택한다.
두 통계는 뫼비우스 반전으로 연결되며, 단순 영점에서

<a id="equation-10-4"></a>

$$
(\mu*c)(n)=\Lambda(n)\log n,\qquad
\frac{a_\rho^{(2)}}{a_\rho^{(\psi)}}=\zeta'(\rho).
\tag{10.4}
$$

공통 컨볼루션 필터의 응답이 복소 주파수 $\nu_\rho=\gamma-i(\beta_\rho-\tfrac12)$에서 정의되고 영이 아니면 계수비가 보존된다. 조화 역산은 $\nu_\rho$를 맞추어 $\beta_\rho$도 추정할 수 있다. 임계선 밖 후보는 별도의 영점 검증으로 확인한다.

높은 영점의 후속 실험에서는 목표 높이로 주파수를 이동하고 해당 대역의 진폭을 보존하는 필터를 사용한다. 절단오차와 경계효과를 평가하며 FFT와 조화 역산(필터 대각화; Main 등)을 비교한다. 성능은 영점 위치를 입력하지 않고 얻은 추정 오차와 누락 수로 측정한다.

<a id="section-10-5"></a>

### 10.5. RH 아래 극한분포

**명제 10.2.** RH 아래
$$
E_2(y)=e^{-y/2}\bigl(C_2(e^y)-e^yP_2(y)\bigr)
$$
는 $B^2$-거의 주기함수이며 극한분포 $\nu$를 가진다. 추가로 단순 영점의 서로 다른 양의 높이가 $\mathbb Q$ 위에서 선형독립(LI)이면

<a id="equation-10-5"></a>

$$
\widehat\nu(\xi)=\prod_{\substack{\gamma>0\\\rho\text{ simple}}}
J_0\!\left(\frac{2|\zeta'(\rho)|\,\xi}{|\rho|}\right).
\tag{10.5}
$$

**증명.** 부록 C는 (10.3)의 정규화한 나머지 $\mathcal E(y,T)=e^{-y/2}R_2(e^y,T)$에 대해
$$
\lim_{Y\to\infty}\frac1Y\int_{\log2}^{Y}|\mathcal E(y,e^Y)|^2\,dy=0
$$
을 준다. 볼록성 상계와 Cauchy 추정으로 $\zeta'(\tfrac12+i\gamma)\ll_\varepsilon\gamma^{1/4+\varepsilon}$이므로 $r_\gamma=-2\zeta'(\rho)/\rho\ll_\varepsilon\gamma^{-3/4+\varepsilon}$이다. 영점 개수 $N(T+1)-N(T)\ll\log T$와 함께 Akbary–Ng–Shahabi의 따름정리 1.3(a)를 적용하면 $B^2$-거의 주기성과 극한분포의 존재가 따른다. LI 아래 곱 공식은 같은 논문의 정리 1.9이다. $\square$

중복 영점은 계수가 $0$이므로 단순영점 가정이 필요 없다. 단순 영점에서 $a_\rho^{(M)}=-(a_\rho^{(\psi)})^2/a_\rho^{(2)}$이지만, 그 제곱합의 제어에는 여전히 역미분값의 추정이 필요하다. 이 결과로 Ng의 $J_{-1}$ 가정을 제거할 수는 없다.

<a id="appendix-b"></a>

## 부록 B. 수치 확인

**방법.** 제곱 인수 없는 수와 $\omega$를 체로 구하고 약수 쌍을 직접 합했다($x\le10^7$).
$T_2$와 $c_n$은 각각 $5\cdot10^7$, $3\cdot10^7$까지 계산했다.
영점과 $\zeta'(\rho)$, $\zeta(2\rho)$는 Euler–Maclaurin 공식으로 계산했으며, 높이 $1000$까지의 영점 649개를 사용했다.

<a id="appendix-b-2"></a>

### B.2. 뫼비우스 경우 (§8)

**서로소 해제 $(-1,1)$, (8.3):** $A^\natural L/x$의 값이다($10^7$).

| $F$ | 실측 | 1차 $-\tfrac6{\pi^2}F'(0)$ | 2차 |
|:--|--:|--:|--:|
| $u$ | −0.6514 | −0.6079 | −0.6456 |
| $u^2$ | −0.0501 | 0 | −0.0435 |
| $e^u$ | −0.6772 | −0.6079 | −0.6674 |

남은 오차에 $L^2$을 곱한 값은 $10^6$과 $10^7$에서 각각 −1.52→−1.50, −1.70→−1.71, −2.53→−2.54이다.

<a id="appendix-b-3"></a>

### B.3. 베타 급수의 계수 (§9, 아크사인)

$G$를 다중원 $|s_i-1|=0.3$ 위에서 계산하고 2차원 FFT로 $G_{ab}$를 구했다. 검산으로 $G_{00}=0.6079271$이다.

| $j$ | 개별 $\max\lvert G_{ab}\rvert$ | $m_{j,0}$ | $m_{j,1}$ | $m_{j,2}$ | $m_{j,2}/m_{j-1,2}$ |
|--:|--:|--:|--:|--:|--:|
| 8 | 5.3 | −0.9947 | 0.8949 | −68.89 | −1.981 |
| 10 | 14.5 | −0.9987 | 0.9116 | −273.26 | −1.995 |
| 12 | 42.8 | −0.9997 | 0.9168 | −1094.77 | −2.003 |
| 14 | 134.1 | −0.9999 | 0.9183 | −4401.24 | −2.005 |

- **$m_{j,0}$:** $s=0$의 극에서 예측되는 $(-1)^{j+1}$로 간다.
- **$m_{j,2}$:** 비가 $-2$로 가므로 $s=\tfrac12$의 단순극이 보인다. $m_{14,2}/(-2)^{14}=-0.2686$이 (9.4)의 $-0.2738$로 다가간다. 느리게 다가가는 것은 $s=\tfrac13$의 다음 층 때문으로 보인다.

<a id="appendix-b-4"></a>

### B.4. 영점 층 (§10)

**알려진 영점으로 잔차 재현.** (10.3)을 적분한 $\int_0^x C_2(t)\,dt$의 잔차를 $x^{3/2}$로 나눈 것과, 알려진 영점 $N$개의 합 $-2\mathrm{Re}\sum\zeta'(\rho)x^{\rho+1}/(\rho(\rho+1))$를 $x\in[10^6,3\cdot10^7]$에서 비교했다.

| 입력 영점 수 | 1 | 10 | 79 | 649 |
|:--|--:|--:|--:|--:|
| 상관 | 0.594 | 0.928 | 0.991 | 0.995 |

잔차의 표준편차는 0.0087, 차이의 표준편차는 0.0009이다. 평활화하지 않은 합에서는 상관이 0.126, 0.383, 0.654, 0.837이다.

**스펙트럼.** 잔차를 $\log x$에 대해 푸리에 변환한 앞의 세 봉우리는 다음과 같다.

| 봉우리 | 14.20 | 21.08 | 25.04 |
|:--|--:|--:|--:|
| $\gamma$ | 14.13 | 21.02 | 25.01 |

**아크사인 설정의 $T_2$.** $x(AL^2+BL+C)+c\sqrt x$로 피팅한 편향은 $c=4.335$로, 예측 $-3\zeta(\tfrac12)=4.381$과 맞는다. 선도 계수 $A=0.30394$도 $1/(2\zeta(2))=0.30396$과 맞는다. 진동은 $\zeta(2s)$ 층과 $s=\tfrac13$ 항이 섞여, 영점 649개로 상관 0.70까지만 나온다.

<a id="appendix-c"></a>

## 부록 C. $C_2$의 절단오차: RH 아래 평균제곱 추정

$Q(s)=\zeta(s)(\zeta'/\zeta)'(s)$, $a=\tfrac14$로 둔다. 소인수분해에서
$$
0\le c_n=\frac12\sum_{p^h\parallel n}h(h+1)(\log p)^2\le(\log n)^2.
$$
따라서 $2\le x\le T$, $c=1+1/\log x$에서 유효 Perron 공식의 오차는
<a id="equation-c-1"></a>

$$
O\!\left(\frac{x\log^3(2x)}T+\log^2(2x)\right)
\tag{C.1}
$$
이다. $n\asymp x$에서는 $\min(1,x/(T|n-x|))$를 합하고, 나머지는 절대수렴으로 처리한다. 정수 끝점의 반가중치 차이도 이 오차에 포함된다.

RH 아래 $\Re s=a$는 모든 비자명 영점과 거리가 적어도 $1/4$이다. 로그미분의 부분분수 전개와 $N(t+1)-N(t)\ll\log(t+2)$에서 $(\zeta'/\zeta)'(a+it)\ll\log(|t|+2)$를 얻는다. 볼록성 상계를 곱하면, 고정된 $0<\varepsilon<1/8$에 대해
<a id="equation-c-2"></a>

$$
\Phi(t):=\frac{Q(a+it)}{a+it}
\ll_\varepsilon(1+|t|)^{-5/8+\varepsilon},
\qquad \Phi\in L^2(\mathbb R).
\tag{C.2}
$$

$X=e^Y$라 두고 모든 영점 높이에서 $\gg1/\log X$ 떨어진 $T\in[X,X+1]$를 고른다. 높이 $\pm T$에서는 $a\le\sigma\le c$에 걸쳐 $(\zeta'/\zeta)'(\sigma\pm iT)\ll\log^3T$이다. Perron 경로를 $\Re s=a$로 옮길 때 두 수평선의 적분은 볼록성 상계로
$$
\ll_\varepsilon T^\varepsilon\left(\frac xT+\frac{x^{1/4}}{T^{5/8}}\right).
$$
지나는 극은 $s=1$과 단순 영점뿐이며 그 유수가 (10.3)의 주항과 영점 합을 준다. 왼쪽 수직선 적분은 $x^{1/4}V_T(\log x)$이고
<a id="equation-c-3"></a>

$$
V_T(y)=\frac1{2\pi}\int_{-T}^{T}\Phi(t)e^{ity}\,dt,
\qquad
\int_{\mathbb R}|V_T(y)|^2\,dy
=\frac1{2\pi}\int_{-T}^{T}|\Phi(t)|^2\,dt\ll1
\tag{C.3}
$$
이다. 마지막 식은 Plancherel 항등식이다. (C.1)과 수평선 추정을 합치면 $\log2\le y\le Y$에서
$$
\mathcal E(y,T)=e^{-y/4}V_T(y)
+O_\varepsilon\!\left((1+y)^2e^{-y/2}+(1+Y)^3e^{-(1/2-\varepsilon)Y}\right).
$$
따라서
<a id="equation-c-4"></a>

$$
\frac1Y\int_{\log2}^{Y}|\mathcal E(y,T)|^2\,dy
\ll_\varepsilon\frac1Y+(1+Y)^6e^{-(1-2\varepsilon)Y}\longrightarrow0.
\tag{C.4}
$$
$T$를 $X=e^Y$로 바꿀 때 추가되는 영점은 $O(\log X)$개이고, 각 정규화한 계수는 $O_\varepsilon(X^{-3/4+\varepsilon})$이다. 이 차이도 균등하게 $0$으로 가므로 §10.5의 평균제곱 조건이 성립한다.

<a id="references"></a>

## 참고 문헌

- H. L. Montgomery and R. C. Vaughan, [Multiplicative Number Theory I, §6.1](https://personal.science.psu.edu/rcv4/personal/Publications/MNTI/10.0_pp_168_198_The_Prime_Number_Theorem.pdf), 정리 6.6–6.7: 영점 없는 영역과 복소 로그 추정.
- G. Tenenbaum, *Introduction to Analytic and Probabilistic Number Theory*, II.5: Selberg–Delange 방법.
- NIST DLMF, [§5.9](https://dlmf.nist.gov/5.9): Hankel 적분. [§13.6](https://dlmf.nist.gov/13.6): Bessel 함수. [§25.16](https://dlmf.nist.gov/25.16): $\psi$의 명시적 공식.
- N. Ng, [The distribution of the summatory function of the Möbius function](https://arxiv.org/abs/math/0310381), Proc. London Math. Soc. (2004): $M(x)$의 극한분포와 $J_{-1}(T)$ 가설.
- A. Akbary, N. Ng, M. Shahabi, [Limiting distributions of the classical error terms of prime number theory](https://arxiv.org/pdf/1306.1657), Q. J. Math.: 따름정리 1.3(a)의 극한분포 기준과 정리 1.9의 LI 아래 곱 공식.
- J. Main, V. A. Mandelshtam, G. Wunner, H. S. Taylor, [Harmonic inversion as a general method for periodic orbit quantization](https://arxiv.org/abs/chao-dyn/9709009) (1998): 필터 대각화를 이용한 영점 복원.
- J. Main, P. A. Dando, Dž. Belkić, H. S. Taylor, [Decimation and harmonic inversion of periodic orbit signals](https://arxiv.org/abs/chao-dyn/9912031): 대역 필터와 조화 역산을 분리한 영점 복원.
