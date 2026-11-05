# M과 ψ: 코사인 핵과 사인 핵

Mertens 함수 $M(x)=\sum_{n\le x}\mu(n)$과 Chebyshev 함수 $\psi(x)=\sum_{n\le x}\Lambda(n)$을 같은 이산핵 틀에서 비교한다. $M$은 지수 $(-\tfrac12,-\tfrac12)$에서, $\psi$는 $(-\tfrac12,\tfrac12)$에서 나온다. 소수마다 전자는 코사인 동전을, 후자는 사인 동전을 가지며, 두 함수의 극점·영점 성분은 모두 아크사인 핵에서 만들어진다. §4는 $\psi$ 쪽 소수핵으로 영점 측도를 정확히 쓴다. 절·명제·수식 번호는 이 문서 안에서 독립적이다.

기호는 [주정리 §8](main-theorem.md#section-8)의 $d_z$, 이산핵 $\mathscr K_n^{z,w}$, 지수합 $\widetilde S_x^{z,w}(\xi)$를 따르고, 베타핵 $\mathcal K_{z,w}$와 그 변환 $\mathcal B_{z,w}(\xi)=\langle\mathcal K_{z,w},e^{i\xi u}\rangle$는 [주정리 §4](main-theorem.md#section-4)와 [주정리 (5.4)](main-theorem.md#equation-5-4)를 따른다. $L=\log x$, $\alpha=\xi/L$이고, $a(u)=1/(\pi\sqrt{u(1-u)})$는 아크사인 밀도, $\varepsilon=d_0$는 디리클레 합성곱의 항등원이다.

| 결과 | 내용 | 범위 |
|:--|:--|:--|
| [명제 1.1](#section-1) | $d_{-1/2}*d_{1/2}=\varepsilon$의 로그 모멘트와 $\Lambda$, $\Lambda_2$ | 정확한 항등식 |
| [명제 2.1](#section-2) | $\psi$의 극·영점 국소 핵, 1차에서는 정확 | **[국소 계산]** |
| [표 (2.4)](#equation-2-4) | $M$과 $\psi$의 네 국소 핵을 아크사인 각도로 쓴 꼬리 | **[국소 계산]** |
| [명제 2.2](#proposition-2-2) | $M$과 $\psi$의 영점 성분을 선분 위 아크사인 가중 적분으로 쓴 표현, $M$의 $1/L$ 전개 (2.7) | (2.5)는 정확한 항등식(단순 영점), (2.7)은 점근 전개 |
| [(3.2)](#section-3) | 사인 쪽의 주파수 판: $\arg\zeta$의 극한분포 | Bohr–Jessen, 가정 없음 |
| [명제 4.1](#section-4) | 대역 제한 시험함수에서 영점합이 $P$ 아래 소수핵의 유한합과 같음 | Guinand–Weil, 가정 없음 |

- [§1. 두 분해와 소수별 동전](#section-1)
- [§2. 국소 핵](#section-2)
- [§3. M과 ψ의 비교](#section-3)
- [§4. 소수핵으로 쓴 영점 측도](#section-4)
- [검산](#checks)
- [참고 문헌](#references)

<a id="section-1"></a>

## 1. 두 분해와 소수별 동전

$\zeta^{-1}=\zeta^{-1/2}\zeta^{-1/2}$와 $1=\zeta^{-1/2}\zeta^{1/2}$에서

<a id="equation-1-1"></a>

$$
d_{-1/2}*d_{-1/2}=\mu,\qquad d_{-1/2}*d_{1/2}=\varepsilon
\tag{1.1}
$$

이다. $\mathfrak M_x(\xi)=\widetilde S_x^{-1/2,-1/2}(\xi)$([주정리 (8.16)](main-theorem.md#equation-8-16))와 함께

$$
\mathfrak P_x(\xi):=\widetilde S_x^{-1/2,1/2}(\xi)=\sum_{dm\le x}d_{-1/2}(d)\,d_{1/2}(m)\,e^{i\xi\log d/L}
$$

로 둔다. [주정리 (8.15)](main-theorem.md#equation-8-15)에서 Perron 피적분함수는 $\zeta(s-i\alpha)^{-1/2}\zeta(s)^{1/2}x^s/s$이다. $\mathfrak P_x(0)=1$이므로 정보는 모두 $\xi$ 방향, 즉 약수의 위치에 있다.

**명제 1.1.** 모든 $n\ge1$에서

<a id="equation-1-2"></a>

$$
\sum_{dm=n}d_{-1/2}(d)d_{1/2}(m)\log d=-\tfrac12\Lambda(n),
\qquad
\sum_{dm=n}d_{-1/2}(d)d_{1/2}(m)\log^2d=\tfrac14\Lambda_2(n)-\tfrac34\Lambda(n)\log n
\tag{1.2}
$$

이다. $\Lambda_2=\Lambda\log+\Lambda*\Lambda$는 Selberg의 일반화 폰 망골트 함수이다. 따라서

<a id="equation-1-3"></a>

$$
\mathfrak M_x(0)=M(x),
\qquad
\mathfrak P_x'(0)=-\frac{i\,\psi(x)}{2L},
\qquad
\mathfrak P_x''(0)=-\frac1{L^2}\sum_{n\le x}\Big(\tfrac14\Lambda_2(n)-\tfrac34\Lambda(n)\log n\Big).
\tag{1.3}
$$

**증명.** 왼쪽은 $(-\partial_{s_1})^k\zeta(s_1)^{-1/2}\big|_{s_1=s}\cdot\zeta(s)^{1/2}$의 디리클레 계수이다. $k=1$이면 $\tfrac12\zeta'/\zeta$이고, $k=2$이면 $\tfrac34(\zeta'/\zeta)^2-\tfrac12\zeta''/\zeta=\tfrac14(\zeta'/\zeta)^2-\tfrac12(\zeta'/\zeta)'$이다. $-\zeta'/\zeta=\sum\Lambda(n)n^{-s}$를 넣는다. $\square$

$\log n$ 좌표에서는 $\widetilde A_x^{-1/2,1/2}(u)=-\tfrac12\sum_{2\le n\le x}\Lambda(n)/\log n$이다. 이는 [주정리 정리 2](main-theorem.md#section-8-3)의 주항 $(x/L)\langle\mathcal K_{-1/2,1/2},u\rangle=-x/(2L)$과 맞는다. $M$의 대응식은 [주정리 (8.13)](main-theorem.md#equation-8-13)이다.

**소수별 동전.** 제곱 인수 없는 $n$에서는 $d_z(d)=z^{\omega(d)}$이므로

<a id="equation-1-4"></a>

$$
\sum_{d\mid n}d_z(d)\,d_w(n/d)\,d^{i\alpha}
=n^{i\alpha/2}\prod_{p\mid n}\big[(z+w)\cos\theta_p+i(z-w)\sin\theta_p\big],
\qquad\theta_p=\tfrac\alpha2\log p
\tag{1.4}
$$

이다. $z=w$이면 코사인, $z=-w$이면 사인만 남는다. DDT의 $(\tfrac12,\tfrac12)$는 $\cos\theta_p$, $M$의 $(-\tfrac12,-\tfrac12)$는 $-\cos\theta_p$, $\psi$의 $(-\tfrac12,\tfrac12)$는 $-i\sin\theta_p$를 준다. $\alpha=\xi/L$이고 $p\le x$이면 $\theta_p=\tfrac\xi2\cdot\log p/L$이다. 따라서 $0\le\xi\le2\pi$에서 $\arccos(\cos\theta_p)$는 접힘 없이 약수 좌표 $\log p/L$을 돌려준다.

<a id="section-2"></a>

## 2. 국소 핵

$g_{x,\alpha}(s)=\zeta(s-i\alpha)^{-1/2}\zeta(s)^{1/2}x^s/s$로 둔다. 쌍 $\{1,1+i\alpha\}$와 각 영점 쌍 $\{\rho,\rho+i\alpha\}$에서 두 지수의 합이 $0$이다. 그러므로 이 쌍들을 잇는 선분을 제거한 영역에서 $g_{x,\alpha}$는 일가이고, 오른쪽의 Euler 가지에 맞춘다. $\alpha=0$이면 $g_{x,0}=x^s/s$이다. [전역 공식 (1.13)](global-formula.md#equation-1-13)처럼 각 선분만 감싸는 작은 폐곡선 $\Gamma_1$, $\Gamma_\rho$를 잡고

$$
\mathcal P^\psi_1(x,\xi)=\frac1{2\pi i}\int_{\Gamma_1}g_{x,\alpha}(s)\,ds,
\qquad
\mathcal I^\psi_\rho(x,\xi)=\frac1{2\pi i}\int_{\Gamma_\rho}g_{x,\alpha}(s)\,ds
$$

로 둔다.

**명제 2.1 (ψ의 국소 핵).** $R>0$을 고정하고 $|\xi|\le R$이라 하자. $\rho$의 중복도가 $m$이면

<a id="equation-2-1"></a>

$$
\mathcal I^\psi_\rho(x,\xi)=\frac{x^\rho}{\rho L}\,\mathcal B_{m/2,-m/2}(\xi)+O_{\rho,R}\!\left(\frac{x^{\Re\rho}}{L^2}\right),
\qquad
\mathcal P^\psi_1(x,\xi)=\frac xL\,\mathcal B_{-1/2,1/2}(\xi)+O_R\!\left(\frac x{L^2}\right)
\tag{2.1}
$$

이다. $\xi=0$에서의 1차 미분은 오차 없이

<a id="equation-2-2"></a>

$$
\partial_\xi\mathcal I^\psi_\rho(x,0)=\frac{im}{2L}\cdot\frac{x^\rho}{\rho},
\qquad
\partial_\xi\mathcal P^\psi_1(x,0)=-\frac{i}{2L}\,x
\tag{2.2}
$$

이다.

**증명.** $s=\rho+q/L$로 두고 $\zeta(\rho+h)=b_\rho h^m(1+O_\rho(h))$를 쓰면 $\zeta(s)/\zeta(s-i\alpha)=q^m(q-i\xi)^{-m}(1+O_{\rho,R}(L^{-1}))$이다. 계수 $b_\rho$는 약분되고, 큰 $|q|$에서 $1$에 가까운 가지가 Euler 가지와 이어진다. $x^s=x^\rho e^q$, $ds=dq/L$, $1/s=1/\rho+O_\rho(L^{-1})$이므로 [주정리 (5.4)](main-theorem.md#equation-5-4)를 $P_\xi(q)=(q-i\xi)^{-m/2}q^{m/2}$에 적용하면 첫 식을 얻는다. 극에서는 $\zeta(1+h)=h^{-1}(1+O(h))$이므로 $P_\xi(q)=(q-i\xi)^{1/2}q^{-1/2}$이다. (2.2)는 $\partial_\alpha g_{x,\alpha}\big|_{\alpha=0}=\tfrac i2(\zeta'/\zeta)(s)\,x^s/s$의 유수이다. $\square$

(2.2)는 고전 명시공식 $\psi(x)=x-\sum_\rho x^\rho/\rho-\cdots$의 항과 같다. 실제로 $\mathfrak P_x'(0)=-i\psi(x)/(2L)$에 극은 $x$를, 각 영점은 $-m\,x^\rho/\rho$를 준다. 따라서 $\psi$에서 새로 얻는 정보는 1차를 넘는 $\xi$ 의존성, 즉 핵의 모양이다. 영점을 무한히 합하는 문제는 [전역 공식 §1.5](global-formula.md#section-1-5)의 주의와 같다.

**네 핵의 아크사인 표현.** [주정리 (7.1)](main-theorem.md#equation-7-1)의 $\mathcal K_{-1/2,1/2}=[(1-u)a]'$와 반사 $R_*$에서 $\mathcal K_{1/2,-1/2}=-(ua)'$이다. [주정리 (4.6)](main-theorem.md#equation-4-6)과 $\langle a,(2u-1)e^{i\xi u}\rangle=ie^{i\xi/2}J_1(\xi/2)$를 쓰면

<a id="equation-2-3"></a>

$$
\begin{aligned}
\mathcal B_{1/2,-1/2}(\xi)&=i\xi\,\langle a,u\,e^{i\xi u}\rangle=\tfrac{i\xi}2e^{i\xi/2}\big(J_0+iJ_1\big)(\tfrac\xi2),\\
\mathcal B_{-1/2,1/2}(\xi)&=-i\xi\,\langle a,(1-u)\,e^{i\xi u}\rangle=-\tfrac{i\xi}2e^{i\xi/2}\big(J_0-iJ_1\big)(\tfrac\xi2),\\
\mathcal B_{-1/2,-1/2}(\xi)&=-\tfrac12\big(\mathcal B_{-1/2,1/2}+\mathcal B_{1/2,-1/2}\big)(\xi)=\tfrac\xi2e^{i\xi/2}J_1(\tfrac\xi2)
\end{aligned}
\tag{2.3}
$$

이다. $\psi$의 두 핵은 아크사인 밀도에 $u$ 또는 $1-u$를 곱해 미분한 것이고, $M$의 극 핵은 그 둘의 평균에 $-1$을 곱한 것이다.

$\vartheta=\arccos(2v-1)\in(0,\pi)$로 두면 아크사인 분포는 $\vartheta$에 대해 균등하다. [주정리 (7.5)](main-theorem.md#equation-7-5)의 분포 꼬리값 $T_{z,w}(v)=\langle\mathbf 1_{(v,\infty)}\mathcal K_{z,w},\chi\rangle$과 반사를 쓰면, 단순 영점에서 네 성분은 표 (2.4)와 같다.

<a id="equation-2-4"></a>

**(2.4) 네 국소 핵**

| 성분 | 크기 | 핵 | 꼬리 $T(v)$ |
|:--|:--|:--|:--|
| $M$의 영점 | $x^\rho/(\rho\zeta'(\rho))$ | $\mathcal K_{1/2,1/2}=a$ | $\vartheta/\pi$ |
| $M$의 극 | $x/L^2$ | $\mathcal K_{-1/2,-1/2}=[(u-\tfrac12)a]'$ | $-\cot\vartheta/\pi$ |
| $\psi$의 극 | $x/L$ | $\mathcal K_{-1/2,1/2}=[(1-u)a]'$ | $-\tan(\vartheta/2)/\pi$ |
| $\psi$의 영점 | $x^\rho/(\rho L)$ | $\mathcal K_{1/2,-1/2}=-(ua)'$ | $\cot(\vartheta/2)/\pi$ |

$M$의 두 줄은 [전역 공식 (1.18)–(1.19)](global-formula.md#equation-1-18)이다. 두 영점 줄은 [명제 2.2](#proposition-2-2)의 정확한 선분 표현에서 $R\equiv1$로 둔 것이기도 하다. $T$는 분포의 꼬리값이며, 산술 꼬리합의 점근식을 주지는 않는다.

**2차 모멘트.** $\langle\mathcal K_{-1/2,1/2},u^2\rangle=-\tfrac14$이므로 (1.3)과 극 성분에서 $\sum_{n\le x}(\tfrac14\Lambda_2-\tfrac34\Lambda\log n)$의 주항은 $-\tfrac14xL$이다. 이는 Selberg 공식 $\sum_{n\le x}\Lambda_2(n)=2xL+O(x)$와 $\sum_{n\le x}\Lambda(n)\log n=xL+O(x)$의 조합과 맞는다. 영점 쪽 2차 계수는 $\langle\mathcal K_{1/2,-1/2},u^2\rangle=\tfrac34$이며, 1차와 달리 $O(x^{\Re\rho}/L^2)$의 보정을 가진다.

<a id="proposition-2-2"></a>

**선분 표현.** 명제 2.1과 표 (2.4)의 영점 핵은 $L\to\infty$에서의 국소 근사이다. 단순 영점에서는 윤곽을 선분 위로 줄여 정확한 식을 얻는다. $M$의 Perron 피적분함수 $h_{x,\alpha}(s)=\zeta(s-i\alpha)^{-1/2}\zeta(s)^{-1/2}x^s/s$([주정리 (8.15)](main-theorem.md#equation-8-15))도 쌍 $\{\rho,\rho+i\alpha\}$에서 지수 합이 정수 $-1$이다. 그러므로 $g_{x,\alpha}$와 같이 선분을 제거한 영역에서 일가인 가지를 오른쪽의 Euler 가지에 맞추고

$$
\mathcal I^M_\rho(x,\xi)=\frac1{2\pi i}\int_{\Gamma_\rho}h_{x,\alpha}(s)\,ds
$$

로 둔다. $\alpha=0$이면 $h_{x,0}=x^s/(s\zeta(s))$이다.

**명제 2.2 (영점 성분의 선분 표현).** $\rho$를 단순 영점이라 하고, 닫힌 선분 $[\rho-i|\alpha|,\rho+i|\alpha|]$ 위에 $\rho$ 외의 영점이 없다고 하자. $f(s)=\zeta(s)/(s-\rho)$, $s_u=\rho+i\alpha u$로 쓰면

<a id="equation-2-5"></a>

$$
\mathcal I^M_\rho(x,\xi)=\frac{x^\rho}{\rho\,\zeta'(\rho)}\int_0^1a(u)\,e^{i\xi u}R_\alpha(u)\,du,
\qquad
\mathcal I^\psi_\rho(x,\xi)=\frac{x^\rho}{\rho}\cdot\frac{i\xi}L\int_0^1u\,a(u)\,e^{i\xi u}R^\psi_\alpha(u)\,du
\tag{2.5}
$$

이다. 여기서

<a id="equation-2-6"></a>

$$
R_\alpha(u)=\frac{\rho}{s_u}\cdot\frac{\zeta'(\rho)}{\sqrt{f(s_u)\,f(s_u-i\alpha)}},
\qquad
R^\psi_\alpha(u)=\frac{\rho}{s_u}\sqrt{\frac{f(s_u)}{f(s_u-i\alpha)}}
\tag{2.6}
$$

이고, 제곱근은 $(\alpha,u)$에 대해 연속이며 $\alpha=0$에서 각각 $\zeta'(\rho)$, $1$이 되는 가지이다. 따라서 $R_0=R^\psi_0\equiv1$이다.

$x^{s_u}=x^\rho e^{i\xi u}$이므로 선분의 매개변수 $u$가 곧 약수 좌표이다. $M$의 아크사인 밀도는 선분 양 끝에 놓인 지수 $-\tfrac12$의 두 분기점에서 나온다. $\psi$에서는 $s=\rho$ 쪽의 지수가 $+\tfrac12$로 바뀌어 무게가 $u\,a(u)=\pi^{-1}\sqrt{u/(1-u)}$가 된다. $R\equiv1$로 두면 (2.5)는 표 (2.4)의 두 영점 줄 $\frac{x^\rho}{\rho\zeta'(\rho)}\mathcal B_{1/2,1/2}(\xi)$와 $\frac{x^\rho}{\rho L}\mathcal B_{1/2,-1/2}(\xi)$가 된다([주정리 (4.6)](main-theorem.md#equation-4-6), (2.3)). (2.5)는 $L$에 대한 극한식이 아니며, $\rho_1$에서는 $|\alpha|<\gamma_2-\gamma_1\approx6.89$이면 성립한다.

**증명.** 선분 근방에서 $\zeta(s)\zeta(s-i\alpha)=(s-\rho)(s-\rho-i\alpha)f(s)f(s-i\alpha)$이고, 가정에 의해 $f(s)f(s-i\alpha)$는 선분 근방에서 영이 아니다. $Q(s)=[(s-\rho)(s-\rho-i\alpha)]^{-1/2}$를 선분 밖에서 일가이고 무한대에서 $1/(s-\rho)$처럼 행동하는 가지로 잡는다. $\Gamma_\rho$ 위에서 $h_{x,\alpha}$와 $Q\,[f(s)f(s-i\alpha)]^{-1/2}x^s/s$는 제곱이 같은 일가함수이므로 부호만 다를 수 있다. $\Gamma_\rho$에서 선분 오른쪽의 한 점을 잡으면 두 값은 모두 $\alpha$에 대해 연속이고, $\alpha=0$에서는 둘 다 $x^s/(s\zeta(s))$이다. 따라서 부호는 $+$이다.

선분 근방에서 해석적인 $\Psi$에 대해 $\Gamma_\rho$를 선분으로 줄인다. 끝점의 특이성은 $|s-\rho|^{-1/2}$, $|s-\rho-i\alpha|^{-1/2}$ 정도이므로 끝점을 감싸는 작은 원의 기여는 사라진다. 선분 위에서는 $(s_u-\rho)(s_u-\rho-i\alpha)=\alpha^2u(1-u)>0$이므로 양쪽 경계값은 $\pm\big(|\alpha|\sqrt{u(1-u)}\big)^{-1}$이다. 따라서 $\Psi$와 무관한 상수 $c$에 대해 $\frac1{2\pi i}\int_{\Gamma_\rho}Q\Psi\,ds=c\int_0^1\Psi(s_u)\,du/\sqrt{u(1-u)}$이다. $\Psi\equiv1$이면 $\Gamma_\rho$를 무한히 키울 수 있어 왼쪽이 $1$이므로 $c=1/\pi$이고,

$$
\frac1{2\pi i}\int_{\Gamma_\rho}Q(s)\Psi(s)\,ds=\int_0^1a(u)\,\Psi(s_u)\,du
$$

이다. $\Psi(s)=[f(s)f(s-i\alpha)]^{-1/2}x^s/s$를 넣으면 (2.5)의 첫 식을 얻는다.

$\psi$에서는 $Q^\psi(s)=(s-\rho)^{1/2}(s-\rho-i\alpha)^{-1/2}$를 무한대에서 $1$에 가까운 가지로 잡고, 부호는 앞과 같이 $\alpha=0$에서 정한다. 선분 위의 경계값은 $\pm i\sqrt{u/(1-u)}$이므로 같은 논증으로 $\frac1{2\pi i}\int_{\Gamma_\rho}Q^\psi\Psi\,ds=c^\psi\int_0^1\sqrt{u/(1-u)}\,\Psi(s_u)\,du$이다. 무한대에서 $Q^\psi=1+\tfrac{i\alpha}2(s-\rho)^{-1}+O(|s|^{-2})$이므로 $\Psi\equiv1$일 때 왼쪽은 $i\alpha/2$이고, $c^\psi=i\alpha/\pi$이다. 즉 $\frac1{2\pi i}\int_{\Gamma_\rho}Q^\psi\Psi\,ds=i\alpha\int_0^1u\,a(u)\Psi(s_u)\,du$이다. $\Psi(s)=[f(s)/f(s-i\alpha)]^{1/2}x^s/s$를 넣고 $\alpha=\xi/L$을 쓰면 둘째 식이다. $\square$

**$1/L$ 전개.** $f'(\rho)/f(\rho)=\zeta''(\rho)/(2\zeta'(\rho))=:\kappa_\rho$이므로 (2.6)에서 $u$에 균등하게

$$
R_\alpha(u)=1-i\alpha\Big[\frac u\rho+\Big(u-\frac12\Big)\kappa_\rho\Big]+O_\rho(\alpha^2)
$$

이다. $-i\xi\langle T,e^{i\xi u}\rangle=\langle T',e^{i\xi u}\rangle$와 $(ua)'=-\mathcal K_{1/2,-1/2}$, $\big((u-\tfrac12)a\big)'=\mathcal K_{-1/2,-1/2}$를 쓰면 $|\xi|\le R$에서

<a id="equation-2-7"></a>

$$
\mathcal I^M_\rho(x,\xi)=\frac{x^\rho}{\rho\,\zeta'(\rho)}\Big[\mathcal B_{1/2,1/2}(\xi)+\frac1L\Big(\kappa_\rho\,\mathcal B_{-1/2,-1/2}(\xi)-\frac{\mathcal B_{1/2,-1/2}(\xi)}{\rho}\Big)\Big]+O_{\rho,R}\!\left(\frac{x^{\Re\rho}}{L^2}\right)
\tag{2.7}
$$

이다. 보정항의 두 모양은 $M$의 극 핵과 $\psi$의 영점 핵이다. $\rho=\tfrac12+i\gamma$가 임계선 위에 있으면 $\zeta(\tfrac12+it)=e^{-i\vartheta_{\rm RS}(t)}Z(t)$에서

$$
\kappa_\rho=-\vartheta_{\rm RS}'(\gamma)-\frac i2\,\frac{Z''(\gamma)}{Z'(\gamma)},
\qquad
\Re\kappa_\rho=-\vartheta_{\rm RS}'(\gamma)=-\tfrac12\log\frac{\gamma}{2\pi}+O(\gamma^{-2})
$$

이다. $\vartheta_{\rm RS}$는 Riemann–Siegel 세타함수로 표 (2.4)의 각도 $\vartheta$와 다르며, $Z$는 Hardy 함수이다. 따라서 보정의 크기는 대략 $\log\gamma/L$이며, 높은 영점일수록 아크사인 모양이 나타나는 데 더 큰 $L$이 필요하다. $\rho_1$에서는 $\kappa_{\rho_1}=-0.405274\ldots-0.082158\ldots i$이다.

<a id="section-3"></a>

## 3. M과 ψ의 비교

- **동전.** $M$은 코사인(대칭, $z=w$), $\psi$는 사인(반대칭, $z=-w$)이다.
- **진폭.** $M$의 영점 진폭 $1/(\rho\zeta'(\rho))$는 미분값에 의존하고, 중복 영점에서는 $L^{m-1}$로 커진다([전역 공식 (1.17)](global-formula.md#equation-1-17)). $\psi$에서는 지수 합이 $0$이라 $b_\rho$가 약분되므로 진폭이 늘 $1/\rho$이다. 중복도는 핵 $\mathcal K_{m/2,-m/2}$의 모양과 중심 $\langle\mathcal K_{m/2,-m/2},u\rangle=m/2$에만 나타난다.
- **위치.** $M$의 영점 성분은 아크사인 꼴로 약수 좌표 전체에 퍼진다. $\psi$의 성분은 끝점에 몰린다. 영점 성분은 $\vartheta\to0$, 즉 $u\to1$인 $d_{-1/2}$ 쪽에, 극 성분은 $\vartheta\to\pi$, 즉 $u\to0$인 $d_{1/2}$ 쪽에 쏠린다. $\Lambda$가 소수 거듭제곱에만 놓인다는 사실이 약수 좌표에서는 이렇게 보인다.
- **규모.** $M$의 극 성분 $x/L^2$은 $F=1,u$에서 사라지고 $\xi=0$에서는 정확히 $0$이므로([전역 공식 §1.5](global-formula.md#section-1-5)), $M(x)$의 진동은 영점 성분이 맡는다. $\psi$의 극 성분 $x/L$은 1차 모멘트에서 주항 $x$가 된다.

**주파수 쪽.** 두 변수 변환 [주정리 (8.5)](main-theorem.md#equation-8-5)에 $s=\sigma-i\alpha$, $t=\sigma+i\alpha$를 넣으면, 실수 $z$와 $\sigma>1$에서

<a id="equation-3-1"></a>

$$
\sum_{d,m\ge1}\frac{d_z(d)d_z(m)}{(dm)^\sigma}\Big(\frac dm\Big)^{i\alpha}=|\zeta(\sigma+i\alpha)|^{2z},
\qquad
\sum_{d,m\ge1}\frac{d_{-z}(d)d_z(m)}{(dm)^\sigma}\Big(\frac dm\Big)^{i\alpha}=e^{2iz\arg\zeta(\sigma+i\alpha)}
\tag{3.1}
$$

이다. 코사인 쪽은 크기를, 사인 쪽은 위상을 준다. 첫 식은 [Ng 쌍대 (2.1)](ng-duality.md#equation-2-1)이고 그 극한분포가 [Ng 쌍대 정리 1](ng-duality.md#section-2-2)이다. 둘째 식의 짝은 Bohr–Jessen 정리의 허수부이다. $\sigma>\tfrac12$에서 $\alpha\in[0,A]$를 균등하게 고르면, $A\to\infty$일 때 $\arg\zeta(\sigma+i\alpha)$(오른쪽에서 수평으로 이어 붙인 가지)는 독립 균등 위상 $\theta_p$에 대한 $\sum_p\arctan\big(p^{-\sigma}\sin\theta_p/(1-p^{-\sigma}\cos\theta_p)\big)$의 분포로 수렴하고, 그 특성함수는

<a id="equation-3-2"></a>

$$
\prod_p{}_2F_1\!\left(-\tfrac\eta2,\tfrac\eta2;1;p^{-2\sigma}\right)
=\sum_{d\ge1}\frac{d_{-\eta/2}(d)\,d_{\eta/2}(d)}{d^{2\sigma}}
\tag{3.2}
$$

이다. 오른쪽은 [Ng 쌍대 (2.4)](ng-duality.md#equation-2-4)와 같은 중심 원자 합이고, 각 인수의 선도항은 아크사인 인수 $J_0(\eta p^{-\sigma})$이다.

<a id="section-4"></a>

## 4. 소수핵으로 쓴 영점 측도

§2에서 $\psi$의 영점 성분은 진폭에 $\zeta'(\rho)$가 없고 중복도만 남았다. 같은 이유로 $\psi$ 쪽의 소수별 양은 영점 측도를 준다. 이 절에서는 $\theta_p=t\log p$로 쓴다.

[Ng 쌍대 (2.3)](ng-duality.md#equation-2-3)의 $z=\tfrac12$ 동전 $Y_p(\sigma,t)=-\log\lvert1-p^{-\sigma}e^{-i\theta_p}\rvert=\sum_kp^{-k\sigma}\cos(k\theta_p)/k$에 대해, $\sigma>1$에서

<a id="equation-4-1"></a>

$$
\sum_pY_p=\log|\zeta(\sigma+it)|,
\qquad
\sum_p\partial_\sigma Y_p=\Re\frac{\zeta'}{\zeta}(\sigma+it)=\partial_t\arg\zeta(\sigma+it)
\tag{4.1}
$$

이다. 중복도 $m$인 영점 $\rho=\beta+i\gamma$ 근처에서 $\Re(\zeta'/\zeta)\approx m(\sigma-\beta)/\big((\sigma-\beta)^2+(t-\gamma)^2\big)$이므로, $\sigma\downarrow\beta$에서 넓이 $m\pi$인 Lorentz 봉우리가 생긴다. Riemann–von Mangoldt 공식 $N(t)=\theta(t)/\pi+1+S(t)$와 함께 보면 임계선 위 영점의 측도는 형식적으로

<a id="equation-4-2"></a>

$$
dN(t)\;\approx\;\Big(\frac{\theta'(t)}\pi+\frac1\pi\sum_p\partial_\sigma Y_p(\tfrac12,t)\Big)dt,
\qquad
\theta'(t)=\tfrac12\Re\frac{\Gamma'}{\Gamma}\big(\tfrac14+\tfrac{it}2\big)-\tfrac12\log\pi
\tag{4.2}
$$

이다. $\theta$는 Riemann–Siegel 세타 함수이다. 오른쪽 소수합은 $\sigma=\tfrac12$에서 수렴하지 않으므로, (4.2)는 아래 (4.3)의 형식적 판이다.

**σ 방향의 계단.** [Ng 쌍대 §4.3](ng-duality.md#section-4-3)의 미분·누적 단계에 한 단을 더하면 다음과 같다. 마지막 열은 $\sigma=\tfrac12$에서 소수별 진폭의 제곱합이다.

| 단계 | 소수별 양 | 합 | 성질 | 진폭 제곱합 |
|:--|:--|:--|:--|:--|
| $\int_{1/2}^\infty d\sigma$ | $\sum_k\cos(k\theta_p)/(k^2p^{k/2}\log p)$ | $\pi S_1(t)$ | 가정 없는 극한분포 ([Ng 쌍대 정리 2](ng-duality.md#section-3-2)) | $\sum_p1/(p\log^2p)<\infty$ |
| 그대로 | $Y_p$ | $\log\lvert\zeta\rvert$ | $\sigma=\tfrac12$에서 퇴화, Selberg 중심극한정리 | $\sum_p1/p=\infty$ |
| $\partial_\sigma$ | $-\log p\sum_kp^{-k\sigma}\cos(k\theta_p)$ | $\Re\zeta'/\zeta$ | 영점 측도, 시험함수로만 의미가 있음 | $\sum_p\log^2p/p=\infty$ |

**명제 4.1 (Guinand–Weil, 대역 제한판).** $P\ge2$이고 $g$가 $[-\log P,\log P]$에 지지된 짝인 $C^2$ 함수이며 $h(t)=\int g(u)e^{iut}du$라 하자. 그러면

<a id="equation-4-3"></a>

$$
\boxed{
\sum_\rho h(\gamma_\rho)
=h(\tfrac i2)+h(-\tfrac i2)
+\frac1\pi\int_{\mathbb R}h(t)\,\theta'(t)\,dt
-2\sum_{n<P}\frac{\Lambda(n)}{\sqrt n}\,g(\log n)
}
\tag{4.3}
$$

이다. 왼쪽은 모든 비자명 영점 $\rho=\tfrac12+i\gamma_\rho$를 중복도와 함께 더하고, 임계선 밖 영점은 복소수 $\gamma_\rho$로 들어간다. RH를 쓰지 않는다.

$\int h(t)\cos(t\log n)\,dt=2\pi g(\log n)$이므로 마지막 항은 $\frac1\pi\int h(t)\sum_{p^k<P}\partial_\sigma\big[p^{-k\sigma}\cos(k\theta_p)/k\big]_{\sigma=1/2}dt$이다. 즉 시험함수의 대역을 $\log P$로 제한하면 (4.2)의 소수합이 정확히 $P$ 아래에서 끊기고, 매끄럽게 한 영점 측도가 유한개 소수핵의 합과 같아진다. $h(\pm\tfrac i2)$는 완비 제타함수의 극 $s=0,1$에서 오며, $s=1$ 쪽이 §2의 $\psi$ 극 성분에 해당한다. $\theta'$ 항은 영점 밀도의 매끄러운 부분이다.

**영점 찾기.** $h_{t_0}(t)=h(t-t_0)+h(t+t_0)$, 즉 $g$ 대신 $2g(u)\cos(t_0u)$를 넣으면

<a id="equation-4-4"></a>

$$
\mathcal Z_P(t_0):=\sum_\rho\big[h(\gamma_\rho-t_0)+h(\gamma_\rho+t_0)\big]
\tag{4.4}
$$

를 (4.3)의 오른쪽, 곧 $P$ 아래의 소수만으로 계산할 수 있다. $\mathcal Z_P$는 영점마다 높이 $h(0)$인 봉우리를 갖는다. 봉우리의 폭은 $1/\log P$에 비례하므로, 해상도는 $\log P$에 반비례해서만 좋아진다. 아래 시험함수 $g(u)=\cos^4(\pi u/2\log P)$에서는 봉우리 주엽의 반폭(첫 영점까지의 거리)이 $3\pi/\log P$이고, 간격이 이보다 꽤 작은 영점 쌍은 한 봉우리로 합쳐진다. 간격 $\delta$인 쌍을 가르려면 대략 $\log P\gtrsim3\pi/\delta$가 필요하다 **[휴리스틱]**.

**[수치]** $g(u)=\cos^4(\pi u/2\log P)$로 둔다.

- (4.3): $P=10^3$, $t_0=0,20,50$에서 높이 $400$ 아래 영점의 합과 소수 쪽이 $10^{-12}$ 안에서 일치한다.
- (4.4): Γ 항을 $4g(0)\theta'(t_0)$로 근사하고 봉우리를 찾았다.
  - $P=10^4$이면 높이 $100$ 아래 영점 29개를 모두 찾고, 위치 오차는 최대 $0.0044$이다.
  - $P=10^5$이면 높이 $150$ 아래 52개를 모두 찾고, 최대 오차는 $0.0105$이다.
  - $P=10^5$, 높이 $300$ 아래 137개 중에서는 간격 $0.61$인 쌍 $294.965$, $295.573$이 한 봉우리($295.27$)로 합쳐진다. 반폭 $3\pi/\log P=0.82$보다 좁은 쌍이다. 나머지 135개는 오차 $0.071$ 안에서 찾는다.

**M 쪽.** 코사인 쪽에서는 [Ng 쌍대 명제 5.1](ng-duality.md#section-5-2)이 영점의 미분값을 준다. $\sigma\downarrow\tfrac12$에서 $\lvert1/(s\zeta(s))\rvert^2$의 봉우리 넓이가 $1/\lvert\rho\zeta'(\rho)\rvert^2$이다. $\sigma>1$에서는 $\lvert\zeta\rvert^{-2}=\prod_p(1-2p^{-\sigma}\cos\theta_p+p^{-2\sigma})$이지만, $\sigma\le1$에서는 이 곱을 점마다 쓸 수 없고 (4.3) 같은 유한 정확식도 없다. 차이는 덧셈꼴과 곱셈꼴에 있다. $\psi$는 $\log\zeta$의 미분이므로 소수에 대해 덧셈꼴이고, 명시공식의 반대쪽이 소수 거듭제곱에만 놓인다. $M$은 $1/\zeta$이므로 반대쪽이 $\mu(n)$, 즉 모든 제곱 인수 없는 정수이다. 그래서 영점의 위치와 중복도는 $\psi$ 쪽 소수핵으로 정확히 쓸 수 있고, 미분값은 $M$ 쪽에서 극한으로만 읽는다. 이는 [영점 기여 §4](zero-contributions.md#section-4)의 역할 분담, 즉 영점은 $\psi$에서 찾고 미분값은 따로 추정하는 방식과 같다.

<a id="checks"></a>

## 검산

[검산 스크립트](../scripts/check_mertens_psi.py)는 numpy와 scipy를 사용하며 다음을 확인한다.

- (1.1)–(1.2): $n\le2\times10^4$에서 기계 정밀도로 성립한다.
- (1.4): 제곱 인수 없는 $n$과 여러 $(z,w)$에서 소수별 동전으로의 분해를 확인한다.
- (2.3)–(2.4): 모멘트 (4.2)로 계산한 $\mathcal B_{z,w}$와 Bessel 표현, 아크사인 평균 표현, 꼬리값의 $\vartheta$ 표현이 $10^{-14}$ 안에서 일치한다.
- (2.1), (2.4): 첫 영점 $\rho_1$과 $s=1$ 주위의 윤곽 적분을 계산한다. $\xi=1,4$에서 네 성분의 상대오차는 $L=25$에서 $400$까지 $L$이 두 배가 될 때마다 반으로 줄어 $O(1/L)$을 따른다. (2.2)의 1차 미분은 $L$과 무관하게 $10^{-8}$ 안에서 정확하다.
- (3.2): $\sigma=2$, $\alpha\in[0,3000]$에서 $e^{i\eta\arg\zeta}$의 경험적 평균과 ${}_2F_1$ 곱의 차이는 $\eta=2,5,10$에서 $2\times10^{-5}$–$10^{-4}$이다. $J_0$ 선도항만 쓰면 차이가 $10^{-3}$–$6\times10^{-3}$으로 커져 소수 거듭제곱 보정이 보인다.
- (4.3)–(4.4): §4의 [수치] 항목. 영점은 Hardy $Z$ 함수의 부호 변화로 높이 $400$까지 구한다.

[선분 표현 검산 스크립트](../scripts/check_mertens_segment.py)도 numpy와 scipy만 사용하며, $\zeta$는 Euler–Maclaurin 합으로 계산한다.

- (2.5)–(2.6): $\rho_1$ 둘레에서 Euler 가지를 따라 계산한 윤곽 적분과 선분 위의 Gauss–Chebyshev 적분을 비교한다. $\alpha=\xi/L$이 $\pm5$까지일 때 $M$과 $\psi$ 모두 상대오차 $3\times10^{-10}$ 안에서 일치한다. $L=5$와 같이 국소 근사가 맞지 않는 범위에서도 일치하므로 (2.5)가 극한식이 아님을 확인한다.
- (2.7): $\xi=1,4$, $L=25$–$400$에서 나머지에 $L^2$을 곱한 값이 거의 일정하여 $O(L^{-2})$을 따른다. $\Re\kappa_{\rho_1}=-\vartheta_{\rm RS}'(\gamma_1)$도 $10^{-13}$ 안에서 성립한다.

이 수치 확인은 증명을 대체하지 않으며, 국소 핵은 하나의 영점에서만 확인했다.

    python -B scripts/check_mertens_psi.py
    python -B scripts/check_mertens_segment.py

<a id="references"></a>

## 참고 문헌

- A. Selberg, An elementary proof of the prime-number theorem, Ann. of Math. (2) 50 (1949), 305–313: $\Lambda_2=\Lambda\log+\Lambda*\Lambda$와 대칭 공식 $\sum_{n\le x}\Lambda_2(n)=2x\log x+O(x)$.
- H. Davenport, *Multiplicative Number Theory*, 3rd ed., GTM 74 (2000), 17장: $\psi(x)$의 명시공식.
- H. Bohr, B. Jessen, Über die Werteverteilung der Riemannschen Zetafunktion, Acta Math. 54 (1930), 1–35; 58 (1932), 1–55: $\sigma>\tfrac12$에서 $\log\zeta(\sigma+it)$의 극한분포.
- A. P. Guinand, A summation formula in the theory of prime numbers, Proc. London Math. Soc. (2) 50 (1948), 107–119: 영점합과 소수합의 합산 공식.
- A. Weil, Sur les « formules explicites » de la théorie des nombres premiers, Comm. Sém. Math. Univ. Lund (1952), 252–265: 명시공식의 시험함수 형태.
- H. Iwaniec, E. Kowalski, *Analytic Number Theory*, AMS Colloq. Publ. 53 (2004), §5.5: (4.3)의 형태와 시험함수 조건.
- E. C. Titchmarsh, *The Theory of the Riemann Zeta-Function*, 2nd ed. (1986), §9.3: Riemann–von Mangoldt 공식 $N(t)=\theta(t)/\pi+1+S(t)$.
