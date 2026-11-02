# 약수 구조를 통한 FTA의 영점 운송

**2026-11-03 (전역 공식과 가우스 평활화)**

여기서 FTA는 산술의 기본정리, 즉 소인수분해의 유일성이다. [주정리 문서 §7](main-theorem.md#section-7)의 이산핵을 소수별 유한 합성곱으로 쓰고, 그 생성함수를 Hadamard 영점 곱으로 옮긴 뒤 Perron 역변환으로 되돌린다.

이 문서에서 구성하는 전역 공식은 고전적인 Euler 곱과 Hadamard 곱을 이산핵의 변환에 적용한 것이다. 소수와 영점을 일대일로 대응시키는 주장은 하지 않는다. 영점 자료에는 위치와 중복도, 감마 인자와 정규화 상수를 함께 포함한다.

- [§11.1. 소수별 약수 선택과 FTA](#section-11-1)
- [§11.2. 정규화한 영점 곱](#section-11-2)
- [§11.3. 전역 재구성과 절단오차](#section-11-3)
- [§11.4. 유한 높이의 특이점별 분해](#section-11-4)
- [§11.5. 국소 베타핵과 전역 공식의 관계](#section-11-5)
- [§11.6. 복원의 범위와 후속 과제](#section-11-6)
- [§12.1. 로그 변수의 가우스 평활화](#section-12-1)
- [§12.2. 두 제타 인자에 공통인 절단 높이](#section-12-2)
- [§12.3. 평활화된 유한 가법 공식](#section-12-3)
- [§12.4. 로그폭 선택과 작은 나머지](#section-12-4)
- [§12.5. 베타핵의 모양과 평활화 제거의 한계](#section-12-5)

정리 3–4의 영점 곱 및 재구성과 §12의 평활화된 유한 가법 공식에는 RH나 단순영점 가정을 쓰지 않는다. §12에서는 정확한 극점 기여를 유지하면서 나머지를 $O(x^{-1})$로 줄인다. 평활화하지 않은 합의 **무한 가법 전개**와 그 나머지는 별도의 문제로 구분한다.

<a id="section-11"></a>

## 11. 전역 공식

<a id="section-11-1"></a>

### 11.1. 소수별 약수 선택과 FTA

$n=\prod_{p\mid n}p^{a_p}\ge2$라 하자. 각 소수 거듭제곱에 로그 약수의 유한 복소측도

$$
\lambda_{p,a}^{z,w}
=\sum_{j=0}^{a}
\frac{(z)_j(w)_{a-j}}{j!(a-j)!}\,\delta_{j\log p}
$$

를 붙인다. $T_c(v)=cv$, $*$를 실수축 위 측도의 가법 합성곱이라 쓰면

<a id="equation-11-1"></a>

$$
\boxed{
\mathscr K_n^{z,w}
=(T_{1/\log n})_*
\left(\mathop{*}_{p\mid n}\lambda_{p,a_p}^{z,w}\right).
}
\tag{11.1}
$$

**증명.** FTA에 의해 약수 $d\mid n$은 지수들의 유일한 선택 $(j_p)_{p\mid n}$, $0\le j_p\le a_p$와 대응한다. 합성곱의 원자 위치는 $\sum_pj_p\log p=\log d$, 계수는 $d_z(d)d_w(n/d)$이다. 마지막으로 $\log n$으로 나누면 이산핵의 정의를 얻는다. $\square$

따라서 FTA가 주는 곱셈적 구조는 로그 약수 좌표의 합성곱 구조로 보존된다. 모든 $d,m$에 대해 생성함수를 만들면 $\Re s,\Re t>1$에서

<a id="equation-11-2"></a>

$$
\sum_{d,m\ge1}\frac{d_z(d)d_w(m)}{d^sm^t}
=\prod_p(1-p^{-s})^{-z}(1-p^{-t})^{-w}
=\zeta(s)^z\zeta(t)^w.
\tag{11.2}
$$

이 단계까지는 소수별 원자 선택의 정확한 재표현이다. 다음 단계에서 오른쪽을 영점 자료로 바꾼다.

<a id="section-11-2"></a>

### 11.2. 정규화한 영점 곱

$\rho$는 서로 다른 비자명 영점을, $m_\rho$는 그 중복도를 뜻한다. 양·음의 허수부를 모두 포함한다. $\gamma_E$를 오일러 상수라 쓰고

$$
b_0=\log(2\pi)-1-\frac{\gamma_E}{2}
$$

로 둔다. 고전적인 Hadamard 곱은

<a id="equation-11-3"></a>

$$
\zeta(s)
=\frac{e^{b_0s}}{2(s-1)\Gamma(1+s/2)}
\prod_\rho
\left[\left(1-\frac{s}{\rho}\right)e^{s/\rho}\right]^{m_\rho}.
\tag{11.3}
$$

이는 완성 제타함수의 정준 곱과 동치이며, 감마 인자가 자명한 영점의 역할을 포함한다. 식과 정규화는 [DLMF 25.2.12](https://dlmf.nist.gov/25.2.E12)를 따른다.

분수 거듭제곱의 가지를 고정하기 위해 기준점을 $2$로 잡는다. 반평면 $\Re v>1$에서

<a id="equation-11-4"></a>

$$
\begin{aligned}
\Psi_\rho(v)
&=\int_2^v\left(\frac1{q-\rho}+\frac1\rho\right)dq\\
&=\log(v-\rho)-\log(2-\rho)+\frac{v-2}{\rho},\\
\ell(v)
&=\log(v-1)+\log\Gamma(1+v/2)-b_0(v-2).
\end{aligned}
\tag{11.4}
$$

$v-\rho$와 $v-1$은 오른쪽 반평면에 있으므로 이 로그는 주가지로 정할 수 있다. $\log\Gamma$는 양의 실축에서 실수값을 갖는 해석적 가지이다. 따라서 $\Psi_\rho(2)=\ell(2)=0$이고, 적분은 경로에 무관하다.

**정리 3 (이산핵 생성함수의 전역 영점 곱).** $\Re s,\Re t>1$와 모든 고정된 $z,w\in\mathbb C$에 대해

<a id="equation-11-5"></a>

$$
\boxed{
\zeta(s)^z\zeta(t)^w
=\zeta(2)^{z+w}e^{-z\ell(s)-w\ell(t)}
\prod_\rho
\exp\!\left[m_\rho\{z\Psi_\rho(s)+w\Psi_\rho(t)\}\right].
}
\tag{11.5}
$$

오른쪽 곱의 로그는 이 반평면들의 콤팩트 집합에서 절대·균등 수렴한다. 왼쪽의 거듭제곱은 Euler 곱의 로그를 사용한다.

**증명.** 영점 개수의 고전적 상계 $N(V)\ll V\log(2V)$로부터 $\sum_\rho m_\rho|\rho|^{-2}<\infty$이다. 고정된 콤팩트 집합에서

$$
\frac1{q-\rho}+\frac1\rho=\frac{q}{\rho(q-\rho)}
=O(|\rho|^{-2})
$$

이므로 $\sum_\rho m_\rho\Psi_\rho(v)$가 정상수렴한다. (11.3)의 로그미분을 $2$에서 $v$까지 적분하면

<a id="equation-11-6"></a>

$$
\log\zeta(v)=\log\zeta(2)-\ell(v)
+\sum_\rho m_\rho\Psi_\rho(v).
\tag{11.6}
$$

를 얻는다. 양쪽은 $v=2$에서 일치하는 해석적 로그이므로 가지의 모호성이 없다. $s,t$에서 각각 $z,w$를 곱하고 지수함수를 취하면 (11.5)가 따른다. $\square$

$N(V)$는 양의 높이의 영점을 중복도와 함께 세는 함수이다. 필요한 상계는 Riemann–von Mangoldt 공식보다 약하다. 예를 들어 [Hasanalizade–Shen–Wong, *Counting zeros of the Riemann zeta function*](https://arxiv.org/abs/2107.06506)의 영점 개수 추정으로도 충분하다.

<a id="section-11-3"></a>

### 11.3. 전역 재구성과 절단오차

이제 $z=w=-1/2$에 집중한다. $x\ge3$은 정수가 아니며, $L=\log x$, $\alpha=\xi/L$로 둔다.

$$
\mathfrak M_x(\xi)
=\sum_{dm\le x}d_{-1/2}(d)d_{-1/2}(m)e^{i\xi\log d/L},
\qquad
F_\alpha(s)=\zeta(s-i\alpha)^{-1/2}\zeta(s)^{-1/2}.
$$

이는 [식 (7.16)](main-theorem.md#equation-7-16)의 변환이다. $\mathfrak M_x(0)=M(x)$이며, $\xi$를 보존하면 로그 약수의 이산측도 $\mathscr Q_x$가 보존된다.

정리 3을 이동 대각선에 적용하여

<a id="equation-11-7"></a>

$$
\begin{aligned}
\mathcal G_\alpha(s)
&=\frac1{\zeta(2)}
\exp\!\left[\frac{\ell(s-i\alpha)+\ell(s)}2\right],\\
\mathcal Z_{\rho,\alpha}(s)
&=\exp\!\left[-\frac{m_\rho}{2}
\{\Psi_\rho(s-i\alpha)+\Psi_\rho(s)\}\right],\\
F_\alpha(s)&=\mathcal G_\alpha(s)\prod_\rho\mathcal Z_{\rho,\alpha}(s)
\qquad(\Re s>1)
\end{aligned}
\tag{11.7}
$$

로 쓴다. $\mathcal G_\alpha$에는 극점·감마 인자·정규화 상수가 들어가고, 각 $\mathcal Z_{\rho,\alpha}$에는 하나의 비자명 영점과 그 중복도가 들어간다.

Perron 적분의 높이 $U$와 영점 곱의 절단 높이 $V$는 서로 다른 매개변수로 둔다.

<a id="equation-11-8"></a>

$$
\begin{aligned}
F_{\alpha,V}(s)
&=\mathcal G_\alpha(s)
\prod_{|\Im\rho|\le V}\mathcal Z_{\rho,\alpha}(s),\\
\mathfrak M_{x;U,V}(\xi)
&=\frac1{2\pi i}\int_{c-iU}^{c+iU}
\frac{x^s}{s}F_{\alpha,V}(s)\,ds,
\qquad c=1+\frac1L.
\end{aligned}
\tag{11.8}
$$

**정리 4 (유한 영점 자료에 의한 재구성).** $U\ge2$, $H=U+|\alpha|+3$, $V\ge2H$라 하자. 어떤 절대상수 $C>0$에 대해

$$
\varepsilon_{U,V}=C\,H^2\frac{\log(2V)}V,
\qquad
\delta_x=\operatorname{dist}(x,\mathbb Z)>0
$$

로 두면

<a id="equation-11-9"></a>

$$
\boxed{
\left|\mathfrak M_x(\xi)-\mathfrak M_{x;U,V}(\xi)\right|
\ll
\frac{x\log(2x)}U
+\min\left(1,\frac{x}{U\delta_x}\right)
+xL\log(2+U)\left(e^{\varepsilon_{U,V}}-1\right).
}
\tag{11.9}
$$

상수는 절대상수이며, 영점은 임계선 밖의 영점도 포함하여 모두 세어야 한다. 이 상계는 정확한 영점 위치로 정의한 (11.8)에 대한 것으로, 수치 영점과 수치 적분 자체의 오차는 별도로 다룬다.

**증명.** 먼저 $|\Im s|\le U$, $\Re s=c$이고 $v=s$ 또는 $s-i\alpha$라 하자. $2$에서 $v$까지의 선분에서는 $|q|\le H$이다. $|\Im\rho|>V\ge2H$이면

$$
|\Psi_\rho(v)|
\le\frac{2H\,|v-2|}{|\rho|^2}
\ll\frac{H^2}{|\rho|^2}.
$$

영점 개수 상계에 부분합산을 적용하면

$$
\sum_{|\Im\rho|>V}\frac{m_\rho}{|\rho|^2}
\ll\frac{\log(2V)}V.
$$

따라서 곱의 꼬리 로그를 직접 사용하여

<a id="equation-11-10"></a>

$$
\left|\log\frac{F_\alpha(s)}{F_{\alpha,V}(s)}\right|
\le\varepsilon_{U,V},
\qquad
|F_\alpha(s)-F_{\alpha,V}(s)|
\le |F_\alpha(s)|(e^{\varepsilon_{U,V}}-1)
\tag{11.10}
$$

을 얻는다. 이 로그는 꼬리 합으로 정한 로그이다.

다음으로 [식 (7.8)](main-theorem.md#equation-7-8)을 더 날카롭게 쓰면

<a id="equation-11-11"></a>

$$
\left|
\sum_{d\mid n}d_{-1/2}(d)d_{-1/2}(n/d)d^{i\alpha}
\right|
\le d_1(n)=1.
\tag{11.11}
$$

그러므로 $\Re s=c$에서 $|F_\alpha(s)|\le\zeta(c)\ll L$이다. (11.10)을 적분하면 $x^c=ex$와 $\int_{-U}^U|c+it|^{-1}dt\ll\log(2+U)$에 의해 (11.9)의 마지막 항을 얻는다.

남은 부분은 $F_\alpha$ 자체의 Perron 절단오차이다. (11.11)에 의해

$$
E_{\mathrm P}(x,U)
\ll\sum_{n\ge1}
\left(\frac xn\right)^c
\min\left(1,\frac1{U|\log(x/n)|}\right)
\ll\frac{x\log(2x)}U
+\min\left(1,\frac{x}{U\delta_x}\right).
$$

$n\notin[x/2,2x]$는 절대수렴으로 처리한다. 중앙에서는 $|\log(x/n)|\asymp|x-n|/x$를 쓰고, 가장 가까운 정수 항을 분리한 뒤 나머지에 조화합 상계를 적용한다. 두 오차를 합하면 정리를 얻는다. $\square$

**따름정리 4.1 (전역 영점 곱을 통한 정확한 이산합).** 고정된 비정수 $x\ge3$와 실수 $\xi$에 대해

<a id="equation-11-12"></a>

$$
\boxed{
\mathfrak M_x(\xi)
=\lim_{U\to\infty}
\frac1{2\pi i}\int_{c-iU}^{c+iU}
\frac{x^s}{s}\,
\mathcal G_{\xi/L}(s)
\prod_{|\Im\rho|\le U^3}\mathcal Z_{\rho,\xi/L}(s)\,ds.
}
\tag{11.12}
$$

실제로 $V=U^3$이면 고정된 $x,\xi$에서 $\varepsilon_{U,V}\ll\log U/U$이고, (11.9)의 모든 항이 $0$으로 간다. 임의의 동시 극한을 가정하지 않고, 적분 높이와 영점 절단의 관계를 지정한 것이다.

같은 구조는 임의의 고정된 복소 $z,w$에도 성립한다. 정리 3에서 $s$를 $s-i\alpha$, $t$를 $s$로 두고 [명제 7.3](main-theorem.md#section-7-5)을 적용하면, $V=U^3$인 극한으로 $\widetilde S_x^{z,w}(\xi)$를 얻는다. 고정된 $x$에서는 계수의 절대수렴으로 Perron 꼬리가 사라지고, 유한 적분에서의 곱 교체 오차는 $O_{x,z,w,\xi}((\log U)^2/U)$이다. (11.9)의 간단한 균등 상계는 음의 대칭 지수에서 계수의 절댓값이 $1$ 이하라는 점을 사용한다.

<a id="section-11-4"></a>

### 11.4. 유한 높이의 특이점별 분해

(11.12)는 모든 영점의 정보를 사용하는 곱 공식이다. 영점별 기여를 더하는 표현은 별도로 경로를 변형하여 얻는다.

$|\alpha|\le A$, $T>|\alpha|+2$라 하자. 직사각형 $-1\le\Re s\le c$, $|\Im s|\le T$의 경계가 특이점을 지나지 않고, 각 쌍 $\{\rho,\rho+i\alpha\}$이 모두 안쪽 또는 모두 바깥쪽에 있도록 $T$를 택한다. 예를 들어 경계 높이 $\pm T$와 모든 영점 높이의 거리가 $|\alpha|$보다 크면 충분하다. 고정된 유한 $T$를 영점 높이와 다르게 택하면, 충분히 작은 $|\alpha|$에서 이 조건이 성립한다.

직사각형 안에서 $1$과 $1+i\alpha$, 그리고 각 $\rho$와 $\rho+i\alpha$를 선분으로 연결한다. 겹치는 영점 선분은 하나의 연결 성분 $C$로 묶는다. $\mathcal C_T$를 이 영점 성분들의 유한 집합이라 한다. 각 성분에는 원래 영점의 완전한 쌍들이 들어 있다.

$F_\alpha(s)^2=1/(\zeta(s-i\alpha)\zeta(s))$의 영점·극점 차수는 각 성분을 합쳐 세면 짝수이다. 따라서 이 선분들을 제거한 영역에서는 오른쪽 경계의 Euler 가지와 일치하는 제곱근 $F_\alpha$를 선택할 수 있다. $\alpha=0$에서는 $F_0=1/\zeta$로 해석한다.

$f_{x,\alpha}(s)=x^sF_\alpha(s)/s$로 두고, 각 성분을 감싸는 양의 방향의 작은 폐곡선을 택하여

<a id="equation-11-13"></a>

$$
\mathcal P_1(x,\xi)
=\frac1{2\pi i}\int_{\Gamma_1}f_{x,\alpha}(s)\,ds,
\qquad
\mathcal I_C(x,\xi)
=\frac1{2\pi i}\int_{\Gamma_C}f_{x,\alpha}(s)\,ds
\tag{11.13}
$$

로 둔다. $\Gamma_1$은 $[1,1+i\alpha]$만 감싸고, $\Gamma_C$는 영점 성분 $C$만 감싼다. 단순성을 가정하지 않으며, 서로 가까운 영점들을 임의로 분리하지 않는다.

**명제 11.1 (유한 높이의 가법 공식).** 위 조건에서

<a id="equation-11-14"></a>

$$
\boxed{
\mathfrak M_x(\xi)
=\mathcal P_1(x,\xi)
+\sum_{C\in\mathcal C_T}\mathcal I_C(x,\xi)
+F_\alpha(0)
+\mathcal R_T(x,\xi)
+O(E_{\mathrm P}(x,T)),
}
\tag{11.14}
$$

여기서 나머지 윤곽 적분은

<a id="equation-11-15"></a>

$$
\mathcal R_T(x,\xi)
=\frac1{2\pi i}
\left(
\int_{c-iT}^{-1-iT}
+\int_{-1-iT}^{-1+iT}
+\int_{-1+iT}^{c+iT}
\right)
f_{x,\alpha}(s)\,ds.
\tag{11.15}
$$

이다. $F_\alpha(0)$은 위에서 이어온 가지의 값으로, $1/s$가 만드는 유수이다.

**증명.** 유효 Perron 공식의 오른쪽 수직선 적분에서 (11.15)의 세 변을 빼면 직사각형의 양의 방향 경계가 된다. 절단선들을 감싸는 내부 경계를 제거한 영역에 Cauchy 정리를 적용한다. 내부 폐곡선의 기여가 (11.13)이고, 남은 고립극 $s=0$의 유수가 $F_\alpha(0)$이다. 직선 $\Re s=-1$의 오른쪽에는 자명한 영점이 없으므로 이들이 목록을 완성한다. $\square$

경계값

$$
\mathcal H_T(\alpha)
=\sup_{\substack{-1\le\sigma\le c\\\epsilon=\pm1}}
|F_\alpha(\sigma+i\epsilon T)|
$$

를 쓰면 나머지는

<a id="equation-11-16"></a>

$$
|\mathcal R_T(x,\xi)|
\ll_A x^{-1}
+\frac{x}{T\log x}\mathcal H_T(\alpha)
\tag{11.16}
$$

로 제어된다. 수평선에서는 $|s|\ge T$와 $\int_{-1}^c x^\sigma d\sigma\le x^c/L$을 사용한다. 왼쪽 수직선에서는 함수방정식과 $\Re(1-s)=2$의 절대수렴으로 $|\zeta(-1+it)|^{-1}\ll(1+|t|)^{-3/2}$이다. $|\alpha|\le A$에서 $|F_\alpha(-1+it)|/|-1+it|$의 적분은 균등 유계이므로 첫 항을 얻는다.

여기서는 $\mathcal H_T(\alpha)$를 작은 거듭제곱으로 가정하지 않았다. 따라서 (11.14)는 정량적인 경계 나머지를 포함한 유한 공식이며, 그 나머지가 영점 신호보다 작다는 주장은 추가 추정이 필요하다.

<a id="section-11-5"></a>

### 11.5. 국소 베타핵과 전역 공식의 관계

고정된 영점 $\rho$를 택하고, $\alpha=\xi/L\to0$일 때 그 쌍이 다른 특이점과 분리되도록 한다. $\mathcal I_\rho$는 해당 성분의 (11.13)이다.

**명제 11.2 (영점의 국소 반사핵).** $\rho$의 중복도가 $m$이고 $b_\rho=\zeta^{(m)}(\rho)/m!$라 하자. 고정된 $R$에 대해 $|\xi|\le R$에서, $\alpha=0$의 $1/\zeta$에 연속적으로 맞춘 가지로

<a id="equation-11-17"></a>

$$
\mathcal I_\rho(x,\xi)
=\frac{x^\rho L^{m-1}}{\rho b_\rho}\,
\mathcal B_{m/2,m/2}(\xi)
+O_{\rho,R}\!\left(x^{\Re\rho}L^{m-2}\right).
\tag{11.17}
$$

특히 단순 영점에서는

<a id="equation-11-18"></a>

$$
\boxed{
\mathcal I_\rho(x,\xi)
=\frac{x^\rho}{\rho\zeta'(\rho)}
e^{i\xi/2}J_0(\xi/2)
+O_{\rho,R}\!\left(\frac{x^{\Re\rho}}L\right).
}
\tag{11.18}
$$

**증명.** $\zeta(\rho+h)=b_\rho h^m(1+O_\rho(h))$를 두 인자에 적용하고 $s=\rho+q/L$로 바꾼다. 두 분기점을 함께 감싸는 유계한 $q$ 경로에서

$$
F_{\xi/L}(\rho+q/L)
=\frac{L^m}{b_\rho}
(q-i\xi)^{-m/2}q^{-m/2}
\left(1+O_{\rho,R}(L^{-1})\right).
$$

또한 $x^s=x^\rho e^q$, $ds=dq/L$, $1/s=1/\rho+O_\rho(L^{-1})$이다. [주정리 문서 (4.4)](main-theorem.md#equation-4-4)의 역라플라스 항등식으로 (11.17)을 얻는다. 두 지수의 합이 정수 $m$이므로 함께 감싸는 폐곡선은 그 항등식의 경로로 변형할 수 있다. $m=1$에서 $\mathcal B_{1/2,1/2}(\xi)=e^{i\xi/2}J_0(\xi/2)$를 쓰면 (11.18)이다. $\square$

극점 쌍에서는 같은 방법으로

<a id="equation-11-19"></a>

$$
\mathcal P_1(x,\xi)
=\frac{x}{L^2}\mathcal B_{-1/2,-1/2}(\xi)
+O_R\!\left(\frac{x}{L^3}\right)
\tag{11.19}
$$

를 얻는다. 그러나 전역 가법 공식 (11.14)에서는 $\mathcal P_1$을 **정확한 윤곽 적분으로 유지**한다. 유한 차수의 로그 전개 나머지는 일반적으로 $x^{1/2}$보다 크기 때문이다.

$\xi=0$은 즉시 검산할 수 있다. $\mathcal P_1(x,0)=0$, $F_0(0)=1/\zeta(0)=-2$이고, 단순 영점의 정확한 폐곡선 적분은 $x^\rho/(\rho\zeta'(\rho))$이다. 중복 영점에서는 $x^s/(s\zeta(s))$의 고차 유수가 된다. 따라서 (11.14)는 고전적인 메르텐스 함수의 절단 명시적 공식으로 돌아간다.

고정된 영점에서의 (11.17)을 무한히 합하는 것은 허용되지 않는다. $T$가 커지면 오차상수의 영점 의존성, 영점 간격, 겹치는 절단선 성분을 함께 제어해야 한다. 모든 단순 영점의 선도 모양이 같은 아크사인 핵이라는 사실도, 그 선도항들만으로 약수별 원자를 복원할 수 있다는 뜻은 아니다.

<a id="section-11-6"></a>

### 11.6. 복원의 범위와 후속 과제

현재 구성의 정보 흐름은 다음과 같다.

$$
\{p^{a_p}\}_{p\mid n}
\ \longrightarrow\
\mathop{*}_{p\mid n}\lambda_{p,a_p}^{z,w}
\ \longrightarrow\
\mathscr K_n^{z,w}
\ \longrightarrow\
\zeta(s)^z\zeta(t)^w
\ \longrightarrow\
\{\rho,m_\rho;\ell,\zeta(2)\}.
$$

정리 3의 정규화한 곱과 따름정리 4.1의 역변환을 사용하면 $\mathscr Q_x$를 되찾을 수 있다. 각 $x$에서 모든 실수 $\xi$의 푸리에 변환을 사용하고, 이어 $x$가 정수 $n$을 지나는 약한 도약을 취하면 $\mathscr K_n$을 얻는다. 이 복원은 [주정리 문서 §7.5](main-theorem.md#section-7-5)의 좌표 변경과 도약 항등식에 따른다.

음의 대칭 지수에서는 $(-1/2)_a\ne0$이므로 모든 약수 원자의 계수 $d_{-1/2}(d)d_{-1/2}(n/d)$가 영이 아니다. 서로 다른 약수는 서로 다른 위치에 놓인다. 따라서 복원된 핵의 지지집합에서 정확히

<a id="equation-11-20"></a>

$$
\{d\in\mathbb N:d\mid n\}
=\left\{n^u:u\in\operatorname{supp}\mathscr K_n^{-1/2,-1/2}\right\}
\tag{11.20}
$$

를 얻는다. 이 약수 집합에서 소수 약수와 각 소수의 최대 지수를 읽으면 원래의 소인수분해가 결정된다. 즉 전체 변환과 역변환을 통해 약수 구조를 잃지 않는다는 의미에서 FTA를 영점 표현까지 운송한다.

이는 전체 영점 자료와 정규화 자료가 전체 이산핵 계열을 결정한다는 진술이다. 한정된 주파수·유한 영점 자료만으로 모든 원자를 정확히 복원한다는 진술은 아니다. 유한 절단의 오차는 (11.9)로 분리된다.

현재 확보한 결과와 다음 과제는 다음과 같이 구분한다.

| 수준 | 결과 | 남는 조건 또는 과제 |
|:--|:--|:--|
| 정확한 전역 곱 | (11.5), (11.12) | 모든 영점과 중복도를 포함하며, RH 불필요 |
| 유한 영점 자료 | (11.9) | 영점 위치의 수치 오차·적분 오차는 별도 |
| 유한 가법 공식 | (11.14) | 경계가 완전한 특이점 쌍들을 분리하는 조건 |
| 영점별 국소핵 | (11.17)–(11.18) | 고정된 영점, 유계한 $\xi$ |
| 평활화된 유한 가법 공식 | [따름정리 5.1](#section-12-4) | 로그폭 $h=L^{-\kappa}$, 정확한 성분 적분, 나머지 $O(x^{-1})$ |
| 평활화 없는 무한 가법 전개 | 후속 과제 | 평활화 제거, 성분의 합산, 역미분값 및 주파수 균등성 |

아래 §12는 (11.14)에 가우스 평활화를 적용하고 절단 높이를 선택하여, 정확한 극점 기여와 유한 영점 성분합 밖의 나머지를 제어한다. 영점 성분들은 정확한 적분으로 유지하므로 개별 역미분값의 상계를 요구하지 않는다.

<a id="section-12"></a>

## 12. 가우스 평활화와 작은 경계 나머지

<a id="section-12-1"></a>

### 12.1. 로그 변수의 가우스 평활화

$x\ge3$, $L=\log x$, $\xi\in\mathbb R$, $\alpha=\xi/L$로 두고, $\Phi$를 표준 정규분포의 누적분포함수라 하자. $h>0$에 대하여

<a id="equation-12-1"></a>

$$
\begin{aligned}
b_\alpha(n)
&=\sum_{d\mid n}d_{-1/2}(d)d_{-1/2}(n/d)e^{i\alpha\log d},\\
\mathfrak M_{x,h}(\xi)
&=\sum_{n\ge1}b_\alpha(n)
  \Phi\!\left(\frac{\log(x/n)}h\right)
\end{aligned}
\tag{12.1}
$$

로 정의한다. [계수 상계 (11.11)](#equation-11-11)에 의해 $|b_\alpha(n)|\le1$이고, 가우스 꼬리 때문에 이 합은 절대수렴한다.

평활화는 $n\le x$라는 경계에 적용한다. $\alpha=\xi/\log x$는 중심 $x$에서 고정하며, 평활화 중에 변하는 절단점의 로그로 다시 나누지 않는다. 따라서 $n>x$인 항도 들어가고, 약수의 좌표 $\log d/L$가 $1$을 넘을 수도 있다.

가우스 푸리에 역변환으로, 임의의 $c>0$에서

$$
\Phi(y/h)=\frac1{2\pi i}\int_{c-i\infty}^{c+i\infty}
       \frac{e^{sy+h^2s^2/2}}s\,ds
$$

이다. 실제로 양변을 $y$로 미분하면 정규밀도가 되고, $y\to-\infty$에서 모두 $0$이다. 이를 $c>1$에서 절대수렴하는 Dirichlet 급수에 적용하면 정확히

<a id="equation-12-2"></a>

$$
\boxed{
\mathfrak M_{x,h}(\xi)
=\frac1{2\pi i}\int_{c-i\infty}^{c+i\infty}
\frac{x^s e^{h^2s^2/2}}s F_\alpha(s)\,ds
}
\tag{12.2}
$$

를 얻는다. 이제 수직 적분은 절대수렴하며, 정수 $x$에서도 별도의 끝점 규약이 필요 없다. $n=x$의 가중치는 정의상 $1/2$이다.

<a id="section-12-2"></a>

### 12.2. 두 제타 인자에 공통인 절단 높이

**보조정리 12.1 (영점 쌍을 가르지 않는 높이).**
절대상수 $\eta,C>0$가 존재하여, 충분히 큰 $Y$마다 $T\in[Y,Y+1]$를 다음과 같이 고를 수 있다.

<a id="equation-12-3"></a>

$$
\begin{aligned}
\operatorname{dist}\bigl(T,\{\Im\rho\}\bigr)
&\ge\frac{\eta}{\log(Y+3)},\\
H_T(\alpha)
&\le \exp\!\bigl(C\log^2(Y+3)\bigr)
\qquad
\left(|\alpha|\le\frac{\eta}{2\log(Y+3)}\right).
\end{aligned}
\tag{12.3}
$$

여기서 $H_T$는 [(11.16)](#equation-11-16)의 경계 최댓값이며 $1<c\le2$에 균등하다. 이 하나의 $T$는 표시한 범위의 모든 $\alpha$에 공통으로 사용할 수 있다. 특히 $\rho,\rho+i\alpha$의 쌍은 직사각형 안에 함께 들어가거나 함께 밖에 놓인다.

**증명.** 길이가 고정된 높이 구간의 영점 개수는 중복도를 포함해 $O(\log Y)$이다. $[Y-1,Y+2]$의 각 영점 높이에서 반지름 $\eta/\log(Y+3)$인 구간을 빼자. $\eta$를 충분히 작게 고르면 제거되는 총길이가 $1/2$보다 작으므로 첫 조건을 만족하는 $T$가 남는다. 켤레대칭에 의해 $-T$에서도 같은 거리가 확보된다.

$t=\pm T$ 또는 $t=\pm T-\alpha$이면 모든 영점 높이와의 거리가 $\eta/(2\log(Y+3))$ 이상이다. 표준 국소식

$$
\frac{\zeta'}{\zeta}(\sigma+it)
=\sum_{|\Im\rho-t|<1}\frac{m_\rho}{\sigma+it-\rho}
 +O(\log|t|),\qquad -1\le\sigma\le2
$$

에 영점 개수와 거리 상계를 넣으면 $\zeta'/\zeta\ll\log^2(Y+3)$이다. 이 식을 $\sigma=2$에서 수평으로 적분하고 $|1/\zeta(2+it)|\le\zeta(2)$를 쓰면
$|1/\zeta(\sigma+it)|\le\exp(C\log^2(Y+3))$를 얻는다. 두 인자의 절댓값을 곱하여 제곱근을 취하면 둘째 조건이다. 사용한 영점 개수와 로그미분 식은 [Kedlaya, 보조정리 9.4와 9.6](https://kskedlaya.org/ant/chap-von-mangoldt.html)에 정리되어 있다. $\square$

<a id="section-12-3"></a>

### 12.3. 평활화된 유한 가법 공식

$c=1+1/L$, $0<h\le1$, $|\alpha|\le A$라 하자. $T>|\alpha|+2$는 [§11.4](#section-11-4)의 완전한 쌍 조건을 만족한다고 가정한다. 절단선의 연결성분 $\mathcal C_T$와 양의 방향 윤곽 $\Gamma_1,\Gamma_C$도 그 절의 정의를 따른다. 가까운 영점들의 절단선이 겹치면 하나의 성분으로 취급한다.

<a id="equation-12-4"></a>

$$
\begin{aligned}
f_{x,\alpha,h}(s)&=\frac{x^s e^{h^2s^2/2}}sF_\alpha(s),\\
\mathcal P_{1,h}(x,\xi)&=\frac1{2\pi i}\int_{\Gamma_1}f_{x,\alpha,h}(s)\,ds,\\
\mathcal I_{C,h}(x,\xi)&=\frac1{2\pi i}\int_{\Gamma_C}f_{x,\alpha,h}(s)\,ds.
\end{aligned}
\tag{12.4}
$$

**정리 5 (평활화된 유한 가법 공식).**
위 조건에서

<a id="equation-12-5"></a>

$$
\mathfrak M_{x,h}(\xi)
=\mathcal P_{1,h}(x,\xi)
 +\sum_{C\in\mathcal C_T}\mathcal I_{C,h}(x,\xi)
 +F_\alpha(0)+\mathcal E_{x,h,T}(\xi)
\tag{12.5}
$$

이고,

<a id="equation-12-6"></a>

$$
\boxed{
|\mathcal E_{x,h,T}(\xi)|
\ll_A x^{-1}
 +\frac{x e^{-h^2T^2/2}}{TL}H_T(\alpha)
 +\frac{xL e^{-h^2T^2/2}}{h^2T^2}.
}
\tag{12.6}
$$

$\mathcal P_{1,h}$와 $\mathcal I_{C,h}$는 근사 베타핵이 아닌 정확한 윤곽 적분이다. $F_\alpha(0)$ 역시 오른쪽 반평면의 가지를 연속시킨 값이며, $\alpha=0$일 때 $-2$이다.

**증명.** (12.2)를 높이 $T$에서 자른 뒤, [(11.14)](#equation-11-14)와 같은 윤곽 변형을 적용한다. 가우스 인자는 전해석적이고 $s=0$에서 $1$이므로 원점의 유수는 여전히 $F_\alpha(0)$이다. 나머지는 오른쪽 수직선의 두 꼬리와

$$
c-iT\ \longrightarrow\ -1-iT\
\longrightarrow\ -1+iT\ \longrightarrow\ c+iT
$$

의 세 변 적분이다. 각 항을 다음과 같이 평가한다.

- 왼쪽 변에서는 함수방정식으로 $|F_\alpha(-1+it)|\ll_A(1+|t|)^{-3/2}$이다. 또한 $|e^{h^2(-1+it)^2/2}|\le e^{1/2}$이므로 $1/s$를 포함한 적분은 $O_A(x^{-1})$이다.
- 수평변에서는 $|e^{h^2s^2/2}|\le e^2e^{-h^2T^2/2}$, $|s|\ge T$이고, $\int_{-1}^c x^\sigma\,d\sigma\ll x/L$이다. 따라서 두 수평변은 (12.6)의 둘째 항으로 제어된다.
- 오른쪽 꼬리에서는 $|F_\alpha(c+it)|\le\zeta(c)\ll L$과
  $\int_T^\infty e^{-h^2t^2/2}\,dt/t\le e^{-h^2T^2/2}/(h^2T^2)$를 쓴다. $x^c=ex$이므로 셋째 항을 얻는다.

날카로운 Perron 절단에서 생기던 $x\log x/T$ 항이나 정수와의 거리 항은 여기에는 없다. $\square$

<a id="section-12-4"></a>

### 12.4. 로그폭 선택과 작은 나머지

**따름정리 5.1 (로그폭 평활화에서 $O(x^{-1})$ 나머지).**
$R>0$, $\kappa\ge0$를 고정한다. 충분히 큰 $x$에 대하여

<a id="equation-12-7"></a>

$$
h=L^{-\kappa},\qquad
Y=\frac{3\sqrt L}{h}=3L^{\kappa+1/2},\qquad
T_x\in[Y,Y+1]
\tag{12.7}
$$

인 높이 $T_x$를 골라, 모든 $|\xi|\le R$에서 균등하게

<a id="equation-12-8"></a>

$$
\boxed{
\mathfrak M_{x,L^{-\kappa}}(\xi)
=\mathcal P_{1,L^{-\kappa}}(x,\xi)
 +F_{\xi/L}(0)
 +\sum_{C\in\mathcal C_{T_x}}
        \mathcal I_{C,L^{-\kappa}}(x,\xi)
 +O_{\kappa,R}(x^{-1})
}
\tag{12.8}
$$

을 얻는다. 영점은 모든 비자명 영점과 중복도를 포함하며, RH나 단순영점 가정은 필요 없다.

**증명.** $\log(Y+3)=O_\kappa(\log L)$이므로
$R/L\le\eta/(2\log(Y+3))$가 충분히 큰 $x$에서 성립한다. 보조정리 12.1로 $|\xi|\le R$ 전체에 공통인 $T_x$를 고른다. 그러면

$$
e^{-h^2T_x^2/2}\le x^{-9/2},\qquad
H_{T_x}(\xi/L)
\le\exp\!\bigl(O_\kappa((\log L)^2)\bigr)=x^{o(1)}.
$$

따라서 (12.6)의 수평변 항은 $O_{\kappa,R}(x^{-3})$이고, 오른쪽 꼬리는 $O(x^{-7/2})$이다. 왼쪽 변의 $O_{\kappa,R}(x^{-1})$를 더하면 결론이다. $\square$

예를 들어 $\kappa=1/2$이면 $h=1/\sqrt{\log x}$이고, $T_x$는 $[3\log x,3\log x+1]$에서 고른다. 이 공식에 들어가는 영점 수는 중복도를 세어
$O_\kappa(L^{\kappa+1/2}\log L)$이다. 정확한 극점 기여를 먼저 빼면, 평활화된 약수 통계의 나머지 부분을 이 유한 성분합과 원점 항으로 $O(x^{-1})$ 오차 안에서 표현한다.

이는 $x\to\infty$에서의 해석적 오차이다. 유한 정밀도의 영점 위치와 윤곽 적분을 사용한 계산 오차까지 $O(x^{-1})$이라고 주장하는 것은 아니다.

<a id="section-12-5"></a>

### 12.5. 베타핵의 모양과 평활화 제거의 한계

고정된 영점 $\rho$의 중복도를 $m$, $b_\rho=\zeta^{(m)}(\rho)/m!$라 하자. 그 영점 쌍이 고립되어 있을 때 [(11.17)](#equation-11-17)의 같은 국소 계산으로, $0<h\le1$, $|\xi|\le R$에서

<a id="equation-12-9"></a>

$$
\mathcal I_{\rho,h}(x,\xi)
=e^{h^2\rho^2/2}
  \frac{x^\rho L^{m-1}}{\rho b_\rho}
  \mathcal B_{m/2,m/2}(\xi)
 +O_{\rho,R}\!\left(x^{\Re\rho}L^{m-2}\right)
\tag{12.9}
$$

이다. $q=L(s-\rho)$로 두면
$e^{h^2s^2/2}=e^{h^2\rho^2/2}(1+O_{\rho,R}(L^{-1}))$이므로 이 오차는 $h$에도 균등하다. 단순 영점의 베타핵은 그대로 $e^{i\xi/2}J_0(\xi/2)$이며, $\rho=\beta+i\gamma$에 붙는 추가 인자의 절댓값은 $e^{h^2(\beta^2-\gamma^2)/2}$이다.

극점 쌍에서는

<a id="equation-12-10"></a>

$$
\mathcal P_{1,h}(x,\xi)
=e^{h^2/2}\frac{x}{L^2}\mathcal B_{-1/2,-1/2}(\xi)
 +O_R(x/L^3).
\tag{12.10}
$$

(12.8)의 $O(x^{-1})$를 보존하려면 이 선도항으로 $\mathcal P_{1,h}$를 대체하지 않는다. 또한 (12.9)의 영점별 오차를 $T_x$까지 합하는 데 필요한 균등성은 아직 주장하지 않는다. 가까운 영점의 분리는 역미분값 문제를 다시 가져오므로 (12.8)에서는 성분 적분을 유지한다.

$\xi=0$에서는 $\mathcal P_{1,h}(x,0)=0$이고, 단순 영점의 정확한 항은
$e^{h^2\rho^2/2}x^\rho/(\rho\zeta'(\rho))$이다. 따라서 이 경우는 가우스 평활화한 메르텐스 함수의 유한 명시적 공식으로 환원된다.

평활화를 제거하는 비용도 별도로 계산할 수 있다. $Z$가 표준 정규변수이면
$\Phi(\log(x/n)/h)=\mathbb P(n\le xe^{hZ})$이다. 절단점이 변해도 $b_\alpha(n)$는 고정하고 $|b_\alpha(n)|\le1$을 사용하면

<a id="equation-12-11"></a>

$$
\begin{aligned}
|\mathfrak M_{x,h}(\xi)-\mathfrak M_x(\xi)|
&\le\mathbb E\left|\lfloor xe^{hZ}\rfloor-\lfloor x\rfloor\right|\\
&\le x e^{h^2/2}\bigl(2\Phi(h)-1\bigr)+1
\ll xh+1,\qquad 0<h\le1.
\end{aligned}
\tag{12.11}
$$

$h=L^{-\kappa}$에서 이 상계는 $x^{1/2}$보다 크다. 따라서 (12.8)의 작은 나머지를 원래의 날카로운 합에 그대로 옮길 수는 없다. 다음에 해결할 문제는 평활화를 제거할 때의 상쇄, 또는 더 좁은 평활폭에서 영점 쌍을 함께 감싸는 윤곽과 그 경계 추정이다. 현재 단계에서 얻은 것은 **평활화된 약수 통계의 정확한 극점 기여와 유한 영점 성분합을 분리하는 무조건적 공식**이다.

<a id="checks"></a>

## 검산

[검산 스크립트](../scripts/check_global_formula.py)는 다음을 확인한다.

- 여러 정수와 유리수 지수에서 (11.1)의 소수별 합성곱과 직접 약수합을 정확한 유리수 연산으로 비교한다.
- 정규화점과 가지가 분명한 $\Re s>1$에서 유한 영점 곱을 직접 계산한 $F_\alpha(s)$와 비교한다.
- $\alpha=0$에서 직사각형 경계 적분과 $s=0$ 및 첫 영점 쌍의 유수를 비교하여 (11.14)–(11.15)의 부호를 검산한다.

[평활화 검산 스크립트](../scripts/check_smoothed_formula.py)는 $\alpha\ne0$에서 (12.1)의 직접 약수합과 (12.2)의 절단 적분을 가우스 꼬리 상계 안에서 비교한다. 또한 극점 쌍과 양·음의 첫 영점 쌍을 각각 감싼 적분으로 (12.5)의 가지와 부호를 확인한다.

이 수치 확인은 증명의 대체가 아니며, 사용한 낮은 영점 목록에서의 계산이다. 정리들의 영점 집합은 임계선에 국한되지 않는다.

<a id="references"></a>

## 참고 문헌

- NIST DLMF, [25.2.11](https://dlmf.nist.gov/25.2.E11), [25.2.12](https://dlmf.nist.gov/25.2.E12): Euler 곱과 Hadamard 곱의 정규화.
- NIST DLMF, [§25.4](https://dlmf.nist.gov/25.4), [§25.10](https://dlmf.nist.gov/25.10): 함수방정식과 영점의 위치.
- E. Hasanalizade, Q. Shen, P.-J. Wong, [Counting zeros of the Riemann zeta function](https://arxiv.org/abs/2107.06506), J. Number Theory 235 (2022), 219–241: 영점 개수 상계.
- N. Ng, [The distribution of the summatory function of the Möbius function](https://arxiv.org/abs/math/0310381): 메르텐스 함수의 절단 공식과 역미분값 모멘트의 역할.
- K. S. Kedlaya, [Notes on analytic number theory, Chapter 9](https://kskedlaya.org/ant/chap-von-mangoldt.html), 보조정리 9.4와 9.6: 짧은 높이 구간의 영점 개수와 로그미분의 국소식.
