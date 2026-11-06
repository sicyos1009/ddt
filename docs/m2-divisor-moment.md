# 양의 약수분포의 이차 모멘트

고정된 정수 $n\ge2$에서 $(1/2,1/2)$ 가중치를 정의하고, 약수의 정규화 로그 좌표에 대한 이차 모멘트 $M_2(n)$을 가중 로그제곱의 약수합으로 표현한다. 이어서 누적 잔차 $\mathcal E_2$에서 제거하는 전체 주항을 명시하고, 무한차 모멘트의 끝점합에서 얻는 잔차 $E$를 정의한다. 마지막으로 $\mathcal E_2$가 $x^{\theta+\varepsilon}$ 이하인 것과 $\Re s>\theta$에 $\zeta$의 **단순영점**이 없는 것이 동치임을 증명하고(정리 8.1), $E$의 제곱근 수준 상계를 $\mathcal E_2$로 전달하는 명제를 그 따름정리로 얻는다(따름정리 8.2). 점근 전개는 정확한 정의와 증명 뒤에 모아 둔다. 절·정리·수식 번호는 이 문서 안에서 독립적이다.

| 결과 | 내용 | 범위 |
|:--|:--|:--|
| [정리 8.1](#section-8-1) | $\mathcal E_2(x)=O(x^{\theta+\varepsilon})\ \forall\varepsilon\iff\Re s>\theta$에 단순영점 없음 | 가정 없음 |
| [따름정리 8.2](#section-8-1) | $E(x)=O(x^{1/2+\eta})\ \forall\eta\Rightarrow\mathcal E_2(x)=O(x^{1/2+\varepsilon})\ \forall\varepsilon$ | 전달 명제 |
| [따름정리 8.3](#section-8-1) | $\mathcal E_2$의 성장 지수 $=\Theta_{\rm simple}\ge\tfrac12$, $\Delta_2=\Omega_\pm(x^{1/2-\varepsilon})$ | 가정 없음 |
| [§8.8](#section-8-8) | 고정 차수 $k$의 필요조건과 충분조건 | 부분 결과 |

- [§1. 약수의 가중치](#section-1)
- [§2. 로그 좌표와 약수쌍의 대칭성](#section-2)
- [§3. 분산에서 이차 모멘트의 약수합으로](#section-3)
- [§4. 주항 제거와 영점 진동의 관계](#section-4)
- [§5. 이차 모멘트의 전체 주항](#section-5)
- [§6. 무한차 모멘트의 끝점합과 잔차 E](#section-6)
- [§7. 잔차 E의 정확한 멜린 변환](#section-7)
- [§8. 𝓔₂의 지수와 단순영점](#section-8)
- [§9. 주항의 점근 전개](#section-9)
- [참고 문헌](#references)

기호 $B,H,\mathcal P,E$는 [짝수 중복도 영점 문서](even-zeros.md#notation)와 같은 정의이다. 비자명 영점 $\rho$의 중복도는 $m_\rho$로 쓴다.

<a id="section-1"></a>

## 1. 약수의 가중치

$f=d_{1/2}$라 쓰자. $f$는 곱셈적 함수이며 소수 거듭제곱에서

$$
f(1)=1,\qquad
f(p^a)=\frac{\binom{2a}{a}}{4^a}
\quad(a\ge0)
$$

로 정의한다. 여기서 곱셈적이라는 것은 $\gcd(r,s)=1$일 때 $f(rs)=f(r)f(s)$라는 뜻이다.

약수 $d\mid n$에 부여하는 가중치는

$$
\boxed{w_n(d)=f(d)f(n/d).}
$$

이 가중치는 양수이고, 약수와 짝약수에 대해 대칭이다.

$$
w_n(d)>0,\qquad w_n(d)=w_n(n/d).
$$

$c_a=\binom{2a}{a}/4^a$라 놓으면

$$
\left(\sum_{a\ge0}c_a z^a\right)^2
=\bigl((1-z)^{-1/2}\bigr)^2
=\frac1{1-z}
$$

이므로 계수 비교를 통해 $\sum_{k=0}^{a}c_kc_{a-k}=1$을 얻는다. 이를 각 소수 거듭제곱에 적용하면

$$
\boxed{\sum_{d\mid n}w_n(d)=1.}
$$

따라서 $w_n$은 $n$의 약수 집합 위의 확률 가중치다. 이는 [주정리 §8.4](main-theorem.md#section-8-4)의 이산핵 $\mathscr K_n^{1/2,1/2}$와 같다.

<a id="section-2"></a>

## 2. 로그 좌표와 약수쌍의 대칭성

$\log$는 자연로그로 쓰고

$$
L_n=\log n,\qquad
u_d=\log_n d=\frac{\log d}{L_n}
$$

라 놓자. $n\ge2$이므로 $L_n>0$이고, 약수 좌표는 $0\le u_d\le1$을 만족한다. 짝약수 $m=n/d$의 좌표는

$$
u_m=1-u_d
$$

이다. 가중치의 대칭성과 정규화를 함께 사용하면

$$
\begin{aligned}
\sum_{d\mid n}w_n(d)u_d
&=\frac12\sum_{d\mid n}w_n(d)(u_d+u_{n/d})\\
&=\frac12\sum_{d\mid n}w_n(d)
=\frac12.
\end{aligned}
$$

확률변수 $D_n$을 $\mathbb P(D_n=d)=w_n(d)$로 정의하고

$$
X_n=\log D_n,\qquad U_n=\frac{X_n}{L_n}
$$

이라 두면, 위 항등식은

$$
\mathbb E[U_n]=\frac12,\qquad
\mathbb E[X_n]=\frac{L_n}{2}
$$

를 뜻한다. 여기서 기대값은 고정된 $n$의 약수에 대해 취한다.

<a id="section-3"></a>

## 3. 분산에서 이차 모멘트의 약수합으로

평균 $\mathbb E[X_n]=L_n/2$를 사용하면 분산은

$$
\begin{aligned}
V(n):=\operatorname{Var}(X_n)
&=\sum_{d\mid n}w_n(d)
\left(\log d-\frac{L_n}{2}\right)^2\\
&=\frac14\sum_{d\mid n}w_n(d)
\left(\log\frac{d}{n/d}\right)^2.
\end{aligned}
$$

같은 양을 로그의 곱으로 쓰면

$$
V(n)=\frac12\sum_{d\mid n}w_n(d)\log d\,
\log\frac{d}{n/d}.
$$

마지막 등식은 합에서 $d$와 $n/d$를 교환한 식을 원래 식과 평균하여 얻는다.

로그의 밑을 $n$으로 바꾸면

$$
\log\frac{d}{n/d}=L_n(2u_d-1)
$$

이므로

$$
\frac{V(n)}{L_n^2}
=\frac14\sum_{d\mid n}w_n(d)(2u_d-1)^2.
$$

제곱을 전개하고 앞서 구한 두 합을 대입하면

$$
\begin{aligned}
\frac{V(n)}{L_n^2}
&=\sum_{d\mid n}w_n(d)u_d^2
-\sum_{d\mid n}w_n(d)u_d
+\frac14\sum_{d\mid n}w_n(d)\\
&=\sum_{d\mid n}w_n(d)u_d^2-\frac14.
\end{aligned}
$$

따라서 $U_n$의 이차 모멘트를 $M_2(n)=\mathbb E[U_n^2]$로 정의하면, 원하는 가중 로그제곱의 약수합은 다음과 같다.

<a id="equation-3-1"></a>

$$
\boxed{
\begin{aligned}
M_2(n)
&=\sum_{d\mid n}w_n(d)\bigl(\log_n d\bigr)^2\\
&=\frac1{(\log n)^2}
\sum_{d\mid n}d_{1/2}(d)d_{1/2}(n/d)(\log d)^2\\
&=\frac14+\frac{V(n)}{(\log n)^2}.
\end{aligned}
}
\tag{3.1}
$$

<a id="section-4"></a>

## 4. 주항 제거와 영점 진동의 관계

양의 $(1/2,1/2)$ 약수분포와 그 이차 모멘트는

$$
\mu_n=\sum_{d\mid n}w_n(d)\delta_{u_d},\qquad
M_2(n)=\int_0^1u^2\,d\mu_n(u)
$$

로 쓴다. $\delta_{u_d}$는 약수 좌표 $u_d$에 놓인 단위 질량이다. $M_2(n)$ 자체는 원래의 양의 분포에 대한 이차 모멘트다.

주항을 제거한 진동을 분석할 때는 먼저 정수에 대한 부분합

$$
T_2(x)=\sum_{2\le n\le x}M_2(n)
$$

을 취한다. 이 부분합의 비진동 주항 전체를 $\mathcal P_2(x)$라 하면, 분석 대상은

$$
\boxed{\mathcal E_2(x)=T_2(x)-\mathcal P_2(x)}
$$

이다. 여기서 전체 주항에는 아크사인 선도항 $3x/8$뿐 아니라 뒤따르는 로그 보정항들도 포함한다.

이 관점에서 $\mathcal E_2(x)$는 **주항을 제거한 누적 약수분포의 이차 모멘트에 해당하는 잔차**다. 원래 분포의 이차 모멘트 $M_2(n)$와 구분하여, 영점의 진동과 연결할 대상은 이 잔차로 잡는다. 영점 성분과 절단오차 등을 구분하는 명시적 공식은 그다음 분석 단계에 해당한다.

이차 모멘트에 나타나는 영점 성분에는 단순 영점의 기여가 남고, 중복 영점의 해당 기여는 짝수든 홀수든 소거된다. 따라서 이 모멘트가 모든 영점의 기여를 검출한다고 해석하지 않는다. 정확한 진술은 [정리 8.1](#section-8-1)이다. 제곱 인수 없는 모형에서 같은 현상은 [영점 기여 §2–3](zero-contributions.md#section-2)에 정리되어 있다.

<a id="section-5"></a>

## 5. 이차 모멘트의 전체 주항

<a id="section-5-1"></a>

### 5.1. 정확한 디리클레급수와 극의 주부분

$\Re s>1$에서 양의 실수축의 제곱근으로 가지를 정하고

$$
F(s):=\sum_{n\ge1}\frac{f(n)}{n^s}=\sqrt{\zeta(s)}
$$

라 쓰자. 로그 분모를 없앤 이차 모멘트의 계수를

$$
a_2(n):=(\log n)^2M_2(n)
=\sum_{d\mid n}f(d)f(n/d)(\log d)^2,
\qquad a_2(1):=0
$$

로 정의하면, 절대수렴 영역에서 항별 미분과 합성곱에 의해

<a id="equation-5-1"></a>

$$
\boxed{
\begin{aligned}
F_2(s)&:=\sum_{n\ge1}\frac{a_2(n)}{n^s}
=F(s)F''(s)\\
&=\frac12\zeta''(s)-\frac14\frac{\zeta'(s)^2}{\zeta(s)}.
\end{aligned}
}
\tag{5.1}
$$

스틸체스 상수의 부호는

$$
\zeta(1+h)=\frac1h+\sum_{j\ge0}\frac{(-1)^j\gamma_j}{j!}h^j,
\qquad \gamma_0=\gamma
$$

로 고정한다. 직접 미분하여 얻는 $s=1$의 주부분은

<a id="equation-5-2"></a>

$$
F_2(1+h)
=\frac{3}{4h^3}+\frac{\gamma}{4h^2}
-\frac{\gamma^2+3\gamma_1}{4h}+O(1).
\tag{5.2}
$$

여기서 $O(1)$은 $h=0$ 근방의 정칙 부분을 뜻한다. 정수 부분합에 대한 점근식을 사용한 단계는 아니다.

$c:=\gamma^2+3\gamma_1$이라 두면, 이 극에서 오는 로그 가중 부분합의 주항은 정확히

<a id="equation-5-3"></a>

$$
\boxed{
\begin{aligned}
Q_2(x)
&:=\operatorname*{Res}_{s=1}\frac{F_2(s)x^s}{s}\\
&=x\left[
\frac38(\log x)^2
+\frac{\gamma-3}{4}\log x
+\frac{3-\gamma-c}{4}
\right].
\end{aligned}
}
\tag{5.3}
$$

<a id="section-5-2"></a>

### 5.2. 로그 분모를 복원한 전체 주항과 잔차

$x\ge2$에서

$$
U_2(x):=\sum_{2\le n\le x}a_2(n),
\qquad \Delta_2(x):=U_2(x)-Q_2(x)
$$

라 놓는다. 부분합 공식은 정확히

$$
T_2(x)=\frac{U_2(x)}{(\log x)^2}
+2\int_2^x\frac{U_2(t)}{t(\log t)^3}\,dt
$$

를 준다. 따라서 이 문서에서 전체 주항의 상수항까지 고정하는 규약은

<a id="equation-5-4"></a>

$$
\boxed{
\mathcal P_2(x):=
\frac{Q_2(x)}{(\log x)^2}
+2\int_2^x\frac{Q_2(t)}{t(\log t)^3}\,dt
}
\tag{5.4}
$$

이다. 이에 대응하는 잔차는

<a id="equation-5-5"></a>

$$
\boxed{
\begin{aligned}
\mathcal E_2(x)
&:=T_2(x)-\mathcal P_2(x)\\
&=\frac{\Delta_2(x)}{(\log x)^2}
+2\int_2^x\frac{\Delta_2(t)}{t(\log t)^3}\,dt.
\end{aligned}
}
\tag{5.5}
$$

전체 주항을 더 간단한 적분으로 쓰려면

$$
I_j(x):=\int_2^x\frac{dt}{(\log t)^j}
$$

라 놓으면 된다. $Q_2'(x)/(\log x)^2=3/8+\gamma/(4\log x)-c/(4(\log x)^2)$이므로

<a id="equation-5-6"></a>

$$
\boxed{
\mathcal P_2(x)
=\frac38x+\frac\gamma4 I_1(x)-\frac c4 I_2(x)+\kappa_2,
\qquad
\kappa_2:=\frac{Q_2(2)}{(\log2)^2}-\frac34.
}
\tag{5.6}
$$

이 식은 $3x/8$과 모든 로그 보정항을 하나의 함수로 정한다. $\kappa_2$를 다른 고정 상수로 바꾸면 $\mathcal E_2$도 상수만큼 달라지며, 아래의 지수 상계에는 영향을 주지 않는다.

$\mathcal P_2'=Q_2'/(\log x)^2$와 $T_2$의 정의에서 $[2,\infty)$ 위의 측도로 $dU_2=(\log t)^2dT_2$, $dQ_2=(\log t)^2d\mathcal P_2$이다. 이를 부분적분하면 (5.5)의 역공식

<a id="equation-5-7"></a>

$$
\Delta_2(x)=(\log x)^2\mathcal E_2(x)-2\int_2^x\mathcal E_2(t)\,\frac{\log t}{t}\,dt\qquad(x\ge2)
\tag{5.7}
$$

를 얻는다. 경계항은 $\Delta_2(2)=a_2(2)-Q_2(2)=(\log2)^2\mathcal E_2(2)$이므로 정확히 상쇄된다.

<a id="section-6"></a>

## 6. 무한차 모멘트의 끝점합과 잔차 $E$

<a id="section-6-1"></a>

### 6.1. 고정된 $n$과 $N$에서 먼저 취하는 극한

정수 $k\ge0$에 대해

$$
M_k(n):=\sum_{d\mid n}w_n(d)
\left(\frac{\log d}{\log n}\right)^k,
\qquad
T_k(N):=\sum_{n=2}^{N}M_k(n)
$$

라 정의한다. $k=0$에서는 모든 좌표의 $0$차 거듭제곱을 $1$로 정한다. 고정된 $n\ge2$에서 $d<n$이면 $\log d/\log n<1$이고, $d=n$이면 그 값이 $1$이므로

$$
\boxed{
\lim_{k\to\infty}M_k(n)=w_n(n)=f(n)=d_{1/2}(n).
}
$$

따라서 고정된 정수 $N\ge2$에서는 유한합과 극한을 교환할 수 있어

$$
\boxed{
\lim_{k\to\infty}T_k(N)=\sum_{n=2}^{N}f(n).
}
$$

아크사인 분포의 $k$차 모멘트를

$$
c_k:=\frac1\pi\int_0^1\frac{u^k}{\sqrt{u(1-u)}}\,du
=\frac{\binom{2k}{k}}{4^k}
$$

라 하면 $c_k\to0$이다. 그러므로 아크사인 선도항만 제거한 합에도

$$
\boxed{
\lim_{k\to\infty}\bigl[T_k(N)-(N-1)c_k\bigr]
=\sum_{n=2}^{N}f(n)
}
$$

가 성립한다. 이 극한에서는 먼저 $N$을 고정한다. $k=k(N)$인 동시 극한이나, 각 유한 차수의 전체 주항을 먼저 제거한 잔차의 극한과는 구분한다.

<a id="section-6-2"></a>

### 6.2. 끝점합의 전체 주항을 실수 적분으로 고정

끝점합을

$$
A(x):=\sum_{n\le x}f(n)
$$

라 놓고, 이 문서에서는 다음 함수를 그 전체 주항으로 선택한다.

$$
\boxed{
\mathcal P(x)
:=1+\frac1\pi\int_0^1
\sqrt{-\zeta(u)}\,\frac{x^u-1}{u}\,du,
\qquad x\ge1.
}
$$

$0<u<1$에서 $\zeta(u)<0$이므로 양의 실수 제곱근을 취한다. $u=0$에서는 $(x^u-1)/u\to\log x$이고, $u=1$에서는 적분 가능한 $(1-u)^{-1/2}$형 특이점이 생긴다. 따라서 이 적분은 수렴하며 $\mathcal P(1)=1$이다.

무한차 모멘트에서 출발한 잔차를

$$
\boxed{E(x):=A(x)-\mathcal P(x)}
$$

로 정의한다. 이 $E$는 [짝수 중복도 영점 문서 (1.3)](even-zeros.md#equation-1-3)의 $E$와 같다. 특히 정수 $N\ge2$에서

$$
\boxed{
E(N)=
\lim_{k\to\infty}\bigl[T_k(N)-(N-1)c_k\bigr]
-\bigl[\mathcal P(N)-1\bigr].
}
$$

즉 이 $E$는 무한차 극한으로 끝점합을 얻은 뒤, 그 끝점합의 전체 주항을 제거한 잔차다. $E=\lim_{k\to\infty}\mathcal E_k$라는 교환 명제를 여기서 가정하지 않는다.

<a id="section-6-3"></a>

### 6.3. 정확한 산술 표현

$g=d_{-1/2}$라 쓰면 $f=\mathbf1*g$이고, 소수 거듭제곱에서

$$
g(1)=1,\qquad
g(p^a)=-\frac1{2a-1}\frac{\binom{2a}{a}}{4^a}
\quad(a\ge1)
$$

이다. 따라서 $\{y\}=y-\lfloor y\rfloor$를 사용하면

$$
\boxed{
\begin{aligned}
E(x)
&=\sum_{d\le x}g(d)\left\lfloor\frac xd\right\rfloor-\mathcal P(x)\\
&=x\sum_{d\le x}\frac{g(d)}d
-\sum_{d\le x}g(d)\left\{\frac xd\right\}
-\mathcal P(x).
\end{aligned}
}
$$

이는 유한합 항등식이다. $E$의 상계를 산술적으로 조사할 때는 가중 부분합과 소수부분 합 사이의 상쇄도 함께 보아야 한다.

<a id="section-7"></a>

## 7. 잔차 $E$의 정확한 멜린 변환

주항의 밀도와 그 변환을

$$
p(x):=\mathcal P'(x)
=\frac1\pi\int_0^1\sqrt{-\zeta(u)}\,x^{u-1}\,du,
$$

$$
B(s):=\frac1\pi\int_0^1\frac{\sqrt{-\zeta(u)}}{s-u}\,du
$$

로 정의한다. $B$는 $\mathbb C\setminus[0,1]$에서 정칙이다. $\Re s>1$에서 절대수렴과 부분적분에 의해

$$
s\int_1^\infty A(x)x^{-s-1}\,dx=F(s),
\qquad
s\int_1^\infty\mathcal P(x)x^{-s-1}\,dx=1+B(s)
$$

가 성립한다. 따라서

<a id="equation-7-1"></a>

$$
\boxed{
H(s):=s\int_1^\infty E(x)x^{-s-1}\,dx
=F(s)-1-B(s),\qquad \Re s>1.
}
\tag{7.1}
$$

이 주항이 $s=1$의 특이 부분 전체를 제거한다는 사실은 국소적으로 확인할 수 있다. $h=0$ 근방에서

$$
\psi(h):=\sqrt{h\zeta(1+h)},\qquad \psi(0)=1
$$

로 정칙인 가지를 택하면

$$
F(1+h)=h^{-1/2}\psi(h),
\qquad
\sqrt{-\zeta(1-t)}=t^{-1/2}\psi(-t)
$$

이다. $B(1+h)$의 적분을 $0<t<\delta$와 나머지로 나누고, 앞부분에서 $\psi(-t)-\psi(h)$가 $t+h$로 나누어짐을 사용하면

$$
B(1+h)=h^{-1/2}\psi(h)+J(h)
$$

를 얻는다. 여기서 $J$는 $h=0$ 근방에서 정칙이다. 그러므로 $F-1-B$도 $s=1$ 근방에서 정칙하게 이어진다. 이 과정에서는 실제 영점의 위치나 영점합을 사용하지 않는다. 짝수 중복도 조건 아래 $H$가 반평면 전체로 연장된다는 것은 [짝수 중복도 영점 보조정리 4.1](even-zeros.md#section-4)이고, 그 역과 $E$의 정확한 성장 지수는 같은 문서의 [따름정리 3](even-zeros.md#section-8), [정리 2](even-zeros.md#section-9)이다.

<a id="section-8"></a>

## 8. $\mathcal E_2$의 지수와 단순영점

이전 판의 이 절은 $E$의 상계에서 $F$의 정칙 연장을 얻고, 코시 미분공식으로 $F''$를 제어했다. 그 논증은 옳지만 필요 이상의 가정을 쓴다. $F_2$는 처음부터 유리형이고(§8.2), 그 극은 단순영점에만 있으며, 평균 크기는 $F$의 정칙성 없이 제어된다(§8.3–8.4). 이 절은 그 결과를 동치 정리로 정리한다.

<a id="section-8-1"></a>

### 8.1. 주정리

**정리 8.1.** $\tfrac12\le\theta<1$에서 다음은 동치이다.

- (a) $\Re\rho>\theta$인 $\zeta$의 단순영점이 없다.
- (b) 모든 $\varepsilon>0$에 대해 $\Delta_2(x)=O_\varepsilon(x^{\theta+\varepsilon})$.
- (c) 모든 $\varepsilon>0$에 대해 $\mathcal E_2(x)=O_\varepsilon(x^{\theta+\varepsilon})$.

(a) 아래에서는 고정된 $0<\delta<\tfrac12$에 대해 $\Delta_2(x)\ll_{\theta,\delta}x^{\theta+\delta}\log x$, $\mathcal E_2(x)\ll_{\theta,\delta}x^{\theta+\delta}/\log x$이다.

중복도 $2$ 이상의 영점은 짝수든 홀수든 조건에 들어가지 않는다. (a)$\Rightarrow$(b)의 증명은 $|\zeta|$의 하계도, $F$의 정칙 연장도 쓰지 않는다.

**따름정리 8.2 (전달).** 모든 $\eta>0$에 대해 $E(x)=O_\eta(x^{1/2+\eta})$이면, 모든 $\varepsilon>0$에 대해 $\mathcal E_2(x)=O_\varepsilon(x^{1/2+\varepsilon})$이다. 특히 $E(x)=O(\sqrt x(\log x)^A)$와 같은 고정된 로그 손실의 제곱근 상계는 이 가정을 충족한다.

**증명.** [짝수 중복도 따름정리 3](even-zeros.md#section-8)에 의해 가정은 $(\mathrm{EM}_{1/2})$, 즉 $\Re\rho>\tfrac12$인 영점이 모두 짝수 중복도라는 것과 동치이다. 그러면 그 영역에 단순영점이 없으므로 정리 8.1의 (a)$\Rightarrow$(c)를 $\theta=\tfrac12$로 적용한다. $\square$

이 전달은 한 방향이다. $\Re s>\tfrac12$에 단순영점은 없고 삼중 이상의 홀수 중복도 영점이 있다면 $E$는 크지만([짝수 중복도 정리 2](even-zeros.md#section-9)) $\mathcal E_2$는 작다. 또한 이 명제는 $E$의 제곱근 상계 자체를 증명한 것이 아니며, 로그 손실 없는 $O(\sqrt x)$의 보존은 포함하지 않는다.

**따름정리 8.3 (가정 없음).** $\Theta_{\rm simple}:=\sup\{\Re\rho:\ m_\rho=1\}$로 두자.

1. $\inf\{\alpha:\ \mathcal E_2(x)=O(x^\alpha)\}=\Theta_{\rm simple}$이고 $\Theta_{\rm simple}\ge\tfrac12$이다.
2. $\rho$가 단순영점이고 $\beta=\Re\rho$이면, 모든 $0<\alpha<\beta$에 대해 $\Delta_2(x)>x^\alpha$인 $x$와 $\Delta_2(x)<-x^\alpha$인 $x$가 각각 임의로 큰 값에서 존재하고, $\mathcal E_2(x)\ne O(x^\alpha)$이다. 특히 $\Delta_2(x)=\Omega_\pm(x^{1/2-\varepsilon})$, $\mathcal E_2(x)=\Omega(x^{1/2-\varepsilon})$이다.

$E$의 성장 지수는 홀수 중복도 영점의 상한 $\Theta_{\rm odd}$이고([짝수 중복도 정리 2](even-zeros.md#section-9)) $\Theta_{\rm simple}\le\Theta_{\rm odd}$이므로, $\mathcal E_2$의 지수는 $E$의 지수를 넘지 않는다. 따름정리 8.2는 이 부등식의 $\Theta_{\rm odd}\le\tfrac12$인 경우이다.

<a id="section-8-2"></a>

### 8.2. $F_2$의 유리형 구조

**보조정리 8.4.** $F_2=\tfrac12\zeta''-\tfrac14\zeta'^2/\zeta$는 $\mathbb C$에서 유리형이다. $\Re s>0$에서 그 극은 정확히 다음과 같다.

- $s=1$: 3차 극, 주부분 (5.2).
- 각 단순영점 $\rho$: 단순극, 유수 $-\tfrac14\zeta'(\rho)$.

중복도 $m\ge2$인 영점에서 $F_2$는 정칙이다.

**증명.** 유리형인 것은 (5.1)의 둘째 식에서 바로 나온다. $\Re s>0$에서 $\zeta''$의 극은 $s=1$뿐이고, $\zeta'^2/\zeta$의 특이점은 $s=1$과 $\zeta$의 영점뿐이다. $\zeta$는 실수 구간 $(0,1)$에서 음수이므로 이 영점들은 비자명 영점이다. 중복도 $m$인 영점에서 $\zeta=c(s-\rho)^m(1+O(s-\rho))$, $c=\zeta^{(m)}(\rho)/m!\ne0$이면

$$
\frac{\zeta'(s)^2}{\zeta(s)}=m^2c\,(s-\rho)^{m-2}\bigl(1+O(s-\rho)\bigr)
$$

이다. $m\ge2$이면 정칙이고, $m=1$이면 유수 $c=\zeta'(\rho)$의 단순극이다. $\square$

유수 $-\tfrac14$는 [주정리 §7](main-theorem.md#section-7)의 $\langle\mathcal K_{-1/2,-1/2},u^2\rangle=-\tfrac14$와 같다. 이는 [영점 기여 (1.1)](zero-contributions.md#equation-1-1)을 $H\equiv1$인 $d_z$ 모형에 적용한 국소 핵과 일치한다. 이전 판 §8.2의 "부호 변화가 $F$와 $F''$에 동시에 작용한다"는 논증은 이 보조정리로 대체된다.

<a id="section-8-3"></a>

### 8.3. $\Delta|g|$ 항등식

**보조정리 8.5.** $\Omega\subset\mathbb C$가 열린 집합이고 $g$가 $\Omega$에서 정칙이며 항등적으로 $0$이 아니면, $|g'|^2/|g|\in L^1_{\rm loc}(\Omega)$이고 분포의 의미에서

<a id="equation-8-1"></a>

$$
\Delta|g|=\frac{|g'|^2}{|g|}.
\tag{8.1}
$$

즉 모든 $\chi\in C_c^\infty(\Omega)$에 대해 $\iint|g|\,\Delta\chi\,dA=\iint\chi\,|g'|^2|g|^{-1}\,dA$이다.

**증명.** $g\ne0$인 곳에서 $v=g\bar g$로 두면 $\partial v=g'\bar g$, $\bar\partial v=g\overline{g'}$, $\partial\bar\partial v=|g'|^2$이다. 따라서

$$
\Delta v^{1/2}=4\partial\bar\partial v^{1/2}
=4\Big(\frac{\partial\bar\partial v}{2v^{1/2}}-\frac{\partial v\,\bar\partial v}{4v^{3/2}}\Big)
=\frac{2|g'|^2}{|g|}-\frac{|g'|^2}{|g|}=\frac{|g'|^2}{|g|}
$$

이다. $m$중 영점 근처에서는 $|g'|^2/|g|\asymp|s-\rho|^{m-2}$이므로 평면에서 국소 적분 가능하다. $\operatorname{supp}\chi$ 안의 유한 개 영점 주위에 반지름 $\varepsilon$의 원판을 빼고 그린 공식을 쓰면, 경계항은

$$
\oint_{|s-\rho|=\varepsilon}\bigl(|g|\,\partial_n\chi-\chi\,\partial_n|g|\bigr)
=O(\varepsilon^m\cdot\varepsilon)+O(\varepsilon^{m-1}\cdot\varepsilon)\to0
$$

이다. 둘째 항에는 $g\ne0$인 곳에서 $|\nabla|g||=|g'|$를 썼다. $\varepsilon\to0$에서 지배수렴정리로 (8.1)을 얻는다. $\square$

특히 $|g|$는 열조화이고 그 Riesz 측도는 $|g'|^2|g|^{-1}dA$이며, 영점에 원자가 없다.

<a id="section-8-4"></a>

### 8.4. 수직선상의 평균 크기

고정된 $\sigma>\tfrac12$에서 성립하는 무조건적 평균제곱 상계

$$
\int_T^{2T}|\zeta(\sigma+it)|^2\,dt\ll_\sigma T
\qquad(T\ge2)
$$

를 사용한다. 이 상계는 $\sigma>\tfrac12$의 고정된 콤팩트 구간에서 균일하다. 표준 평균제곱 공식과 출처는 참고 문헌 [1]의 §1, 식 (1.1)에 정리되어 있다. 코시–슈바르츠 부등식으로 같은 범위에서 $\int_T^{2T}|\zeta(\sigma+it)|\,dt\ll T$이다.

**보조정리 8.6.** $\tfrac12\le\theta<1$, $0<\delta<\tfrac12$이고 정리 8.1의 (a)를 가정한다. 그러면 $\theta+\delta\le\sigma\le2$, $T\ge3$에서 균등하게

<a id="equation-8-2"></a>

$$
\int_T^{2T}|F_2(\sigma+it)|\,dt\ll_{\theta,\delta}T.
\tag{8.2}
$$

**증명.** $r=\delta/4$로 둔다. 아래의 모든 원판과 띠는 $\{\Re s\ge\theta+\delta/2,\ \Im s\ge1\}$ 안에 있고, 거기서 $\zeta$는 정칙이며 위의 평균값이 균일하게 성립한다.

*$\zeta''$ 항.* 코시 공식으로 $|\zeta''(s)|\le\frac{2}{r^2}\cdot\frac1{2\pi}\int_0^{2\pi}|\zeta(s+re^{i\phi})|\,d\phi$이다. $t$에 대해 적분하면 $\int_T^{2T}|\zeta''(\sigma+it)|\,dt\ll r^{-2}T$이다.

*$\zeta'^2/\zeta$ 항.* (a)와 보조정리 8.4에 의해 $h=\zeta'^2/\zeta$는 $\{\Re s>\theta,\ \Im s>1\}$에서 정칙이므로 $|h|$는 열조화이다. 따라서 $|h(\sigma+it)|\le(\pi r^2)^{-1}\iint_{D(\sigma+it,r)}|h|\,dA$이고, $t\in[T,2T]$에 대해 적분하면

$$
\int_T^{2T}|h(\sigma+it)|\,dt\le\frac{2}{\pi r}\iint_S|h|\,dA,\qquad
S=[\sigma-r,\sigma+r]\times[T-r,2T+r]
$$

이다. $\chi=\chi_1(\Re s)\chi_2(\Im s)$를 $S$ 위에서 $1$, $S'=[\sigma-2r,\sigma+2r]\times[T-r-1,2T+r+1]$ 밖에서 $0$이고 $|\chi_1''|\ll r^{-2}$, $|\chi_2''|\ll1$인 비음 함수로 잡는다. $|h|=|\zeta'|^2/|\zeta|$이므로 보조정리 8.5에 의해

$$
\iint_S|h|\,dA\le\iint\chi\,\frac{|\zeta'|^2}{|\zeta|}\,dA=\iint|\zeta|\,\Delta\chi\,dA
\ll r^{-2}\int_{\sigma-2r}^{\sigma+2r}\int_{T-r-1}^{2T+r+1}|\zeta(\sigma'+it)|\,dt\,d\sigma'\ll\frac Tr
$$

이다. 그러므로 $\int_T^{2T}|h|\ll T/r^2$이다. 두 항을 합하면 (8.2)이다. $\square$

이 증명은 $|\zeta|$의 하계를 쓰지 않는다. $\zeta'^2/\zeta$의 평균을 $|\zeta|$ 자체의 평균으로 바꾸는 것은 (8.1)이며, 열조화성에는 $h$의 정칙성, 즉 단순영점이 없다는 것만 쓰인다.

<a id="section-8-5"></a>

### 8.5. 절단 Perron 공식과 로그 가중 잔차

(a)를 가정하고 $0<\delta<\tfrac12$를 고정한다. $0\le M_2(n)\le1$이므로 $0\le a_2(n)\le(\log n)^2$이다. 먼저 $x$가 반정수인 경우를 고려하고, $c_x=1+1/\log x$에서 높이 $T\in[x^2,2x^2]$의 절단 Perron 공식을 적용한다. 계수 상계와 $F_2(c_x)\ll(\log x)^3$으로부터 절단오차는

$$
O\!\left(\frac{x(\log x)^3}{T}\right)=O\!\left(\frac{(\log x)^3}{x}\right)
$$

로 제어된다.

적분선을 $\sigma_0=\theta+\delta$로 옮긴다. 직사각형 $\sigma_0\le\Re s\le c_x$, $|\Im s|\le T$ 안에서 $F_2(s)x^s/s$의 특이점은 (a)와 보조정리 8.4에 의해 $s=1$의 3차 극뿐이고, 그 유수가 $Q_2(x)$이다. 보조정리 8.6을 $\sigma\in[\sigma_0,c_x]$에서 적분하면 $\int_{x^2}^{2x^2}\int_{\sigma_0}^{c_x}|F_2(\sigma+it)|\,d\sigma\,dt\ll x^2$이므로

$$
\int_{\sigma_0}^{c_x}|F_2(\sigma+iT)|\,d\sigma\ll_{\theta,\delta}1
$$

인 높이 $T\in[x^2,2x^2]$를 택할 수 있다. 계수가 실수이므로 $F_2(\bar s)=\overline{F_2(s)}$이고 아래쪽 수평변도 같은 추정을 만족한다. 따라서 수평변의 적분은 $O(x^{c_x}/T)=O(1/x)$이다.

왼쪽 수직변에서 $|t|\le3$인 부분은 $F_2$가 그 선분에서 연속이므로 $O(x^{\sigma_0})$이다. 높이 $3$ 아래에는 비자명 영점이 없고 $\sigma_0<1$이기 때문이다. 나머지는 이진 구간별로 보조정리 8.6을 적용하여

$$
\int_3^T\frac{|F_2(\sigma_0+it)|}{t}\,dt\ll_{\theta,\delta}\log T
$$

를 얻는다. 그러므로

<a id="equation-8-3"></a>

$$
\boxed{
\Delta_2(x)=U_2(x)-Q_2(x)
\ll_{\theta,\delta}x^{\theta+\delta}\log x.
}
\tag{8.3}
$$

일반 실수 $x$에서는 $\lfloor x\rfloor+1/2$로 옮겨 부분합을 그대로 유지할 수 있다. 이때 $Q_2$의 변화는 $O((\log x)^2)$이므로 같은 상계가 성립한다.

<a id="section-8-6"></a>

### 8.6. 로그 정규화 복원

(5.5)에 (8.3)을 넣으면 $x\ge3$에서

$$
\begin{aligned}
\mathcal E_2(x)
&=\frac{\Delta_2(x)}{(\log x)^2}
+2\int_2^x\frac{\Delta_2(t)}{t(\log t)^3}\,dt\\
&\ll_{\theta,\delta}\frac{x^{\theta+\delta}}{\log x}
+\int_2^x\frac{t^{\theta+\delta-1}}{(\log t)^2}\,dt
\ll_{\theta,\delta}\frac{x^{\theta+\delta}}{\log x}.
\end{aligned}
$$

$\delta$를 임의로 작게 고를 수 있으므로 (a)$\Rightarrow$(b), (a)$\Rightarrow$(c)가 증명되었다.

<a id="section-8-7"></a>

### 8.7. 필요조건과 진동

**(c)$\Rightarrow$(b).** (5.7)에서 $\mathcal E_2(x)=O(x^{\theta+\varepsilon})$이면 $\Delta_2(x)=O(x^{\theta+\varepsilon}(\log x)^2)$이다.

**(b)$\Rightarrow$(a).** $1\le x<2$에서는 $U_2=0$이므로 $\Delta_2=-Q_2$로 둔다. $\Re s>1$에서

$$
s\int_1^\infty\Delta_2(x)x^{-s-1}\,dx=F_2(s)-R_2(s),\qquad
R_2(s):=s\int_1^\infty Q_2(x)x^{-s-1}\,dx
$$

이다. $\int_1^\infty x^{-s}(\log x)^jdx=j!/(s-1)^{j+1}$이므로 $R_2$는 $s=1$에서만 극을 갖는 유리함수이다. 더 정확히는 $Q_2(x)=\frac1{2\pi i}\oint_{|w-1|=\epsilon}F_2(w)x^ww^{-1}dw$를 대입하고 리만 구 위의 유수 정리를 쓰면 $R_2(s)=P(s)-P(0)$이다. 여기서 $P$는 $F_2$의 $s=1$에서의 주부분이다. 따라서 $F_2-R_2$는 $s=1$에서 정칙이다.

(b)를 가정하면 좌변은 $\Re s>\theta$에서 정칙이므로 $F_2=R_2+(\Re s>\theta$에서 정칙$)$이다. $F_2$는 $\Re s>\theta$에서 $s=1$ 외의 극을 갖지 않고, 보조정리 8.4에 의해 (a)가 성립한다. $\square$

**따름정리 8.3의 증명.** 단순영점 $\rho$, $\beta=\Re\rho$, $0<\alpha<\beta$를 고정하고, 큰 $x$에서 $\Delta_2(x)\le x^\alpha$라 하자. $g(x)=x^\alpha-\Delta_2(x)$는 결국 비음이고

$$
\int_1^\infty g(x)x^{-s-1}dx=\frac1{s-\alpha}-\frac{F_2(s)-R_2(s)}s\qquad(\Re s>1)
$$

이다. 보조정리 8.4에 의해 $F_2-R_2$는 양의 실수 반직선의 한 근방에서 정칙이다. $\Re s>0$의 극 가운데 $s=1$은 상쇄되었고 나머지는 실수가 아닌 단순영점이기 때문이다. $\Delta_2(x)\ll x(\log x)^2$이므로 수렴 가로좌표 $\sigma_c$는 유한하고, Landau의 정리 [5]에 의해 $\sigma_c\le\alpha$이다. 따라서 $F_2-R_2$는 $\Re s>\alpha$로 정칙하게 연장되고, 이는 $\rho$와 모순이다. $g=x^\alpha+\Delta_2$에 대해서도 같다. 또한 $\mathcal E_2(x)=O(x^\alpha)$이면 (5.7)에서 $\Delta_2(x)=O(x^\alpha(\log x)^2)$이고 위와 같은 모순이 생긴다.

Levinson의 방법으로 얻는 임계선 위 영점은 단순영점이며, 그 비율은 양수이다 [3, 4]. 따라서 $\beta=\tfrac12$인 단순영점이 존재하여 $\Omega$ 결과와 $\Theta_{\rm simple}\ge\tfrac12$를 얻는다. 지수의 상계는 $\Theta_{\rm simple}<1$이면 상한의 정의에 의해 (a)가 $\theta=\Theta_{\rm simple}$에서 성립하므로 정리 8.1이고, $\Theta_{\rm simple}=1$이면 자명한 상계 $\mathcal E_2(x)\ll x$이다. $\square$

<a id="section-8-8"></a>

### 8.8. 고정 차수 $k$

고정된 $k\ge2$에 대해

$$
F_k(s):=\sum_{n\ge2}\frac{(\log n)^kM_k(n)}{n^s}
=(-1)^kF(s)F^{(k)}(s)
$$

이고, $s=1$의 극에서 얻은 주항 $Q_k(x)=\operatorname{Res}_{s=1}F_k(s)x^s/s$를 (5.4)와 같은 방식으로 로그 정규화에 맞게 복원하여 $\mathcal P_k$, $\mathcal E_k$를 정한다.

**보조정리 8.7.** $F_k$는 $\mathbb C$에서 유리형이다. $\Re s>0$에서 그 극은 $s=1$($k+1$차)과, 중복도 $m<k$인 **홀수** 중복도 영점($k-m$차)뿐이다.

**증명.** 국소적으로 $\phi=\tfrac12\log\zeta$로 두면 Faà di Bruno 공식에 의해 $F^{(k)}=F\cdot Y_k(\phi',\dots,\phi^{(k)})$이고 $Y_k$는 다항식이다. $\phi^{(j)}=\tfrac12(\zeta'/\zeta)^{(j-1)}$은 일가 유리형이므로 $F_k=(-1)^k\zeta\,Y_k$는 유리형이다. 짝수 중복도 영점에서는 $F$가 정칙이므로 $F_k$도 정칙이다. 홀수 $m$이면 $F=c^{1/2}(s-\rho)^{m/2}(1+O(s-\rho))$에서

$$
F(s)F^{(k)}(s)=c\,\Big(\frac m2\Big)^{\underline k}(s-\rho)^{m-k}\bigl(1+O(s-\rho)\bigr)
$$

이고, $m/2$가 정수가 아니므로 하강계승 $(m/2)^{\underline k}$는 $0$이 아니다. $s=1$에서는 $F=(s-1)^{-1/2}\times$(정칙, $1$에서 $1$)이므로 $k+1$차 극이다. $\square$

이로부터 두 방향의 결과를 얻는다.

- **필요조건 (가정 없음).** 모든 $\varepsilon>0$에 대해 $\mathcal E_k(x)=O(x^{\theta+\varepsilon})$이면, $\Re\rho>\theta$에 중복도 $m<k$인 홀수 중복도 영점이 없다. 증명은 §8.7과 같다. $d\Delta_k=(\log t)^kd\mathcal E_k$를 쓰고, $R_k$는 $s=1$에서만 극을 갖는 유리함수이다.
- **충분조건.** $(\mathrm{EM}_\theta)$, 즉 $\Re\rho>\theta$인 영점이 모두 짝수 중복도이면 모든 $\varepsilon>0$에 대해 $\mathcal E_k(x)=O(x^{\theta+\varepsilon})$이다. 이것이 이전 판의 논증이다. $F$가 $\mathcal H_\theta$의 절단 밖에서 정칙이므로 반지름 $r$의 원 위 코시 공식과 코시–슈바르츠 부등식으로 $|F^{(k)}(s)|^2\le(k!)^2r^{-2k}\frac1{2\pi}\int_0^{2\pi}|F(s+re^{i\phi})|^2d\phi$이다. $|F|^2=|\zeta|$이므로 $\int_T^{2T}|F^{(k)}|^2\ll T$이고, $\int_T^{2T}|F_k|\le(\int|F|^2)^{1/2}(\int|F^{(k)}|^2)^{1/2}\ll T$이다. 나머지는 §8.5–8.6과 같다.

$k=2$에서는 두 조건이 "단순영점 없음"으로 일치한다(정리 8.1). $k\ge3$에서 두 조건 사이, 즉 중복도 $k$ 이상의 홀수 중복도 영점이 있는 경우를 닫는 평균 추정은 여기서 다루지 않는다. $F_k$에는 $\zeta'^k/\zeta^{k-1}$ 같은 항이 있어 보조정리 8.5를 그대로 적용할 수 없다.

<a id="section-9"></a>

## 9. 주항의 점근 전개

앞 절들의 정의와 증명은 전체 주항을 정확한 함수로 유지한 상태에서 진행했다. 크기를 비교하기 위해 마지막에만 이를 전개하면, $L=\log x$에 대해

<a id="equation-9-1"></a>

$$
\boxed{
\mathcal P_2(x)
=\frac38x+\frac\gamma4\frac{x}{L}
+\frac{\gamma-\gamma^2-3\gamma_1}{4}\frac{x}{L^2}
+O\!\left(\frac{x}{L^3}\right).
}
\tag{9.1}
$$

앞의 두 항은 [주정리 (8.12)](main-theorem.md#equation-8-12)에서 $F=u^2$로 얻는 $3/8+\gamma/(4L)$과 일치한다. 끝점합의 주항은 실수 적분의 $u=1$ 근방을 전개하여

<a id="equation-9-2"></a>

$$
\boxed{
\mathcal P(x)
=\frac{x}{\sqrt{\pi L}}
\left[1+\frac{2-\gamma}{4L}+O\!\left(L^{-2}\right)\right]
}
\tag{9.2}
$$

를 얻는다. 잔차 $E$와 $\mathcal E_2$의 정의에는 이 유한 전개가 아니라 각각 §6과 §5의 전체 주항을 사용한다.

<a id="references"></a>

## 참고 문헌

1. Kohji Matsumoto, *Recent Developments in the Mean Square Theory of the Riemann Zeta and Other Zeta-Functions*, in *Number Theory* (2000), pp. 241–286. 특히 §1, 식 (1.1)의 $\sigma>1/2$ 평균제곱 공식. [나고야대학교 공개본](https://nagoya.repo.nii.ac.jp/records/18330), [DOI](https://doi.org/10.1007/978-3-0348-7023-8_14).
2. $s=1$ 분지점의 적분을 일반화된 약수합의 전체 주항으로 고정하는 배경: Hideki Nakaya, *On the Generalized Divisor Problem in Arithmetic Progressions*, *Sci. Rep. Kanazawa Univ.* 37(1) (1992), pp. 23–47, 특히 서론 p. 25. [논문](https://scirep.w3.kanazawa-u.ac.jp/articles/37-01-002.pdf). 본문 §6의 실수 적분은 하한과 상수항을 명시한 이 문서의 규약이다.
3. N. Levinson, *More than one third of zeros of Riemann's zeta-function are on $\sigma=1/2$*, Adv. Math. 13 (1974), 383–436; D. R. Heath-Brown, *Simple zeros of the Riemann zeta-function on the critical line*, Bull. London Math. Soc. 11 (1979), 17–18: Levinson의 방법으로 얻는 임계선 위 영점의 단순성.
4. J. B. Conrey, *More than two fifths of the zeros of the Riemann zeta function are on the critical line*, J. Reine Angew. Math. 399 (1989), 1–26.
5. H. L. Montgomery, R. C. Vaughan, *Multiplicative Number Theory I*, Cambridge (2007), §15.1: 비음 함수의 멜린 변환에 대한 Landau의 정리.
