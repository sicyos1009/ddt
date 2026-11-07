# 양의 약수분포의 이차 모멘트와 잔차 $\mathcal E_2$

$(1/2,1/2)$ 약수분포의 이차 모멘트 $M_2$와 전체 주항을 제거한 잔차 $\mathcal E_2$를 다룬다. 아래의 정확식, 단순극 상쇄와 로그 전개는 RH나 임계선 밖 단순영점의 부재를 가정하지 않는다. $\mathcal E_{2,\rho}$는 고정된 단순영점 하나의 기여이며, 전체 잔차의 성장 지수는 실제 단순영점들의 상한 $\Theta_{\rm simple}$로 표현한다.

| 결과 | 내용 |
|:--|:--|
| [§5](#section-5) | 전체 주항과 $\mathcal E_2$의 정확한 정규화 |
| [§8.2–8.3](#section-8-2) | 단순영점 기여의 정확식과 $\sigma,r$의 완전 상쇄 |
| [§8.4–8.5](#section-8-4) | 단순극을 보존하는 평균 추정과 임의 유한 차수의 로그 전개 |
| [정리 8.1](#section-8-1) | $\mathcal E_2$의 성장 지수 $=\Theta_{\rm simple}$ |
| [명제 5.1](#section-5-3) | $a_2(n)$의 비산술 부분 $b(n)=\sum_{d\mid n}\Lambda(d)\log d$의 점별 영점 표현 |
| [명제 10.1](#section-10-2) | $4\Delta_2=x\,\mathcal G+O(\sqrt x\log x)$, $\mathcal G$는 $\psi-x$의 두 이산–연속 차 |
| [명제 10.2](#section-10-3) | 영점 하나가 $\mathcal G$에 주는 기여 $\zeta'(\rho)x^{\rho-1}$, 중복 영점의 소멸 **[국소 계산]** |
| [따름정리 10.3–10.4](#section-10-4) | $x\mathcal G$의 성장 판정과 $\Theta_{\rm simple}$, 첫 번째 차 $R=\Psi_1-x\log x+x$의 판정과 RH |

- [§1–3. 약수분포와 이차 모멘트](#section-1)
- [§4–5. 부분합, 전체 주항, $b(n)$의 영점 표현](#section-4)
- [§6. 관련 잔차 $E$의 기호](#section-6)
- [§7. $\mathcal E_2$의 멜린 변환](#section-7)
- [§8. 단순영점의 기여와 전역 성장](#section-8)
- [§9. 전체 주항의 점근 전개](#section-9)
- [§10. ψ를 통한 이산–연속 분해](#section-10)
- [참고 문헌](#references)

자연로그를 쓰며, 비자명 영점 $\rho=\beta+i\gamma_\rho$의 중복도는 $m_\rho$로 쓴다. 절과 수식 번호는 이 문서 안에서 독립적이다.

<a id="section-1"></a>

## 1. 약수의 가중치

$f=d_{1/2}$는 곱셈적이고

$$
f(p^a)=\frac{\binom{2a}{a}}{4^a},\qquad
w_n(d):=f(d)f(n/d)\quad(d\mid n)
$$

이다. 생성함수 $(\sum_{a\ge0}f(p^a)z^a)^2=(1-z)^{-1}$에서

$$
w_n(d)>0,\qquad w_n(d)=w_n(n/d),\qquad
\sum_{d\mid n}w_n(d)=1
$$

을 얻는다. 이는 [주정리 §8.4](main-theorem.md#section-8-4)의 이산핵 $\mathscr K_n^{1/2,1/2}$이다.

<a id="section-2"></a>

## 2. 로그 좌표와 대칭성

$n\ge2$에서 $u_d=\log d/\log n$라 놓으면 $u_{n/d}=1-u_d$이므로

$$
\sum_{d\mid n}w_n(d)u_d=\frac12,\qquad
\mathbb E[\log D_n]=\frac12\log n,
\quad \mathbb P(D_n=d)=w_n(d).
$$

<a id="section-3"></a>

## 3. 이차 모멘트의 정확한 약수합

<a id="equation-3-1"></a>

$$
\boxed{
M_2(n):=\sum_{d\mid n}w_n(d)u_d^2
=\frac14+\frac{\operatorname{Var}(\log D_n)}{(\log n)^2},
\qquad \frac14\le M_2(n)\le\frac12.
}
\tag{3.1}
$$

$a_2(n):=(\log n)^2M_2(n)$, $a_2(1):=0$로 두면 또 다른 정확식은

$$
4a_2(n)=(\log n)^2+\sum_{d\mid n}\Lambda(d)\log d.
$$

여기서 $\Lambda$는 폰 망골트 함수이다. 이 식은 §5의 $4F_2=\zeta''+\zeta(\log\zeta)''$에서 계수를 비교하면 얻는다.

<a id="section-4"></a>

## 4. 부분합과 잔차

$$
T_2(x):=\sum_{2\le n\le x}M_2(n),\qquad
U_2(x):=\sum_{2\le n\le x}a_2(n),\qquad
\mathcal E_2(x):=T_2(x)-\mathcal P_2(x).
$$

전체 주항 $\mathcal P_2$에는 $3x/8$과 모든 로그 보정항을 포함한다. 그 상수항까지 §5의 하한 $2$ 규약으로 고정한다.

<a id="section-5"></a>

## 5. 전체 주항과 로그 정규화

<a id="section-5-1"></a>

### 5.1. 디리클레급수와 $s=1$의 극

$\Re s>1$에서 $F(s)=\sum f(n)n^{-s}=\sqrt{\zeta(s)}$이고
<a id="equation-5-1"></a>

$$
\boxed{
F_2(s):=\sum_{n\ge1}\frac{a_2(n)}{n^s}
=F(s)F''(s)
=\frac12\zeta''(s)-\frac14\frac{\zeta'(s)^2}{\zeta(s)}.
}
\tag{5.1}
$$

$\zeta(1+h)=h^{-1}+\gamma-\gamma_1h+\cdots$와 $c:=\gamma^2+3\gamma_1$를 사용하면
<a id="equation-5-2"></a>

$$
F_2(1+h)=\frac3{4h^3}+\frac{\gamma}{4h^2}-\frac{c}{4h}+O(1).
\tag{5.2}
$$

따라서
<a id="equation-5-3"></a>

$$
\boxed{
Q_2(x):=\operatorname*{Res}_{s=1}F_2(s)\frac{x^s}{s}
=x\left[\frac38(\log x)^2+\frac{\gamma-3}{4}\log x
+\frac{3-\gamma-c}{4}\right].
}
\tag{5.3}
$$

<a id="section-5-2"></a>

### 5.2. 전체 주항과 잔차의 정확식

$x\ge2$에서

$$
\mathscr L f(x):=\frac{f(x)}{(\log x)^2}
+2\int_2^x\frac{f(t)}{t(\log t)^3}\,dt,\qquad
\Delta_2(x):=U_2(x)-Q_2(x)
$$

로 정의한다. 부분합 공식으로 $T_2=\mathscr L U_2$이므로
<a id="equation-5-4"></a>

$$
\boxed{\mathcal P_2:=\mathscr L Q_2}
\tag{5.4}
$$

<a id="equation-5-5"></a>

$$
\boxed{\mathcal E_2=\mathscr L\Delta_2}
\tag{5.5}
$$

는 정확한 항등식이다. $I_j(x):=\int_2^x(\log t)^{-j}\,dt$라 쓰면
<a id="equation-5-6"></a>

$$
\boxed{
\mathcal P_2(x)=\frac38x+\frac\gamma4I_1(x)-\frac c4I_2(x)+\kappa_2,
\qquad \kappa_2=\frac{Q_2(2)}{(\log2)^2}-\frac34.
}
\tag{5.6}
$$

역공식도 정확히
<a id="equation-5-7"></a>

$$
\boxed{
\Delta_2(x)=(\log x)^2\mathcal E_2(x)
-2\int_2^x\mathcal E_2(t)\frac{\log t}{t}\,dt
}
\tag{5.7}
$$

이다. 하한에서는 $\Delta_2(2)=(\log2)^2\mathcal E_2(2)$이다.

<a id="section-5-3"></a>

### 5.3. $b(n)$의 점별 영점 표현

§3의 항등식에서 $(\log n)^2$을 뺀 부분을
<a id="equation-5-8"></a>

$$
b(n):=\sum_{d\mid n}\Lambda(d)\log d
=\frac12\sum_{p^a\parallel n}a(a+1)(\log p)^2,
\qquad
4a_2(n)=(\log n)^2+b(n)
\tag{5.8}
$$

로 쓴다. 두 번째 표현은 $\Lambda(p^k)\log p^k=k(\log p)^2$를 $k=1,\dots,a$에 대해 더한 것이다. $b(n)$은 [극한분포 §1](limiting-distribution.md#section-1)의 $c_n$과 같은 수이다. $\sum_n\Lambda(n)n^{-s}=-\zeta'/\zeta$를 미분하면 $\sum_n\Lambda(n)\log n\,n^{-s}=(\zeta'/\zeta)'$이고, 여기에 $\zeta$를 곱하면 $\Re s>1$에서

$$
\widehat b(s):=\sum_{n\ge1}\frac{b(n)}{n^s}
=\zeta(s)\Bigl(\frac{\zeta'}{\zeta}\Bigr)'(s)
=\zeta''(s)-\frac{\zeta'(s)^2}{\zeta(s)},
\qquad
4F_2=\zeta''+\widehat b
$$

이다. 두 번째 등식은 (5.1)과 같다.

**주항.** $\zeta''(1+h)=2h^{-3}+O(1)$을 (5.2)의 $4F_2$에서 빼면 $\widehat b(1+h)=h^{-3}+\gamma h^{-2}-ch^{-1}+O(1)$이다. $x^{1+h}/(1+h)=x\bigl[1+(L-1)h+(\tfrac12L^2-L+1)h^2+O(h^3)\bigr]$이므로
<a id="equation-5-9"></a>

$$
\boxed{
Q_b(x):=\operatorname*{Res}_{s=1}\widehat b(s)\frac{x^s}{s}
=x\left[\frac12L^2+(\gamma-1)L+1-\gamma-c\right].
}
\tag{5.9}
$$

같은 전개로 $\operatorname{Res}_{s=1}\zeta''(s)x^s/s=x(L^2-2L+2)$이고, 두 유수의 합 $4Q_2(x)=x(L^2-2L+2)+Q_b(x)$는 (5.3)과 일치한다. [극한분포 §1](limiting-distribution.md#section-1)의 주항은 $xP_2(\log x)=Q_b(x)$이다.

**경로 이동.** $S_b(x):=\sum_{n\le x}b(n)$으로 두고, 영점의 높이가 아닌 $T>0$과 정수 $K\ge0$을 고정한다. 직사각형 $[-2K-1,2]\times[-T,T]$의 경계에서 $\widehat b(s)x^s/s$를 반시계 방향으로 적분하고, 오른쪽 변의 적분을 $I_{\rm right}$, 나머지 세 변의 적분의 합을 $I_{\rm rest}$라 하자. 경로 안의 유수는 다음과 같다.

- $s=1$: $Q_b(x)$.
- $s=0$: $\zeta(0)=-\tfrac12\ne0$이므로 $\widehat b$는 정칙이고, 유수는 상수 $\widehat b(0)$이다.
- 비자명 영점 $\rho$: §8.2의 계산에서 $\zeta'^2/\zeta$는 $\rho$에서 $(s-\rho)^{m_\rho-2}$의 차수를 가지므로 $m_\rho\ge2$이면 $\widehat b$가 정칙이고, $m_\rho=1$이면 $\operatorname{Res}_{s=\rho}\widehat b=-\zeta'(\rho)$이다. 임계띠의 실수축에는 영점이 없으므로($0<\sigma<1$에서 $\zeta(\sigma)<0$) $\rho$와 $\bar\rho$를 묶어 $\gamma_\rho>0$만 센다.
- 자명 영점 $s=-2k$ ($1\le k\le K$): $\zeta''$는 정칙이고 $-\zeta'^2/\zeta$의 유수는 $-\zeta'(-2k)$이므로 유수는 $\zeta'(-2k)x^{-2k}/(2k)$이다.

따라서 $\mathcal R_{T,K}(x):=\bigl(S_b(x)-I_{\rm right}\bigr)-I_{\rm rest}$로 정의하면, $x\ge1$에서 정확히
<a id="equation-5-10"></a>

$$
S_b(x)=Q_b(x)+\widehat b(0)
-2\Re\sum_{\substack{0<\gamma_\rho\le T\\ m_\rho=1}}\frac{\zeta'(\rho)}{\rho}x^\rho
+\sum_{k=1}^{K}\frac{\zeta'(-2k)}{2k}x^{-2k}
+\mathcal R_{T,K}(x)
\tag{5.10}
$$

이다. $S_b(x)-I_{\rm right}$는 높이 $T$에서 자른 Perron 공식의 절단오차이다. 정수 $x$에서 Perron 적분은 $S_b(x)-\tfrac12b(x)$를 근사하므로, 끝점의 $\tfrac12b(x)$도 이 절단오차에 들어간다.

**명제 5.1 (가정 없음).** 정수 $n\ge2$에서
<a id="equation-5-11"></a>

$$
\boxed{
\begin{aligned}
b(n)
&=Q_b(n)-Q_b(n-1)
-2\Re\sum_{\substack{0<\gamma_\rho\le T\\ m_\rho=1}}
\frac{\zeta'(\rho)}{\rho}\bigl[n^{\rho}-(n-1)^{\rho}\bigr]\\
&\quad+\sum_{k=1}^{K}\frac{\zeta'(-2k)}{2k}\bigl[n^{-2k}-(n-1)^{-2k}\bigr]
+\mathcal R_{T,K}(n)-\mathcal R_{T,K}(n-1).
\end{aligned}
}
\tag{5.11}
$$

**증명.** (5.10)을 $x=n$과 $x=n-1$에 쓰고 빼면 $b(n)=S_b(n)-S_b(n-1)$이고 상수 $\widehat b(0)$은 소거된다. $n=2$일 때는 $S_b(1)=0$이다. $\square$

영점의 위치나 단순성은 가정하지 않았다. 단순영점은 임계선 밖에 있어도 그대로 들어가고, 중복 영점은 항을 남기지 않는다. §8.2의 $\Delta_{2,\rho}=c_\rho x^\rho/\rho$, $c_\rho=-\zeta'(\rho)/4$에서

$$
-\frac{\zeta'(\rho)}{\rho}x^\rho=4\Delta_{2,\rho}(x)
$$

이므로, (5.10)의 영점항은 $\Delta_{2,\rho}+\Delta_{2,\bar\rho}$를 4배한 것이고 (5.11)의 영점항은 그 차분이다. (5.8)의 $(\log n)^2$ 부분은 $\zeta''$에서 오며 영점항을 주지 않으므로, $a_2(n)$의 영점 성분은 모두 $b(n)$에서 나온다. (5.10)은 [극한분포 §1](limiting-distribution.md#section-1)의 명시적 공식에서 자명 영점의 $K$개 항과 상수 $\widehat b(0)$을 분리한 형태이다. 이 절에서는 $\mathcal R_{T,K}$의 크기를 추정하지 않으며, 절단오차의 평균제곱 추정은 [극한분포 §3](limiting-distribution.md#section-3)에 있다.

<a id="section-6"></a>

## 6. 관련 잔차 $E$의 기호

[짝수 중복도 영점 문서](even-zeros.md#notation)와의 기호 대응은 다음과 같다.

$$
A(x)=\sum_{n\le x}f(n),\quad
\mathcal P(x)=1+\frac1\pi\int_0^1
\sqrt{-\zeta(u)}\,\frac{x^u-1}{u}\,du,\quad
E(x)=A(x)-\mathcal P(x),
$$

$$
B(s)=\frac1\pi\int_0^1\frac{\sqrt{-\zeta(u)}}{s-u}\,du,\qquad
H(s)=F(s)-1-B(s).
$$

$E$의 정칙 연장과 성장에 관한 논의는 그 문서에 둔다. 아래에서는 $\mathcal E_2$를 다룬다.

<a id="section-7"></a>

## 7. $\mathcal E_2$의 멜린 변환

$1\le x<2$에는 $U_2(x)=0$, $\Delta_2(x)=-Q_2(x)$로 둔다. $\Re s>1$에서
<a id="equation-7-1"></a>

$$
\boxed{
s\int_1^\infty\Delta_2(x)x^{-s-1}\,dx
=F_2(s)-R_2(s),\qquad
R_2(s):=s\int_1^\infty Q_2(x)x^{-s-1}\,dx.
}
\tag{7.1}
$$

$P_1(s)$을 (5.2)의 주부분이라 하면 $R_2(s)=P_1(s)-P_1(0)$이다. 따라서 $F_2-R_2$는 $s=1$에서 정칙이고, $\Re s>0$에서 남는 극은 §8.2의 단순영점뿐이다.

<a id="section-8"></a>

## 8. 단순영점의 기여와 전역 성장

<a id="section-8-1"></a>

### 8.1. 가정 없는 성장 지수

**정리 8.1.** 실제 단순영점들의 상한을

$$
\Theta_{\rm simple}:=\sup\{\Re\rho:m_\rho=1\}
$$

이라 하면

$$
\boxed{
\inf\{\alpha:\mathcal E_2(x)=O(x^\alpha)\}
=\inf\{\alpha:\Delta_2(x)=O(x^\alpha)\}
=\Theta_{\rm simple},\qquad
\frac12\le\Theta_{\rm simple}\le1.
}
$$

동치로, 각 $\tfrac12\le\theta<1$에 대해

$$
\mathcal E_2(x)=O_\varepsilon(x^{\theta+\varepsilon})
\ \text{모든 }\varepsilon>0
\quad\Longleftrightarrow\quad
\Re\rho>\theta\text{인 단순영점이 없음}.
$$

고정된 단순영점 $\rho$, $\beta=\Re\rho$와 $0<\alpha<\beta$에 대해 $\mathcal E_2(x)\ne O(x^\alpha)$이다. 증명은 §8.7에 둔다.

<a id="section-8-2"></a>

### 8.2. 단순영점 기여의 정확식

$F_2$는 유리형이다. $\zeta(s)=a(s-\rho)^m(1+O(s-\rho))$이면

$$
\frac{\zeta'(s)^2}{\zeta(s)}
=m^2a(s-\rho)^{m-2}(1+O(s-\rho)).
$$

따라서 $m\ge2$에서는 $F_2$가 정칙이고, $m=1$에서는

$$
c_\rho:=\operatorname*{Res}_{s=\rho}F_2(s)
=-\frac14\zeta'(\rho),\qquad
\Delta_{2,\rho}(x):=\frac{c_\rho}{\rho}x^\rho.
$$

계수 $-1/4$는 [주정리 §7](main-theorem.md#section-7)의 국소 핵 $\langle\mathcal K_{-1/2,-1/2},u^2\rangle$와 일치한다.

<a id="simple-zero-exact-contribution"></a>

$L=\log x$, $\ell=\log2$라 쓰고 $\mathcal E_{2,\rho}:=\mathscr L\Delta_{2,\rho}$로 정의하면

$$
\boxed{
\mathcal E_{2,\rho}(x)
=\frac{\zeta'(\rho)}4
\left[
\frac{x^\rho}{L}-\frac{2^\rho}{\ell}
-\rho\bigl(\operatorname{Ei}(\rho L)-\operatorname{Ei}(\rho\ell)\bigr)
-\frac{2^\rho}{\rho\ell^2}
\right],\qquad x\ge2.
}
$$

여기서 $x^\rho=e^{\rho\log x}$이고 $\operatorname{Ei}$는 주가지이다. 이 식은

$$
\mathcal E_{2,\rho}(2)=-\frac{\zeta'(\rho)2^\rho}{4\rho\ell^2},
\qquad
\mathcal E_{2,\rho}'(x)=-\frac{\zeta'(\rho)}4
\frac{x^{\rho-1}}{(\log x)^2}
$$

를 만족한다. $\beta>1/2$인 단순영점에도 그대로 적용된다.

<a id="section-8-3"></a>

### 8.3. 평균항과 단순극항의 완전 상쇄

고정된 단순영점에 대해

$$
G_\rho(s):=\frac{c_\rho}{s-\rho}+\frac{c_\rho}{\rho}
=\frac{c_\rho s}{\rho(s-\rho)}
$$

를 떼어내면 $F_2-G_\rho$는 $\rho$에서 정칙이다. $r>0$인 원판 평균과 보정항은 정확히

$$
A_{r,\rho}(s)=\frac{c_\rho}{\rho}
+c_\rho\begin{cases}
\overline{s-\rho}/r^2,&|s-\rho|<r,\\
1/(s-\rho),&|s-\rho|\ge r,
\end{cases}
$$

$$
\mathfrak S_{r,\rho}(s)
=c_\rho\left(\frac1{s-\rho}-\frac{\overline{s-\rho}}{r^2}\right)
\mathbf1_{\{|s-\rho|<r\}},\qquad
G_\rho=A_{r,\rho}+\mathfrak S_{r,\rho}
$$

이다. $\sigma>0$, $\sigma\ne\beta$에서 대칭 극한으로

$$
P_{\rho,r,\sigma}(x):=\frac1{2\pi i}\int_{(\sigma)}
\mathfrak S_{r,\rho}(s)\frac{x^s}{s}\,ds,\qquad
B_{\rho,r,\sigma}(x):=\frac1{2\pi i}\int_{(\sigma)}
A_{r,\rho}(s)\frac{x^s}{s}\,ds
$$

를 정의한다. $G_\rho(s)x^s/s=c_\rho x^s/(\rho(s-\rho))$에 유수정리를 적용하면

$$
P_{\rho,r,\sigma}+B_{\rho,r,\sigma}
=\mathbf1_{\{\beta<\sigma\}}\Delta_{2,\rho}.
$$

따라서 경로 이동에서 얻은 유수까지 합친 정확식은

$$
\boxed{
\mathbf1_{\{\beta>\sigma\}}\Delta_{2,\rho}
+P_{\rho,r,\sigma}+B_{\rho,r,\sigma}
=\Delta_{2,\rho},
\qquad
\mathbf1_{\{\beta>\sigma\}}\mathcal E_{2,\rho}
+\mathscr L P_{\rho,r,\sigma}+\mathscr L B_{\rho,r,\sigma}
=\mathcal E_{2,\rho}.
}
$$

$P_{\rho,r,\sigma}$만 추정할 때 생기는 잔여항은 $B_{\rho,r,\sigma}$와 정확히 상쇄된다. 최종 기여에는 $\sigma,r$가 남지 않는다.

<a id="section-8-4"></a>
<a id="simple-zero-mean-decomposition"></a>

### 8.4. 단순극을 보존하는 무조건적 평균 추정

단순영점 근방의 $|F_2(\rho+w)|\sim|\zeta'(\rho)|/(4|w|)$는 면적 적분 가능하지만, 영점을 지나는 수직선의 절댓값 적분은 발산한다. 면적 평균에서 수직선 추정으로 넘어갈 때 단순극 보정항을 보존한다.

$0<r<1/4$, $\sigma-2r>1/2$, $\sigma\le2$, $s=\sigma+it$, $t\ge3$에서

$$
\mathfrak A_rF_2(s):=\frac1{\pi r^2}\iint_{|w-s|<r}F_2(w)\,dA(w),
\qquad
\mathfrak S_r(s):=\sum_{\substack{m_\rho=1\\|\rho-s|<r}}
\frac{c_\rho}{s-\rho}\left(1-\frac{|s-\rho|^2}{r^2}\right).
$$

원 위의 유수정리를 반지름에 대해 평균하면

$$
\boxed{
F_2=\mathfrak A_rF_2+\mathfrak S_r,\qquad
\int_T^{2T}|F_2(\sigma+it)-\mathfrak S_r(\sigma+it)|\,dt
\ll_{\sigma,r}T.
}
$$

첫 등식은 극 밖에서 성립하고, 차는 극에서 연속으로 연장한다.

**평균 추정의 증명.** 정칙 함수 $g$에 대해 분포 항등식

$$
\Delta|g|=\frac{|g'|^2}{|g|}
$$

이 성립한다. 영점 밖에서 직접 미분하고, $m$중 영점 주위의 반지름 $\varepsilon$인 원을 제거하면 경계항이 $O(\varepsilon^m)\to0$이므로 영점에서도 성립한다. 표준 평균제곱 공식 [1]과 코시–슈바르츠로 $\int_T^{2T}|\zeta(\sigma'+it)|\,dt\ll_{\sigma'}T$이다.

$S=[\sigma-r,\sigma+r]\times[T-r,2T+r]$를 덮는 매끄러운 절단함수 $\chi$를 $\Re s>1/2$ 안에서 잡으면

$$
\iint_S\frac{|\zeta'|^2}{|\zeta|}\,dA
\le\iint|\zeta|\,\Delta\chi\,dA\ll_{\sigma,r}T.
$$

$\zeta''$도 코시 공식으로 같은 면적 추정을 만족한다. 따라서

$$
\int_T^{2T}|\mathfrak A_rF_2(\sigma+it)|\,dt
\le\frac2{\pi r}\iint_S|F_2|\,dA\ll_{\sigma,r}T.
$$

단순영점의 부재는 쓰지 않았다. $\square$

<a id="section-8-5"></a>

### 8.5. 로그 전개와 정확한 잔여항

고정된 단순영점 $\rho$에 대해 $x\to\infty$일 때

$$
\boxed{
\mathcal E_{2,\rho}(x)
=-\frac{\zeta'(\rho)}{4\rho}\frac{x^\rho}{L^2}
\left[
1+\frac2{\rho L}+\frac6{(\rho L)^2}
+\frac{24}{(\rho L)^3}+\frac{120}{(\rho L)^4}
+\frac{720}{(\rho L)^5}+O_\rho(L^{-6})
\right].
}
$$

일반계수는 $(k+1)!$이다. §8.2의 미분식에 대입하면 연속한 계수 사이의 비가 $(k+2)/\rho$임을 얻는다.

임의의 정수 $N\ge0$에 대한 유한 전개의 정확식은

$$
\boxed{
\mathcal E_{2,\rho}(x)
=-\frac{\zeta'(\rho)x^\rho}{4}
\sum_{k=0}^{N}\frac{(k+1)!}{\rho^{k+1}L^{k+2}}
+C_{\rho,N}+R_{\rho,N}(x),
}
$$

$$
C_{\rho,N}
=\frac{\zeta'(\rho)2^\rho}{4}
\sum_{k=1}^{N}\frac{(k+1)!}{\rho^{k+1}\ell^{k+2}},
$$

$$
R_{\rho,N}(x)
=\frac{(-1)^N(N+2)!}{4}\rho\zeta'(\rho)
\left[\Gamma(-N-2,-\rho L)-\Gamma(-N-2,-\rho\ell)\right].
$$

$\Gamma(a,z)$는 상부 불완전 감마함수의 주가지이며, 빈 합은 $0$이다. 반복 부분적분 또는 직접 미분과 $x=2$의 경계값으로 확인할 수 있다. 고정된 $\rho,N$에서

$$
C_{\rho,N}+R_{\rho,N}(x)
=O_{\rho,N}\!\left(\frac{x^\beta}{L^{N+3}}\right).
$$

계수가 팩토리얼로 증가하므로 고정된 $x$에서 무한합으로 보내는 수렴급수는 아니다.

<a id="section-8-6"></a>

### 8.6. 전체 잔차와의 연결

유한한 단순영점 집합 $\mathcal Z$에 대해

$$
\Delta_{2,\mathcal Z}^{\rm reg}
:=\Delta_2-\sum_{\rho\in\mathcal Z}\Delta_{2,\rho},
\qquad
\mathcal E_{2,\mathcal Z}^{\rm reg}
:=\mathscr L\Delta_{2,\mathcal Z}^{\rm reg}
$$

라 놓으면 선형성에 의해

$$
\boxed{
\mathcal E_2
=\sum_{\rho\in\mathcal Z}\mathcal E_{2,\rho}
+\mathcal E_{2,\mathcal Z}^{\rm reg}.
}
$$

§8.3의 상쇄는 유한합에도 그대로 적용된다. 무한 영점합으로 넘기려면 합산 규약과 수렴 제어가 필요하며, 전체 $\mathcal E_2$의 상계에는 정칙 부분의 추정도 포함한다.

<a id="section-8-7"></a>

### 8.7. 정리 8.1의 증명

**하계.** $\mathcal E_2(x)=O(x^\alpha)$, $0<\alpha<\beta$이면 (5.7)에서 $\Delta_2(x)=O(x^\alpha(\log x)^2)$이다. 그러면 (7.1)의 적분이 $\Re s>\alpha$에서 정칙이므로 $F_2-R_2$도 정칙이어야 한다. 이는 $\rho$에서의 유수 $-\zeta'(\rho)/4\ne0$과 모순이다. $\Delta_2(x)=O(x^\alpha)$일 때도 같다.

**상계.** $\Theta_{\rm simple}<1$이면 $\Theta_{\rm simple}<\sigma_0<1$을 고정한다. 충분히 작은 $r$을 택하면 $\sigma\in[\sigma_0,2]$의 원판에는 단순영점이 없으므로, §8.4의 정확한 분해에서 $\mathfrak S_r=0$이고 $\int_T^{2T}|F_2(\sigma+it)|\,dt\ll T$이다. 이는 실제 영점 상한의 정의를 사용한 것이다.

반정수 $x$에서 $c_x=1+1/\log x$의 절단 Perron 적분을 $\sigma_0$로 옮긴다. $a_2(n)\ll(\log n)^2$이므로 절단오차는 $O(x(\log x)^3/T)$이다. 위 평균 추정을 $\sigma$에 대해서도 적분하면 $T\in[x^2,2x^2]$ 중 수평변 적분이 $O(1/x)$인 높이를 고를 수 있다. 왼쪽 수직변은 이진 구간별 평균 추정으로 $O(x^{\sigma_0}\log T)$이며, 경로 안의 유수는 $Q_2(x)$뿐이다. 따라서

$$
\Delta_2(x)\ll_{\sigma_0}x^{\sigma_0}\log x,\qquad
\mathcal E_2(x)\ll_{\sigma_0}\frac{x^{\sigma_0}}{\log x}.
$$

일반 $x$는 $\lfloor x\rfloor+1/2$로 옮기면 된다. $Q_2$의 변화는 $O((\log x)^2)$이다. $\sigma_0\downarrow\Theta_{\rm simple}$로 보내면 원하는 성장 지수를 얻는다. $\Theta_{\rm simple}=1$인 경우에는 $M_2(n)\le1/2$와 (5.6)에서 $\mathcal E_2\ll x$, $\Delta_2\ll x(\log x)^2$이다.

임계선 위에는 단순영점이 존재하므로 [2, 3] $\Theta_{\rm simple}\ge1/2$이다. (5.5), (5.7)은 $\Delta_2$와 $\mathcal E_2$의 성장 지수가 같음을 준다. $\square$

더 강하게, 각 단순영점 $\rho$와 $0<\alpha<\beta$에 대해 $\Delta_2(x)>x^\alpha$와 $\Delta_2(x)<-x^\alpha$가 각각 임의로 큰 $x$에서 성립한다. 그렇지 않으면 $x^\alpha\mp\Delta_2(x)$가 결국 비음이 되고, 그 멜린 변환에 대한 Landau의 정리 [4]가 $\rho$에서의 극과 모순된다.

<a id="section-9"></a>

## 9. 전체 주항의 점근 전개

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

앞의 두 항은 [주정리 (8.12)](main-theorem.md#equation-8-12)와 일치한다. $\mathcal E_2$의 정의에는 (5.4)의 전체 주항을 사용한다.

<a id="section-10"></a>

## 10. ψ를 통한 이산–연속 분해

$\Delta_2$를 소수 쪽 오차 $\Delta_\psi(x):=\psi(x)-x$, $\psi(x)=\sum_{n\le x}\Lambda(n)$로 직접 쓴다. 입력은 §5.3의 $b=(\Lambda\log)*1$ 하나이고, 이산–연속 차가 두 번 들어간다. $y=\sqrt x$, $\Psi_1(t):=\sum_{n\le t}\Lambda(n)\log n$으로 둔다.

<a id="section-10-1"></a>

### 10.1. 두 이산–연속 차

첫 번째 차는 소수 로그 질량 $\Lambda(n)\log n$을 연속 질량 $\log t\,dt$와 비교한다. $t\ge1$에서
<a id="equation-10-1"></a>

$$
\mathcal C(t):=\frac1t\left[\Psi_1(t)-\int_1^t\log u\,du\right],
\qquad
(\mathcal Wf)(t):=\frac1t\left[f(t)\log t-\int_1^t\frac{f(u)}u\,du\right]
\tag{10.1}
$$

로 두면
<a id="equation-10-2"></a>

$$
\boxed{\mathcal C=\mathcal W\Delta_\psi.}
\tag{10.2}
$$

**증명.** Stieltjes 부분적분으로 $\Psi_1(t)=\int_{1}^{t}\log u\,d\psi(u)=\psi(t)\log t-\int_1^t\psi(u)u^{-1}\,du$이다. 같은 식을 $\psi$ 대신 $u\mapsto u$에 적용하면 $\int_1^t\log u\,du=t\log t-\int_1^t du$이고, 두 식을 빼면 된다. $\square$

상수함수에 대해서는 $\mathcal W1=0$이다. 또 $\int_1^t\log u\,du=t\log t-t+1$이므로 $\Psi_1(t)=t\,\mathcal C(t)+t\log t-t+1$이다.

두 번째 차는 첫 번째 오차를 약수 방향으로 평균한다.
<a id="equation-10-3"></a>

$$
(\mathcal Df)(x):=\sum_{m\le\sqrt x}\frac{f(x/m)}m-\int_{\sqrt x}^{\infty}\frac{f(t)}t\,dt.
\tag{10.3}
$$

$x=e^{\tau}$로 쓰면 $1/m$ 가중 이산합과 로그 좌표의 연속 적분을 비교하는 연산이다. 두 차는 같은 꼴이지만, 첫 번째는 $\Lambda(n)\log n$ 질량을, 두 번째는 $1/m$ 가중치를 쓴다.

**수렴.** 고전적 소수정리 $\Delta_\psi(t)\ll t\exp(-\delta\sqrt{\log t})$ ($\delta>0$) [4, 6장]에서, (10.1)의 첫 항은 $\Delta_\psi(t)\log t/t$이고 둘째 항은 $\int_1^t$를 $\sqrt t$에서 나누어 추정하면 $\mathcal C(t)\ll\exp(-\delta'\sqrt{\log t})$이다. 따라서
<a id="equation-10-4"></a>

$$
\kappa_\Lambda:=\int_1^\infty\frac{\mathcal C(t)}t\,dt
\tag{10.4}
$$

는 절대수렴하고, $(\mathcal D\mathcal C)(x)$의 적분도 수렴한다.

<a id="section-10-2"></a>

### 10.2. 합성식

**명제 10.1 (가정 없음).** $\mathcal G:=\mathcal D\mathcal C=\mathcal D\mathcal W\Delta_\psi$로 두면 $x\ge2$에서
<a id="equation-10-5"></a>

$$
\boxed{
4\Delta_2(x)=x\,\mathcal G(x)+O(\sqrt x\log x),
\qquad
\mathcal E_2(x)=\frac14\mathscr L(\mathrm{id}\cdot\mathcal G)(x)+O\!\left(\frac{\sqrt x}{\log x}\right).
}
\tag{10.5}
$$

여기서 $\mathscr L(\mathrm{id}\cdot\mathcal G)(x)=x\mathcal G(x)/L^2+2\int_2^x\mathcal G(t)(\log t)^{-3}\,dt$이다. 또
<a id="equation-10-6"></a>

$$
\boxed{\kappa_\Lambda=-(\gamma^2+2\gamma_1).}
\tag{10.6}
$$

**증명.** (i) *로그 제곱 부분.* (5.8)과 $4Q_2=x(L^2-2L+2)+Q_b$에서

$$
4\Delta_2(x)=\Bigl[\sum_{2\le n\le x}(\log n)^2-x(L^2-2L+2)\Bigr]+\bigl[S_b(x)-Q_b(x)\bigr].
$$

$\int_1^x(\log t)^2\,dt=x(L^2-2L+2)-2$이고 $(\log t)^2$은 단조증가하므로 첫 괄호는 $O(L^2)$이다.

(ii) *쌍곡선법.* $S_b(x)=\sum_{dm\le x}\Lambda(d)\log d$에 분할점 $y=\sqrt x$의 디리클레 쌍곡선법을 쓰면

$$
S_b(x)=\sum_{m\le y}\Psi_1(x/m)+\sum_{d\le y}\Lambda(d)\log d\Bigl\lfloor\frac xd\Bigr\rfloor-\Psi_1(y)\lfloor y\rfloor
$$

이다. 바닥함수를 풀 때 생기는 $\sum_{d\le y}\Lambda(d)\log d\,\{x/d\}$와 $\{y\}\Psi_1(y)$는 Chebyshev 상계로 둘 다 $O(\Psi_1(y))=O(\sqrt x\log x)$이다. 남는 $d$ 쪽 합은 Stieltjes 적분으로 정확히

$$
x\sum_{d\le y}\frac{\Lambda(d)\log d}d=\frac xy\Psi_1(y)+x\int_1^y\frac{\Psi_1(t)}{t^2}\,dt
$$

이고, $x/y=y$이므로 첫 항이 $y\Psi_1(y)$와 상쇄된다. 따라서
<a id="equation-10-7"></a>

$$
S_b(x)=\sum_{m\le y}\Psi_1(x/m)+x\int_1^y\frac{\Psi_1(t)}{t^2}\,dt+O(\sqrt x\log x).
\tag{10.7}
$$

(iii) *이산–연속 분리.* (10.7)에 $\Psi_1(t)=t\,\mathcal C(t)+t\log t-t+1$을 넣고 $\int_1^y\mathcal C(t)\,dt/t=\kappa_\Lambda-\int_y^\infty\mathcal C(t)\,dt/t$를 쓰면, $\mathcal C$가 들어간 항은 정확히 $x\,\mathcal G(x)+x\kappa_\Lambda$로 모인다. 나머지는

$$
P(x):=\sum_{m\le y}\Bigl[\frac xm\log\frac xm-\frac xm+1\Bigr]+x\int_1^y\frac{t\log t-t+1}{t^2}\,dt
$$

이다. $\sum_{m\le y}1/m=\log y+\gamma+O(1/y)$, $\sum_{m\le y}\log m/m=\tfrac12(\log y)^2+\gamma_1+O(\log y/y)$와 $\log y=L/2$에서

$$
P(x)=x\Bigl[\frac12L^2+(\gamma-1)L+1-\gamma-\gamma_1\Bigr]+O(\sqrt x\log x)
$$

이다. 따라서 $S_b(x)=x\,\mathcal G(x)+x\kappa_\Lambda+P(x)+O(\sqrt x\log x)$이다.

(iv) *상수 $\kappa_\Lambda$.* $\Re s>1$에서 $\int_1^\infty\Psi_1(t)t^{-s-1}\,dt=s^{-1}(\zeta'/\zeta)'(s)$이고, $t\log t-t+1$의 같은 변환은 $(s-1)^{-2}-(s-1)^{-1}+s^{-1}$이므로

$$
\int_1^\infty\mathcal C(t)\,t^{-s}\,dt=\frac1s\Bigl(\frac{\zeta'}{\zeta}\Bigr)'(s)-\frac1{(s-1)^2}+\frac1{s-1}-\frac1s
$$

이다. 왼쪽은 (10.4)의 절대수렴과 지배수렴 정리로 $s\to1^+$에서 $\kappa_\Lambda$로 간다. $\log\zeta(1+h)=-\log h+\gamma h-(\gamma_1+\tfrac12\gamma^2)h^2+O(h^3)$에서 $(\zeta'/\zeta)'(1+h)=h^{-2}-(\gamma^2+2\gamma_1)+O(h)$이고, 오른쪽에 넣으면 $h^{-2}$와 $h^{-1}$ 항이 상쇄되어 (10.6)을 얻는다.

(v) *주항의 상쇄.* (10.6)과 $c=\gamma^2+3\gamma_1$에서 $\kappa_\Lambda+1-\gamma-\gamma_1=1-\gamma-c$이므로 $x\kappa_\Lambda+P(x)=Q_b(x)+O(\sqrt x\log x)$이다. (iii)과 합치면 $S_b(x)-Q_b(x)=x\,\mathcal G(x)+O(\sqrt x\log x)$이고, (i)과 합쳐 (10.5)의 첫 식을 얻는다.

(vi) *$\mathscr L$의 적용.* (5.5)와 선형성에서 $\mathcal E_2=\tfrac14\mathscr L(\mathrm{id}\cdot\mathcal G)+\tfrac14\mathscr L r$이고, $|r(t)|\ll\sqrt t\log t$이면

$$
\mathscr Lr(x)\ll\frac{\sqrt x}{L}+\int_2^x\frac{dt}{\sqrt t(\log t)^2}\ll\frac{\sqrt x}{L}.
\qquad\square
$$

오차 $O(\sqrt x\log x)$는 (ii)의 바닥함수에서만 나온다. 나머지 단계는 정확하거나 그보다 작은 오차를 준다. 이 오차는 임계선 위 영점 기여의 크기 $x^{1/2}$보다 크므로, (10.5)는 $1/2$보다 큰 성장 지수를 정확히 옮기지만 $\sqrt x$ 규모의 진동을 분해하지는 않는다(§10.4).

<a id="section-10-3"></a>

### 10.3. 영점 하나의 기여 [국소 계산]

$\psi$의 명시적 공식 [4, 12장]에서 중복도 $m_\rho$인 비자명 영점은 $\Delta_\psi$에 $-m_\rho t^\rho/\rho$를 준다. 따라서 이 절에서는 $t\mapsto t^\rho$ ($t\ge1$)에 $\mathcal D\mathcal W$를 적용한다.

**명제 10.2.** $\rho=\beta+i\gamma_\rho$, $0<\beta<1$을 고정하면 $x\to\infty$일 때
<a id="equation-10-8"></a>

$$
\boxed{
\mathcal D\mathcal W[t^\rho](x)
=x^{\rho-1}\Bigl(\zeta(\rho)L+\zeta'(\rho)-\frac{\zeta(\rho)}{\rho}\Bigr)
+O_\rho\bigl(x^{\beta/2-1}L\bigr).
}
\tag{10.8}
$$

특히 $\rho$가 영점이면 $\mathcal D\mathcal W[t^\rho](x)=\zeta'(\rho)x^{\rho-1}+O_\rho(x^{\beta/2-1}L)$이다.

**증명.** $\int_1^tu^{\rho-1}\,du=(t^\rho-1)/\rho$에서

$$
\mathcal W[t^\rho](t)=t^{\rho-1}\log t-\frac{t^{\rho-1}-t^{-1}}{\rho}.
$$

$\sigma>0$, $s\ne1$에서 정확히

$$
\sum_{m\le y}m^{-s}=\zeta(s)+\frac{y^{1-s}}{1-s}+E(s,y),
\qquad
E(s,y)=s\int_y^\infty\{t\}t^{-s-1}\,dt-\{y\}y^{-s}
$$

이다. 이를 $s$로 미분하면 $E$와 $\partial_sE$는 고정된 $\rho$에서 $O_\rho(y^{-\beta}\log y)$이다. $a=1-\rho$라 쓰면 $\Re a>0$이고, 위 식과 그 미분에서

$$
\sum_{m\le y}\frac{(x/m)^{\rho-1}\log(x/m)}m
=x^{\rho-1}\bigl(\zeta(\rho)L+\zeta'(\rho)\bigr)
+x^{\rho-1}y^{a}\Bigl(\frac{L-\log y}{a}+\frac1{a^2}\Bigr)
+O_\rho(x^{\beta-1}y^{-\beta}L)
$$

이다. $y^2=x$이므로 $x^{\rho-1}y^a=y^{-a}$이고 $L-\log y=\log y$이다. 한편

$$
\int_y^\infty t^{\rho-2}\log t\,dt=y^{-a}\Bigl(\frac{\log y}{a}+\frac1{a^2}\Bigr)
$$

이므로 경계항이 정확히 상쇄되어 $\mathcal D[t^{\rho-1}\log t](x)=x^{\rho-1}(\zeta(\rho)L+\zeta'(\rho))+O(x^{\beta/2-1}L)$이다. 같은 방법으로 $\mathcal D[t^{\rho-1}](x)=x^{\rho-1}\zeta(\rho)+O(x^{\beta/2-1})$이고, $\mathcal D[t^{-1}](x)=(\lfloor y\rfloor-y)/x=O(x^{-1})$이다. 세 식을 합치면 (10.8)이다. $\square$

따라서 중복도 $m_\rho$인 영점의 $x\mathcal G$에 대한 국소 기여는

$$
-m_\rho\frac{\zeta'(\rho)}{\rho}x^\rho+O_\rho(x^{\beta/2}L)
=\begin{cases}
4\Delta_{2,\rho}(x)+O_\rho(x^{\beta/2}L),&m_\rho=1,\\
O_\rho(x^{\beta/2}L),&m_\rho\ge2
\end{cases}
$$

이다. 이는 §8.2와 명제 10.1에 맞는다. 오차 $x^{\beta/2}L$은 (10.5)의 오차 규모 안에 있다. $\mathcal D$의 이산합은 영점항 $x^{\rho-1}\log x$를 $\zeta(\rho)L+\zeta'(\rho)$로 바꾼다. 앞의 인자는 모든 영점에서, 뒤의 인자는 중복 영점에서 사라진다. $\Delta_\psi$의 성장 지수는 모든 영점의 실수부 상한이지만 $x\mathcal G$는 단순영점만 보는 이유가 이 상쇄이다. $|\mathcal C|$를 항별로 추정하면 이 상쇄를 잃으므로, 성장률을 추정할 때는 합과 적분을 함께 유지해야 한다. 무한 영점합에서 $\mathcal D\mathcal W$와 합의 순서를 바꾸는 문제는 다루지 않는다.

<a id="section-10-4"></a>

### 10.4. 성장 판정

**따름정리 10.3 (가정 없음).** 각 $\tfrac12\le\theta<1$에 대해

$$
x\,\mathcal G(x)=O_\varepsilon(x^{\theta+\varepsilon})\ \text{모든 }\varepsilon>0
\quad\Longleftrightarrow\quad
\Re\rho>\theta\text{인 단순영점이 없음}.
$$

$\Theta_{\rm simple}>\tfrac12$이면 $\inf\{\alpha:x\mathcal G(x)=O(x^\alpha)\}=\Theta_{\rm simple}$이다.

**증명.** $\theta\ge\tfrac12$이므로 명제 10.1의 오차 $O(\sqrt x\log x)$는 $O(x^{\theta+\varepsilon})$에 들어간다. 따라서 $x\mathcal G=O(x^{\theta+\varepsilon})$과 $\Delta_2=O(x^{\theta+\varepsilon})$은 동치이고, 정리 8.1의 동치를 $\Delta_2$에 대해 쓰면 된다($\Delta_2$와 $\mathcal E_2$의 성장 지수는 같다). 둘째 주장은 $\alpha>\tfrac12$에서 같은 논법을 쓴 것이다. $\square$

두 번째 차 $\mathcal D$를 거치기 전의 첫 번째 오차 $\mathcal C$로 같은 판정을 하면, 이번에는 모든 영점이 대상이 된다.

**따름정리 10.4 (가정 없음).**
<a id="equation-10-9"></a>

$$
R(x):=\Psi_1(x)-x\log x+x=x\,\mathcal C(x)+1
\tag{10.9}
$$

로 두면 각 $\tfrac12\le\theta<1$에 대해

$$
R(x)=O_\varepsilon(x^{\theta+\varepsilon})\ \text{모든 }\varepsilon>0
\quad\Longleftrightarrow\quad
\Re\rho>\theta\text{인 비자명 영점이 없음}
$$

이다. 영점의 중복도는 따지지 않는다. 특히

$$
\boxed{R(x)=O_\varepsilon(x^{1/2+\varepsilon})\ \text{모든 }\varepsilon>0\quad\Longleftrightarrow\quad\text{RH}.}
$$

**증명.** (⇐) $\Re\rho>\theta$인 영점이 없으면 명시적 공식 [4, 12장]과 $N(T+1)-N(T)\ll\log T$에서 $\Delta_\psi(x)\ll x^\theta(\log x)^2$이다. (10.2)를 풀어 쓰면

$$
R(x)-1=x\,\mathcal C(x)=\Delta_\psi(x)\log x-\int_1^x\frac{\Delta_\psi(t)}t\,dt
$$

이므로 $R(x)\ll x^\theta(\log x)^3$이다.

(⇒) $R(x)=O(x^{\theta+\varepsilon})$이면 $\mathcal C(t)=(R(t)-1)/t=O(t^{\theta+\varepsilon-1})$이므로 $\int_1^\infty\mathcal C(t)t^{-s}\,dt$는 $\Re s>\theta+\varepsilon$에서 절대수렴하고 정칙이다. §10.2 (iv)에서 이 적분은 $\Re s>1$에서

$$
\frac1s\Bigl(\frac{\zeta'}{\zeta}\Bigr)'(s)-\frac1{(s-1)^2}+\frac1{s-1}-\frac1s
$$

와 같다. 따라서 이 함수는 $\Re s>\theta+\varepsilon$으로 정칙하게 연장된다. $s=1$에서는 $(\zeta'/\zeta)'(1+h)=h^{-2}+O(1)$이므로 극이 상쇄된다. 한편 $m$중 영점 $\rho$ 근방에서는 $\log\zeta(s)=m\log(s-\rho)+(\text{정칙})$이므로

$$
\Bigl(\frac{\zeta'}{\zeta}\Bigr)'(s)=(\log\zeta(s))''=-\frac{m}{(s-\rho)^2}+O(1)
$$

이고, $m\ge1$이면 이중극이 남는다. $1/s$는 $\rho$에서 $0$이 아니므로 $\Re\rho>\theta+\varepsilon$인 영점은 없다. $\varepsilon\to0$으로 보내면 된다. $\square$

따름정리 10.3과 10.4를 나란히 놓으면 두 이산–연속 차의 역할이 갈린다. 첫 번째 차까지의 $\mathcal C$는 $(\log\zeta)''$의 이중극 $-m_\rho/(s-\rho)^2$ 때문에 중복도와 관계없이 모든 영점을 보고, 판정 대상은 $\Theta:=\sup\Re\rho$이다($\theta=\tfrac12$이면 RH). 두 번째 차 $\mathcal D$까지 거친 $\mathcal G$는 명제 10.2의 상쇄 때문에 단순영점만 보고, 판정 대상은 $\Theta_{\rm simple}\le\Theta$이다. 이는 [영점 기여 §4](zero-contributions.md#section-4)에서 $\psi-x$가 중복 영점도 계수 $-m_\rho/\rho$로 포착한다는 관찰과 맞는다.

<a id="references"></a>

## 참고 문헌

1. Kohji Matsumoto, *Recent Developments in the Mean Square Theory of the Riemann Zeta and Other Zeta-Functions*, in *Number Theory* (2000), pp. 241–286, §1, 식 (1.1). [공개본](https://nagoya.repo.nii.ac.jp/records/18330), [DOI](https://doi.org/10.1007/978-3-0348-7023-8_14).
2. N. Levinson, *More than one third of zeros of Riemann's zeta-function are on $\sigma=1/2$*, Adv. Math. 13 (1974), 383–436; D. R. Heath-Brown, *Simple zeros of the Riemann zeta-function on the critical line*, Bull. London Math. Soc. 11 (1979), 17–18.
3. J. B. Conrey, *More than two fifths of the zeros of the Riemann zeta function are on the critical line*, J. Reine Angew. Math. 399 (1989), 1–26.
4. H. L. Montgomery, R. C. Vaughan, *Multiplicative Number Theory I*, Cambridge (2007): 6장(소수정리의 오차항), 12장(명시적 공식), §15.1(Landau의 정리).
