# 양의 약수분포의 이차 모멘트와 잔차 $\mathcal E_2$

$(1/2,1/2)$ 약수분포의 이차 모멘트 $M_2$와 전체 주항을 제거한 잔차 $\mathcal E_2$를 다룬다. 아래의 정확식, 단순극 상쇄와 로그 전개는 RH나 임계선 밖 단순영점의 부재를 가정하지 않는다. $\mathcal E_{2,\rho}$는 고정된 단순영점 하나의 기여이며, 전체 잔차의 성장 지수는 실제 단순영점들의 상한 $\Theta_{\rm simple}$로 표현한다.

| 결과 | 내용 |
|:--|:--|
| [§5](#section-5) | 전체 주항과 $\mathcal E_2$의 정확한 정규화 |
| [§8.2–8.3](#section-8-2) | 단순영점 기여의 정확식과 $\sigma,r$의 완전 상쇄 |
| [§8.4–8.5](#section-8-4) | 단순극을 보존하는 평균 추정과 임의 유한 차수의 로그 전개 |
| [정리 8.1](#section-8-1) | $\mathcal E_2$의 성장 지수 $=\Theta_{\rm simple}$ |

- [§1–3. 약수분포와 이차 모멘트](#section-1)
- [§4–5. 부분합과 전체 주항](#section-4)
- [§6. 관련 잔차 $E$의 기호](#section-6)
- [§7. $\mathcal E_2$의 멜린 변환](#section-7)
- [§8. 단순영점의 기여와 전역 성장](#section-8)
- [§9. 전체 주항의 점근 전개](#section-9)
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

<a id="references"></a>

## 참고 문헌

1. Kohji Matsumoto, *Recent Developments in the Mean Square Theory of the Riemann Zeta and Other Zeta-Functions*, in *Number Theory* (2000), pp. 241–286, §1, 식 (1.1). [공개본](https://nagoya.repo.nii.ac.jp/records/18330), [DOI](https://doi.org/10.1007/978-3-0348-7023-8_14).
2. N. Levinson, *More than one third of zeros of Riemann's zeta-function are on $\sigma=1/2$*, Adv. Math. 13 (1974), 383–436; D. R. Heath-Brown, *Simple zeros of the Riemann zeta-function on the critical line*, Bull. London Math. Soc. 11 (1979), 17–18.
3. J. B. Conrey, *More than two fifths of the zeros of the Riemann zeta function are on the critical line*, J. Reine Angew. Math. 399 (1989), 1–26.
4. H. L. Montgomery, R. C. Vaughan, *Multiplicative Number Theory I*, Cambridge (2007), §15.1: Landau의 정리.
