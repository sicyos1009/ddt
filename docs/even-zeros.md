# 짝수 중복도 영점과 $d_{1/2}$ 오차항

$\zeta(s)^{1/2}$의 계수 $f=d_{1/2}$의 부분합에서 $s=1$을 포함하는 실수축 분기절단의 기여 전체를 실수 적분으로 빼고, 남은 잔차 $E$를 다룬다. 주정리의 지수로는 $(z,w)=(\tfrac12,0)$, 즉 $r=\tfrac12$인 한쪽 성분의 총질량이다([주정리 §8](main-theorem.md#section-8)). 오른쪽 반평면의 영점이 모두 짝수 중복도이면 $E$가 제곱근 수준으로 작고(정리 1), 역도 성립하며(따름정리 3), $E$의 성장 지수는 정확히 홀수 중복도 영점의 실수부 상한이다(정리 2). 영점 중복도에 관한 조건 자체는 증명하지 않는다. 절·정리·수식 번호는 이 문서 안에서 독립적이다.

| 결과 | 내용 | 범위 |
|:--|:--|:--|
| [정리 1](#section-1) | $(\mathrm{EM}_\theta)$ 아래 $E(x)\ll x^\theta\log x$, $\theta=\tfrac12$이면 $x^{1/2}(\log x)^{5/4}$ | 조건부 |
| [정리 2](#section-1) | $E$의 성장 지수 $=\Theta_{\rm odd}$, 홀수 영점마다 $\Omega_\pm$ | 가정 없음 |
| [따름정리 3](#section-1) | $(\mathrm{EM}_\theta)\iff E(x)=O(x^{\theta+\eta})\ \forall\eta>0$ | 가정 없음 |
| [(10.1)](#equation-10-1) | 홀수 영점 하나의 국소 기여와 $1/\Gamma(-m/2)$ | **[국소 계산]** |

- [§1. 정의와 주정리](#section-1)
- [§2. 연속 주항과 멜린 변환](#section-2)
- [§3. 짝수 중복도와 해석적 제곱근](#section-3)
- [§4. 실수축 분기절단의 정확한 상쇄](#section-4)
- [§5. 균등한 평균제곱 추정](#section-5)
- [§6. 절단 Perron 공식과 수평변](#section-6)
- [§7. 정리 1의 증명](#section-7)
- [§8. 역방향: 따름정리 3](#section-8)
- [§9. 정확한 지수와 진동: 정리 2](#section-9)
- [§10. 홀수 영점의 국소 기여](#section-10)
- [참고 문헌](#references)

<a id="notation"></a>

**기호와 전제.** 자연로그를 쓴다. $\mathcal H_\theta=\{s:\Re s>\theta\}$이다. $G(s)=(s-1)\zeta(s)$로 두고 $G(1)=1$로 제거가능 특이점을 채운다. 비자명 영점 $\rho$의 중복도를 $m_\rho$로 쓴다. 아래의 $B,H,\mathcal P,E$는 [이차 모멘트 문서](m2-divisor-moment.md#section-6)와 같은 정의이다. 이 문서의 이전 판(PDF 노트)과는 $I(x)=p(x)$, $A(s)=B(s)$, $R(s)=H(s)$, $D(x)=E(x)$로 대응한다.

<a id="section-1"></a>

## 1. 정의와 주정리

$\Re s>1$에서 양의 실수축의 값이 양수가 되는 가지를 택하고

<a id="equation-1-1"></a>

$$
F(s)=\zeta(s)^{1/2}=\sum_{n\ge1}\frac{f(n)}{n^s},\qquad
f(p^a)=\frac{1}{4^a}\binom{2a}{a}
\tag{1.1}
$$

로 $f=d_{1/2}$를 정한다. 오일러 곱과 $(1-v)^{-1/2}$의 이항급수에서 $f$는 곱셈적이고 $0<f(n)\le1$, $f(1)=1$이다.

실수 구간 $0\le u<1$에서 $\zeta(u)<0$이다. 실제로 $\zeta(0)=-\tfrac12$이고, $0<u<1$에서는 교대급수 $\eta(u)=\sum(-1)^{n-1}n^{-u}>0$과 $\zeta(u)=\eta(u)/(1-2^{1-u})$, $1-2^{1-u}<0$에서 따른다. 양의 실수 제곱근으로

<a id="equation-1-2"></a>

$$
w(u)=\frac{\sqrt{-\zeta(u)}}{\pi}\quad(0\le u<1),\qquad
w(u)=\frac{1+O(1-u)}{\pi\sqrt{1-u}}\quad(u\uparrow1),\qquad
W:=\int_0^1w(u)\,du<\infty
\tag{1.2}
$$

로 둔다. 둘째 식은 $s=1$의 극 $\zeta(u)=1/(u-1)+O(1)$에서 나온다. 끝점합과 그 전체 주항을

<a id="equation-1-3"></a>

$$
A(x)=\sum_{n\le x}f(n),\qquad
\mathcal P(x)=1+\int_0^1w(u)\,\frac{x^u-1}{u}\,du,\qquad
E(x)=A(x)-\mathcal P(x)\qquad(x\ge1)
\tag{1.3}
$$

로 정의한다. $u=0$에서 $(x^u-1)/u$는 연속 연장값 $\log x$로 읽는다. 주항의 밀도는 $p(x)=\mathcal P'(x)=\int_0^1w(u)x^{u-1}\,du$이며 양의 감소함수이고 $p(1)=W$이다.

**가정.** $\tfrac12\le\theta<1$을 고정한다.

<a id="equation-1-4"></a>

$$
(\mathrm{EM}_\theta)\qquad \Re\rho>\theta\text{인 }\zeta\text{의 모든 영점은 짝수 중복도이다.}
\tag{1.4}
$$

경계선 $\Re\rho=\theta$ 위 영점의 중복도는 제한하지 않는다. 또한

$$
\Theta_{\rm odd}:=\sup\{\Re\rho:\ \zeta(\rho)=0,\ m_\rho\text{ 홀수}\}
$$

로 둔다. §9에서 보이듯 $\Theta_{\rm odd}\ge\tfrac12$이다.

**정리 1 (조건부 상계).** $(\mathrm{EM}_\theta)$를 가정하면 $x\ge3$에서

<a id="equation-1-5"></a>

$$
E(x)\ll_\theta
\begin{cases}
x^\theta\log x, & \tfrac12<\theta<1,\\[2pt]
x^{1/2}(\log x)^{5/4}, & \theta=\tfrac12.
\end{cases}
\tag{1.5}
$$

**정리 2 (정확한 지수와 진동).** 가정 없이 다음이 성립한다.

1. $\inf\{\alpha:\ E(x)=O(x^\alpha)\}=\Theta_{\rm odd}$.
2. $\rho$가 홀수 중복도 영점이고 $\beta=\Re\rho$이면, 모든 $0<\alpha<\beta$에 대해 $E(x)>x^\alpha$인 $x$와 $E(x)<-x^\alpha$인 $x$가 각각 임의로 큰 값에서 존재한다. 특히 $E(x)=\Omega_\pm(x^{1/2-\varepsilon})$이다.

**따름정리 3 (동치).** $\tfrac12\le\theta<1$에서 다음은 동치이다.

- (a) $(\mathrm{EM}_\theta)$.
- (b) 모든 $\eta>0$에 대해 $E(x)=O_{\theta,\eta}(x^{\theta+\eta})$.

정리 1은 두 점근 전개를 형식적으로 빼서 얻는 결과가 아니다. $s=1$을 포함하는 실수축 분기절단을 변환 수준에서 정확히 제거하고(§4), 왼쪽으로 옮긴 경로의 적분까지 평가한다(§6–7). 로그 지수의 최적성은 주장하지 않는다.

<a id="section-2"></a>

## 2. 연속 주항과 멜린 변환

**보조정리 2.1 (정수 표본과 반정수 이동).** 정수 $N\ge2$에 대해

<a id="equation-2-1"></a>

$$
\sum_{n=2}^N\bigl(f(n)-p(n)\bigr)=E(N)+c_N,\qquad
c_N:=\int_1^Np(t)\,dt-\sum_{n=2}^Np(n)\in[0,W].
\tag{2.1}
$$

또한 모든 $x\ge1$에서 $|E(x)-E(\lfloor x\rfloor+\tfrac12)|\le W/2$이다.

**증명.** $p$가 감소하므로 $p(n)\le\int_{n-1}^np\le p(n-1)$이고, $n=2,\dots,N$에 대해 합하면 $0\le c_N\le p(1)-p(N)\le W$이다. $A(N)-1=\sum_{n=2}^Nf(n)$, $\mathcal P(N)-1=\int_1^Np$이므로 (2.1)이 따른다. $x$와 $\lfloor x\rfloor+\tfrac12$ 사이에는 정수가 없어 $A$가 같고, $|\mathcal P(x)-\mathcal P(\lfloor x\rfloor+\tfrac12)|\le\tfrac12\sup p=W/2$이다. $\square$

따라서 정수에서 $p(n)$을 빼고 유계인 수열을 더 빼는 이전 판의 정의는 $E$와 $O(1)$만 다르며, 아래의 모든 결론은 그대로 옮겨진다. 특히 $c_N$을 (2.1)의 값으로 잡으면 정확히 $E(N)$이다.

$s\notin[0,1]$에서 코시형 변환

<a id="equation-2-2"></a>

$$
B(s)=\int_0^1\frac{w(u)}{s-u}\,du,\qquad
|B(s)|\le\frac{W}{|\Im s|}
\tag{2.2}
$$

는 $\mathbb C\setminus[0,1]$에서 정칙이다. 부등식은 $|s-u|\ge|\Im s|$에서 나온다.

**보조정리 2.2 (멜린 변환).** $\Re s>1$에서

<a id="equation-2-3"></a>

$$
s\int_1^\infty E(x)x^{-s-1}\,dx=H(s):=F(s)-1-B(s).
\tag{2.3}
$$

**증명.** 절대수렴 영역에서 $s\int_1^\infty A(x)x^{-s-1}dx=F(s)$이다. $\mathcal P(x)-1=\int_1^xp$이고 $\mathcal P(x)=O(x)$이므로 부분적분과 토넬리 정리로

$$
s\int_1^\infty\mathcal P(x)x^{-s-1}dx
=1+\int_1^\infty p(x)x^{-s}dx
=1+\int_0^1w(u)\int_1^\infty x^{u-1-s}\,dx\,du=1+B(s).
\qquad\square
$$

<a id="section-3"></a>

## 3. 짝수 중복도와 해석적 제곱근

**보조정리 3.1 (반평면의 가지).** $(\mathrm{EM}_\theta)$ 아래 $G$는 $\mathcal H_\theta$에서 정칙인 제곱근 $\varphi$를 가지며, $\varphi(s)>0$ $(s\in(\theta,\infty))$, $\varphi(1)=1$로 정규화할 수 있다. 그러면

<a id="equation-3-1"></a>

$$
F(s)=\frac{\varphi(s)}{(s-1)^{1/2}}\qquad\bigl(s\in\mathcal H_\theta\setminus(\theta,1]\bigr)
\tag{3.1}
$$

는 (1.1)의 해석적 연속이다. $(s-1)^{1/2}$는 $s>1$에서 양수인 주가지이다. 또한 $\theta<u<1$에서

<a id="equation-3-2"></a>

$$
w(u)=\frac{\varphi(u)}{\pi\sqrt{1-u}}.
\tag{3.2}
$$

**증명.** $G$는 $\mathcal H_\theta$에서 정칙이고, 그 영점은 $\Re\rho>\theta$인 비자명 영점과 같다. $\zeta$는 $\Re s\ge1$과 실수 구간 $(0,1)$에서 영이 아니고 $G(1)=1$이기 때문이다. 가정에 의해 이 영점들은 모두 짝수 중복도이다.

*정의와 경로 무관성.* $Z_G$를 $G$의 영점 집합이라 하자. $s\in\mathcal H_\theta\setminus Z_G$와 $2$에서 $s$로 가는 $\mathcal H_\theta\setminus Z_G$ 안의 경로 $\gamma$에 대해

$$
\varphi(s)=\sqrt{G(2)}\exp\Big(\frac12\int_\gamma\frac{G'}{G}\Big)
$$

로 둔다. 두 경로의 차이는 $\mathcal H_\theta\setminus Z_G$ 안의 닫힌 경로 $\Gamma$이다. $\mathcal H_\theta$가 단일연결이므로 $\Gamma$는 $\mathcal H_\theta$에서 영호몰로그이고, 편각원리의 일반형에 의해

$$
\frac12\int_\Gamma\frac{G'}{G}=\pi i\sum_{\rho}n(\Gamma,\rho)\,m_\rho\in2\pi i\,\mathbb Z
$$

이다. $n(\Gamma,\rho)$는 회전수이며 영이 아닌 항은 유한 개이고, 각 $m_\rho$가 짝수이므로 마지막 소속이 성립한다. 따라서 $\varphi$는 잘 정의된 정칙함수이고, $(\varphi^2/G)'=0$과 $s=2$의 값에서 $\varphi^2=G$이다.

*영점에서의 연장.* $2k$중 영점 $\rho$의 작은 원판 $D$에서 $G=(s-\rho)^{2k}g$, $g\ne0$이다. $D$는 볼록이므로 $g$의 정칙 제곱근 $\sqrt g$가 있다. $D\setminus\{\rho\}$에서 $(\varphi/(s-\rho)^k)^2=g$이므로 연결 집합 위의 연속함수 $\varphi/((s-\rho)^k\sqrt g)$는 $\pm1$ 중 한 값이다. 따라서 $\varphi=\pm(s-\rho)^k\sqrt g$로 $\rho$까지 정칙하게 연장된다.

*실수축에서의 양수성.* $\varphi^*(s):=\overline{\varphi(\bar s)}$도 $\mathcal H_\theta$에서 정칙이고, $G$가 실수축에서 실수이므로 $(\varphi^*)^2=G$, $\varphi^*(2)=\varphi(2)>0$이다. 두 제곱근의 비 $\varphi^*/\varphi$는 연결 집합 $\mathcal H_\theta\setminus Z_G$에서 $\pm1$ 값의 연속함수이므로 $\varphi^*=\varphi$이다. 따라서 $\varphi$는 $(\theta,\infty)$에서 실수이다. 이 구간에서 $G(u)=(u-1)\zeta(u)>0$이므로 $\varphi\ne0$이고, 연결 구간에서 부호가 일정하므로 $\varphi>0$, $\varphi(1)=\sqrt{G(1)}=1$이다.

*연속과 (3.2).* $\Re s>1$에서 (3.1)의 우변은 제곱이 $\zeta$이고 $s>1$에서 양수이므로 $F$와 같다. $\theta<u<1$에서 $\varphi(u)=\sqrt{G(u)}=\sqrt{1-u}\,\sqrt{-\zeta(u)}$이므로 (3.2)가 성립한다. $\square$

<a id="section-4"></a>

## 4. 실수축 분기절단의 정확한 상쇄

경곗값 이론 없이 상쇄를 보이기 위해 다음 함수를 쓴다.

<a id="equation-4-1"></a>

$$
\Phi(\zeta_*)=\int_0^1\frac{dt}{1+\zeta_*t^2}.
\tag{4.1}
$$

$1+\zeta_*t^2=0$인 $t\in[0,1]$이 있는 것은 $\zeta_*\in(-\infty,-1]$일 때뿐이므로 $\Phi$는 $\mathbb C\setminus(-\infty,-1]$에서 정칙이다. $y>0$이면 $\Phi(y^2)=\arctan(y)/y$이다.

**보조정리 4.1 (상쇄된 변환의 연장).** $(\mathrm{EM}_\theta)$ 아래 $H=F-1-B$는 $\mathcal H_\theta$ 전체로 정칙하게 연장된다. 특히 $s=1$의 특이점과 $(\theta,1)$의 분기 불연속이 모두 제거된다.

**증명.** $b\in(\theta,1)$을 고정하고 $\ell=1-b$, $B=B_{0,b}+B_{b,1}$로 적분 구간을 나눈다. $B_{0,b}$는 $\mathbb C\setminus[0,b]$에서 정칙이다. $s\in\mathcal H_\theta\setminus(\theta,1]$에서

<a id="equation-4-2"></a>

$$
F(s)-B_{b,1}(s)
=\frac{2\varphi(s)}{\pi\sqrt\ell}\,\Phi\!\Big(\frac{s-1}{\ell}\Big)
+\frac1\pi\int_b^1\frac{\varphi(s)-\varphi(u)}{(s-u)\sqrt{1-u}}\,du
\tag{4.2}
$$

임을 보인다.

*실수 $s>1$에서의 확인.* $z=s-1>0$으로 두면 (3.2)에서

$$
B_{b,1}(s)=\frac{\varphi(s)}\pi\int_b^1\frac{du}{(s-u)\sqrt{1-u}}
-\frac1\pi\int_b^1\frac{\varphi(s)-\varphi(u)}{(s-u)\sqrt{1-u}}\,du
$$

이다. $v=\sqrt{1-u}$로 치환하면

$$
\int_b^1\frac{du}{(s-u)\sqrt{1-u}}
=2\int_0^{\sqrt\ell}\frac{dv}{z+v^2}
=\frac{2}{\sqrt z}\arctan\sqrt{\frac\ell z}
=\frac{\pi}{\sqrt z}-\frac{2}{\sqrt z}\arctan\sqrt{\frac z\ell}
=\frac{\pi}{\sqrt z}-\frac{2}{\sqrt\ell}\,\Phi\!\Big(\frac z\ell\Big)
$$

이다. 셋째 등호는 양수 $y$에서 $\arctan(1/y)+\arctan y=\pi/2$이다. $F(s)=\varphi(s)/\sqrt z$를 빼면 (4.2)를 얻는다.

*양변의 정칙성.* 좌변은 $\mathcal H_\theta\setminus(\theta,1]$에서 정칙이다. 우변 첫 항은 $(s-1)/\ell\in(-\infty,-1]$, 즉 $s\in(-\infty,b]$에서만 특이할 수 있으므로 $\mathcal H_\theta\setminus(\theta,b]$에서 정칙이다. 둘째 항의 피적분함수는

$$
\frac{\varphi(s)-\varphi(u)}{s-u}=\int_0^1\varphi'(u+\tau(s-u))\,d\tau
$$

로 $s$에 관해 정칙이고 $(s,u)$에 관해 연속이다. 콤팩트 집합 $K\subset\mathcal H_\theta$에 대해 $K\cup[b,1]$의 볼록포는 $\mathcal H_\theta$의 콤팩트 부분집합이므로 이 몫은 $\sup|\varphi'|$로 균등 유계이고, $(1-u)^{-1/2}$는 적분 가능하다. 따라서 둘째 항은 $\mathcal H_\theta$ 전체에서 정칙이다.

$\mathcal H_\theta\setminus(\theta,1]$은 연결이고 (4.2)의 양변은 $(1,\infty)$에서 일치하므로, 항등정리에 의해 이 집합 전체에서 일치한다. 그러므로 $H=F-1-B_{0,b}-B_{b,1}$은 $\mathcal H_\theta\setminus(\theta,b]$로 정칙하게 연장된다. 서로 다른 $b$에서의 연장은 실수축 밖에서 모두 $H$와 같으므로 서로 일치하고, $b\in(\theta,1)$이 임의이므로 $H$는 $\bigcup_b(\mathcal H_\theta\setminus(\theta,b])=\mathcal H_\theta$에서 정칙이다. $\square$

이 상쇄는 $x/(\log x)^j$ 꼴의 항을 유한 개 빼는 근사가 아니라 변환 $H$의 정확한 해석적 연속이다. 짝수 중복도 영점 $\rho$에서도 $H$는 정칙이므로, 작은 원 $C_\rho$에 대해 $\oint_{C_\rho}H(s)x^s s^{-1}ds=0$이다. 이는 경로적분에 관한 등식이며 $H(\rho)=0$이라는 주장이 아니다.

**보조정리 4.2 (실수축 근방, 가정 없음).** $0<\theta_0<1$이면 어떤 $\delta\in(0,\theta_0/4)$에 대해 $H$는 $N_\delta=\{s:\operatorname{dist}(s,[\theta_0,2])<\delta\}$로 정칙하게 연장된다. 따라서 $H$는 양의 실수 반직선 $(0,\infty)$의 한 근방에서 정칙하다.

**증명.** $G$는 콤팩트 구간 $K_0=[\theta_0/2,2]$에서 양수이다. 연속성과 콤팩트성으로 $\delta\in(0,\theta_0/4)$를 골라 $K_0$의 $\delta$-근방 $N'$에서 $G\ne0$이 되게 한다. $N'$은 볼록이므로 $G$는 $N'$에서 정칙인 제곱근 $\tilde\varphi$를 갖고, 보조정리 3.1과 같은 반사 논증으로 $\tilde\varphi>0$ $(K_0$ 위$)$이다. 이제 보조정리 4.1의 증명을 $\mathcal H_\theta$ 대신 볼록 집합 $N'$에서, $b\in(\theta_0/2,\theta_0-\delta)$로 되풀이한다. (4.2)의 항등식은 연결 집합 $N'\setminus(-\infty,1]$에서 성립하고, 둘째 항의 정칙성에는 $N'$의 볼록성만 쓰인다. 따라서 $H$는 $N'\setminus(-\infty,b]\supset N_\delta$로 연장된다. 이 연장은 볼록 집합 $N'\cap\{\Re s>1\}$에서 원래의 $H$와 같다. $\theta_0\downarrow0$으로 보내고 $\Re s>1$을 더하면 마지막 주장이 따른다. $\square$

**따름정리 4.3 (낮은 높이의 균등 유계).** $(\mathrm{EM}_\theta)$ 아래

$$
\sup\Bigl\{|H(s)|:\ \theta<\Re s\le\tfrac{1+\theta}2,\ |\Im s|\le2\Bigr\}<\infty .
$$

**증명.** 보조정리 4.2를 $\theta_0=\theta$로 적용해 $\delta$를 얻는다. 두 연장은 볼록 집합 $N_\delta\cap\mathcal H_\theta$에서 일치한다. $|\Im s|\le\delta/2$인 부분은 $N_\delta$의 콤팩트 부분집합 $\{\theta\le\Re s\le\tfrac{1+\theta}2,\ |\Im s|\le\delta/2\}$에 들어가므로 $H$가 유계이다. 나머지 부분에서는 $F^2=\zeta$에서 $|F|=|\zeta|^{1/2}$이고, 이 집합은 $s=1$에서 떨어진 콤팩트 집합 안에 있으므로 $|F|$가 유계이며, $|B|\le2W/\delta$이다. $\square$

이 따름정리는 이전 판 §7에서 "$a\downarrow\theta$에서 상수가 유계"라고만 적었던 부분을 대신한다. $\{\theta<\Re s\}$는 $\mathcal H_\theta$에서 콤팩트가 아니므로, 실수축 근방에서 $H$를 $\Re s=\theta$ 너머로 연장하는 단계가 필요하다.

<a id="section-5"></a>

## 5. 균등한 평균제곱 추정

**보조정리 5.1.** $U\ge2$, $\tfrac12\le\sigma\le\tfrac32$에서 균등하게

<a id="equation-5-1"></a>

$$
\int_U^{2U}|\zeta(\sigma+it)|^2\,dt\ll U\log U,
\qquad
\int_U^{2U}|\zeta(\sigma+it)|^{1/2}\,dt\ll U(\log U)^{1/4}.
\tag{5.1}
$$

$\sigma\ge\theta_1>\tfrac12$이면 두 우변을 $O_{\theta_1}(U)$로 바꿀 수 있다. 또한 $|t|\ge2$, $\tfrac12\le\sigma\le\tfrac32$에서 $|\zeta(\sigma+it)|\ll|t|^{1/2}$이다.

**증명.** $Y=\lceil4U\rceil$로 둔다. Arias de Reyna의 근사 [2]는 $\sigma\ge0$, $0<t\le Y$에서

$$
\zeta(s)=\sum_{n\le Y}n^{-s}+\frac{Y^{1-s}}{s-1}+\vartheta\,\frac{29}{14}\,Y^{-\sigma},\qquad|\vartheta|\le1
$$

를 준다. $U\le t\le2U\le Y$에서 $|Y^{1-s}/(s-1)|\le Y^{1-\sigma}/U\ll U^{-\sigma}$이므로 $\zeta(\sigma+it)=\sum_{n\le Y}n^{-\sigma-it}+O(U^{-\sigma})$이고, 상수는 $\sigma$에 균등하다. 디리클레 다항식의 평균값 정리 [3]

$$
\int_0^T\Big|\sum_{n\le Y}a_nn^{-it}\Big|^2dt=\sum_{n\le Y}|a_n|^2\bigl(T+O(n)\bigr)
$$

를 계수 $a_nn^{-iU}$에 적용하면

$$
\int_U^{2U}\Big|\sum_{n\le Y}n^{-\sigma-it}\Big|^2dt=U\sum_{n\le Y}n^{-2\sigma}+O\Big(\sum_{n\le Y}n^{1-2\sigma}\Big)
\le U(1+\log Y)+O(Y)
$$

이다. $\sigma\ge\theta_1>\tfrac12$이면 첫 합은 $\zeta(2\theta_1)$ 이하이다. 오차 항의 제곱 적분은 $O(U^{1-2\sigma})=O(1)$이고 $|a+b|^2\le2|a|^2+2|b|^2$이므로 첫 추정이 따른다. 둘째 추정은 횔더 부등식 $\int|\zeta|^{1/2}\le U^{3/4}(\int|\zeta|^2)^{1/4}$이다. 점별 상계는 $\sum_{n\le Y}n^{-1/2}\le2\sqrt Y$에서 나오고, 음의 $t$는 켤레대칭으로 처리한다. $\square$

<a id="section-6"></a>

## 6. 절단 Perron 공식과 수평변

충분히 큰 정수 $N$에 대해 $x=N+\tfrac12$로 두고

<a id="equation-6-1"></a>

$$
c=1+\frac1{\log x},\qquad a=\theta+\frac1{\log x},\qquad T=x
\tag{6.1}
$$

를 택한다. 그러면 $\theta<a<c\le\tfrac32$, $x^c=ex$, $x^a=ex^\theta$이다.

**보조정리 6.1.** 위 매개변수에서 가정 없이

<a id="equation-6-2"></a>

$$
E(x)=\frac1{2\pi i}\int_{c-ix}^{c+ix}H(s)\frac{x^s}s\,ds+O(\log x).
\tag{6.2}
$$

**증명.** $F(s)-1=\sum_{n\ge2}f(n)n^{-s}$에 높이 $T$의 절단 Perron 공식을 적용하면 오차는

$$
\ll\sum_{n\ge2}f(n)\Big(\frac xn\Big)^c\min\Big(1,\frac1{x|\log(x/n)|}\Big)
$$

이다. $n\notin[x/2,2x]$에서는 $|\log(x/n)|\ge\log2$이고 $f\le1$이므로 이 부분은 $\ll x^{c-1}\zeta(c)\ll\log x$이다. $x/2\le n\le2x$에서는 $(x/n)^c\ll1$, $|\log(x/n)|\gg|x-n|/x$이며 반정수 $x$에서 $|x-n|\ge\tfrac12$이므로 이 부분도 $\ll\sum_n\min(1,|x-n|^{-1})\ll\log x$이다. 따라서 $A(x)-1$은 $(F-1)x^s/s$의 절단 적분과 $O(\log x)$만큼 다르다.

$B(s)/s$의 수직선 적분은 $|B(s)/s|\le W/(|\Im s||s|)$로 절대수렴한다. 적분 순서를 바꾸고, 각 $u$에서 $x^s/(s(s-u))$의 유수를 $s=0,u$에서 모으면

$$
\frac1{2\pi i}\int_{c-i\infty}^{c+i\infty}B(s)\frac{x^s}{s}\,ds=\int_0^1w(u)\frac{x^u-1}{u}\,du=\mathcal P(x)-1
$$

이다. $|t|\ge x$의 꼬리는 (2.2)에 의해 $\ll Wx^c\int_x^\infty t^{-2}dt=O(1)$이다. 두 식을 빼면 (6.2)이다. $\square$

$(\mathrm{EM}_\theta)$ 아래 보조정리 4.1에 의해 $H(s)x^s/s$는 직사각형 $a\le\Re s\le c$, $|\Im s|\le x$의 근방에서 정칙이다. $a>0$이므로 $s=0$은 밖에 있다. 따라서 수직선을 $c$에서 $a$로 옮길 수 있다. 수평변에서는 보조정리 5.1과 (2.2)로

$$
|H(\sigma\pm ix)|\le|\zeta(\sigma\pm ix)|^{1/2}+1+\frac Wx\ll x^{1/4}
$$

이고, 길이는 $c-a\le1$, $|s|\ge x$, $|x^s|\le ex$이므로 두 수평변의 기여는 $O(x^{1/4})$이다. 그러므로

<a id="equation-6-3"></a>

$$
E(x)=\frac1{2\pi i}\int_{a-ix}^{a+ix}H(s)\frac{x^s}{s}\,ds+O(x^{1/4}+\log x).
\tag{6.3}
$$

<a id="section-7"></a>

## 7. 정리 1의 증명

$|t|\le2$의 부분은 따름정리 4.3에 의해 $O_\theta(x^a)=O_\theta(x^\theta)$이다. $|t|\ge2$에서는 (2.2)로

<a id="equation-7-1"></a>

$$
|H(a+it)|\le|\zeta(a+it)|^{1/2}+1+\frac{W}{|t|}.
\tag{7.1}
$$

$[2,x]$를 $[U,2U]$ 꼴의 $O(\log x)$개 구간으로 덮는다. 마지막 불완전 구간은 완전한 구간으로 늘려 상계한다. 각 구간에서 보조정리 5.1을 $\sigma=a$에 적용하면

<a id="equation-7-2"></a>

$$
\int_U^{2U}\frac{|H(a+it)|}{|a+it|}\,dt
\le\frac1U\int_U^{2U}|H(a+it)|\,dt\ll_\theta
\begin{cases}1,&\theta>\tfrac12,\\(\log U)^{1/4},&\theta=\tfrac12\end{cases}
\tag{7.2}
$$

이다. $\theta>\tfrac12$이면 $a\ge\theta$이므로 개선된 평균값을 쓴다. 음의 $t$도 같다. 합하면

<a id="equation-7-3"></a>

$$
\Big|\frac1{2\pi i}\int_{a-ix}^{a+ix}H(s)\frac{x^s}s\,ds\Big|
\ll_\theta x^a
\begin{cases}\log x,&\theta>\tfrac12,\\(\log x)^{5/4},&\theta=\tfrac12\end{cases}
\tag{7.3}
$$

이다. 둘째 경우의 지수는 $\sum_{j\le O(\log x)}j^{1/4}\ll(\log x)^{5/4}$에서 나온다. $x^a=ex^\theta$이고 $\theta\ge\tfrac12$이므로 (6.3)의 $O(x^{1/4}+\log x)$는 흡수된다. 일반 실수 $x$는 보조정리 2.1로 $\lfloor x\rfloor+\tfrac12$에 옮긴다. $\square$

<a id="section-8"></a>

## 8. 역방향: 따름정리 3

**보조정리 8.1.** $0\le\alpha<1$이고 $H$가 $\{\Re s>\alpha\}$로 정칙하게 연장되면, $\Re\rho>\alpha$인 $\zeta$의 모든 영점은 짝수 중복도이다.

**증명.** $\Omega_\alpha=\{\Re s>\alpha\}\setminus[0,1]$은 연결이고, $B$가 $\mathbb C\setminus[0,1]$에서 정칙이므로 $\tilde F=1+B+H$는 $\Omega_\alpha$에서 정칙이다. $\Re s>1$에서 $\tilde F=F$이므로 항등정리에 의해 $\Omega_\alpha$에서 $\tilde F^2=\zeta$이다. $\Re\rho>\alpha$인 영점은 실수가 아니므로 $\rho\in\Omega_\alpha$이고, $\operatorname{ord}_\rho\zeta=2\operatorname{ord}_\rho\tilde F$는 짝수이다. $\square$

**따름정리 3의 증명.** (a)$\Rightarrow$(b)는 정리 1이다. (b)를 가정하면 각 $\eta>0$에 대해 $s\int_1^\infty E(x)x^{-s-1}dx$가 $\Re s>\theta+\eta$의 콤팩트 집합에서 균등하게 절대수렴하여 정칙함수를 정의하고, (2.3)에 의해 이는 $H$의 연장이다. $\eta$가 임의이므로 $H$는 $\Re s>\theta$로 연장되고, 보조정리 8.1이 (a)를 준다. $\square$

이전 판은 정수 표본 $\sum(f(n)-p(n))$의 디리클레 급수 $\sqrt{\zeta(s)}-J(s)$, $J(s)=1+\int_0^1w(u)[\zeta(s+1-u)-1]\,du$를 썼다. $J(s)-1-B(s)=\int_0^1w(u)\bigl[\zeta(s+1-u)-\tfrac1{s-u}-1\bigr]du$가 전체함수이므로(괄호 안이 $s$의 전체함수) 두 변환의 특이점은 같다. 연속 주항 $\mathcal P$를 쓰면 이 단계가 필요 없다.

<a id="section-9"></a>

## 9. 정확한 지수와 진동: 정리 2

**Hardy의 정리.** $Z(t)=e^{i\vartheta(t)}\zeta(\tfrac12+it)$는 실함수이고 부호를 무한히 자주 바꾼다 [4]. 부호가 바뀌는 점은 홀수 차수 영점이므로 임계선 위에 홀수 중복도 영점이 무한히 많고, $\Theta_{\rm odd}\ge\tfrac12$이다.

**Landau의 정리 [5].** $g$가 $[1,\infty)$의 유한 구간마다 유계 가측이고 충분히 큰 $x$에서 $g(x)\ge0$이라 하자. $\int_1^\infty g(x)x^{-s-1}dx$의 수렴 가로좌표 $\sigma_c$가 유한하면, 이 적분이 정의하는 함수는 실수점 $s=\sigma_c$에서 정칙이 아니다.

**정리 2의 증명.** (2)부터 보인다. 어떤 $0<\alpha<\beta$와 $x_0$에 대해 $x\ge x_0$에서 $E(x)\le x^\alpha$라 하자. $g(x)=x^\alpha-E(x)$로 두면 $\Re s>1$에서

$$
\int_1^\infty g(x)x^{-s-1}dx=\frac1{s-\alpha}-\frac{H(s)}s
$$

이다. 보조정리 4.2에 의해 우변은 $(\alpha,\infty)$의 각 실수점 근방으로 정칙하게 연장되고, 그 연장은 볼록 집합 $N_\delta\cap\{\Re s>\sigma_c\}$에서 적분과 일치한다. $E(x)=O(x)$이므로 $\sigma_c\le1$이고, Landau의 정리에 의해 $\sigma_c\le\alpha$이다. 따라서 적분은 $\Re s>\alpha$에서 정칙이고, $H=s/(s-\alpha)-s\int g\,x^{-s-1}$이 $\Re s>\alpha$로 연장된다. 보조정리 8.1에 의해 $\Re\rho>\alpha$인 영점은 모두 짝수인데, 이는 $\rho$와 모순이다. $g=x^\alpha+E$에 대해서도 같다. Hardy의 정리가 주는 $\beta=\tfrac12$의 영점에 적용하면 $\Omega_\pm(x^{1/2-\varepsilon})$를 얻는다.

(1)에서, $\Theta_{\rm odd}<1$이면 상한의 정의에 의해 $(\mathrm{EM}_{\Theta_{\rm odd}})$가 성립하므로 정리 1이 $E(x)\ll x^{\Theta_{\rm odd}}(\log x)^{5/4}$를 준다. $\Theta_{\rm odd}=1$이면 자명한 상계 $E(x)\ll x$를 쓴다. 반대 부등식은 (2)이다. $\square$

따라서 $E$의 크기를 정하는 것은 영점의 위치가 아니라 **홀수 중복도 영점**의 위치이다. 짝수 중복도 영점은 $F$를 정칙으로 남기므로 어떤 크기에서도 보이지 않는다. 이는 $M(x)$의 지수가 모든 영점의 실수부 상한 $\Theta$로 정해지는 것과 대조된다.

<a id="section-10"></a>

## 10. 홀수 영점의 국소 기여

**[국소 계산]** $\rho$가 중복도 $m$인 영점이고 $b_\rho=\zeta^{(m)}(\rho)/m!$라 하자. $\rho$에서 왼쪽 수평으로 절단을 두고, [영점 기여 (1.1)](zero-contributions.md#equation-1-1)의 계산을 $H\equiv1$인 $d_z$ 모형([주정리 명제 8.2](main-theorem.md#section-8-2))과 $(z,w)=(\tfrac12,0)$, $\xi=0$에 적용한다. [주정리 (4.2)](main-theorem.md#equation-4-2)에서 $\langle\mathcal K_{-m/2,0},1\rangle=1/\Gamma(-m/2)$이므로, $E(x)$에 대한 $\rho$의 기여는

<a id="equation-10-1"></a>

$$
\frac{b_\rho^{1/2}}{\Gamma(-m/2)}\,\frac{x^\rho}{\rho}\,(\log x)^{-1-m/2}
+O_\rho\!\left(x^{\Re\rho}(\log x)^{-2-m/2}\right)
\tag{10.1}
$$

이다. 가지 $b_\rho^{1/2}$는 $\Re s>1$에서의 연속으로 정해진다. $m$이 짝수이면 $1/\Gamma(-m/2)=0$이어서 보조정리 4.1과 맞고, $m=1$이면 $\Gamma(-\tfrac12)=-2\sqrt\pi$이다. 이 식은 고정된 $\rho$의 기여이며 영점합의 수렴은 다루지 않는다.

같은 계산을 대칭 지수 $(\tfrac12,\tfrac12)$에 적용하면 국소 핵은 $\mathcal K_{-m/2,-m/2}$이고, $m$이 짝수이면 두 인수의 지지가 $\{0\},\{1\}$로 분리되어 핵이 $0$이다([주정리 §4](main-theorem.md#section-4)). 이 경우의 이차 모멘트는 [이차 모멘트 문서](m2-divisor-moment.md)에서 다룬다.

<a id="references"></a>

## 참고 문헌

1. NIST DLMF, [§25.2](https://dlmf.nist.gov/25.2), [§25.10](https://dlmf.nist.gov/25.10): 제타 함수의 해석적 연속, 극, 교대급수 표현과 영점의 기본 성질.
2. J. Arias de Reyna, [On the approximation of the zeta function by Dirichlet polynomials](https://arxiv.org/abs/2406.16667), arXiv:2406.16667 (2024): $\sigma\ge0$, $0<t\le x$에서 오차 $\tfrac{29}{14}x^{-\sigma}$의 근사. 보조정리 5.1에 사용.
3. H. L. Montgomery, R. C. Vaughan, *Hilbert's inequality*, J. London Math. Soc. (2) 8 (1974), 73–82: 디리클레 다항식의 평균값 정리. 보조정리 5.1에 사용.
4. G. H. Hardy, *Sur les zéros de la fonction $\zeta(s)$ de Riemann*, C. R. Acad. Sci. Paris 158 (1914), 1012–1014; E. C. Titchmarsh, *The Theory of the Riemann Zeta-Function*, 2nd ed., Ch. X: 임계선 위 홀수 차수 영점의 무한성.
5. H. L. Montgomery, R. C. Vaughan, *Multiplicative Number Theory I*, Cambridge (2007), §15.1: 비음 함수의 멜린 변환에 대한 Landau의 정리.
