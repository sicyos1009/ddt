# 일반화된 베타 핵과 가중 약수합

**2026-11-01 (개정판)**

고정된 복소 지수와 충분히 매끄러운 시험함수에 대해, 제곱 인수 없는 정수의 가중 약수합을 이차 항까지 전개한다. 유효 Perron 공식과 제타함수의 고전적 영점 없는 영역을 사용한다.

**[개요]**는 증명의 뼈대만 제시한 부분이며, **RH**는 리만 가설을 뜻한다.

## 1. 정리

$\omega(n)$은 서로 다른 소인수의 개수, $\mu$는 뫼비우스 함수이다. $z,w\in\mathbb C$, $F:[0,1]\to\mathbb C$에 대해

$$
A_x^{z,w}(F)=
\sum_{\substack{2\le n\le x\\\mu^2(n)=1}}
\sum_{d\mid n}z^{\omega(d)}w^{\omega(n/d)}
F\!\left(\frac{\log d}{\log n}\right),
\qquad 0^0:=1
$$

로 둔다. 이하 $r=z+w$, $\beta=\Re r$, $L=\log x$이며, $L^\lambda=\exp(\lambda\log L)$이다. $\gamma$는 오일러 상수이다.

$$
C(r)=\prod_p\left(1+\frac rp\right)\left(1-\frac1p\right)^r,
\qquad
\mathcal K_{z,w}(u)=g_z(u)g_w(1-u),
\tag{1.1}
$$

여기서 $g_\alpha(t)=t_+^{\alpha-1}/\Gamma(\alpha)$는 §3에서 정의하는 전체해석적 분포 가족이다. §2의 $H=H_{z,w}$에 대해 $H_s=\partial_sH(1,1)$, $H_t=\partial_tH(1,1)$로 쓴다. $H(1,1)=C(r)$이다.

**정리 1 (이차 베타 법칙).** 고정된 $z,w$에 대해 충분히 큰 유한 정수 $N=N(z,w)$가 존재하여, $F\in C^N([0,1])$이면 $x\to\infty$일 때

$$
\boxed{
A_x^{z,w}(F)
=xL^{r-1}\Big\langle C(r)\,\mathcal K_{z,w}+\frac1L\,\mathcal K^{(1)}_{z,w},\ F\Big\rangle
+O_{F,z,w}(xL^{\beta-3}),
}
\tag{1.2}
$$

$$
\mathcal K^{(1)}_{z,w}
=\bigl(\gamma zC+H_s\bigr)\mathcal K_{z-1,w}
+\bigl(\gamma wC+H_t\bigr)\mathcal K_{z,w-1}
-(r-1)\,C\,\mathcal K_{z,w},
\qquad C=C(r).
\tag{1.3}
$$

(2.4)의 $C\eta$를 전체함수로 연장하고 $C\kappa:=\gamma C+C\eta$로 읽으면, 모든 $r\in\mathbb C$에서

$$
\mathcal K^{(1)}_{z,w}
=(C\kappa)(r)\bigl(z\,\mathcal K_{z-1,w}+w\,\mathcal K_{z,w-1}\bigr)-(r-1)C(r)\mathcal K_{z,w},
\qquad
\kappa(r)=\gamma+\sum_p\frac{(r+1)\log p}{(p-1)(p+r)}.
\tag{1.4}
$$

증명은 §5와 부록 A에 있다.

특히 $A_x^{z,w}(F)=C(r)xL^{r-1}\langle\mathcal K_{z,w},F\rangle+O_{F,z,w}(xL^{\beta-2})$이므로, 정규화한 합은 분포의 의미에서 $C(r)\mathcal K_{z,w}$로 수렴한다.

**예.**

- **아크사인 $(\tfrac12,\tfrac12)$.** $a(u)=1/(\pi\sqrt{u(1-u)})$라 쓰면, §6의 핵 항등식에서
$$
A_x^{1/2,1/2}(F)=\frac{x}{\zeta(2)}\Big[\int_0^1a(u)F(u)\,du+\frac{\kappa(1)}L\int_0^1\big(u-\tfrac12\big)a(u)F'(u)\,du\Big]+O\!\left(\frac x{L^2}\right),
\tag{1.5}
$$
$$
\kappa(1)=\gamma-2\frac{\zeta'(2)}{\zeta(2)}=1.7171378\ldots
$$
이차 보정은 질량을 양 끝점 쪽으로 옮긴다(B.1).
- **총질량 $F=1$.** $\langle\mathcal K_{a,b},1\rangle=1/\Gamma(a+b)$에서
$$
\langle\mathcal K^{(1)}_{z,w},1\rangle=\frac{\Phi'(1)}{\Gamma(r-1)},\qquad
\Phi(s)=\frac{((s-1)\zeta(s))^rH(s,s)}{s}.
$$
이는 $\sum_{n\le x}\mu^2(n)r^{\omega(n)}$에 대한 Selberg–Delange 전개의 둘째 계수와 같다.

양의 실수 $z,w$에서는 지시함수에도 1차 법칙이 성립한다(§5.3).

## 2. 오일러 곱

$\mu^2(dm)=1$이면 각 소수는 $d,m$ 중 한쪽에만 속한다. 따라서 $\Re s,\Re t>1$에서

$$
\begin{aligned}
D_{z,w}(s,t)
&=\sum_{\substack{d,m\ge1\\\mu^2(dm)=1}}
\frac{z^{\omega(d)}w^{\omega(m)}}{d^s m^t}\\
&=\prod_p(1+zp^{-s}+wp^{-t})\\
&=\zeta(s)^z\zeta(t)^wH_{z,w}(s,t),
\end{aligned}
\tag{2.1}
$$

$$
H_{z,w}(s,t)=
\prod_p(1+zp^{-s}+wp^{-t})
(1-p^{-s})^z(1-p^{-t})^w.
\tag{2.2}
$$

거듭제곱의 가지는 오일러 곱의 로그로 정한다. 국소 인수의 일차항이 상쇄되어

$$
H_p(s,t)-1
=O_{z,w}\!\left(p^{-2\Re s}+p^{-\Re(s+t)}+p^{-2\Re t}\right).
$$

그러므로 고정된 $\sigma_0>1/2$에서 $H$와 필요한 유한 차수 편미분들은 $\Re s,\Re t\ge\sigma_0$에서 균등 유계이다. 미분으로 생기는 $(\log p)^j$도 합산 가능하다. 이 논증은 $H$의 영점도 허용한다.

$$
H_{z,w}(1,1)=C(r),\qquad
C(0)=C(-1)=1,\quad
C(1)=\frac1{\zeta(2)},\quad C(-2)=0.
\tag{2.3}
$$

특히 $z=w=1/2$일 때 $x$로 정규화한 총질량은 $1/\zeta(2)$이다. 제곱 인수 없는 정수의 개수로 정규화하면 질량 $1$의 아크사인 법칙이 된다.

국소 인수를 $s$로 미분하면, 우선 $C(r)\ne0$에서

$$
H_s=C(r)\,z\,\eta(r),\qquad H_t=C(r)\,w\,\eta(r),\qquad
\eta(r)=\sum_p\frac{(r+1)\log p}{(p-1)(p+r)},
\tag{2.4}
$$

이다. $C$는 $r=-p$에서 단순영점, $\eta$는 같은 점에서 단순극을 가진다. 그 곱의 연장값은
$$
(C\eta)(-p)=-\frac{\log p}{p}\left(1-\frac1p\right)^{-p}
\prod_{q\ne p}\left(1-\frac pq\right)\left(1-\frac1q\right)^{-p}.
$$
국소 곱의 미분도 $H_s=zC\eta$, $H_t=wC\eta$의 같은 값을 주므로 (2.4)는 모든 $r$로 연장된다. $C\eta$와 $C\kappa=\gamma C+C\eta$는 전체함수이다.

### 2.1. 서로소 조건을 푼 변형

제곱 인수 없는 $d,m$의 모든 쌍을 서로소 조건 없이 세면

$$
D^\natural_{z,w}(s,t)=\sum_{d,m\ge1}\frac{\mu^2(d)\mu^2(m)z^{\omega(d)}w^{\omega(m)}}{d^sm^t}
=\prod_p(1+zp^{-s})(1+wp^{-t})
=\zeta(s)^z\zeta(t)^wH^\natural(s,t),
$$

$$
H^\natural(s,t)=H_z(s)H_w(t),\qquad H_\alpha(s)=\prod_p(1+\alpha p^{-s})(1-p^{-s})^\alpha
\tag{2.5}
$$

이다. 국소 인수의 차이는 $H^\natural_p/H_p-1=zwp^{-s-t}/(1+zp^{-s}+wp^{-t})=O(p^{-\Re(s+t)})$로, 위의 $H_p-1$ 상계 안에 들어간다. 계수도 같은 $\tau_{k_0}$로 지배된다. 따라서

$$
A_x^{\natural}(F)=\sum_{2\le dm\le x}\mu^2(d)\mu^2(m)z^{\omega(d)}w^{\omega(m)}F\!\left(\frac{\log d}{\log dm}\right)
$$

에 대해 정리 1이 그대로 성립하고, $C(r),H_s,H_t$가 각각 $C(z)C(w)$, $C(z)C(w)z\eta(z)$, $C(z)C(w)w\eta(w)$로 바뀐다. 위와 같이 곱을 연장하면 모든 $z,w$에서

$$
\mathcal K^{(1)\natural}_{z,w}=C(z)C(w)\Big[z\kappa(z)\mathcal K_{z-1,w}+w\kappa(w)\mathcal K_{z,w-1}-(r-1)\mathcal K_{z,w}\Big].
\tag{2.6}
$$

### 2.2. $H$의 2차 구조

$e_1=z(z+1)/2$, $e_2=w(w+1)/2$로 두면

$$
H(s,t)=\zeta(2s)^{-e_1}\,\zeta(s+t)^{-zw}\,\zeta(2t)^{-e_2}\,\widetilde H(s,t),
\tag{2.7}
$$

$$
\widetilde H_p=H_p\,(1-p^{-2s})^{-e_1}(1-p^{-s-t})^{-zw}(1-p^{-2t})^{-e_2}=1+O\!\left(p^{-3\min(\Re s,\Re t)}\right).
$$

$X=p^{-s}$, $Y=p^{-t}$로 두고 로그를 전개하면 이차항은 $-\tfrac12(zX+wY)^2-\tfrac z2X^2-\tfrac w2Y^2+e_1X^2+zwXY+e_2Y^2=0$으로 정확히 상쇄된다. 따라서 $\widetilde H$는 $\Re s,\Re t>\tfrac13$에서 해석적이다. 서로소 변형에서는 교차 인수 $\zeta(s+t)^{-zw}$가 없다. 즉 서로소 조건은 정확히 이 인수에 해당한다.

실축 근처($|\Im s|,|\Im t|<7$)의 $\Re>\tfrac13$에서 $H$의 특이점은 $s=\tfrac12$, $s+t=1$, $t=\tfrac12$뿐이다. $\zeta(2s)$, $\zeta(s+t)$의 영점은 높이 $14$ 이상에서만 나오기 때문이다. 대각선에서는 세 인수가 합쳐져

$$
H(s,s)=\zeta(2s)^{-r(r+1)/2}\,\widetilde H(s,s)
\tag{2.8}
$$

가 된다. $r=1$이면 $H(s,s)=\prod_p(1-p^{-2s})=1/\zeta(2s)$이고, $r=0,-1$이면 지수가 $0$이다.

## 3. 분포 핵과 모멘트

$\Re\alpha>0$에서 $g_\alpha=t_+^{\alpha-1}/\Gamma(\alpha)$는 국소 적분 가능 함수이다. 비음의 정수 $N$에 대해

$$
\langle g_\alpha,\varphi\rangle
=\frac{(-1)^N}{\Gamma(\alpha+N)}
\int_0^\infty t^{\alpha+N-1}\varphi^{(N)}(t)\,dt,
\qquad \Re\alpha>-N
\tag{3.1}
$$

로 연장한다. $\varphi\in C_c^\infty(\mathbb R)$이며, 겹치는 반평면에서 부분적분으로 일치하므로 전체해석적 분포 가족을 얻는다. 특히

$$
g_0=\delta_0,\qquad g_{-1}=\delta_0',
\qquad g_\alpha'=g_{\alpha-1}.
$$

$g_z(u)$와 $g_w(1-u)$의 특이 지지집합은 각각 $\{0\}$, $\{1\}$이다. 각 끝점에서 다른 인수가 매끄러우므로 곱 $\mathcal K_{z,w}$가 정의되며, 이는 $[0,1]$에 지지된 전체해석적 분포 가족이다. $F\in C^\infty([0,1])$의 작용은 매끄러운 콤팩트 지지 연장으로 정의한다. 연장들은 $[0,1]$에서 모든 제트가 같으므로 결과는 선택에 무관하다.

베타적분을 수렴 영역에서 계산한 뒤 연장하면

$$
\boxed{
\langle\mathcal K_{z,w},u^k\rangle
=\frac{(z)_k}{\Gamma(r+k)}
\qquad(k\ge0).
}
\tag{3.2}
$$

$(z)_k$는 상승계승이며 $(z)_0=1$, $1/\Gamma$는 전체함수로 읽는다. 특히 총질량은 $1/\Gamma(r)$이다. 질량이 영이어도 핵은 비영일 수 있다. 한 지수가 $0$이면

$$
\mathcal K_{0,w}=\frac{\delta_0}{\Gamma(w)},
\qquad
\mathcal K_{z,0}=\frac{\delta_1}{\Gamma(z)}.
$$

두 지수가 모두 비양의 정수이면 두 인수의 지지가 분리되어 핵 자체가 영이다.

**핵 항등식.** $tg_{\alpha-1}(t)=(\alpha-1)g_\alpha(t)$와 $g_\alpha'=g_{\alpha-1}$에서

$$
u\,\mathcal K_{z-1,w}=(z-1)\mathcal K_{z,w},\qquad
(1-u)\,\mathcal K_{z,w-1}=(w-1)\mathcal K_{z,w},\qquad
\partial_u\mathcal K_{z,w}=\mathcal K_{z-1,w}-\mathcal K_{z,w-1}.
\tag{3.3}
$$

이로부터 $(u\mathcal K_{z,w})'=(r-1)\mathcal K_{z,w}-\mathcal K_{z,w-1}$이고, 따라서

$$
\langle\mathcal K_{z,w},uF'\rangle=\langle\mathcal K_{z,w-1},F\rangle-(r-1)\langle\mathcal K_{z,w},F\rangle.
\tag{3.4}
$$

또한 $u^a(1-u)^b\mathcal K_{z,w}=(z)_a(w)_b\,\mathcal K_{z+a,w+b}$이다.

**로그 가중치.** 정수 $a,b\ge0$에 대해 위 항등식과 정리 1의 부분합산으로

$$
\begin{aligned}
&\sum_{\substack{2\le n\le x\\\mu^2(n)=1}}\sum_{d\mid n}
z^{\omega(d)}w^{\omega(n/d)}(\log d)^a
\Big(\log\frac nd\Big)^bF\!\left(\frac{\log d}{\log n}\right)\\
&\qquad=C(r)(z)_a(w)_b\,xL^{r-1+a+b}
\langle\mathcal K_{z+a,w+b},F\rangle
+O_{F,z,w,a,b}(xL^{\beta-2+a+b}).
\end{aligned}
\tag{3.5}
$$

실제로 가중치는 $(\log n)^{a+b}u^a(1-u)^bF(u)$이다. 아크사인에 $\log d\log(n/d)$를 곱하면 핵이 $\tfrac14\mathcal K_{3/2,3/2}=\pi^{-1}\sqrt{u(1-u)}$로 바뀌어, 정규화 후 반원 법칙을 얻는다.

아크사인 핵의 푸리에 변환은 Bessel 함수 $J_0$로 표현된다.

$$
\langle\mathcal K_{1/2,1/2},e^{i\xi u}\rangle
=e^{i\xi/2}J_0\!\left(\frac\xi2\right).
\tag{3.6}
$$

## 4. 이동 인수와 역라플라스 변환

실수 $\xi$에 대해

$$
\begin{aligned}
b_{\xi,L}(n)
&=\mu^2(n)\sum_{d\mid n}
z^{\omega(d)}w^{\omega(n/d)}e^{i\xi\log d/L},\\
S_x(\xi)&=\sum_{n\le x}b_{\xi,L}(n)
\end{aligned}
\tag{4.1}
$$

로 둔다. $n=1$ 항은 $1$이다. 대응하는 디리클레 급수는

$$
\mathfrak D_{\xi,L}(s)
=D_{z,w}\!\left(s-\frac{i\xi}{L},s\right).
\tag{4.2}
$$

(2.1)의 두 제타 인수는 $s=1$, $s=1+i\xi/L$에 특이점 후보를 갖는다. 후보는 합류하거나 상쇄될 수 있다. 예컨대 $r=0,\xi=0$이면 정확히 $D_{z,w}(s,s)=1$, $S_x(0)=1$이다.

$s=1+q/L$로 확대하면 주인수는 $L^rP_\xi(q)$이며

$$
P_\xi(q)=(q-i\xi)^{-z}q^{-w}.
$$

두 거듭제곱은 각각 왼쪽 수평 반직선을 절단선으로 하는 주가지를 쓴다. $Y>|\xi|+1$에 대해

$$
\mathcal U_Y:
(-\infty-iY)\to(1-iY)\to(1+iY)\to(-\infty+iY)
\tag{4.3}
$$

로 경로를 정한다. 이 경로는 두 후보에서 거리가 적어도 $1$이고 $\Re q\le1$이다.

$$
\mathcal B_{z,w}(\xi)
:=\frac1{2\pi i}\int_{\mathcal U_Y}e^qP_\xi(q)\,dq
=\langle\mathcal K_{z,w},e^{i\xi u}\rangle.
\tag{4.4}
$$

**증명.** 먼저 $\Re z,\Re w>1$에서, $\Re q>0$인 수직선 위의 라플라스 변환

$$
\mathcal L(e^{i\xi t}g_z(t))(q)=(q-i\xi)^{-z},
\qquad \mathcal L(g_w)(q)=q^{-w}
$$

에 합성곱 정리를 적용한다. 역변환을 $t=1$에서 평가하면

$$
\frac1{\Gamma(z)\Gamma(w)}
\int_0^1e^{i\xi u}u^{z-1}(1-u)^{w-1}\,du
$$

를 얻는다. 수직선을 $\mathcal U_Y$로 변형할 때, 왼쪽 보조 수직선 $\Re q=-R$, $|\Im q|\ge Y$의 적분은 $O_{z,w,\xi}(e^{-R}R^{1-\beta})$로 사라진다. 높이를 무한대로 보내는 수평 연결 적분도 사라진다.

$\mathcal U_Y$의 적분은 모든 $z,w$에서 절대수렴하고, 매개변수의 콤팩트 집합에서 균등 수렴한다. 경로 높이의 독립성은 Cauchy 정리로, (4.4)의 전 매개변수 확장은 전체해석적 연장으로 따른다.

(4.4)는 모든 $z,w$에 대해 성립하므로, 특히 $(q-i\xi)P_\xi(q)$, $qP_\xi(q)$의 적분은 각각 $\mathcal B_{z-1,w}(\xi)$, $\mathcal B_{z,w-1}(\xi)$이다. 이차 법칙의 보정 핵은 이렇게 나온다.

콤팩트 지지 분포의 유한 차수 성질에 의해, 각 고정된 $j$와 $(a,b)\in\{(0,0),(1,0),(0,1)\}$에 대해

$$
|\mathcal B_{z-a,w-b}^{(j)}(\xi)|
\ll_{z,w,j}(1+|\xi|)^{M_j}
\tag{4.5}
$$

인 $M_j\ge0$가 존재한다. 경로 표현을 미분할 때는 $\xi$의 작은 근방을 감싸는 하나의 $Y$를 고정한다.

## 5. 산술 점근식의 증명

$H_s,H_t$는 §1과 같고, $C=C(r)$이다. 이 절에서

$$
G_1=\gamma zC+H_s,\qquad G_2=\gamma wC+H_t-C
\tag{5.0}
$$

로 둔다. $G_1,G_2$는 부록 A.3의 국소 인수 $G$를 $(s_1,s_2)=(s-i\alpha,s)$에 대해 $(1,1)$에서 편미분한 값이다. $G_2$의 $-C$는 Perron 인자 $1/s$에서 온다.

### 5.1. 이동 균등 추정

**보조정리 1.** 고정된 $z,w$에 대해 어떤 $B\ge0$가 존재하여 $j=0,1$, $|\xi|\le L^{1/4}$에서 균등하게

$$
\partial_\xi^jS_x(\xi)
=xL^{r-1}\,\partial_\xi^j\Big[C\mathcal B_{z,w}(\xi)+\frac1L\big(G_1\mathcal B_{z-1,w}(\xi)+G_2\mathcal B_{z,w-1}(\xi)\big)\Big]
+O_{z,w}\!\left(xL^{\beta-3}(1+|\xi|)^B\right).
\tag{5.1}
$$

증명은 부록 A에 있다.

### 5.2. 정규화 변수의 변경

$\Theta_x(\xi):=A_x^{z,w}(e^{i\xi u})$로 둔다. 각 항에 대해

$$
\frac{d}{dt}e^{i\xi\log d/\log t}
=-\frac{\xi}{t\log t}\,
\partial_\xi e^{i\xi\log d/\log t}
$$

를 $dm\le t\le x$에서 적분하고 합하면

$$
\Theta_x(\xi)
=S_x(\xi)-1+
\xi\int_2^x\frac{\partial_\xi S_t(\xi)}{t\log t}\,dt.
\tag{5.2}
$$

여기서 $S_t$는 정규화 $\log t$를 쓴다. $k_0\ge\max(2,|z|+|w|)$인 정수를 고정하고, $\tau_{k_0}$를 $\zeta(s)^{k_0}$의 계수라 하자. 절댓값 합은

$$
V(t):=\sum_{n\le t}\mu^2(n)(|z|+|w|)^{\omega(n)}
\le\sum_{n\le t}\tau_{k_0}(n)
\ll_{z,w}t(\log(2t))^{k_0-1}.
$$

$|\partial_\xi S_t|\le V(t)$이므로 작은 구간의 기여는

$$
\left|\xi\int_2^{\sqrt x}
\frac{\partial_\xi S_t(\xi)}{t\log t}\,dt\right|
\ll_{z,w}|\xi|\sqrt x\,L^{k_0}.
$$

$|\xi|\le L^{1/8}$, $t\ge\sqrt x$이면 충분히 큰 $x$에서 $|\xi|\le(\log t)^{1/4}$이다. 이 구간에는 (5.1)을 $j=1$로 적용한다. 주항 $Ct(\log t)^{r-1}\mathcal B'(\xi)$의 기여는

$$
\xi\,C\mathcal B_{z,w}'(\xi)\int_{\sqrt x}^x(\log t)^{r-2}\,dt
=\xi\,C\mathcal B_{z,w}'(\xi)\,xL^{r-2}\bigl(1+O(L^{-1})\bigr)
$$

이고, 나머지는 (4.5)에 의해 $O(xL^{\beta-3}(1+|\xi|)^{M})$이다. $S_x(\xi)$에 (5.1)을 $j=0$으로 대입하면, 어떤 $M_*\ge0$에 대해 $|\xi|\le L^{1/8}$에서

$$
\Theta_x(\xi)
=xL^{r-1}\Big[C\mathcal B_{z,w}+\frac1L\big(G_1\mathcal B_{z-1,w}+G_2\mathcal B_{z,w-1}+C\,\xi\mathcal B_{z,w}'\big)\Big](\xi)
+O_{z,w}\!\left(xL^{\beta-3}(1+|\xi|)^{M_*}\right).
\tag{5.3}
$$

### 5.3. 푸리에 적분

$F$의 연장 $\widetilde F\in C_c^N(\mathbb R)$를 $\|\widetilde F\|_{C^N}\ll_N\|F\|_{C^N}$이 되도록 고정한다. 푸리에 규약은

$$
\widetilde F(u)=\frac1{2\pi}
\int_{\mathbb R}\widehat{\widetilde F}(\xi)e^{i\xi u}\,d\xi
$$

이며 $|\widehat{\widetilde F}(\xi)|\ll_N\|F\|_{C^N}(1+|\xi|)^{-N}$이다. $N>M_*+1$로 잡고 $R=L^{1/8}$로 나누면 저주파의 오차 적분은 (5.3)에 의해 $O_{F,z,w}(xL^{\beta-3})$이다. 고주파에서는

$$
\begin{aligned}
\int_{|\xi|>R}|\widehat{\widetilde F}(\xi)\Theta_x(\xi)|\,d\xi
&\ll_{F,z,w,N}xL^{k_0-1}R^{1-N},\\
xL^{\beta-1}\int_{|\xi|>R}
|\widehat{\widetilde F}(\xi)|\,(1+|\xi|)\max_{a,b}|\mathcal B_{z-a,w-b}(\xi)|\,d\xi
&\ll_{F,z,w,N}xL^{\beta-1}R^{M_0+2-N}.
\end{aligned}
$$

여기서 $M_0$는 필요한 핵과 그 1차 주파수 미분을 모두 지배하도록 키운다. $N\ge M_0+18$, $N\ge1+8(k_0+2-\beta)$도 만족시키면 두 꼬리도 $O(xL^{\beta-3})$ 이하이다. 마지막으로 $\xi\mathcal B'_{z,w}(\xi)=-\langle(u\mathcal K_{z,w})',e^{i\xi u}\rangle$이므로

$$
\frac1{2\pi}\int_{\mathbb R}\widehat{\widetilde F}(\xi)\Big[\cdots\Big]d\xi
=C\langle\mathcal K_{z,w},F\rangle+\frac1L\Big(G_1\langle\mathcal K_{z-1,w},F\rangle+G_2\langle\mathcal K_{z,w-1},F\rangle+C\langle\mathcal K_{z,w},uF'\rangle\Big).
$$

(3.4)를 쓰면 $G_2\langle\mathcal K_{z,w-1},F\rangle+C\langle\mathcal K_{z,w},uF'\rangle=(G_2+C)\langle\mathcal K_{z,w-1},F\rangle-(r-1)C\langle\mathcal K_{z,w},F\rangle$이고, $G_2+C=\gamma wC+H_t$이므로 (1.2)–(1.3)이 따른다. 분포와 적분의 교환은 (4.5)와 위의 감소 추정으로 정당화된다. 오차 상수는 $\|F\|_{C^N}$에 선형으로 의존한다.

**양의 지수와 지시함수.** $z,w>0$, $0<v<1$이면 가중치가 양수이므로 $\mathbf1_{(v,1]}$를 폭 $\varepsilon$의 매끄러운 함수로 위아래에서 끼울 수 있다. 정규화한 오차는 $O(\varepsilon+L^{-1}\varepsilon^{-N})$이므로, $\delta=1/(N+1)$에 대해
$$
A_x^{z,w}(\mathbf1_{(v,1]})
=C(r)xL^{r-1}\int_v^1\frac{u^{z-1}(1-u)^{w-1}}{\Gamma(z)\Gamma(w)}\,du
+O_{z,w,v}(xL^{r-1-\delta}).
$$
부호가 섞이는 §6의 두 경우에는 이 양의 가중치 논증을 적용할 수 없다.

## 6. 특수핵과 끝점 항의 확인

$$
a(u)=\frac{\mathbf1_{(0,1)}(u)}{\pi\sqrt{u(1-u)}},
\qquad
h_1(u)=(1-u)a(u),\quad
h_2(u)=\left(u-\frac12\right)a(u)
$$

로 둔다. 도함수는 $\mathbb R$ 위에서 영으로 연장한 함수의 분포적 도함수이다.

**아크사인 핵**

$$
\mathcal K_{1/2,1/2}=a,\qquad
\langle\mathcal K_{1/2,1/2},F\rangle=\int_0^1a(u)F(u)\,du.
$$

**한쪽 지수가 음수인 경우**

$$
\mathcal K_{-1/2,1/2}=h_1',
\tag{6.1}
$$

$$
\langle\mathcal K_{-1/2,1/2},F\rangle
=-\frac1\pi\int_0^1\sqrt{\frac{1-u}{u}}\,F'(u)\,du.
\tag{6.2}
$$

**두 지수가 음수인 경우**

$$
\mathcal K_{-1/2,-1/2}=h_2',
\tag{6.3}
$$

$$
\langle\mathcal K_{-1/2,-1/2},F\rangle
=-\frac1\pi\int_0^1
\frac{u-1/2}{\sqrt{u(1-u)}}\,F'(u)\,du.
\tag{6.4}
$$

두 작용 적분은 보통의 수렴 적분이다. 반면 내부 $0<u<1$의 식은

$$
\mathcal K_{-1/2,1/2}(u)
=-\frac1{2\pi u^{3/2}(1-u)^{1/2}},
$$

$$
\mathcal K_{-1/2,-1/2}(u)
=\frac1{4\pi[u(1-u)]^{3/2}}
$$

이며, 각각 $0$ 또는 양 끝에서 적분 불가능하므로 전체 구간의 밀도로 사용할 수 없다.

**증명.** (3.3)은 분포 항등식이므로

$$
\begin{aligned}
[(1-u)a]'&=\mathcal K_{-1/2,1/2},\\
[(u-\tfrac12)a]'&=-\tfrac12\bigl(\mathcal K_{-1/2,1/2}
+\mathcal K_{1/2,-1/2}\bigr)=\mathcal K_{-1/2,-1/2}.
\end{aligned}
$$

따라서 끝점 항을 포함해 (6.1), (6.3)이 성립한다. $\square$

| $(z,w)$ | $\langle\mathcal K,1\rangle$ | $\langle\mathcal K,u\rangle$ | $\langle\mathcal K,u^2\rangle$ | $\langle\mathcal K,u^3\rangle$ |
|:--|--:|--:|--:|--:|
| $(1/2,1/2)$ | $1$ | $1/2$ | $3/8$ | $5/16$ |
| $(-1/2,1/2)$ | $0$ | $-1/2$ | $-1/4$ | $-3/16$ |
| $(-1/2,-1/2)$ | $0$ | $0$ | $-1/4$ | $-3/8$ |

세 경우의 산술 상수는 각각 $1/\zeta(2),1,1$이고 기준 규모는 $x,x/L,x/L^2$이다. 마지막 경우의 산술적 의미는 §7.4에 있다.

**꼬리값.** $0<v<1$이고 $h$가 $[0,1]$에 지지된 적분 가능 함수이며 $v$ 근방에서 매끄러우면
$\mathbf1_{(v,\infty)}h'=(\mathbf1_{(v,\infty)}h)'-h(v)\delta_v$를 만족한다.
따라서 $T_{z,w}(v):=\langle\mathbf1_{(v,\infty)}\mathcal K_{z,w},\chi\rangle$라 두면
($\chi\in C_c^\infty(\mathbb R)$는 $[0,1]$ 근방에서 $1$),

$$
T_{-1/2,1/2}(v)=-\frac1\pi\sqrt{\frac{1-v}{v}},\qquad
T_{-1/2,-1/2}(v)=\frac{1-2v}{2\pi\sqrt{v(1-v)}}.
\tag{6.5}
$$

이는 분포의 꼬리값이다. 산술 꼬리합의 점근식을 얻으려면 $u=v$ 근방의 합을 추가로 제어해야 하며, 정리 1의 매끄러운 시험함수에 대한 수렴만으로는 충분하지 않다.

## 7. 뫼비우스 경우

### 7.1. 정수 지수

**명제 7.1.** $z,w\in\{0,-1,-2,\dots\}$, $F\in C^\infty([0,1])$이면 모든 $k$에 대해 $A_x^{z,w}(F)=O_{F,z,w,k}(xL^{-k})$이다. $F$가 다항식이면 $O(xe^{-c\sqrt L})$이다.

**증명.** 지수가 정수이므로 $\mathfrak D_{\xi,L}(s)=\zeta(s-i\alpha)^z\zeta(s)^wH$에는 분기점이 없고, 영점 없는 영역 $1-\delta\le\Re s\le c$, $|\Im s|\le T$에서 해석적이다. (A.5)에서 $|\zeta|^{\pm1}\ll(\log\tau)^{O(1)}$이므로 (A.6)도 그대로 성립한다. 따라서 Perron 경로를 우회 없이 $\Re s=1-\delta$로 옮길 수 있고, $|\xi|\le L^{1/4}$, $j=0,1$에서 $\partial_\xi^jS_x(\xi)\ll xe^{-c\sqrt L}$이다. (5.2)에서 $|\xi|\le L^{1/8}$에 대해 $\Theta_x(\xi)\ll xe^{-c\sqrt L}(1+|\xi|)+|\xi|\sqrt xL^{k_0}$이다. §5.3의 고주파 꼬리는 $xL^{k_0-1-(N-1)/8}$이고 $N$은 임의로 크게 잡을 수 있다. $F=u^k$이면 $\sum_{n\le t}\sum_{d\mid n}(\text{가중치})(\log d)^k$의 디리클레 급수 $(-\partial_{s_1})^kD(s_1,s)|_{s_1=s}$도 같은 영역에서 해석적이므로, $(\log t)^{-k}$에 대한 부분합산으로 $xe^{-c\sqrt L}$을 얻는다. $\square$

### 7.2. $(-1,-1)$: 뫼비우스의 곱

제곱 인수 없는 $n=dm$에서 $d,m$은 서로소이므로 쌍 가중치는 $(-1)^{\omega(d)}(-1)^{\omega(m)}=\mu(d)\mu(m)=\mu(n)$이고

$$
A_x^{-1,-1}(F)=\sum_{2\le n\le x}\mu(n)\sum_{d\mid n}F\!\left(\frac{\log d}{\log n}\right).
$$

$F=1$이면 $\sum_{2\le n\le x}\mu(n)\tau(n)$이다. 주항은 서로 독립인 두 이유로 사라진다.

- **산술:** $C(-2)=0$이다. $p=2$의 국소 인수 $1+r/p=1-2/2$가 $0$이기 때문이다.
- **해석:** $\mathcal K_{-1,-1}=0$이다. 지수가 정수이기 때문이다.

$D(s,s)=\prod_p(1-2p^{-s})=\zeta(s)^{-2}H(s,s)$는 $s=1$에서 3중 영점을 가진다($\zeta^{-2}$에서 2, $p=2$에서 1).

서로소 조건을 풀면(§2.1) 쌍 가중치는 $\mu(d)\mu(m)$이고, 디리클레 급수는 정확히 $1/(\zeta(s)\zeta(t))$, 즉 $H^\natural\equiv1$이다. 상수는 $C(-1)^2=1$이 되지만 핵은 여전히 $0$이다. 따라서 서로소 조건을 풀어도 명제 7.1의 소멸 추정은 유지된다.

### 7.3. $(-1,1)$: 한쪽만 뫼비우스

쌍 가중치는 $\mu(d)$이다. $C(0)=1$이고, $g_1(1-u)=\mathbf 1_{u<1}$이 $u=0$ 근처에서 $1$이므로 $\mathcal K_{-1,1}=\delta_0'$, $\mathcal K_{-2,1}=\delta_0''$, $\mathcal K_{-1,0}=0$이다. (1.4)에서

$$
A_x^{-1,1}(F)=-\frac xL\Big[F'(0)+\frac{F'(0)+\kappa(0)F''(0)}{L}\Big]+O\!\left(\frac x{L^3}\right),
\qquad \kappa(0)=\gamma+\sum_p\frac{\log p}{p(p-1)}=1.3325822\ldots
\tag{7.1}
$$

뫼비우스 반전은 이 차수까지 $F$의 원점에서의 1차, 2차 테일러 계수만 남긴다.

- **$F=u$:** $\sum_{d\mid n}\mu(d)\log d=-\Lambda(n)$이므로 좌변은 정확히 $-\pi(x)$이다. (7.1)은 $\pi(x)=x/L+x/L^2+O(x/L^3)$, 즉 $\mathrm{li}(x)$의 앞 두 항이다.
- **$F=u^2$:** 제곱 인수 없는 $n$에 대해 $\sum_{d\mid n}\mu(d)(\log d/\log n)^2$는 $\omega(n)\le2$일 때만 영이 아니고, $n=p$에서 $-1$, $n=pq$에서 $2ab$이다($a=\log p/\log n$, $b=\log q/\log n$). 따라서
$$
\sum_{\substack{p<q\\pq\le x}}\frac{2\log p\log q}{\log^2(pq)}=\pi(x)-\frac{2\kappa(0)\,x}{L^2}+O\!\left(\frac x{L^3}\right).
\tag{7.2}
$$

**서로소 해제 변형.** 쌍 가중치 $\mu(d)\mu^2(m)$이면 (2.6)에서 상수는 $C(-1)C(1)=6/\pi^2$이고, $H_{-1}\equiv1$이므로 $\eta(-1)=0$, $\kappa(-1)=\gamma$이다.

$$
A_x^{\natural,-1,1}(F)=-\frac{6}{\pi^2}\frac xL\Big[F'(0)+\frac{F'(0)+\gamma F''(0)}{L}\Big]+O\!\left(\frac x{L^3}\right).
\tag{7.3}
$$

$6/\pi^2$은 $\mu*\mu^2$가 제곱수에만 놓인다는 사실($\sum(\mu*\mu^2)(n)n^{-s}=1/\zeta(2s)$)에서 온다. **[수치]** B.1, B.2.

### 7.4. $(-\tfrac12,-\tfrac12)$

제곱 인수 없는 $n=dm$에서 쌍 가중치는 $\mu(n)/\tau(n)$이다. 따라서 $F=1,u$에 대한 합은 정확히 $M(x)-1$, $(M(x)-1)/2$이며, §6 표의 영인 주계수와 일치한다. 여기서 $M(x)=\sum_{n\le x}\mu(n)$이다.

## 8. 베타 급수의 수렴 반경과 $s=\tfrac12$

### 8.1. 베타 급수

부록 A.3의 국소 인수 $G$를 전개하여

$$
G(s_1,s_2)=\sum_{a,b\ge0}G_{ab}(s_1-1)^a(s_2-1)^b
\tag{8.1}
$$

로 쓴다. $G_{00}=C$, $G_{10}=G_1$, $G_{01}=G_2$이다.

**명제 8.1.** 각 $J\ge1$에 대해, 같은 범위의 $\xi$에서

$$
S_x(\xi)=xL^{r-1}\sum_{a+b<J}G_{ab}L^{-a-b}\mathcal B_{z-a,w-b}(\xi)+O\!\left(xL^{\beta-1-J}(1+|\xi|)^{B_J}\right).
\tag{8.2}
$$

$J=2$가 (5.1)이고, 증명은 부록 A.3의 테일러 전개를 $J$차까지 늘린 것이다. 즉 베타핵 전개는 $G$의 $s=1$에서의 테일러 전개이다.

$j$차 핵들의 $k$차 모멘트를 함께 모으면 $\sum_{a+b=j}G_{ab}\langle\mathcal K_{z-a,w-b},u^k\rangle=m_{j,k}/\Gamma(r-j+k)$이고

$$
m_{j,k}=\sum_{a+b=j}G_{ab}\,(z-a)_k.
\tag{8.3}
$$

$m_{j,0}$은 $G(1+t,1+t)$의 테일러 계수이다. $m_{j,1}$, $m_{j,2}$에는 대각선 위의 $\partial_1G$, $\partial_1^2G$가 들어간다.

### 8.2. 수렴 반경과 핵들 사이의 상쇄

(2.7)에서 가장 가까운 실축 쪽 특이점 후보는 $s_1=\tfrac12$, $s_1+s_2=1$, $s_2=\tfrac12$이다. 소거되지 않는 후보는 해당 방향의 테일러 수렴 반경을 $\tfrac12$로 제한한다. 더 왼쪽에는 고차 소수 거듭제곱의 층과 $1/s_2$의 극이 있다.

대각선에서는 (2.8)에 의해 세 특이점이 $\zeta(2s)^{-r(r+1)/2}$로 합쳐진다.

- **$r(r+1)/2$가 0 이상의 정수인 경우:** 합친 결과가 $s=\tfrac12$에서 해석적이므로 $m_{j,0}$의 반경은 다음 특이점까지 늘어난다. $r=0,\pm1$에서는 $s=0$의 극까지 거리 $1$이다. 즉 같은 차수의 핵들이 총질량에서 서로 지워진다.
- **비상수 모멘트:** $\partial_1^2G$ 등은 대각선에서도 $s=\tfrac12$에서 특이하다. 그래서 상쇄가 불완전하다.

**[개요]** $r\notin\mathbb Z$에서 $1/\Gamma(r-j)$의 계승 성장은, 반경 $\tfrac12$인 모멘트 급수의 최소항이 $j\approx L/2$에 있고 크기가 대략 $x^{1/2}$임을 시사한다. (8.2)는 고정된 $J$의 추정이므로, 이를 최적 절단 오차로 주장하려면 상수와 나머지의 $J$에 대한 균등 추정을 별도로 증명해야 한다.

### 8.3. 아크사인의 세 점 핵

$z=w=\tfrac12$이면 $e_1=e_2=\tfrac38$, $zw=\tfrac14$이고, 지수 합이 $1$이다. 대각선에서 $H(s,s)=1/\zeta(2s)$는 $s=\tfrac12$에서 오히려 영점이므로 $m_{j,0}$에는 $s=\tfrac12$이 보이지 않는다. $m_{j,1}$도 $d\leftrightarrow m$ 대칭 때문에 보이지 않는다. 반면 $m_{j,2}$의 생성함수에는 대각선에서 다음 단순극이 생긴다.

$$
\partial_1^2\big[(2s_1-1)^{3/8}(s_1+s_2-1)^{1/4}\big](2s_2-1)^{3/8}\Big|_{s_1=s_2=s}=-\frac34(2s-1)^{-1}.
$$

따라서

$$
m_{j,2}\sim c\,(-2)^j,\qquad c=-\frac3{16}\cdot2Z(\tfrac12)=-\frac{3}{16}\bigl(-\zeta(\tfrac12)\bigr)=-0.27382\ldots
\tag{8.4}
$$

**[수치]** B.3.

비틀림이 있으면 $s-\tfrac12=q/L$에서 $H\approx\frac2L(q-i\xi)^{3/8}(q-\tfrac{i\xi}2)^{1/4}q^{3/8}$이다($\widetilde H(\tfrac12,\tfrac12)=1$). $s=1$에서 뻗은 두 절단을 위아래로 감싸는 Hankel 고리 위에서는 $\zeta(s_1)^{1/2}\zeta(s_2)^{1/2}\approx\zeta(\tfrac12)$이다. 세 점 핵을

$$
\langle\mathcal T,e^{i\xi u}\rangle=\frac1{2\pi i}\int_{\mathcal U_Y}e^q(q-i\xi)^{3/8}\Big(q-\frac{i\xi}2\Big)^{1/4}q^{3/8}\,dq
\tag{8.5}
$$

로 정의하면 $\mathcal T$는 $[0,1]$에 지지된 분포이고

$$
\langle\mathcal T,1\rangle=0,\qquad\langle\mathcal T,u\rangle=0,\qquad\langle\mathcal T,u^2\rangle=-\frac3{16}.
$$

질량은 $q$의 총 차수가 $1$이라서, 1차 모멘트는 상수항의 Hankel 적분이라서 $0$이다. 2차 모멘트는 $(1-t)^{3/8}(1-t/2)^{1/4}$의 $t^2$ 계수 $-3/32$에서 온다.

**[개요]** 고정된 $\xi$에서 $s=\tfrac12$ 근방의 국소 기여는 $4\zeta(\tfrac12)\,x^{1/2}L^{-2}\langle\mathcal T,e^{i\xi u}\rangle+O_\xi(x^{1/2}L^{-3})$이다. $F=u^2$(S 수준)이면 $-\tfrac34\zeta(\tfrac12)x^{1/2}L^{-2}$이고, 이는 §9.2의 정확한 계산과 일치한다. 다음 층은 소수 세제곱 항의 $s=\tfrac13$이다.

$r\in\mathbb Z$이면 $1$과 $1+i\alpha$의 절단을 두 점 사이의 선분으로 묶을 수 있다. 두 점을 함께 도는 모노드로미가 $e^{-2\pi ir}=1$이므로 $s=\tfrac12$ 근방의 가지가 단일값이고, 그 극한은 $\zeta(\tfrac12)^r$이다. 이는 $r=-1$에도 적용된다. 일반 정수 $r$의 세 점 핵은 (8.5)의 지수들을 $e_1,zw,e_2$로 바꾸어 얻으며, 총차수는 $r(r+1)/2$이다. 비정수 $r$에서는 절단의 도약을 따로 처리해야 한다.

### 8.4. 임계띠와의 관계

실축 근처에서는 무조건적으로 $\Re s=\tfrac12$까지 연장된다. $\zeta$는 $(0,1)$ 위에서, 그리고 $|t|<14$에서 영이 아니기 때문이다. 세 점 핵을 더하면 $\tfrac13$까지 간다. 베타 급수는 실축 근처의 특이점만 보므로, 전체 오차를 막는 $|\Im s|\ge14.13$ 쪽 경로는 핵을 어떻게 재배열해도 건드릴 수 없다.

**명제 8.2.** 모든 $\varepsilon>0$에 대해 $\sum_{n\le x}(\mu*\mu)(n)\ll x^{\theta+\varepsilon}$이면 $\sigma>\theta$에서 $\zeta(s)\ne0$이다. 역으로 RH이면 $\sum_{n\le x}(\mu*\mu)(n)\ll x^{1/2+\varepsilon}$이다.

**증명.** 부분합산으로 $1/\zeta(s)^2=\sum(\mu*\mu)(n)n^{-s}$가 $\sigma>\theta$에서 수렴하고 해석적이다. 그 영역의 영점은 극이 되므로 있을 수 없다. 역은 RH 아래 $\sigma\ge\tfrac12+\varepsilon$에서 $1/\zeta(s)\ll|t|^\varepsilon$인 것과 Perron 공식으로 따른다. $\square$

따라서 뫼비우스 경우의 거듭제곱 오차 추정은 영점 없는 영역과 직결된다.

## 9. 영점 층과 복원

### 9.1. 반사된 핵

**[개요]** 중복도 $m$인 비자명 영점 $\rho=\beta_\rho+i\gamma$에서 $H$가 $(\rho,\rho)$ 근방에 정칙이라고 하자. 이 정칙성은 $\beta_\rho\ge\tfrac12$이면 (2.7)로 확보된다. $b_\rho=\zeta^{(m)}(\rho)/m!$로 두고 가지를 고정하면 $\zeta(s_1)^z\zeta(s_2)^w\approx b_\rho^r(s_1-\rho)^{mz}(s_2-\rho)^{mw}$이다. $s-\rho=q/L$로 옮긴 국소 기여는 고정된 $\rho,\xi$에서

$$
x^\rho L^{-mr-1}\,\frac{b_\rho^rH(\rho,\rho)}{\rho}\,\mathcal B_{-mz,-mw}(\xi)
+O_{\rho,z,w,\xi}\!\left(x^{\beta_\rho}L^{-m\Re r-2}\right)
\tag{9.1}
$$

이다. 반사된 핵은 $\mathcal K_{-mz,-mw}$이며, $m=1$이면 $\mathcal K_{-z,-w}$이다. 이 국소식에는 RH가 필요 없다. RH 아래에서는 $x$의 거듭제곱이 $x^{1/2}$이고, 로그 인자를 포함한 기준 규모는 $x^{1/2}L^{-m\Re r-1}$이다.

### 9.2. 아크사인: 2차 모멘트부터 영점이 보인다

$z=w=\tfrac12$이고 $\Re\rho\ge\tfrac12$인 단순 영점에서는 $H(\rho,\rho)=1/\zeta(2\rho)$이고 핵은 $\mathcal K_{-1/2,-1/2}=h_2'$이다. §6의 표에서 $\langle h_2',1\rangle=\langle h_2',u\rangle=0$, $\langle h_2',u^2\rangle=-\tfrac14$이다. 따라서 이 영점의 기여는 2차 모멘트부터 나타난다. 실제로 제곱 인수 없는 수의 개수의 급수 $\zeta(s)/\zeta(2s)$는 $\rho$에서 특이점이 없다.

이 핵은 아크사인의 이차 보정 (1.5)에도 나타난다. 대칭 핵에서
$z(\mathcal K_{z-1,z}+\mathcal K_{z,z-1})=\frac{z}{z-1}\mathcal K_{z-1,z-1}$ ($z\ne1$)이므로,
보정 핵과 반사된 핵의 지수가 일치하는 경우는 $z-1=-z$, 즉 $z=\tfrac12$이다.

**정확한 확인.** $D$를 $n$의 균등한 무작위 약수라 하면 $\log D=\sum_{p\mid n}\varepsilon_p\log p$이고, $\varepsilon_p$는 독립인 공정한 동전이다. 따라서

$$
2^{-\omega(n)}\sum_{d\mid n}(\log d)^2=\frac14\Big(\log^2n+\sum_{p\mid n}\log^2p\Big).
\tag{9.2}
$$

$T_2(x)=\sum_{n\le x}\mu^2(n)\sum_{p\mid n}\log^2p$의 디리클레 급수는 $\dfrac{\zeta(s)}{\zeta(2s)}\sum_p\dfrac{\log^2p\,p^{-s}}{1+p^{-s}}$이다.

- **영점에서:** $\sum_p\log^2p\,p^{-s}=(\zeta'/\zeta)'(s)-\sum_{k\ge2}k\sum_p\log^2p\,p^{-ks}\approx-(s-\rho)^{-2}$이므로, 단순극과 유수 $-\zeta'(\rho)/\zeta(2\rho)$가 생긴다. 이를 $4L^2$로 나눈 S 수준 $u^2$ 통계의 영점 항은 (9.1)에 $\langle h_2',u^2\rangle=-\tfrac14$를 넣은 것과 정확히 같다.
- **$s=\tfrac12$에서:** $-3\sum_p\log^2p\,p^{-2s}\approx-\tfrac34(s-\tfrac12)^{-2}$이고 $\zeta(s)/\zeta(2s)\approx2\zeta(\tfrac12)(s-\tfrac12)$이다. 그래서 진동 없는 편향 $-3\zeta(\tfrac12)\sqrt x=4.381\sqrt x$가 생긴다. S 수준 $u^2$로는 $-\tfrac34\zeta(\tfrac12)x^{1/2}L^{-2}=4\zeta(\tfrac12)\langle\mathcal T,u^2\rangle x^{1/2}L^{-2}$로, (8.5)의 세 점 핵과 일치한다.
- **더 낮은 층:** $s=\tfrac13$의 $x^{1/3}(\log)$ 항, 그리고 $\zeta(2s)$ 영점에서 오는 진동수 $\gamma/2$의 $x^{1/4}(\log x)^2$ 항이 있다.

**보조 통계.** $c_n=\sum_{d\mid n}\Lambda(d)\log d$로 두면 제곱 인수 없는 $n$에서 $c_n=\sum_{p\mid n}\log^2p$, 즉 $\log D$ 분산의 4배이다. 디리클레 급수 $Q(s)=\zeta(s)(\zeta'/\zeta)'(s)$에는 $s=\tfrac12$ 항도 $\zeta(2s)$ 층도 없다. 절단된 명시적 공식을

$$
\sum_{n\le x}c_n=xP_2(\log x)-2\,\mathrm{Re}\sum_{0<\gamma\le T}\frac{\zeta'(\rho)}{\rho}\,x^\rho+R_2(x,T)
\tag{9.3}
$$

로 쓴다. $xP_2(\log x)=\operatorname{Res}_{s=1}Q(s)x^s/s$이고, 합은 서로 다른 비자명 영점에 걸친다. 중복도 $m\ge2$이면 $Q(s)=O((s-\rho)^{m-2})$로 정칙이고 $\zeta'(\rho)=0$이므로 그 항은 사라진다. $R_2$의 RH 아래 평균제곱 추정은 부록 C에서 증명한다. **[수치]** B.4.

### 9.3. 모멘트의 차수와 영점의 중복도

공정한 동전의 누적량은 $\kappa_2=\tfrac14$, $\kappa_3=0$, $\kappa_4=-\tfrac18$이고, $\kappa_k(\log D)=\kappa_k(\varepsilon)\sum_{p\mid n}(\log p)^k$이다.

**명제 9.1.** $\Re\rho\ge\tfrac12$인 중복도 $m$의 비자명 영점에서, $\sum_n\mu^2(n)\sum_{p\mid n}(\log p)^k\,n^{-s}$는 $k>m$이면 $k-m$차 극을 가지고, $k\le m$이면 정칙이다.

**증명.** 해당 급수는 $\zeta(s)/\zeta(2s)$와 $\sum_p(\log p)^kp^{-s}/(1+p^{-s})$의 곱이다. 두 번째 인수에서 $(-d/ds)^{k-1}(-\zeta'/\zeta)(s)\sim-m(k-1)!(s-\rho)^{-k}$를 빼면 $\rho$ 근방에 정칙이다. 소수 제곱 항도 $\zeta(2\rho)\ne0$에 의해 이 근방으로 연장된다. 첫 번째 인수는 정확히 $m$차 영점이므로 결론이 따른다. $\square$

진폭이 $\zeta'(\rho)$이므로 2차 모멘트는 단순 영점만 본다. 특히 $C_2$의 오차에서 얻는 영점 배제도 단순 영점에 한정된다. 아크사인에서는 홀수 누적량이 $0$이라 이중 영점은 4차부터 보인다.

### 9.4. 영점 높이와 미분값의 복원

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

$$
(\mu*c)(n)=\Lambda(n)\log n,\qquad
\frac{a_\rho^{(2)}}{a_\rho^{(\psi)}}=\zeta'(\rho).
\tag{9.4}
$$

공통 컨볼루션 필터의 응답이 복소 주파수 $\nu_\rho=\gamma-i(\beta_\rho-\tfrac12)$에서 정의되고 영이 아니면 계수비가 보존된다. 조화 역산은 $\nu_\rho$를 맞추어 $\beta_\rho$도 추정할 수 있다. 임계선 밖 후보는 별도의 영점 검증으로 확인한다.

높은 영점의 후속 실험에서는 목표 높이로 주파수를 이동하고 해당 대역의 진폭을 보존하는 필터를 사용한다. 절단오차와 경계효과를 평가하며 FFT와 조화 역산(필터 대각화; Main 등)을 비교한다. 성능은 영점 위치를 입력하지 않고 얻은 추정 오차와 누락 수로 측정한다.

### 9.5. RH 아래 극한분포

**명제 9.2.** RH 아래
$$
E_2(y)=e^{-y/2}\bigl(C_2(e^y)-e^yP_2(y)\bigr)
$$
는 $B^2$-거의 주기함수이며 극한분포 $\nu$를 가진다. 추가로 단순 영점의 서로 다른 양의 높이가 $\mathbb Q$ 위에서 선형독립(LI)이면

$$
\widehat\nu(\xi)=\prod_{\substack{\gamma>0\\\rho\text{ simple}}}
J_0\!\left(\frac{2|\zeta'(\rho)|\,\xi}{|\rho|}\right).
\tag{9.5}
$$

**증명.** 부록 C는 (9.3)의 정규화한 나머지 $\mathcal E(y,T)=e^{-y/2}R_2(e^y,T)$에 대해
$$
\lim_{Y\to\infty}\frac1Y\int_{\log2}^{Y}|\mathcal E(y,e^Y)|^2\,dy=0
$$
을 준다. 볼록성 상계와 Cauchy 추정으로 $\zeta'(\tfrac12+i\gamma)\ll_\varepsilon\gamma^{1/4+\varepsilon}$이므로 $r_\gamma=-2\zeta'(\rho)/\rho\ll_\varepsilon\gamma^{-3/4+\varepsilon}$이다. 영점 개수 $N(T+1)-N(T)\ll\log T$와 함께 Akbary–Ng–Shahabi의 따름정리 1.3(a)를 적용하면 $B^2$-거의 주기성과 극한분포의 존재가 따른다. LI 아래 곱 공식은 같은 논문의 정리 1.9이다. $\square$

중복 영점은 계수가 $0$이므로 단순영점 가정이 필요 없다. 단순 영점에서 $a_\rho^{(M)}=-(a_\rho^{(\psi)})^2/a_\rho^{(2)}$이지만, 그 제곱합의 제어에는 여전히 역미분값의 추정이 필요하다. 이 결과로 Ng의 $J_{-1}$ 가정을 제거할 수는 없다.

## 부록 A. 이동 균등 추정

### A.1. 계수 상계와 Perron 절단

§5의 $k_0$를 쓰고 $j=0,1$로 둔다. 계수는

$$
|\partial_\xi^j b_{\xi,L}(n)|
\le\tau_{k_0}(n)\left(\frac{\log n}{L}\right)^j
\tag{A.1}
$$

로 지배된다. $j=0,n=1$의 마지막 인자는 $1$로 읽는다.

짧은 구간에서는 $X\ge2,\ 1\le h\le X$에 대해

$$
\sum_{\substack{X/2<n\le2X\\|n-X|\le h}}\tau_{k_0}(n)
\ll_{k_0}
h(\log(2X))^{k_0-1}
+X^{1-1/k_0}(\log(2X))^{k_0-2}.
\tag{A.2}
$$

각 순서 있는 인수분해에서 가장 큰 인수를 제거하면 나머지 곱은 $(2X)^{1-1/k_0}$ 이하이다. 따라서 왼쪽은
$\ll_{k_0}\sum_{d\le(2X)^{1-1/k_0}}\tau_{k_0-1}(d)(h/d+1)$이며, 표준 약수합 상계로 (A.2)가 따른다.

$$
c=1+\frac1L,\qquad T=e^{\sqrt L},\qquad
\mathfrak D_j(s):=\partial_\xi^j\mathfrak D_{\xi,L}(s)
$$

로 둔다. 유효 Perron 공식은

$$
\partial_\xi^jS_x(\xi)
=\frac1{2\pi i}\int_{c-iT}^{c+iT}
\mathfrak D_j(s)\frac{x^s}{s}\,ds+O(E_j),
\tag{A.3}
$$

$$
E_j\ll\sum_{n\ge1}\tau_{k_0}(n)
\left(\frac{\log n}{L}\right)^j
\left(\frac xn\right)^c
\min\left(1,\frac1{T|\log(x/n)|}\right)
$$

를 준다. 최소값은 $n=x$에서 $1$로 읽으며, 정수 끝점의 반가중치 차이도 포함한다. 이 오차는 역 Mellin 핵의 꼬리 부분적분과 핵 자체의 균등 상계에서 얻는다.

$n\notin[x/2,2x]$의 기여는
$\sum\tau_{k_0}(n)(\log n)^jn^{-c}\ll L^{k_0+j}$로 처리한다. 중앙에서는 $|n-x|\le x/T$와 이진 구간들에 (A.2)를 적용한다. 그 결과, 어떤 $c_1>0$에 대해

$$
E_j\ll_{z,w}\frac{x}{T}L^{k_0+1}
+x^{1-1/k_0}L^{k_0}
\ll_{z,w}xe^{-c_1\sqrt L}.
\tag{A.4}
$$

### A.2. 영점 없는 영역과 윤곽선

$\tau=|t|+4$로 두자. 고전적 영점 없는 영역의 상수 $a_0>0$를 충분히 작게 잡으면
$\sigma>1-a_0/(2\log\tau)$, $|t|\ge7/8$에서

$$
|\log\zeta(\sigma+it)|\le\log\log\tau+O(1),
\qquad
\left|\frac{\zeta'}{\zeta}(\sigma+it)\right|\ll\log\tau.
\tag{A.5}
$$

로그는 $\sigma>1$의 오일러 곱에서 연속한 가지이다. 같은 영역의 $|t|\le7/8$, $\sigma\le2$에서는
$\log((s-1)\zeta(s))=O(1)$,
$\zeta'/\zeta=-1/(s-1)+O(1)$이다.
이는 [Montgomery–Vaughan, 정리 6.6–6.7](https://personal.science.psu.edu/rcv4/personal/Publications/MNTI/10.0_pp_168_198_The_Prime_Number_Theorem.pdf)의 추정이다.

이제 $|\xi|\le L^{1/4}$에서

$$
\delta=\frac{a_0}{100\log(T+5)},\qquad
\alpha=\frac{\xi}{L},\qquad Y=2+|\xi|
$$

로 둔다. $1-\delta\le\Re s\le c,\ |\Im s|\le T$에서 $s$와 $s-i\alpha$는 위 영점 없는 영역 안에 있다. 두 후보 $1,1+i\alpha$에서 왼쪽 변까지 수평 절단선을 두어 가지를 고정한다.

(A.3)의 경로를 왼쪽 직선 $\Re s=1-\delta$로 옮기되, $|\Im s|\le Y/L$에서는

$$
(1-\delta-iY/L)\to(c-iY/L)
\to(c+iY/L)\to(1-\delta+iY/L)
$$

의 세 변으로 우회한다. 높이 $\pm T$의 수평 연결선도 포함한다. 모든 경로는 후보에서 거리가 적어도 $1/L$이며 절단선을 가로지르지 않는다.

(A.5)와 유계 높이의 추정을 각각 적용하면 $\log\zeta=O(\log L)$, $\zeta'/\zeta=O(L)$이다. $H$와 편미분은 유계이고 주파수 미분에는 $1/L$이 생기므로, 어떤 $K=K(z,w)$에 대해

$$
|\mathfrak D_0(s)|+|\mathfrak D_1(s)|\ll_{z,w}L^K.
\tag{A.6}
$$

왼쪽 직선과 위아래 연결선의 기여는 각각
$O(x^{1-\delta}L^K\log(T+5))$, $O(xL^K/T)$이다.
$\delta L\asymp\sqrt L$이므로 둘 다 $O_{z,w}(xe^{-c_2\sqrt L})$로 흡수된다. 이 선택은 절단 오차 $T^{-1}$와 왼쪽 이동의 이득 $e^{-\delta L}$을 같은 지수 규모로 맞춘 것이다.

### A.3. 국소 전개와 주항의 경로 연장

$q=L(s-1)$로 옮긴 국소 경로는

$$
\mathcal U_{Y,\delta L}:
(-\delta L-iY)\to(1-iY)\to(1+iY)\to(-\delta L+iY).
$$

여기서 $|q|/L+|\alpha|=O(L^{-1/2})$이다. $Z(s)=(s-1)\zeta(s)$, $\log Z(1)=0$으로 두고

$$
G(s_1,s_2)=
\frac{Z(s_1)^zZ(s_2)^wH_{z,w}(s_1,s_2)}{s_2},\qquad (s_1,s_2)=(s-i\alpha,s)
$$

로 정의한다. 이는 $(1,1)$ 근방에서 해석적이며 $G(1,1)=C(r)$이다. $\Re s>1$에서
$\log\zeta(s)=\log Z(s)-\log(s-1)$이므로, 이 국소 가지는 A.2에서 선택한 가지와 일치한다. $Z'(1)=\gamma$이므로

$$
\partial_{s_1}G(1,1)=\gamma zC+H_s=G_1,\qquad
\partial_{s_2}G(1,1)=\gamma wC+H_t-C=G_2.
$$

$G_L(q,\xi)=G(1+(q-i\xi)/L,\,1+q/L)$라 쓰면 2차 테일러 정리로

$$
\Big|G_L-C-G_1\frac{q-i\xi}L-G_2\frac qL\Big|\ll_{z,w}\frac{(|q|+|\xi|)^2}{L^2},
\qquad
\Big|\partial_\xi\Big(G_L+G_1\frac{i\xi}{L}\Big)\Big|\ll_{z,w}\frac{|q|+|\xi|}{L^2}.
\tag{A.7}
$$

Perron 인자 $1/s$까지 포함한 국소 적분은

$$
\frac{xL^{r-1}}{2\pi i}
\int_{\mathcal U_{Y,\delta L}}
e^q\partial_\xi^j\!\left(P_\xi(q)G_L(q,\xi)\right)\,dq.
$$

$\mathfrak D_j$에 먼저 Perron 공식을 적용했으므로 움직이는 경로를 미분할 필요는 없다.

$m=\max(0,-\Re z)+\max(0,-\Re w)$로 두자. 경로 위의 $|q|,|q-i\xi|\ge1$로부터
$|\partial_\xi^jP_\xi(q)|\ll_{z,w}(1+|q|+|\xi|)^m$이다. 수직 부분의 길이는 $2Y$, 수평 부분에는 $e^{\Re q}$의 감소가 있으므로, $G_L^{[1]}=C+G_1(q-i\xi)/L+G_2q/L$로 두면

$$
\int_{\mathcal U_{Y,\delta L}}
\left|e^q\partial_\xi^j\{P_\xi(q)(G_L-G_L^{[1]})\}\right|\,|dq|
\ll_{z,w}\frac{(1+|\xi|)^{m+3}}{L^2}.
\tag{A.8}
$$

따라서 정수 $B\ge m+3$를 택하면

$$
\begin{aligned}
\partial_\xi^jS_x(\xi)
&=\frac{xL^{r-1}}{2\pi i}
\int_{\mathcal U_{Y,\delta L}}e^q\partial_\xi^j\Big[P_\xi(q)\Big(C+G_1\frac{q-i\xi}L+G_2\frac qL\Big)\Big]\,dq\\
&\quad+O_{z,w}\!\left(xL^{\beta-3}(1+|\xi|)^B\right).
\end{aligned}
\tag{A.9}
$$

(A.9)의 주항에서 두 수평 끝을 $-\infty$까지 늘릴 때 추가 적분은

$$
\ll_{z,w}e^{-\delta L}(1+\delta L+Y)^{m+1}
$$

이다. 앞의 $xL^{r-1}$을 곱해도 지수적으로 작아 (A.9)의 오차에 흡수된다. (4.4)와 그 미분판을 $(z,w)$, $(z-1,w)$, $(z,w-1)$에 적용하면

$$
\frac1{2\pi i}\int_{\mathcal U_Y}e^q\,\partial_\xi^j\Big[P_\xi(q)\Big(C+G_1\frac{q-i\xi}L+G_2\frac qL\Big)\Big]dq
=\partial_\xi^j\Big[C\mathcal B_{z,w}+\frac1L\big(G_1\mathcal B_{z-1,w}+G_2\mathcal B_{z,w-1}\big)\Big](\xi)
$$

이고, 이를 대입하면 (5.1)을 얻는다. 추정 과정에서 $\xi$로 나누지 않았으므로 후보의 합류와 상쇄도 포함한다. 테일러 전개를 $J$차까지 늘리면 같은 방법으로 (8.2)를 얻는다.

## 부록 B. 수치 확인

**방법.** 제곱 인수 없는 수와 $\omega$를 체로 구하고 약수 쌍을 직접 합했다($x\le10^7$).
$T_2$와 $c_n$은 각각 $5\cdot10^7$, $3\cdot10^7$까지 계산했다.
영점과 $\zeta'(\rho)$, $\zeta(2\rho)$는 Euler–Maclaurin 공식으로 계산했으며, 높이 $1000$까지의 영점 649개를 사용했다.

### B.1. 이차 법칙 (정리 1)

값은 $A/(C\,xL^{r-1})$, $x=10^7$이다.

| $(z,w)$ | $F$ | 실측 | 1차 | 2차 |
|:--|:--|--:|--:|--:|
| $(\tfrac12,\tfrac12)$ | $u^2$ | 0.39957 | 0.37500 | 0.40163 |
| $(\tfrac12,\tfrac12)$ | $(u-\tfrac12)^4$ | 0.03173 | 0.02344 | 0.03343 |
| $(1,1)$ | $u^2$ | 0.43441 | 0.33333 | 0.43526 |
| $(2,1)$ | $u^2$ | 0.46624 | 0.25000 | 0.44359 |
| $(2,1)$ | $1$ | 0.87273 | 0.50000 | 0.84226 |
| $(-1,1)$ | $u$ | −1.07117 | −1.00000 | −1.06204 |
| $(-1,1)$ | $u^2$ | −0.19026 | 0 | −0.16535 |

남은 오차에 $L^2$을 곱하면 $10^6$과 $10^7$에서 거의 같다: $(\tfrac12,\tfrac12),u^2$은 −0.550→−0.535, $(2,1),1$은 7.91→7.92, $(-1,1),u$는 −2.31→−2.37. 이 계산은 예측한 $L^{-2}$ 오차와 부합한다.

### B.2. 뫼비우스 경우 (§7)

**서로소 해제 $(-1,1)$, (7.3):** $A^\natural L/x$의 값이다($10^7$).

| $F$ | 실측 | 1차 $-\tfrac6{\pi^2}F'(0)$ | 2차 |
|:--|--:|--:|--:|
| $u$ | −0.6514 | −0.6079 | −0.6456 |
| $u^2$ | −0.0501 | 0 | −0.0435 |
| $e^u$ | −0.6772 | −0.6079 | −0.6674 |

남은 오차에 $L^2$을 곱한 값은 $10^6$과 $10^7$에서 각각 −1.52→−1.50, −1.70→−1.71, −2.53→−2.54이다.

### B.3. 베타 급수의 계수 (§8, 아크사인)

$G$를 다중원 $|s_i-1|=0.3$ 위에서 계산하고 2차원 FFT로 $G_{ab}$를 구했다. 검산으로 $G_{00}=0.6079271$이다.

| $j$ | 개별 $\max\lvert G_{ab}\rvert$ | $m_{j,0}$ | $m_{j,1}$ | $m_{j,2}$ | $m_{j,2}/m_{j-1,2}$ |
|--:|--:|--:|--:|--:|--:|
| 8 | 5.3 | −0.9947 | 0.8949 | −68.89 | −1.981 |
| 10 | 14.5 | −0.9987 | 0.9116 | −273.26 | −1.995 |
| 12 | 42.8 | −0.9997 | 0.9168 | −1094.77 | −2.003 |
| 14 | 134.1 | −0.9999 | 0.9183 | −4401.24 | −2.005 |

- **$m_{j,0}$:** $s=0$의 극에서 예측되는 $(-1)^{j+1}$로 간다.
- **$m_{j,2}$:** 비가 $-2$로 가므로 $s=\tfrac12$의 단순극이 보인다. $m_{14,2}/(-2)^{14}=-0.2686$이 (8.4)의 $-0.2738$로 다가간다. 느리게 다가가는 것은 $s=\tfrac13$의 다음 층 때문으로 보인다.

### B.4. 영점 층 (§9)

**알려진 영점으로 잔차 재현.** (9.3)을 적분한 $\int_0^x C_2(t)\,dt$의 잔차를 $x^{3/2}$로 나눈 것과, 알려진 영점 $N$개의 합 $-2\mathrm{Re}\sum\zeta'(\rho)x^{\rho+1}/(\rho(\rho+1))$를 $x\in[10^6,3\cdot10^7]$에서 비교했다.

| 입력 영점 수 | 1 | 10 | 79 | 649 |
|:--|--:|--:|--:|--:|
| 상관 | 0.594 | 0.928 | 0.991 | 0.995 |

잔차의 표준편차는 0.0087, 차이의 표준편차는 0.0009이다. 평활화하지 않은 합에서는 상관이 0.126, 0.383, 0.654, 0.837이다.

**스펙트럼.** 잔차를 $\log x$에 대해 푸리에 변환한 앞의 세 봉우리는 다음과 같다.

| 봉우리 | 14.20 | 21.08 | 25.04 |
|:--|--:|--:|--:|
| $\gamma$ | 14.13 | 21.02 | 25.01 |

**아크사인 설정의 $T_2$.** $x(AL^2+BL+C)+c\sqrt x$로 피팅한 편향은 $c=4.335$로, 예측 $-3\zeta(\tfrac12)=4.381$과 맞는다. 선도 계수 $A=0.30394$도 $1/(2\zeta(2))=0.30396$과 맞는다. 진동은 $\zeta(2s)$ 층과 $s=\tfrac13$ 항이 섞여, 영점 649개로 상관 0.70까지만 나온다.

## 부록 C. $C_2$의 절단오차: RH 아래 평균제곱 추정

$Q(s)=\zeta(s)(\zeta'/\zeta)'(s)$, $a=\tfrac14$로 둔다. 소인수분해에서
$$
0\le c_n=\frac12\sum_{p^h\parallel n}h(h+1)(\log p)^2\le(\log n)^2.
$$
따라서 $2\le x\le T$, $c=1+1/\log x$에서 유효 Perron 공식의 오차는
$$
O\!\left(\frac{x\log^3(2x)}T+\log^2(2x)\right)
\tag{C.1}
$$
이다. $n\asymp x$에서는 $\min(1,x/(T|n-x|))$를 합하고, 나머지는 절대수렴으로 처리한다. 정수 끝점의 반가중치 차이도 이 오차에 포함된다.

RH 아래 $\Re s=a$는 모든 비자명 영점과 거리가 적어도 $1/4$이다. 로그미분의 부분분수 전개와 $N(t+1)-N(t)\ll\log(t+2)$에서 $(\zeta'/\zeta)'(a+it)\ll\log(|t|+2)$를 얻는다. 볼록성 상계를 곱하면, 고정된 $0<\varepsilon<1/8$에 대해
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
지나는 극은 $s=1$과 단순 영점뿐이며 그 유수가 (9.3)의 주항과 영점 합을 준다. 왼쪽 수직선 적분은 $x^{1/4}V_T(\log x)$이고
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
$$
\frac1Y\int_{\log2}^{Y}|\mathcal E(y,T)|^2\,dy
\ll_\varepsilon\frac1Y+(1+Y)^6e^{-(1-2\varepsilon)Y}\longrightarrow0.
\tag{C.4}
$$
$T$를 $X=e^Y$로 바꿀 때 추가되는 영점은 $O(\log X)$개이고, 각 정규화한 계수는 $O_\varepsilon(X^{-3/4+\varepsilon})$이다. 이 차이도 균등하게 $0$으로 가므로 §9.5의 평균제곱 조건이 성립한다.

## 참고 문헌

- H. L. Montgomery and R. C. Vaughan, [Multiplicative Number Theory I, §6.1](https://personal.science.psu.edu/rcv4/personal/Publications/MNTI/10.0_pp_168_198_The_Prime_Number_Theorem.pdf), 정리 6.6–6.7: 영점 없는 영역과 복소 로그 추정.
- G. Tenenbaum, *Introduction to Analytic and Probabilistic Number Theory*, II.5: Selberg–Delange 방법.
- NIST DLMF, [§5.9](https://dlmf.nist.gov/5.9): Hankel 적분. [§13.6](https://dlmf.nist.gov/13.6): Bessel 함수. [§25.16](https://dlmf.nist.gov/25.16): $\psi$의 명시적 공식.
- N. Ng, [The distribution of the summatory function of the Möbius function](https://arxiv.org/abs/math/0310381), Proc. London Math. Soc. (2004): $M(x)$의 극한분포와 $J_{-1}(T)$ 가설.
- A. Akbary, N. Ng, M. Shahabi, [Limiting distributions of the classical error terms of prime number theory](https://arxiv.org/pdf/1306.1657), Q. J. Math.: 따름정리 1.3(a)의 극한분포 기준과 정리 1.9의 LI 아래 곱 공식.
- J. Main, V. A. Mandelshtam, G. Wunner, H. S. Taylor, [Harmonic inversion as a general method for periodic orbit quantization](https://arxiv.org/abs/chao-dyn/9709009) (1998): 필터 대각화를 이용한 영점 복원.
- J. Main, P. A. Dando, Dž. Belkić, H. S. Taylor, [Decimation and harmonic inversion of periodic orbit signals](https://arxiv.org/abs/chao-dyn/9912031): 대역 필터와 조화 역산을 분리한 영점 복원.
