# 주정리: DDT에서 베타핵과 약수별 이산핵으로

**2026-11-03 (DDT 도입 및 이산핵 확장)**, **2026-11-06 (임의 차수 전개, §9)**

**DDT의 평균 법칙 → 일반화된 베타핵 → 약수별 이산핵**의 순서로 전개한다. 정리 1은 제곱 인수 없는 정수의 가중 약수합을, 정리 2는 모든 정수의 $d_z$ 가중 약수합을 이차 항까지 근사한다. 정리 3은 두 정리를 임의의 $J$차로 늘려 오차를 $O(xL^{\beta-1-J})$로 줄인다. 세 정리는 유효 Perron 공식과 고전적 영점 없는 영역을 사용하며 RH를 가정하지 않는다.

주제별 확장은 [후속 연구 안내](further-results.md)에, 정확한 영점 표현과 이산핵 복원은 [전역 공식](global-formula.md)에 있다.

- [§1. DDT에서 베타핵으로](#section-1)
- [§2. 이차 베타 법칙](#section-2)
- [§3. 오일러 곱](#section-3)
- [§4. 분포 핵과 모멘트](#section-4)
- [§5. 이동 인수와 역라플라스 변환](#section-5)
- [§6. 증명](#section-6)
- [§7. 특수핵과 끝점 항](#section-7)
- [§8. 약수별 이산핵으로의 확장](#section-8)
- [§9. 임의 차수 전개](#section-9)
- [부록 A. 이동 균등 추정](#appendix-a)
- [부록 B.1. 주정리의 수치 확인](#appendix-b-1)
- [부록 B.2. J차 전개의 수치 확인](#appendix-b-2)

주제별 수치 기록은 [후속 연구 안내](further-results.md#numerics)에서 찾을 수 있다.

<a id="section-1"></a>

## 1. DDT에서 베타핵으로

$\tau(n)$을 약수의 개수라 하고, $n\ge2$에서 약수를 균등하게 고르는 확률측도를

$$
\nu_n=\frac1{\tau(n)}\sum_{d\mid n}\delta_{\log d/\log n}
$$

로 둔다. Deshouillers–Dress–Tenenbaum(DDT) 정리는 $0\le v\le1$에서 균등하게

<a id="equation-1-1"></a>

$$
\frac1x\sum_{2\le n\le x}\nu_n([0,v])
=\frac2\pi\arcsin\sqrt v+O\!\left((\log x)^{-1/2}\right)
\tag{1.1}
$$

임을 말한다. $n=1$ 항의 생략은 $O(x^{-1})$ 차이만 만든다. 이 고전적 결과의 서술과 베타분포로의 확장은 [Bareikis–Manstavičius, *On the DDT theorem*](https://web.vu.lt/mif/e.manstavicius/wp-content/uploads/2016/08/EM-5.pdf), §1–2를 참조한다.

극한밀도는 $a(u)=1/(\pi\sqrt{u(1-u)})$이다. 여기서 수렴하는 것은 정수들에 대한 평균이다. 개별 $\nu_n$의 수렴을 뜻하지 않는다. 예를 들어 소수 $p$에서는 항상 $\nu_p=(\delta_0+\delta_1)/2$이다.

이 문서는 먼저 제곱 인수 없는 정수에서 이 구조를 살핀다. 이때 $\tau(n)=2^{\omega(n)}$이므로 $(z,w)=(1/2,1/2)$의 가중치는 정확히 $1/\tau(n)$이다. 정리 1은 이 제한된 합의 주항이 $x\langle a,F\rangle/\zeta(2)$임을 보이며, 제곱 인수 없는 정수의 개수로 정규화하면 같은 아크사인 법칙을 얻는다. 이는 모든 정수의 균등 약수를 다루는 (1.1)과 합산 대상이 다르다.

양의 실수 $z,w$로 가중치를 바꾸면 아크사인 핵은

$$
\mathcal K_{z,w}(u)=\frac{u^{z-1}(1-u)^{w-1}}{\Gamma(z)\Gamma(w)}
$$

가 된다. 이 핵의 총질량은 $1/\Gamma(z+w)$이므로, 확률밀도는 $\Gamma(z+w)\mathcal K_{z,w}$이다. §4는 이를 복소 지수까지 분포의 의미로 연장하며, §2–7은 산술 합과 이 연속핵 사이의 이차 점근식을 다룬다.

§8에서는 평균을 내기 전의 약수별 원자들을 $\mathscr K_n^{z,w}$에 보존한다. $d_z$ 가중치를 사용하여 모든 소수 거듭제곱을 포함하되, 제곱 인수 없는 정수에서는 기존 가중치와 일치하게 한다. 이렇게 얻는 이산핵의 평균을 다시 베타핵과 연결하고, 이후 전역 공식을 위한 정확한 Mellin–Fourier 변환을 마련한다.

<a id="section-2"></a>

## 2. 이차 베타 법칙

$\omega(n)$은 서로 다른 소인수의 개수, $\mu$는 뫼비우스 함수이다. $z,w\in\mathbb C$, $F:[0,1]\to\mathbb C$에 대해

$$
A_x^{z,w}(F)=
\sum_{\substack{2\le n\le x\\\mu^2(n)=1}}
\sum_{d\mid n}z^{\omega(d)}w^{\omega(n/d)}
F\!\left(\frac{\log d}{\log n}\right),
\qquad 0^0:=1
$$

로 둔다. 이하 $r=z+w$, $\beta=\Re r$, $L=\log x$이며, $L^\lambda=\exp(\lambda\log L)$이다. $\gamma_E$는 오일러 상수이다.

<a id="equation-2-1"></a>

$$
C(r)=\prod_p\left(1+\frac rp\right)\left(1-\frac1p\right)^r,
\qquad
\mathcal K_{z,w}(u)=g_z(u)g_w(1-u),
\tag{2.1}
$$

여기서 $g_\alpha(t)=t_+^{\alpha-1}/\Gamma(\alpha)$는 §4에서 정의하는 전체해석적 분포 가족이다. §3의 $H=H_{z,w}$에 대해 $H_s=\partial_sH(1,1)$, $H_t=\partial_tH(1,1)$로 쓴다. $H(1,1)=C(r)$이다.

**정리 1 (이차 베타 법칙).** 고정된 $z,w$에 대해 충분히 큰 유한 정수 $N=N(z,w)$가 존재하여, $F\in C^N([0,1])$이면 $x\to\infty$일 때

<a id="equation-2-2"></a>

$$
\boxed{
A_x^{z,w}(F)
=xL^{r-1}\Big\langle C(r)\,\mathcal K_{z,w}+\frac1L\,\mathcal K^{(1)}_{z,w},\ F\Big\rangle
+O_{F,z,w}(xL^{\beta-3}),
}
\tag{2.2}
$$

<a id="equation-2-3"></a>

$$
\mathcal K^{(1)}_{z,w}
=\bigl(\gamma_E zC+H_s\bigr)\mathcal K_{z-1,w}
+\bigl(\gamma_E wC+H_t\bigr)\mathcal K_{z,w-1}
-(r-1)\,C\,\mathcal K_{z,w},
\qquad C=C(r).
\tag{2.3}
$$

(3.4)의 $C\eta$를 전체함수로 연장하고 $C\kappa:=\gamma_E C+C\eta$로 읽으면, 모든 $r\in\mathbb C$에서

<a id="equation-2-4"></a>

$$
\mathcal K^{(1)}_{z,w}
=(C\kappa)(r)\bigl(z\,\mathcal K_{z-1,w}+w\,\mathcal K_{z,w-1}\bigr)-(r-1)C(r)\mathcal K_{z,w},
\qquad
\kappa(r)=\gamma_E+\sum_p\frac{(r+1)\log p}{(p-1)(p+r)}.
\tag{2.4}
$$

증명은 §6과 부록 A에 있다. 오차 $O(xL^{\beta-3})$은 국소 인수와 정규화 변경을 2차에서 자른 결과이다. 일반적으로 다음 항 $xL^{r-3}\langle\mathcal K^{(2)}_{z,w},F\rangle$이 영이 아니므로 이 차수에서는 최선이다. 이를 임의의 $J$차로 늘린 것이 [§9의 정리 3](#section-9-1)이다.

특히 $A_x^{z,w}(F)=C(r)xL^{r-1}\langle\mathcal K_{z,w},F\rangle+O_{F,z,w}(xL^{\beta-2})$이므로, 정규화한 합은 분포의 의미에서 $C(r)\mathcal K_{z,w}$로 수렴한다.

**아크사인 예 $(\tfrac12,\tfrac12)$.** $a(u)=1/(\pi\sqrt{u(1-u)})$라 쓰면, §7의 핵 항등식에서

<a id="equation-2-5"></a>

$$
A_x^{1/2,1/2}(F)=\frac{x}{\zeta(2)}\Big[\int_0^1a(u)F(u)\,du+\frac{\kappa(1)}L\int_0^1\big(u-\tfrac12\big)a(u)F'(u)\,du\Big]+O\!\left(\frac x{L^2}\right),
\tag{2.5}
$$
$$
\kappa(1)=\gamma_E-2\frac{\zeta'(2)}{\zeta(2)}=1.7171378\ldots
$$
이차 보정은 질량을 양 끝점 쪽으로 옮긴다([부록 B.1](#appendix-b-1)).

**총질량 $F=1$.** $\langle\mathcal K_{a,b},1\rangle=1/\Gamma(a+b)$에서
$$
\langle\mathcal K^{(1)}_{z,w},1\rangle=\frac{\Phi'(1)}{\Gamma(r-1)},\qquad
\Phi(s)=\frac{((s-1)\zeta(s))^rH(s,s)}{s}.
$$
이는 $\sum_{n\le x}\mu^2(n)r^{\omega(n)}$에 대한 Selberg–Delange 전개의 둘째 계수와 같다.

양의 실수 $z,w$에서는 지시함수에도 1차 법칙이 성립한다(§6.3).

<a id="section-3"></a>

## 3. 오일러 곱

$\mu^2(dm)=1$이면 각 소수는 $d,m$ 중 한쪽에만 속한다. 따라서 $\Re s,\Re t>1$에서

<a id="equation-3-1"></a>

$$
\begin{aligned}
D_{z,w}(s,t)
&=\sum_{\substack{d,m\ge1\\\mu^2(dm)=1}}
\frac{z^{\omega(d)}w^{\omega(m)}}{d^s m^t}\\
&=\prod_p(1+zp^{-s}+wp^{-t})\\
&=\zeta(s)^z\zeta(t)^wH_{z,w}(s,t),
\end{aligned}
\tag{3.1}
$$

<a id="equation-3-2"></a>

$$
H_{z,w}(s,t)=
\prod_p(1+zp^{-s}+wp^{-t})
(1-p^{-s})^z(1-p^{-t})^w.
\tag{3.2}
$$

거듭제곱의 가지는 오일러 곱의 로그로 정한다. 국소 인수의 일차항이 상쇄되어

$$
H_p(s,t)-1
=O_{z,w}\!\left(p^{-2\Re s}+p^{-\Re(s+t)}+p^{-2\Re t}\right).
$$

그러므로 고정된 $\sigma_0>1/2$에서 $H$와 필요한 유한 차수 편미분들은 $\Re s,\Re t\ge\sigma_0$에서 균등 유계이다. 미분으로 생기는 $(\log p)^j$도 합산 가능하다. 이 논증은 $H$의 영점도 허용한다.

<a id="equation-3-3"></a>

$$
H_{z,w}(1,1)=C(r),\qquad
C(0)=C(-1)=1,\quad
C(1)=\frac1{\zeta(2)},\quad C(-2)=0.
\tag{3.3}
$$

특히 $z=w=1/2$일 때 $x$로 정규화한 총질량은 $1/\zeta(2)$이다. 제곱 인수 없는 정수의 개수로 정규화하면 질량 $1$의 아크사인 법칙이 된다.

국소 인수를 $s$로 미분하면, 우선 $C(r)\ne0$에서

<a id="equation-3-4"></a>

$$
H_s=C(r)\,z\,\eta(r),\qquad H_t=C(r)\,w\,\eta(r),\qquad
\eta(r)=\sum_p\frac{(r+1)\log p}{(p-1)(p+r)},
\tag{3.4}
$$

이다. $C$는 $r=-p$에서 단순영점, $\eta$는 같은 점에서 단순극을 가진다. 그 곱의 연장값은
$$
(C\eta)(-p)=-\frac{\log p}{p}\left(1-\frac1p\right)^{-p}
\prod_{q\ne p}\left(1-\frac pq\right)\left(1-\frac1q\right)^{-p}.
$$
국소 곱의 미분도 $H_s=zC\eta$, $H_t=wC\eta$의 같은 값을 주므로 (3.4)는 모든 $r$로 연장된다. $C\eta$와 $C\kappa=\gamma_E C+C\eta$는 전체함수이다.

<a id="section-3-1"></a>

### 3.1. 서로소 조건을 푼 변형

제곱 인수 없는 $d,m$의 모든 쌍을 서로소 조건 없이 세면

$$
D^\natural_{z,w}(s,t)=\sum_{d,m\ge1}\frac{\mu^2(d)\mu^2(m)z^{\omega(d)}w^{\omega(m)}}{d^sm^t}
=\prod_p(1+zp^{-s})(1+wp^{-t})
=\zeta(s)^z\zeta(t)^wH^\natural(s,t),
$$

<a id="equation-3-5"></a>

$$
H^\natural(s,t)=H_z(s)H_w(t),\qquad H_\alpha(s)=\prod_p(1+\alpha p^{-s})(1-p^{-s})^\alpha
\tag{3.5}
$$

이다. 국소 인수의 차이는 분모 없이

$$
H^\natural_p-H_p
=zwp^{-s-t}(1-p^{-s})^z(1-p^{-t})^w
=O_{z,w}\!\left(p^{-\Re(s+t)}\right)
$$

로 쓸 수 있어 $H_p=0$인 경우도 포함한다. 같은 반평면에서 $H^\natural$와 그 편미분이 유계이고, 계수도 $\tau_{k_0}$로 지배된다. 따라서

$$
A_x^{\natural}(F)=\sum_{2\le dm\le x}\mu^2(d)\mu^2(m)z^{\omega(d)}w^{\omega(m)}F\!\left(\frac{\log d}{\log dm}\right)
$$

에 대해 정리 1이 그대로 성립하고, $C(r),H_s,H_t$가 각각 $C(z)C(w)$, $C(z)C(w)z\eta(z)$, $C(z)C(w)w\eta(w)$로 바뀐다. 위와 같이 곱을 연장하면 모든 $z,w$에서

<a id="equation-3-6"></a>

$$
\mathcal K^{(1)\natural}_{z,w}=C(z)C(w)\Big[z\kappa(z)\mathcal K_{z-1,w}+w\kappa(w)\mathcal K_{z,w-1}-(r-1)\mathcal K_{z,w}\Big].
\tag{3.6}
$$

<a id="section-3-2"></a>

### 3.2. $H$의 2차 구조

$e_1=z(z+1)/2$, $e_2=w(w+1)/2$로 두면

<a id="equation-3-7"></a>

$$
H(s,t)=\zeta(2s)^{-e_1}\,\zeta(s+t)^{-zw}\,\zeta(2t)^{-e_2}\,\widetilde H(s,t),
\tag{3.7}
$$

$$
\widetilde H_p=H_p\,(1-p^{-2s})^{-e_1}(1-p^{-s-t})^{-zw}(1-p^{-2t})^{-e_2}=1+O\!\left(p^{-3\min(\Re s,\Re t)}\right).
$$

$X=p^{-s}$, $Y=p^{-t}$로 두고 충분히 큰 $p$에서 로그를 전개하면 이차항은 $-\tfrac12(zX+wY)^2-\tfrac z2X^2-\tfrac w2Y^2+e_1X^2+zwXY+e_2Y^2=0$으로 상쇄된다. 따라서 $\widetilde H$는 $\Re s,\Re t>\tfrac13$에서 해석적이다. 서로소 조건을 푼 모형에서는 교차 인수 $\zeta(s+t)^{-zw}$가 없다. 이 인수는 두 모형의 **이차항 수준의 차이**를 나타내며, 삼차 이상을 담는 잔여 인자도 달라진다.

실축 근처($|\Im s|,|\Im t|<7$)의 $\Re>\tfrac13$에서 $H$의 특이점은 $s=\tfrac12$, $s+t=1$, $t=\tfrac12$뿐이다. $\zeta(2s)$, $\zeta(s+t)$의 영점은 높이 $14$ 이상에서만 나오기 때문이다. 대각선에서는 세 인수가 합쳐져

<a id="equation-3-8"></a>

$$
H(s,s)=\zeta(2s)^{-r(r+1)/2}\,\widetilde H(s,s)
\tag{3.8}
$$

가 된다. $r=1$이면 $H(s,s)=\prod_p(1-p^{-2s})=1/\zeta(2s)$이고, $r=0,-1$이면 지수가 $0$이다.

<a id="section-4"></a>

## 4. 분포 핵과 모멘트

$\Re\alpha>0$에서 $g_\alpha=t_+^{\alpha-1}/\Gamma(\alpha)$는 국소 적분 가능 함수이다. 비음의 정수 $N$에 대해

<a id="equation-4-1"></a>

$$
\langle g_\alpha,\varphi\rangle
=\frac{(-1)^N}{\Gamma(\alpha+N)}
\int_0^\infty t^{\alpha+N-1}\varphi^{(N)}(t)\,dt,
\qquad \Re\alpha>-N
\tag{4.1}
$$

로 연장한다. $\varphi\in C_c^\infty(\mathbb R)$이며, 겹치는 반평면에서 부분적분으로 일치하므로 전체해석적 분포 가족을 얻는다. 특히

$$
g_0=\delta_0,\qquad g_{-1}=\delta_0',
\qquad g_\alpha'=g_{\alpha-1}.
$$

$g_z(u)$와 $g_w(1-u)$의 특이 지지집합은 각각 $\{0\}$, $\{1\}$이다. 각 끝점에서 다른 인수가 매끄러우므로 곱 $\mathcal K_{z,w}$가 정의되며, 이는 $[0,1]$에 지지된 전체해석적 분포 가족이다. $F\in C^\infty([0,1])$의 작용은 매끄러운 콤팩트 지지 연장으로 정의한다. 연장들은 $[0,1]$에서 모든 제트가 같으므로 결과는 선택에 무관하다.

베타적분을 수렴 영역에서 계산한 뒤 연장하면

<a id="equation-4-2"></a>

$$
\boxed{
\langle\mathcal K_{z,w},u^k\rangle
=\frac{(z)_k}{\Gamma(r+k)}
\qquad(k\ge0).
}
\tag{4.2}
$$

$(z)_k$는 상승계승이며 $(z)_0=1$, $1/\Gamma$는 전체함수로 읽는다. 특히 총질량은 $1/\Gamma(r)$이다. 질량이 영이어도 핵은 비영일 수 있다. 한 지수가 $0$이면

$$
\mathcal K_{0,w}=\frac{\delta_0}{\Gamma(w)},
\qquad
\mathcal K_{z,0}=\frac{\delta_1}{\Gamma(z)}.
$$

두 지수가 모두 비양의 정수이면 두 인수의 지지가 분리되어 핵 자체가 영이다.

**핵 항등식.** $tg_{\alpha-1}(t)=(\alpha-1)g_\alpha(t)$와 $g_\alpha'=g_{\alpha-1}$에서

<a id="equation-4-3"></a>

$$
u\,\mathcal K_{z-1,w}=(z-1)\mathcal K_{z,w},\qquad
(1-u)\,\mathcal K_{z,w-1}=(w-1)\mathcal K_{z,w},\qquad
\partial_u\mathcal K_{z,w}=\mathcal K_{z-1,w}-\mathcal K_{z,w-1}.
\tag{4.3}
$$

이로부터 $(u\mathcal K_{z,w})'=(r-1)\mathcal K_{z,w}-\mathcal K_{z,w-1}$이고, 따라서

<a id="equation-4-4"></a>

$$
\langle\mathcal K_{z,w},uF'\rangle=\langle\mathcal K_{z,w-1},F\rangle-(r-1)\langle\mathcal K_{z,w},F\rangle.
\tag{4.4}
$$

또한 $u^a(1-u)^b\mathcal K_{z,w}=(z)_a(w)_b\,\mathcal K_{z+a,w+b}$이다.

**로그 가중치.** 정수 $a,b\ge0$에 대해 위 항등식과 정리 1의 부분합산으로

<a id="equation-4-5"></a>

$$
\begin{aligned}
&\sum_{\substack{2\le n\le x\\\mu^2(n)=1}}\sum_{d\mid n}
z^{\omega(d)}w^{\omega(n/d)}(\log d)^a
\Big(\log\frac nd\Big)^bF\!\left(\frac{\log d}{\log n}\right)\\
&\qquad=C(r)(z)_a(w)_b\,xL^{r-1+a+b}
\langle\mathcal K_{z+a,w+b},F\rangle
+O_{F,z,w,a,b}(xL^{\beta-2+a+b}).
\end{aligned}
\tag{4.5}
$$

실제로 가중치는 $(\log n)^{a+b}u^a(1-u)^bF(u)$이다. 아크사인에 $\log d\log(n/d)$를 곱하면 핵이 $\tfrac14\mathcal K_{3/2,3/2}=\pi^{-1}\sqrt{u(1-u)}$로 바뀌어, 정규화 후 반원 법칙을 얻는다.

아크사인 핵의 푸리에 변환은 Bessel 함수 $J_0$로 표현된다.

<a id="equation-4-6"></a>

$$
\langle\mathcal K_{1/2,1/2},e^{i\xi u}\rangle
=e^{i\xi/2}J_0\!\left(\frac\xi2\right).
\tag{4.6}
$$

<a id="section-5"></a>

## 5. 이동 인수와 역라플라스 변환

실수 $\xi$에 대해

<a id="equation-5-1"></a>

$$
\begin{aligned}
b_{\xi,L}(n)
&=\mu^2(n)\sum_{d\mid n}
z^{\omega(d)}w^{\omega(n/d)}e^{i\xi\log d/L},\\
S_x(\xi)&=\sum_{n\le x}b_{\xi,L}(n)
\end{aligned}
\tag{5.1}
$$

로 둔다. $n=1$ 항은 $1$이다. 대응하는 디리클레 급수는

<a id="equation-5-2"></a>

$$
\mathfrak D_{\xi,L}(s)
=D_{z,w}\!\left(s-\frac{i\xi}{L},s\right).
\tag{5.2}
$$

(3.1)의 두 제타 인수는 $s=1$, $s=1+i\xi/L$에 특이점 후보를 갖는다. 후보는 합류하거나 상쇄될 수 있다. 예컨대 $r=0,\xi=0$이면 정확히 $D_{z,w}(s,s)=1$, $S_x(0)=1$이다.

$s=1+q/L$로 확대하면 주인수는 $L^rP_\xi(q)$이며

$$
P_\xi(q)=(q-i\xi)^{-z}q^{-w}.
$$

두 거듭제곱은 각각 왼쪽 수평 반직선을 절단선으로 하는 주가지를 쓴다. $Y>|\xi|+1$에 대해

<a id="equation-5-3"></a>

$$
\mathcal U_Y:
(-\infty-iY)\to(1-iY)\to(1+iY)\to(-\infty+iY)
\tag{5.3}
$$

로 경로를 정한다. 이 경로는 두 후보에서 거리가 적어도 $1$이고 $\Re q\le1$이다.

<a id="equation-5-4"></a>

$$
\mathcal B_{z,w}(\xi)
:=\frac1{2\pi i}\int_{\mathcal U_Y}e^qP_\xi(q)\,dq
=\langle\mathcal K_{z,w},e^{i\xi u}\rangle.
\tag{5.4}
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

$\mathcal U_Y$의 적분은 모든 $z,w$에서 절대수렴하고, 매개변수의 콤팩트 집합에서 균등 수렴한다. 경로 높이의 독립성은 Cauchy 정리로, (5.4)의 전 매개변수 확장은 전체해석적 연장으로 따른다.

(5.4)는 모든 $z,w$에 대해 성립하므로, 특히 $(q-i\xi)P_\xi(q)$, $qP_\xi(q)$의 적분은 각각 $\mathcal B_{z-1,w}(\xi)$, $\mathcal B_{z,w-1}(\xi)$이다. 이차 법칙의 보정 핵은 이렇게 나온다.

콤팩트 지지 분포의 유한 차수 성질에 의해, 각 고정된 $j$와 $(a,b)\in\{(0,0),(1,0),(0,1)\}$에 대해

<a id="equation-5-5"></a>

$$
|\mathcal B_{z-a,w-b}^{(j)}(\xi)|
\ll_{z,w,j}(1+|\xi|)^{M_j}
\tag{5.5}
$$

인 $M_j\ge0$가 존재한다. 경로 표현을 미분할 때는 $\xi$의 작은 근방을 감싸는 하나의 $Y$를 고정한다.

<a id="section-6"></a>

## 6. 산술 점근식의 증명

$H_s,H_t$는 §2와 같고, $C=C(r)$이다. 이 절에서

<a id="equation-6-1"></a>

$$
G_1=\gamma_E zC+H_s,\qquad G_2=\gamma_E wC+H_t-C
\tag{6.1}
$$

로 둔다. $G_1,G_2$는 부록 A.3의 국소 인수 $G$를 $(s_1,s_2)=(s-i\alpha,s)$에 대해 $(1,1)$에서 편미분한 값이다. $G_2$의 $-C$는 Perron 인자 $1/s$에서 온다.

<a id="section-6-1"></a>

### 6.1. 이동 균등 추정

**보조정리 1.** 고정된 $z,w$에 대해 어떤 $B\ge0$가 존재하여 $j=0,1$, $|\xi|\le L^{1/4}$에서 균등하게

<a id="equation-6-2"></a>

$$
\partial_\xi^jS_x(\xi)
=xL^{r-1}\,\partial_\xi^j\Big[C\mathcal B_{z,w}(\xi)+\frac1L\big(G_1\mathcal B_{z-1,w}(\xi)+G_2\mathcal B_{z,w-1}(\xi)\big)\Big]
+O_{z,w}\!\left(xL^{\beta-3}(1+|\xi|)^B\right).
\tag{6.2}
$$

증명은 부록 A에 있다.

<a id="section-6-2"></a>

### 6.2. 정규화 변수의 변경

$\Theta_x(\xi):=A_x^{z,w}(e^{i\xi u})$로 둔다. 각 항에 대해

$$
\frac{d}{dt}e^{i\xi\log d/\log t}
=-\frac{\xi}{t\log t}\,
\partial_\xi e^{i\xi\log d/\log t}
$$

를 $dm\le t\le x$에서 적분하고 합하면

<a id="equation-6-3"></a>

$$
\Theta_x(\xi)
=S_x(\xi)-1+
\xi\int_2^x\frac{\partial_\xi S_t(\xi)}{t\log t}\,dt.
\tag{6.3}
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

$|\xi|\le L^{1/8}$, $t\ge\sqrt x$이면 충분히 큰 $x$에서 $|\xi|\le(\log t)^{1/4}$이다. 이 구간에는 (6.2)을 $j=1$로 적용한다. 주항 $Ct(\log t)^{r-1}\mathcal B'(\xi)$의 기여는

$$
\xi\,C\mathcal B_{z,w}'(\xi)\int_{\sqrt x}^x(\log t)^{r-2}\,dt
=\xi\,C\mathcal B_{z,w}'(\xi)\,xL^{r-2}\bigl(1+O(L^{-1})\bigr)
$$

이고, 나머지는 (5.5)에 의해 $O(xL^{\beta-3}(1+|\xi|)^{M})$이다. $S_x(\xi)$에 (6.2)을 $j=0$으로 대입하면, 어떤 $M_*\ge0$에 대해 $|\xi|\le L^{1/8}$에서

<a id="equation-6-4"></a>

$$
\Theta_x(\xi)
=xL^{r-1}\Big[C\mathcal B_{z,w}+\frac1L\big(G_1\mathcal B_{z-1,w}+G_2\mathcal B_{z,w-1}+C\,\xi\mathcal B_{z,w}'\big)\Big](\xi)
+O_{z,w}\!\left(xL^{\beta-3}(1+|\xi|)^{M_*}\right).
\tag{6.4}
$$

<a id="section-6-3"></a>

### 6.3. 푸리에 적분

$F$의 연장 $\widetilde F\in C_c^N(\mathbb R)$를 $\|\widetilde F\|_{C^N}\ll_N\|F\|_{C^N}$이 되도록 고정한다. 푸리에 규약은

$$
\widetilde F(u)=\frac1{2\pi}
\int_{\mathbb R}\widehat{\widetilde F}(\xi)e^{i\xi u}\,d\xi
$$

이며 $|\widehat{\widetilde F}(\xi)|\ll_N\|F\|_{C^N}(1+|\xi|)^{-N}$이다. $N>M_*+1$로 잡고 $R=L^{1/8}$로 나누면 저주파의 오차 적분은 (6.4)에 의해 $O_{F,z,w}(xL^{\beta-3})$이다. 고주파에서는

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

(4.4)를 쓰면 $G_2\langle\mathcal K_{z,w-1},F\rangle+C\langle\mathcal K_{z,w},uF'\rangle=(G_2+C)\langle\mathcal K_{z,w-1},F\rangle-(r-1)C\langle\mathcal K_{z,w},F\rangle$이고, $G_2+C=\gamma_E wC+H_t$이므로 (2.2)–(2.3)이 따른다. 분포와 적분의 교환은 (5.5)와 위의 감소 추정으로 정당화된다. 오차 상수는 $\|F\|_{C^N}$에 선형으로 의존한다.

**양의 지수와 지시함수.** $z,w>0$, $0<v<1$이면 가중치가 양수이므로 $\mathbf1_{(v,1]}$를 폭 $\varepsilon$의 매끄러운 함수로 위아래에서 끼울 수 있다. 정규화한 오차는 $O(\varepsilon+L^{-1}\varepsilon^{-N})$이므로, $\delta=1/(N+1)$에 대해
$$
A_x^{z,w}(\mathbf1_{(v,1]})
=C(r)xL^{r-1}\int_v^1\frac{u^{z-1}(1-u)^{w-1}}{\Gamma(z)\Gamma(w)}\,du
+O_{z,w,v}(xL^{r-1-\delta}).
$$
부호가 섞이는 §7의 두 경우에는 이 양의 가중치 논증을 적용할 수 없다.

<a id="section-7"></a>

## 7. 특수핵과 끝점 항의 확인

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

<a id="equation-7-1"></a>

$$
\mathcal K_{-1/2,1/2}=h_1',
\tag{7.1}
$$

<a id="equation-7-2"></a>

$$
\langle\mathcal K_{-1/2,1/2},F\rangle
=-\frac1\pi\int_0^1\sqrt{\frac{1-u}{u}}\,F'(u)\,du.
\tag{7.2}
$$

**두 지수가 음수인 경우**

<a id="equation-7-3"></a>

$$
\mathcal K_{-1/2,-1/2}=h_2',
\tag{7.3}
$$

<a id="equation-7-4"></a>

$$
\langle\mathcal K_{-1/2,-1/2},F\rangle
=-\frac1\pi\int_0^1
\frac{u-1/2}{\sqrt{u(1-u)}}\,F'(u)\,du.
\tag{7.4}
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

**증명.** (4.3)은 분포 항등식이므로

$$
\begin{aligned}
[(1-u)a]'&=\mathcal K_{-1/2,1/2},\\
[(u-\tfrac12)a]'&=-\tfrac12\bigl(\mathcal K_{-1/2,1/2}
+\mathcal K_{1/2,-1/2}\bigr)=\mathcal K_{-1/2,-1/2}.
\end{aligned}
$$

따라서 끝점 항을 포함해 (7.1), (7.3)이 성립한다. $\square$

| $(z,w)$ | $\langle\mathcal K,1\rangle$ | $\langle\mathcal K,u\rangle$ | $\langle\mathcal K,u^2\rangle$ | $\langle\mathcal K,u^3\rangle$ |
|:--|--:|--:|--:|--:|
| $(1/2,1/2)$ | $1$ | $1/2$ | $3/8$ | $5/16$ |
| $(-1/2,1/2)$ | $0$ | $-1/2$ | $-1/4$ | $-3/16$ |
| $(-1/2,-1/2)$ | $0$ | $0$ | $-1/4$ | $-3/8$ |

세 경우의 산술 상수는 각각 $1/\zeta(2),1,1$이고 기준 규모는 $x,x/L,x/L^2$이다. 마지막 경우의 산술적 의미는 [뫼비우스 경우 §4](mobius-cases.md#section-4)에 있다.

**꼬리값.** $0<v<1$이고 $h$가 $[0,1]$에 지지된 적분 가능 함수이며 $v$ 근방에서 매끄러우면
$\mathbf1_{(v,\infty)}h'=(\mathbf1_{(v,\infty)}h)'-h(v)\delta_v$를 만족한다.
따라서 $T_{z,w}(v):=\langle\mathbf1_{(v,\infty)}\mathcal K_{z,w},\chi\rangle$라 두면
($\chi\in C_c^\infty(\mathbb R)$는 $[0,1]$ 근방에서 $1$),

<a id="equation-7-5"></a>

$$
T_{-1/2,1/2}(v)=-\frac1\pi\sqrt{\frac{1-v}{v}},\qquad
T_{-1/2,-1/2}(v)=\frac{1-2v}{2\pi\sqrt{v(1-v)}}.
\tag{7.5}
$$

이는 분포의 꼬리값이다. 산술 꼬리합의 점근식을 얻으려면 $u=v$ 근방의 합을 추가로 제어해야 하며, 정리 1의 매끄러운 시험함수에 대한 수렴만으로는 충분하지 않다.

<a id="section-8"></a>

## 8. 약수별 이산핵으로의 확장

베타핵 $\mathcal K_{z,w}$는 평균의 연속적인 모양을 기술한다. 이제 각 정수의 약수에 놓인 원자들을 따로 기록한다. 새 모형에는 $\widetilde A,\widetilde D,\widetilde S$를 사용하여 §2–7의 제곱 인수 없는 모형과 구분한다.

<a id="section-8-1"></a>

### 8.1. 일반화된 약수계수와 이산핵

$d_z(n)$을 $\zeta(s)^z$의 디리클레 계수로 정의한다. $\Re s>1$에서 오일러 곱의 로그를 사용하면

<a id="equation-8-1"></a>

$$
\zeta(s)^z=\sum_{n\ge1}\frac{d_z(n)}{n^s},
\qquad
d_z(p^a)=\frac{(z)_a}{a!},
\qquad d_z(1)=1.
\tag{8.1}
$$

$d_z$는 곱셈적이고, 고정된 $n$에서는 $z$의 다항식이다. 특히 $d_0(n)=\mathbf1_{n=1}$, $d_1(n)=1$, $d_{-1}(n)=\mu(n)$이다.

**정의 7.1 (약수별 이산핵).** $n\ge2$에서 $[0,1]$ 위의 유한 복소측도를

<a id="equation-8-2"></a>

$$
\boxed{
\mathscr K_n^{z,w}
=\sum_{d\mid n}d_z(d)d_w(n/d)\,
\delta_{\log d/\log n}
}
\tag{8.2}
$$

로 둔다. 따라서 $\langle\mathscr K_n^{z,w},F\rangle$는 약수마다 $F$를 평가한 정확한 유한합이다. 아래 합산 변환의 $n=1$ 항은 별도로 $1$로 정의하며, $\log1/\log1$은 사용하지 않는다.

**명제 8.1 (총질량과 반사).** $R(u)=1-u$라 하면 모든 $z,w\in\mathbb C$와 $n\ge2$에 대해

<a id="equation-8-3"></a>

$$
\langle\mathscr K_n^{z,w},1\rangle=d_{z+w}(n),
\qquad
R_*\mathscr K_n^{z,w}=\mathscr K_n^{w,z}.
\tag{8.3}
$$

특히 $z=w$이면 $\langle\mathscr K_n^{z,z},u\rangle=d_{2z}(n)/2$이다.

**증명.** 절대수렴 영역에서 $\zeta(s)^z\zeta(s)^w=\zeta(s)^{z+w}$의 계수를 비교하면 $d_z*d_w=d_{z+w}$를 얻는다. 반사식은 약수 $d$를 $n/d$로 바꾸면 따르고, 대칭인 경우 $u+(1-u)=1$을 적용한다. $\square$

소수 거듭제곱에서는 원자가 명시적으로

$$
\mathscr K_{p^a}^{z,w}
=\sum_{j=0}^a
\frac{(z)_j(w)_{a-j}}{j!(a-j)!}\,\delta_{j/a}
\qquad(a\ge1)
$$

에 놓인다. 일반적인 $n=\prod_p p^{a_p}$에서는 각 $0\le j_p\le a_p$의 선택이 $d=\prod_p p^{j_p}$와 대응하며, 원자 위치는 $\sum_p j_p\log p/\log n$이다. 이산핵은 이 소수 거듭제곱별 선택을 모두 보존한다.

<a id="section-8-2"></a>

### 8.2. 기존 모형과의 관계 및 생성함수

제곱 인수 없는 $n$과 그 약수 $d$에 대해서는 $d_z(d)=z^{\omega(d)}$이다. 그러므로

<a id="equation-8-4"></a>

$$
A_x^{z,w}(F)
=\sum_{\substack{2\le n\le x\\\mu^2(n)=1}}
\langle\mathscr K_n^{z,w},F\rangle,
\qquad
\widetilde A_x^{z,w}(F)
:=\sum_{2\le n\le x}\langle\mathscr K_n^{z,w},F\rangle.
\tag{8.4}
$$

새 모형은 기존 가중치의 제곱 인수 없는 부분을 그대로 포함한다. 다만 $\mu^2(n)$ 제한을 없앨 때 소수 거듭제곱의 가중치는 (8.1)로 정한다. §3.1의 서로소 조건 해제 모형은 $d,m$ 각각에 제곱 인수 없는 제한을 유지하므로, 현재의 $d_z$ 모형과 다르다.

**명제 8.2 (정확한 두 변수 변환).** $\Re s,\Re t>1$에서

<a id="equation-8-5"></a>

$$
\boxed{
\widetilde D_{z,w}(s,t)
=\sum_{d,m\ge1}\frac{d_z(d)d_w(m)}{d^sm^t}
=\zeta(s)^z\zeta(t)^w,
\qquad
\widetilde D_{z,w}(s,s)=\zeta(s)^{z+w}.
}
\tag{8.5}
$$

**증명.** $(1-v)^{-z}=\sum_{a\ge0}(z)_a v^a/a!$를 각 소수에 적용한다. $|(z)_a|\le(|z|)_a$이므로 각 급수는 주어진 영역에서 절대수렴한다. 두 급수를 곱하면 첫 항등식을 얻고, 대각선에서는 오일러 곱으로 선택한 로그가 같으므로 두 번째 항등식을 얻는다. $\square$

즉 새 모형에서는 (3.1)의 $H$에 해당하는 해석적 인자가 $1$이다.

<a id="section-8-3"></a>

### 8.3. 이산핵 평균의 이차 베타 법칙

**정리 2 (이산핵 평균의 이차 베타 법칙).** 고정된 $z,w$에 대해 충분히 큰 유한 정수 $N=N(z,w)$가 존재하여, $F\in C^N([0,1])$이면 $x\to\infty$일 때

<a id="equation-8-6"></a>

$$
\boxed{
\widetilde A_x^{z,w}(F)
=xL^{r-1}
\left\langle\mathcal K_{z,w}+\frac1L\widetilde{\mathcal K}^{(1)}_{z,w},F\right\rangle
+O_{F,z,w}(xL^{\beta-3}),
}
\tag{8.6}
$$

<a id="equation-8-7"></a>

$$
\widetilde{\mathcal K}^{(1)}_{z,w}
=\gamma_E z\mathcal K_{z-1,w}
+\gamma_E w\mathcal K_{z,w-1}
-(r-1)\mathcal K_{z,w}.
\tag{8.7}
$$

여기서 $r=z+w$, $\beta=\Re r$, $L=\log x$이며, $\gamma_E$는 §2와 같은 오일러 상수이다.

**증명.** 정리 1의 증명에서 사용한 계수 지배와 국소 전개를 새 모형에 대해 확인한다. $k_0\ge\max(2,|z|+|w|)$인 정수를 택하면

<a id="equation-8-8"></a>

$$
\sum_{d\mid n}|d_z(d)d_w(n/d)|
\le(d_{|z|}*d_{|w|})(n)
=d_{|z|+|w|}(n)
\le\tau_{k_0}(n).
\tag{8.8}
$$

마지막 부등식은 소수 거듭제곱에서 $(|z|+|w|)_a\le(k_0)_a$를 적용한 것이다. 이 지배는 (A.1)의 $j=0,1$ 주파수 미분에도 적용된다.

§5와 같이

<a id="equation-8-9"></a>

$$
\begin{aligned}
\widetilde S_x^{z,w}(\xi)
&:=\sum_{dm\le x}d_z(d)d_w(m)e^{i\xi\log d/L},\\
\sum_{n\ge1}\frac{\sum_{d\mid n}d_z(d)d_w(n/d)e^{i\xi\log d/L}}{n^s}
&=\zeta(s-i\xi/L)^z\zeta(s)^w.
\end{aligned}
\tag{8.9}
$$

로 둔다. (8.8)에 의해 부록 A.1의 Perron 절단 오차가 그대로 성립하고, 부록 A.2에서는 $H\equiv1$을 사용한다. 부록 A.3의 국소 인자는

$$
\widetilde G(s_1,s_2)
=\frac{Z(s_1)^zZ(s_2)^w}{s_2},
\qquad Z(s)=(s-1)\zeta(s)
$$

이며, $(1,1)$에서 값은 $1$, 두 편미분은 각각 $\gamma_E z$, $\gamma_E w-1$이다. 따라서 어떤 $B\ge0$에 대해 $j=0,1$, $|\xi|\le L^{1/4}$에서

<a id="equation-8-10"></a>

$$
\partial_\xi^j\widetilde S_x^{z,w}(\xi)
=xL^{r-1}\partial_\xi^j
\left[
\mathcal B_{z,w}(\xi)
+\frac{\gamma_E z\mathcal B_{z-1,w}(\xi)+(\gamma_E w-1)\mathcal B_{z,w-1}(\xi)}L
\right]
+O_{z,w}\!\left(xL^{\beta-3}(1+|\xi|)^B\right).
\tag{8.10}
$$

분모 $\log x$를 원래 이산핵의 $\log n$으로 바꾸려면 §6.2의 정확한 항등식을 사용한다.

<a id="equation-8-11"></a>

$$
\widetilde\Theta_x^{z,w}(\xi)
:=\widetilde A_x^{z,w}(e^{i\xi u})
=\widetilde S_x^{z,w}(\xi)-1
+\xi\int_2^x
\frac{\partial_\xi\widetilde S_t^{z,w}(\xi)}{t\log t}\,dt.
\tag{8.11}
$$

여기서 $\widetilde S_t$의 정규화는 $\log t$이고, $-1$은 $n=1$ 항을 제거한다. (8.8)은 §6.2의 작은 $t$ 구간과 §6.3의 고주파 꼬리에도 같은 지배를 준다. 따라서 (8.10)을 넣고 $F$의 푸리에 역변환을 적용하면 보정항은

$$
\gamma_E z\langle\mathcal K_{z-1,w},F\rangle
+(\gamma_E w-1)\langle\mathcal K_{z,w-1},F\rangle
+\langle\mathcal K_{z,w},uF'\rangle
$$

이다. (4.4)로 마지막 항을 바꾸면 (8.7)을 얻는다. 필요한 유한 차수 $N$과 오차의 제어는 §6.3과 같다. $\square$

$J$차 판은 [§9의 정리 3](#section-9-1)이며, 이 모형에서는 계수가 Stieltjes 상수만으로 정해진다((9.7)).

특히 $z,w>0$이면 양의 측도들의 정규화한 합이

$$
\frac{\sum_{2\le n\le x}\mathscr K_n^{z,w}}
{\sum_{2\le n\le x}d_r(n)}
\Longrightarrow
\Gamma(r)\mathcal K_{z,w}(u)\,du
$$

로 약수렴한다. 매끄러운 시험함수에서의 수렴과 양성으로 연속함수까지 확장할 수 있다. 이 평균은 정수 $n$을 $d_r(n)$에 비례해 고르고, 그 안에서 $\mathscr K_n^{z,w}/d_r(n)$를 사용하는 가중 평균이다. $r=1$이면 $d_r(n)=1$이므로 정수를 균등하게 고르는 평균이 된다.

<a id="section-8-4"></a>

### 8.4. 양의 아크사인 핵과 음의 뫼비우스 핵

**양의 지수 $(1/2,1/2)$.** 각 $\mathscr K_n^{1/2,1/2}$는 총질량 $1$의 확률측도이다. 하지만 일반적인 정수에서는 DDT의 균등 약수측도 $\nu_n$과 다르다. 예를 들어

$$
\mathscr K_{p^2}^{1/2,1/2}
=\frac38\delta_0+\frac14\delta_{1/2}+\frac38\delta_1,
\qquad
\nu_{p^2}=\frac13(\delta_0+\delta_{1/2}+\delta_1).
$$

두 모형은 평균에서 같은 아크사인 극한을 갖는다. 새 모형에서는 정리 2와 §7의 핵 항등식으로 더 정확히

<a id="equation-8-12"></a>

$$
\widetilde A_x^{1/2,1/2}(F)
=x\left[
\int_0^1a(u)F(u)\,du
+\frac{\gamma_E}{L}\int_0^1\left(u-\frac12\right)a(u)F'(u)\,du
\right]
+O_F\!\left(\frac{x}{L^2}\right).
\tag{8.12}
$$

를 얻는다. 기존 (2.5)와 비교하면 주항의 산술 상수 $1/\zeta(2)$가 $1$로, 보정계수 $\kappa(1)$이 $\gamma_E$로 바뀐다. $F=u^2$에서는 $\widetilde A_x^{1/2,1/2}(u^2)/x=3/8+\gamma_E/(4L)+O(L^{-2})$이다.

**음의 지수 $(-1/2,-1/2)$.** 총질량은 $\mu(n)$이며, 소수 거듭제곱에서는 부호가 섞인 원자들이 상쇄된다. 예를 들어

$$
\mathscr K_{p^2}^{-1/2,-1/2}
=-\frac18\delta_0+\frac14\delta_{1/2}-\frac18\delta_1.
$$

이 측도는 총질량이 $0$이지만 영측도는 아니다. 명제 8.1에서 정확히

<a id="equation-8-13"></a>

$$
\widetilde A_x^{-1/2,-1/2}(1)=M(x)-1,
\qquad
\widetilde A_x^{-1/2,-1/2}(u)=\frac{M(x)-1}{2},
\qquad M(x)=\sum_{n\le x}\mu(n).
\tag{8.13}
$$

이다. 이 총질량과 일차 모멘트는 기존 제곱 인수 없는 모형의 [뫼비우스 경우 §4](mobius-cases.md#section-4)와 같다. 더 높은 모멘트에서는 제곱 인수가 있는 정수의 원자들이 영향을 줄 수 있다.

정리 2의 주항은 $xL^{-2}\langle\mathcal K_{-1/2,-1/2},F\rangle$이다. §7에서 이 핵은 $F=1,u$에 대해 $0$이므로 (8.13)과 양립한다. $F=u^2$에서는 $\langle\mathcal K_{-1/2,-1/2},u^2\rangle=-1/4$이므로

$$
\widetilde A_x^{-1/2,-1/2}(u^2)
=-\frac{x}{4L^2}+O\!\left(\frac{x}{L^3}\right).
$$

총질량의 뫼비우스 상쇄와 약수 위치에 따른 비영인 주항이 하나의 이산핵 안에 공존한다.

<a id="section-8-5"></a>

### 8.5. 두 로그 좌표와 정확한 Mellin–Fourier 변환

원래 이산핵의 좌표는 $\log d/\log n$이고, $\widetilde S_x$의 좌표는 $\log d/\log x$이다. $T_c(u)=cu$라 쓰면 후자의 유한 복소측도는 정확히

<a id="equation-8-14"></a>

$$
\mathscr Q_x^{z,w}
:=\delta_0+\sum_{2\le n\le x}
(T_{\log n/L})_*\mathscr K_n^{z,w},
\qquad
\widetilde S_x^{z,w}(\xi)
=\langle\mathscr Q_x^{z,w},e^{i\xi u}\rangle.
\tag{8.14}
$$

여기서 $\langle(T_c)_*\nu,F\rangle=\langle\nu,F(cu)\rangle$이다. (8.11)은 이 좌표 변경을 합 수준에서 되돌리는 항등식이며, 이차 보정항을 계산할 때도 필요하다.

**명제 8.3 (정확한 Perron 표현).** $x>1$이 정수가 아니고 $L=\log x$, $c>1$, $\xi\in\mathbb R$이면

<a id="equation-8-15"></a>

$$
\boxed{
\widetilde S_x^{z,w}(\xi)
=\lim_{T\to\infty}\frac1{2\pi i}
\int_{c-iT}^{c+iT}
\zeta(s-i\xi/L)^z\zeta(s)^w\,\frac{x^s}{s}\,ds.
}
\tag{8.15}
$$

두 복소 거듭제곱은 $\Re s>1$의 오일러 곱 로그로 정의한다.

**증명.** (8.8)로 디리클레 급수 (8.9)가 $\Re s>1$에서 절대수렴한다. $x$를 고정하고 $L$도 고정한 상태에서 유효 Perron 공식을 적용한다. 정수가 아닌 $x$에서는 절단오차의 각 항이 $T\to\infty$일 때 $0$으로 가며, 절대수렴하는 계수합으로 지배되므로 (8.15)를 얻는다. 정수 $x$에서는 Perron 적분이 끝점에 반가중치를 주므로 해당 항의 절반을 추가해야 한다. $\square$

음의 대칭 지수에서는 $\mathfrak M_x(\xi):=\widetilde S_x^{-1/2,-1/2}(\xi)$로 두면

<a id="equation-8-16"></a>

$$
\mathfrak M_x(0)=M(x),
\qquad
F_{\xi/L}(s)
=\exp\!\left[-\frac12\left\{\log\zeta(s-i\xi/L)+\log\zeta(s)\right\}\right].
\tag{8.16}
$$

따라서 $F_{\xi/L}(s)$를 $1/\sqrt{\zeta(s-i\xi/L)\zeta(s)}$로 쓸 때도 위 로그로 정한 가지를 뜻한다. $\xi=0$이면 정확히 $1/\zeta(s)$가 된다. §5의 제곱 인수 없는 모형 $\mathfrak D_{\xi,L}$와는 구분한다.

고정된 $x$에서 모든 실수 $\xi$의 값을 알면 푸리에 변환의 유일성으로 $\mathscr Q_x$가 결정된다. $x$도 변화시키면 $n\ge2$에서의 약한 좌우 극한의 차이가

$$
\lim_{\varepsilon\downarrow0}
\left(\mathscr Q_{n+\varepsilon}^{z,w}
-\mathscr Q_{n-\varepsilon}^{z,w}\right)
=\mathscr K_n^{z,w}
$$

이므로 정수별 이산핵을 구분할 수 있다. 기존 항들은 좌표가 연속적으로 움직이고, 새로 들어오는 $n$의 항에서 $\log n/\log x\to1$이기 때문이다.

[전역 공식](global-formula.md)은 이 정확한 변환을 Hadamard 영점 곱과 연결하고, 특이점별 윤곽 적분을 분리한다. [전역 공식 §2](global-formula.md#section-2)는 가우스 평활화된 합을, [전역 공식 §3](global-formula.md#section-3)은 영점 쪽 가우스 합산 규칙으로 원래 정수합과 그 도약을 복원하는 공식을 다룬다.

<a id="section-9"></a>

## 9. 임의 차수 전개

**2026-11-06 (오차항 전개)**

정리 1·2의 오차 $O(xL^{\beta-3})$에는 두 종류의 항이 들어 있다. 하나는 국소 인수 $G$의 테일러 전개를 2차에서 자른 나머지(부록 A.3)이다. 다른 하나는 §6.2에서 $\int_{\sqrt x}^x(\log t)^{r-2}dt=xL^{r-2}(1+O(L^{-1}))$로 근사하며 버린 항이다. 반면 Perron 절단과 왼쪽 윤곽선 이동의 손실은 $xe^{-c\sqrt L}$뿐이다(부록 A.1–A.2). 이 절에서는 앞의 두 단계를 모두 $J$차까지 유지한다.

<a id="section-9-1"></a>

### 9.1. J차 베타 법칙

부록 A.3의 $Z(s)=(s-1)\zeta(s)$를 쓰고, Perron 인자 $1/s$를 뺀 국소 인수와 로그 적분을

<a id="equation-9-1"></a>

$$
\widetilde G(s_1,s_2)=Z(s_1)^zZ(s_2)^wH_{z,w}(s_1,s_2)=\sum_{a,b\ge0}\widetilde G_{ab}(s_1-1)^a(s_2-1)^b,
\qquad
\mathrm{Li}_\alpha(x)=\int_2^x(\log t)^\alpha\,dt
\tag{9.1}
$$

로 둔다. 거듭제곱의 가지는 부록 A.3과 같다. $\widetilde G=s_2G$이므로 $\widetilde G_{ab}=G_{ab}+G_{a,b-1}$이고($G_{a,-1}:=0$), $\widetilde G_{00}=C(r)$이다. 각 $\widetilde G_{ab}$는 $z,w$의 전체함수이다. $\mathrm{Li}_{-1}(x)=\mathrm{li}(x)-\mathrm{li}(2)$이다.

**정리 3 ($J$차 베타 법칙).** 고정된 $z,w$와 정수 $J\ge1$에 대해 충분히 큰 유한 정수 $N=N_J(z,w)$가 존재하여, $F\in C^N([0,1])$이면

<a id="equation-9-2"></a>

$$
\boxed{
A_x^{z,w}(F)=\sum_{a+b<J}\widetilde G_{ab}\,\big\langle\mathcal K_{z-a,w-b},F\big\rangle\,\mathrm{Li}_{r-1-a-b}(x)
+O_{F,z,w,J}\!\left(xL^{\beta-1-J}\right).
}
\tag{9.2}
$$

$H\equiv1$로 둔 $\widetilde G$에 대해 이산핵 평균 $\widetilde A_x^{z,w}(F)$에도 같은 식이 성립한다(정리 2의 모형).

$\mathrm{Li}_\alpha$를 부분적분하면 $\mathrm{Li}_\alpha(x)=x\sum_{m<M}(-1)^m\alpha^{\underline m}L^{\alpha-m}+O_\alpha(xL^{\Re\alpha-M})$이다. 여기서 $\alpha^{\underline m}=\alpha(\alpha-1)\cdots(\alpha-m+1)$이다. 따라서 (9.2)는 거듭제곱 형태

<a id="equation-9-3"></a>

$$
A_x^{z,w}(F)=xL^{r-1}\sum_{j<J}L^{-j}\big\langle\mathcal K^{(j)}_{z,w},F\big\rangle+O(xL^{\beta-1-J}),
\qquad
\mathcal K^{(j)}_{z,w}=\sum_{i=0}^{j}(-1)^{j-i}(r-1-i)^{\underline{j-i}}\sum_{a+b=i}\widetilde G_{ab}\,\mathcal K_{z-a,w-b}
\tag{9.3}
$$

와 동치이다. $J=1$은 1차 법칙이다. $\widetilde G_{10}=\gamma_EzC+H_s$, $\widetilde G_{01}=\gamma_EwC+H_t$이므로 $J=2$는 (2.2)–(2.3)이고, $H\equiv1$이면 (8.6)–(8.7)이다.

두 형태의 오차 차수는 같다. 하지만 (9.2)는 $\mathrm{Li}_\alpha$의 발산하는 점근 급수를 펼치지 않으므로 수치적으로 훨씬 정확하다([부록 B.2](#appendix-b-2)). 극단적인 예로 $(z,w)=(-1,1)$, $F(u)=u$이면 $A_x(u)=-\pi(x)$이다. $i\ge1$이면 $\langle\mathcal K_{-1-a,1-b},u\rangle=(-1-a)/\Gamma(1-i)=0$이므로 (9.2)는 모든 $J$에서 $-\mathrm{Li}_{-1}(x)$, 즉 $\mathrm{li}(x)$를 쓰는 소수 정리이다. (9.3)은 그 점근 급수 $-\tfrac xL(1+\tfrac1L+\tfrac2{L^2}+\cdots)$이다.

<a id="section-9-2"></a>

### 9.2. 2차 계수

$\gamma_1$을 $\zeta(s)=\frac1{s-1}+\gamma_E-\gamma_1(s-1)+\cdots$의 Stieltjes 상수($\gamma_1=-0.0728158\ldots$)라 하고, (2.4)의 $\kappa$와 함께

<a id="equation-9-4"></a>

$$
\nu(r)=\sum_p\frac{(\log p)^2}{(p+r)^2},
\qquad
\lambda(r)=\sum_p(\log p)^2\left[\frac{r+1}{(p-1)(p+r)}+\frac1{(p-1)^2}\right],
\qquad
\theta=\gamma_1+\frac{\gamma_E^2}2+\frac\lambda2
\tag{9.4}
$$

로 둔다. 국소 인수를 두 번 미분하면 $\partial_s^2\log H(1,1)=-z\lambda-z^2\nu$, $\partial_s\partial_t\log H(1,1)=-zw\nu$이고, $\log Z(s)=\gamma_E(s-1)-(\gamma_1+\gamma_E^2/2)(s-1)^2+\cdots$이다. 따라서 $C(r)\ne0$에서

<a id="equation-9-5"></a>

$$
\widetilde G_{10}=zC\kappa,\quad \widetilde G_{01}=wC\kappa,
\qquad
\sum_{a+b=2}\widetilde G_{ab}\varepsilon_1^a\varepsilon_2^b
=C\left[\frac{\kappa^2-\nu}2\,(z\varepsilon_1+w\varepsilon_2)^2-\theta\,(z\varepsilon_1^2+w\varepsilon_2^2)\right]
\tag{9.5}
$$

이다. 이를 (9.3)에 넣으면 셋째 핵은

<a id="equation-9-6"></a>

$$
\mathcal K^{(2)}_{z,w}
=(r-1)(r-2)C\,\mathcal K_{z,w}
-(r-2)C\kappa\,\big(z\mathcal K_{z-1,w}+w\mathcal K_{z,w-1}\big)
+\widetilde G_{20}\mathcal K_{z-2,w}+\widetilde G_{11}\mathcal K_{z-1,w-1}+\widetilde G_{02}\mathcal K_{z,w-2}
\tag{9.6}
$$

이다. 이산핵 모형에서는 $\nu=\lambda=0$, $\kappa=\gamma_E$, $C=1$이며 계수가 정확히 분리된다.

<a id="equation-9-7"></a>

$$
\widetilde G_{ab}=c_a(z)c_b(w),\qquad Z(s)^z=\sum_{a\ge0}c_a(z)(s-1)^a,\qquad
c_1(z)=\gamma_Ez,\quad c_2(z)=\frac{z(z-1)}2\gamma_E^2-z\gamma_1.
\tag{9.7}
$$

고차 계수는 [검산 스크립트](../scripts/check_beta_expansion.py)가 $(1,1)$ 둘레의 2변수 Cauchy 적분으로 계산한다. 그 값은 (9.5)와 소수 합의 절단 오차 수준인 약 $10^{-6}$까지 일치한다. 수치값은 다음과 같다.

| $(z,w)$ | $C$ | $\kappa$ | $\nu$ | $\lambda$ | $\widetilde G_{20}$ | $\widetilde G_{11}$ | $\widetilde G_{02}$ |
|:--|--:|--:|--:|--:|--:|--:|--:|
| $(\tfrac12,\tfrac12)$ | 0.607927 | 1.71714 | 0.54190 | 2.99602 | −0.300961 | 0.365768 | −0.300961 |
| $(1,1)$ | 0.286747 | 1.97622 | 0.43936 | 3.48209 | −0.029182 | 0.993895 | −0.029182 |
| $(2,1)$ | 0.114884 | 2.17220 | 0.37508 | 3.88716 | 0.529852 | 0.997971 | 0.015433 |
| $(-1,1)$ | 1 | 1.33258 | 0.74160 | 2.36756 | 1.794643 | −1.034178 | −0.760465 |

<a id="section-9-3"></a>

### 9.3. 증명

**1단계 (J차 이동 균등 추정).** 부록 A.3의 $G=\widetilde G/s_2$는 $(1,1)$ 중심의 고정된 다원판에서 해석적이다. 경로 위에서는 $|q|/L+|\xi|/L\ll L^{-1/2}$이므로, $G_L^{[J]}=\sum_{a+b<J}G_{ab}\big(\tfrac{q-i\xi}L\big)^a\big(\tfrac qL\big)^b$에 대해 Cauchy 추정으로 (A.7)을 다음과 같이 바꿀 수 있다.

$$
\big|G_L-G_L^{[J]}\big|\ll_{z,w,J}\frac{(|q|+|\xi|)^J}{L^J},
\qquad
\big|\partial_\xi\big(G_L-G_L^{[J]}\big)\big|\ll_{z,w,J}\frac{(|q|+|\xi|)^{J-1}}{L^J}.
$$

(A.8)의 우변은 $(1+|\xi|)^{m+J+1}L^{-J}$가 된다. 각 단항식에 (5.4)를 $(z-a,w-b)$로 적용하면, $j=0,1$, $|\xi|\le L^{1/4}$에서

<a id="equation-9-8"></a>

$$
\partial_\xi^jS_x(\xi)=xL^{r-1}\partial_\xi^j\sum_{a+b<J}G_{ab}L^{-a-b}\mathcal B_{z-a,w-b}(\xi)
+O_{z,w,J}\!\left(xL^{\beta-1-J}(1+|\xi|)^{m+J+1}\right)
\tag{9.8}
$$

이다. 이는 [고차 전개 (1.2)](higher-order.md#equation-1-2)에 $j=1$의 경우를 더한 것이다.

**2단계 (정규화 변수의 변경을 버리지 않기).** (6.3)의 적분에서 $t\in[\sqrt x,x]$에 (9.8)을 적용한다. 이때 $\int_{\sqrt x}^x(\log t)^{r-2-a-b}dt$를 근사하지 않고 $\mathrm{Li}_{r-2-a-b}(x)$로 바꾸면 $O(\sqrt x\,L^{O(1)})$만 달라진다. 작은 구간은 §6.2와 같이 $|\xi|\sqrt xL^{k_0}$이다. 나머지 적분은 $|\xi|(1+|\xi|)^{m+J+1}xL^{\beta-2-J}$이다. 따라서 $\mathcal B_{ab}:=\mathcal B_{z-a,w-b}$라 쓰면 어떤 $M_J$에 대해 $|\xi|\le L^{1/8}$에서

<a id="equation-9-9"></a>

$$
\Theta_x(\xi)=\sum_{a+b<J}G_{ab}\Big[xL^{r-1-a-b}\mathcal B_{ab}(\xi)+\mathrm{Li}_{r-2-a-b}(x)\,\xi\mathcal B_{ab}'(\xi)\Big]
+O_{z,w,J}\!\left(xL^{\beta-1-J}(1+|\xi|)^{M_J}\right).
\tag{9.9}
$$

**3단계 (푸리에 적분).** §6.3과 같이 $R=L^{1/8}$로 나눈다. 꼬리 추정의 지수만 바뀌어

$$
N\ge\max\big(M_J+2,\ M_0^{(J)}+2+8J,\ 1+8(k_0+J-\beta)\big)
$$

이면 충분하다. 여기서 $M_0^{(J)}$는 $a+b<J$인 $\mathcal B_{z-a,w-b}$와 그 1차 주파수 미분을 지배한다. $\mathrm{Li}_\alpha(x)\ll_\alpha xL^{\Re\alpha}$이므로 주항의 꼬리도 같은 방식으로 처리된다. $\xi\mathcal B'_{ab}$의 역변환은 §6.3과 같이 $\langle\mathcal K_{z-a,w-b},uF'\rangle$이므로

<a id="equation-9-10"></a>

$$
A_x^{z,w}(F)=\sum_{a+b<J}G_{ab}\Big[xL^{r-1-a-b}\langle\mathcal K_{z-a,w-b},F\rangle+\mathrm{Li}_{r-2-a-b}(x)\langle\mathcal K_{z-a,w-b},uF'\rangle\Big]+O(xL^{\beta-1-J}).
\tag{9.10}
$$

**4단계 (정리).** $\alpha=r-2-a-b$라 하자. (4.4)에서 $\langle\mathcal K_{z-a,w-b},uF'\rangle=\langle\mathcal K_{z-a,w-b-1},F\rangle-(\alpha+1)\langle\mathcal K_{z-a,w-b},F\rangle$이다. $\frac{d}{dt}\,t(\log t)^{\alpha+1}=(\log t)^{\alpha+1}+(\alpha+1)(\log t)^\alpha$에서는 $xL^{\alpha+1}=\mathrm{Li}_{\alpha+1}(x)+(\alpha+1)\mathrm{Li}_\alpha(x)+O(1)$를 얻는다. 따라서 (9.10)의 괄호는

$$
\mathrm{Li}_{r-1-a-b}(x)\langle\mathcal K_{z-a,w-b},F\rangle+\mathrm{Li}_{r-2-a-b}(x)\langle\mathcal K_{z-a,w-b-1},F\rangle+O(1)
$$

이다. 둘째 항의 첨자를 $b\to b-1$로 옮기고 $\widetilde G_{ab}=G_{ab}+G_{a,b-1}$을 쓰면 (9.2)가 된다. 이때 남는 $a+b=J-1$의 항은 $\mathrm{Li}_{r-1-J}(x)\ll xL^{\beta-1-J}$이므로 오차에 들어간다. 이산핵 모형은 정리 2의 증명과 같이 (8.8)–(8.11)과 $H\equiv1$을 쓴다. $\square$

<a id="section-9-4"></a>

### 9.4. 남는 오차

- **모든 차수.** $F\in C^\infty([0,1])$이면 (9.2)는 모든 $J$에서 성립하므로 점근 전개이다. 그러나 수렴하지는 않는다. $\widetilde G$의 수렴 반경은 일반적으로 $\tfrac12$ 이하이고([고차 전개 §2](higher-order.md#section-2)), $\langle\mathcal K_{z-a,w-b},u^k\rangle=(z-a)_k/\Gamma(r-a-b+k)$는 계승으로 자란다. 고정된 $x$에서는 [고차 전개의 휴리스틱](higher-order.md#section-2)처럼 최적 차수가 있다.
- **유한한 경우.** $r$이 정수이고 $F$가 차수 $k$의 다항식이면 $i\ge r+k$에서 $1/\Gamma(r-i+k)=0$이다. 이때 (9.2)는 유한합이 되고 오차는 모든 $A$에 대해 $O_A(xL^{-A})$이다. 예를 들어 $\sum_{n\le x}\mu^2(n)3^{\omega(n)}$은 세 항 $\mathrm{Li}_2,\mathrm{Li}_1,\mathrm{Li}_0$로 끝난다.
- **바닥.** 어떤 $J$로도 지울 수 없는 부분은 영점의 진동이다. $(z,w)=(-1,1)$, $F=u$에서는 남는 오차가 정확히 $\mathrm{li}(x)-\mathrm{li}(2)-\pi(x)$이다. 부록 B.2의 다른 경우에도 $J=3$–$5$에서 상대오차가 $10^{-5}$–$10^{-8}$ 규모의 바닥에서 멈춘다. 이 부분은 [전역 공식](global-formula.md)의 영점 성분이 다룬다.

<a id="appendix-a"></a>

## 부록 A. 이동 균등 추정

<a id="appendix-a-1"></a>

### A.1. 계수 상계와 Perron 절단

§6의 $k_0$를 쓰고 $j=0,1$로 둔다. 계수는

<a id="equation-a-1"></a>

$$
|\partial_\xi^j b_{\xi,L}(n)|
\le\tau_{k_0}(n)\left(\frac{\log n}{L}\right)^j
\tag{A.1}
$$

로 지배된다. $j=0,n=1$의 마지막 인자는 $1$로 읽는다.

짧은 구간에서는 $X\ge2,\ 1\le h\le X$에 대해

<a id="equation-a-2"></a>

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

<a id="equation-a-3"></a>

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

<a id="equation-a-4"></a>

$$
E_j\ll_{z,w}\frac{x}{T}L^{k_0+1}
+x^{1-1/k_0}L^{k_0}
\ll_{z,w}xe^{-c_1\sqrt L}.
\tag{A.4}
$$

<a id="appendix-a-2"></a>

### A.2. 영점 없는 영역과 윤곽선

$\tau=|t|+4$로 두자. 고전적 영점 없는 영역의 상수 $a_0>0$를 충분히 작게 잡으면
$\sigma>1-a_0/(2\log\tau)$, $|t|\ge7/8$에서

<a id="equation-a-5"></a>

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

<a id="equation-a-6"></a>

$$
|\mathfrak D_0(s)|+|\mathfrak D_1(s)|\ll_{z,w}L^K.
\tag{A.6}
$$

왼쪽 직선과 위아래 연결선의 기여는 각각
$O(x^{1-\delta}L^K\log(T+5))$, $O(xL^K/T)$이다.
$\delta L\asymp\sqrt L$이므로 둘 다 $O_{z,w}(xe^{-c_2\sqrt L})$로 흡수된다. 이 선택은 절단 오차 $T^{-1}$와 왼쪽 이동의 이득 $e^{-\delta L}$을 같은 지수 규모로 맞춘 것이다.

<a id="appendix-a-3"></a>

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
$\log\zeta(s)=\log Z(s)-\log(s-1)$이므로, 이 국소 가지는 A.2에서 선택한 가지와 일치한다. $Z'(1)=\gamma_E$이므로

$$
\partial_{s_1}G(1,1)=\gamma_E zC+H_s=G_1,\qquad
\partial_{s_2}G(1,1)=\gamma_E wC+H_t-C=G_2.
$$

$G_L(q,\xi)=G(1+(q-i\xi)/L,\,1+q/L)$라 쓰면 2차 테일러 정리로

<a id="equation-a-7"></a>

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

<a id="equation-a-8"></a>

$$
\int_{\mathcal U_{Y,\delta L}}
\left|e^q\partial_\xi^j\{P_\xi(q)(G_L-G_L^{[1]})\}\right|\,|dq|
\ll_{z,w}\frac{(1+|\xi|)^{m+3}}{L^2}.
\tag{A.8}
$$

따라서 정수 $B\ge m+3$를 택하면

<a id="equation-a-9"></a>

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

이다. 앞의 $xL^{r-1}$을 곱해도 지수적으로 작아 (A.9)의 오차에 흡수된다. (5.4)와 그 미분판을 $(z,w)$, $(z-1,w)$, $(z,w-1)$에 적용하면

$$
\frac1{2\pi i}\int_{\mathcal U_Y}e^q\,\partial_\xi^j\Big[P_\xi(q)\Big(C+G_1\frac{q-i\xi}L+G_2\frac qL\Big)\Big]dq
=\partial_\xi^j\Big[C\mathcal B_{z,w}+\frac1L\big(G_1\mathcal B_{z-1,w}+G_2\mathcal B_{z,w-1}\big)\Big](\xi)
$$

이고, 이를 대입하면 (6.2)을 얻는다. 추정 과정에서 $\xi$로 나누지 않았으므로 후보의 합류와 상쇄도 포함한다. 테일러 전개를 $J$차까지 늘리면 같은 방법으로 [고차 전개 (1.2)](higher-order.md#equation-1-2)와 (9.8)을 얻는다.

<a id="appendix-b"></a>

## 부록 B. 수치 확인

B.1의 표는 초기 원고에 보고된 수치이다. [check_beta_expansion.py](../scripts/check_beta_expansion.py)를 `--x 1e7`로 실행하면 일곱 값과 2차 예측을 모두 다시 계산하며, 다섯 자리까지 일치한다(2026-11-06). [전역 공식의 검산 스크립트](global-formula.md#checks)는 별도의 항등식과 윤곽 적분을 확인한다.

**방법.** 체로 소인수를 구하고, 약수 모멘트 $\sum_{d\mid n}z^{\omega(d)}w^{\omega(n/d)}(\log d)^k$ ($k\le4$)를 소수마다 곱셈적으로 갱신해 정확히 합했다($x\le10^7$).

<a id="appendix-b-1"></a>

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

남은 오차에 $L^2$을 곱하면 $10^6$과 $10^7$에서 거의 같다: $(\tfrac12,\tfrac12),u^2$은 −0.550→−0.535, $(2,1),1$은 7.91→7.92, $(-1,1),u$는 −2.31→−2.37. 이 계산은 예측한 $L^{-2}$ 오차와 부합한다. 이 값들은 (9.6)의 $\langle\mathcal K^{(2)},F\rangle/C$, 즉 −0.464, 7.917, −2에 $O(L^{-1})$ 보정을 더한 것이다.

<a id="appendix-b-2"></a>

### B.2. J차 전개 (정리 3)

$x=10^7$에서 상대오차 $(\text{실측}-\text{예측})/(|C|xL^{r-1})$이다. "정리 1·2"는 (2.2) 또는 (8.6), 즉 (9.3)의 $J=2$이다. 위의 일곱 행은 제곱 인수 없는 모형(정리 1)이고, 아래 세 행은 이산핵 모형(정리 2)이다.

| $(z,w)$ | $F$ | 정리 1·2 | (9.3), $J=6$ | (9.2), $J=3$ | (9.2), $J=4$ | (9.2), $J=6$ |
|:--|:--|--:|--:|--:|--:|--:|
| $(\tfrac12,\tfrac12)$ | $u^2$ | $-2.1\cdot10^{-3}$ | $-1.8\cdot10^{-6}$ | $4.0\cdot10^{-6}$ | $4.0\cdot10^{-6}$ | $4.0\cdot10^{-6}$ |
| $(\tfrac12,\tfrac12)$ | $(u-\tfrac12)^4$ | $-1.7\cdot10^{-3}$ | $-4.2\cdot10^{-6}$ | $7.3\cdot10^{-4}$ | $-1.8\cdot10^{-4}$ | $3.7\cdot10^{-6}$ |
| $(1,1)$ | $u^2$ | $-8.5\cdot10^{-4}$ | $-1.2\cdot10^{-7}$ | $-1.0\cdot10^{-5}$ | $4.1\cdot10^{-7}$ | $4.1\cdot10^{-7}$ |
| $(2,1)$ | $1$ | $3.0\cdot10^{-2}$ | $-8.3\cdot10^{-6}$ | $-8.3\cdot10^{-6}$ | $-8.3\cdot10^{-6}$ | $-8.3\cdot10^{-6}$ |
| $(2,1)$ | $u^2$ | $2.3\cdot10^{-2}$ | $-3.9\cdot10^{-6}$ | $-3.3\cdot10^{-4}$ | $2.1\cdot10^{-5}$ | $-3.9\cdot10^{-6}$ |
| $(-1,1)$ | $u$ | $-9.1\cdot10^{-3}$ | $4.6\cdot10^{-4}$ | $5.5\cdot10^{-4}$ | $5.5\cdot10^{-4}$ | $5.5\cdot10^{-4}$ |
| $(-1,1)$ | $u^2$ | $-2.5\cdot10^{-2}$ | $6.8\cdot10^{-4}$ | $9.0\cdot10^{-4}$ | $9.0\cdot10^{-4}$ | $9.0\cdot10^{-4}$ |
| $(\tfrac12,\tfrac12)$ | $u^2$ | $5.1\cdot10^{-4}$ | $1.2\cdot10^{-6}$ | $-7.0\cdot10^{-8}$ | $-7.0\cdot10^{-8}$ | $-7.0\cdot10^{-8}$ |
| $(\tfrac12,\tfrac12)$ | $(u-\tfrac12)^4$ | $1.7\cdot10^{-4}$ | $3.9\cdot10^{-7}$ | $-1.4\cdot10^{-6}$ | $-8.9\cdot10^{-9}$ | $-4.1\cdot10^{-8}$ |
| $(2,1)$ | $u^2$ | $1.7\cdot10^{-3}$ | $8.4\cdot10^{-7}$ | $3.6\cdot10^{-5}$ | $5.2\cdot10^{-7}$ | $7.8\cdot10^{-7}$ |

- 거듭제곱 형태 (9.3)에서는 한 차수마다 오차가 대략 $1/5$–$1/10$로 줄어든다. $(-1,1)$을 제외하면 $J=6$의 오차는 정리 1·2보다 약 $500$–$7000$배 작다.
- (9.2)는 더 낮은 $J$에서 바닥에 닿는다. 바닥에서는 $J$를 올려도 값이 변하지 않는다. $(2,1),F=1$은 §9.4의 유한한 경우이므로 $J=3$에서 끝난다.
- $(-1,1)$의 바닥은 소수 계량의 진동이다. $F=u$에서 $5.5\cdot10^{-4}=(\mathrm{li}(10^7)-\mathrm{li}(2)-\pi(10^7))/(10^7/L)$이다($\pi(10^7)=664579$). 이 행에서 (9.3)이 더 작아 보이는 것은 점근 급수가 진동항과 우연히 상쇄된 결과이다.
- (9.2)의 $J=2$는 정리 1·2와 같은 차수이지만 상수가 다르며, 이 표의 몇 경우에는 오히려 크다. 개선은 $J\ge3$에서 나온다.

<a id="references"></a>

## 참고 문헌

- J.-M. Deshouillers, F. Dress, G. Tenenbaum, *Lois de répartition des diviseurs, I*, Acta Arith. 34 (1979), 273–285: DDT의 평균 아크사인 법칙.
- G. Bareikis, E. Manstavičius, [On the DDT theorem](https://web.vu.lt/mif/e.manstavicius/wp-content/uploads/2016/08/EM-5.pdf), Acta Arith. 126.2 (2007), 155–168: §1의 DDT 정리 서술과 §2의 베타분포 확장.
- H. L. Montgomery and R. C. Vaughan, [Multiplicative Number Theory I, §6.1](https://personal.science.psu.edu/rcv4/personal/Publications/MNTI/10.0_pp_168_198_The_Prime_Number_Theorem.pdf), 정리 6.6–6.7: 영점 없는 영역과 복소 로그 추정.
- G. Tenenbaum, *Introduction to Analytic and Probabilistic Number Theory*, II.5: Selberg–Delange 방법.
- NIST DLMF, [§5.9](https://dlmf.nist.gov/5.9): Hankel 적분. [§13.6](https://dlmf.nist.gov/13.6): Bessel 함수. [§25.16](https://dlmf.nist.gov/25.16): $\psi$의 명시적 공식.
