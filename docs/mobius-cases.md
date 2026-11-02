# 뫼비우스 가중치와 끝점 항

제곱 인수 없는 정수의 가중 약수합에서 음의 정수·반정수 지수가 만드는 소멸과 끝점 항을 다룬다. 마지막 절은 뫼비우스 합성곱의 거듭제곱 오차와 영점 없는 영역을 연결한다. RH는 리만 가설이며, 사용하는 곳에 따로 표시한다.

- [§1. 정수 지수](#section-1)
- [§2. 두 뫼비우스 가중치](#section-2)
- [§3. 한쪽만 뫼비우스인 경우](#section-3)
- [§4. 음의 반정수 지수](#section-4)
- [§5. 거듭제곱 오차와 영점 없는 영역](#section-5)
- [부록 A. 수치 기록](#appendix-a)

<a id="notation"></a>

**기호와 전제.** $\mu,\omega,\tau,\Lambda$는 각각 뫼비우스 함수, 서로 다른 소인수의 개수, 약수의 개수, 폰 망골트 함수이다. $r=z+w$, $L=\log x$, $\gamma_E$는 오일러 상수이며

$$
A_x^{z,w}(F)=\sum_{\substack{2\le n\le x\\\mu^2(n)=1}}
\sum_{d\mid n}z^{\omega(d)}w^{\omega(n/d)}
F\!\left(\frac{\log d}{\log n}\right),\qquad 0^0=1
$$

로 둔다. $F$는 충분히 매끄러운 시험함수이다. $\Re s,\Re t>1$에서 기본 변환은

$$
D(s,t)=\prod_p(1+zp^{-s}+wp^{-t})
=\zeta(s)^z\zeta(t)^wH(s,t),\qquad C(r)=H(1,1)
$$

이며, $\mathcal K_{z,w}(u)=g_z(u)g_w(1-u)$, $g_z(v)=v_+^{z-1}/\Gamma(z)$는 전체해석적으로 연장한 분포이다. 이차 점근식과 오일러 곱의 성질은 [주정리 §2–4](main-theorem.md#section-2)를 사용한다. $A^\natural$은 $d,m$ 각각의 제곱 인수 없는 제한을 유지하면서 서로소 조건만 푼 합이다.

증명에서 $S_x(\xi)$는 위 합에서 $F(u)=e^{i\xi u}$로 놓은 뒤 분모 $\log n$을 $L$로 바꾸고 $n=1$ 항을 더한 지수합, $\Theta_x(\xi)=A_x^{z,w}(e^{i\xi u})$이다. $\alpha=\xi/L$이며, $\mathfrak D_{\xi,L}(s)=D(s-i\alpha,s)$이다. $k_0$와 절단 경로의 선택은 [주정리 §6 및 부록 A](main-theorem.md#section-6)를 따른다.

<a id="section-1"></a>

## 1. 정수 지수

**명제 1.** $z,w\in\{0,-1,-2,\dots\}$, $F\in C^\infty([0,1])$이면 모든 $k$에 대해 $A_x^{z,w}(F)=O_{F,z,w,k}(xL^{-k})$이다. $F$가 다항식이면 $O(xe^{-c\sqrt L})$이다.

**증명.** 지수가 비양의 정수이므로 $\mathfrak D_{\xi,L}(s)=\zeta(s-i\alpha)^z\zeta(s)^wH$는 [주정리 부록 A.2](main-theorem.md#appendix-a-2)의 영점 없는 영역에서 해석적이다. [주정리 (A.5)](main-theorem.md#equation-a-5)–[주정리 (A.6)](main-theorem.md#equation-a-6)의 추정 아래 Perron 경로를 우회 없이 $\Re s=1-\delta$로 옮기면, $|\xi|\le L^{1/4}$, $j=0,1$에서 $\partial_\xi^jS_x(\xi)\ll xe^{-c\sqrt L}$이다.

[주정리 (6.3)](main-theorem.md#equation-6-3)는 $|\xi|\le L^{1/8}$에서 $\Theta_x(\xi)\ll xe^{-c\sqrt L}(1+|\xi|)+|\xi|\sqrt xL^{k_0}$를 준다. [주정리 §6.3](main-theorem.md#section-6-3)의 고주파 꼬리는 $xL^{k_0-1-(N-1)/8}$이고 $F\in C^\infty$이면 $N$을 임의로 크게 잡을 수 있다. 따라서 모든 로그 거듭제곱의 절약을 얻는다.

$F=u^k$에서는 로그 가중합의 디리클레 급수 $(-\partial_{s_1})^kD(s_1,s)|_{s_1=s}$도 같은 영역에서 해석적이다. 여기에 $(\log t)^{-k}$의 부분합산을 적용하면 $O(xe^{-c\sqrt L})$이다. $\square$

<a id="section-2"></a>

## 2. $(-1,-1)$: 뫼비우스의 곱

제곱 인수 없는 $n=dm$에서 $d,m$은 서로소이므로 쌍 가중치는 $(-1)^{\omega(d)}(-1)^{\omega(m)}=\mu(d)\mu(m)=\mu(n)$이고

$$
A_x^{-1,-1}(F)=\sum_{2\le n\le x}\mu(n)\sum_{d\mid n}F\!\left(\frac{\log d}{\log n}\right).
$$

$F=1$이면 $\sum_{2\le n\le x}\mu(n)\tau(n)$이다. 주항은 서로 독립인 두 이유로 사라진다.

- **산술:** $C(-2)=0$이다. $p=2$의 국소 인수 $1+r/p=1-2/2$가 $0$이기 때문이다.
- **해석:** $\mathcal K_{-1,-1}=0$이다. 지수가 정수이기 때문이다.

$D(s,s)=\prod_p(1-2p^{-s})=\zeta(s)^{-2}H(s,s)$는 $s=1$에서 3중 영점을 가진다($\zeta^{-2}$에서 2, $p=2$에서 1).

서로소 조건을 풀면([주정리 §3.1](main-theorem.md#section-3-1)) 쌍 가중치는 $\mu(d)\mu(m)$이고, 디리클레 급수는 정확히 $1/(\zeta(s)\zeta(t))$, 즉 $H^\natural\equiv1$이다. 상수는 $C(-1)^2=1$이 되지만 핵은 여전히 $0$이다. 따라서 서로소 조건을 풀어도 명제 1의 소멸 추정은 유지된다.

<a id="section-3"></a>

## 3. $(-1,1)$: 한쪽만 뫼비우스

쌍 가중치는 $\mu(d)$이다. $C(0)=1$이고, $g_1(1-u)=\mathbf 1_{u<1}$이 $u=0$ 근처에서 $1$이므로 $\mathcal K_{-1,1}=\delta_0'$, $\mathcal K_{-2,1}=\delta_0''$, $\mathcal K_{-1,0}=0$이다. [주정리 (2.4)](main-theorem.md#equation-2-4)에서

<a id="equation-3-1"></a>

$$
A_x^{-1,1}(F)=-\frac xL\Big[F'(0)+\frac{F'(0)+\kappa(0)F''(0)}{L}\Big]+O\!\left(\frac x{L^3}\right),
\qquad \kappa(0)=\gamma_E+\sum_p\frac{\log p}{p(p-1)}=1.3325822\ldots
\tag{3.1}
$$

뫼비우스 반전은 이 차수까지 $F$의 원점에서의 1차, 2차 테일러 계수만 남긴다.

**$F=u$.** $\sum_{d\mid n}\mu(d)\log d=-\Lambda(n)$이므로 좌변은 정확히 $-\pi(x)$이다. (3.1)은 $\pi(x)=x/L+x/L^2+O(x/L^3)$, 즉 $\mathrm{li}(x)$의 앞 두 항이다.

**$F=u^2$.** 제곱 인수 없는 $n$에 대해 $\sum_{d\mid n}\mu(d)(\log d/\log n)^2$는 $\omega(n)\le2$일 때만 영이 아니고, $n=p$에서 $-1$, $n=pq$에서 $2ab$이다($a=\log p/\log n$, $b=\log q/\log n$). 따라서

<a id="equation-3-2"></a>

$$
\sum_{\substack{p<q\\pq\le x}}\frac{2\log p\log q}{\log^2(pq)}=\pi(x)-\frac{2\kappa(0)\,x}{L^2}+O\!\left(\frac x{L^3}\right).
\tag{3.2}
$$

**서로소 해제 변형.** 쌍 가중치 $\mu(d)\mu^2(m)$이면 [주정리 (3.6)](main-theorem.md#equation-3-6)에서 상수는 $C(-1)C(1)=6/\pi^2$이고, $H_{-1}\equiv1$이므로 $\eta(-1)=0$, $\kappa(-1)=\gamma_E$이다.

<a id="equation-3-3"></a>

$$
A_x^{\natural,-1,1}(F)=-\frac{6}{\pi^2}\frac xL\Big[F'(0)+\frac{F'(0)+\gamma_E F''(0)}{L}\Big]+O\!\left(\frac x{L^3}\right).
\tag{3.3}
$$

$6/\pi^2$은 $\mu*\mu^2$가 제곱수에만 놓인다는 사실($\sum(\mu*\mu^2)(n)n^{-s}=1/\zeta(2s)$)에서 온다. **[수치]** [부록 A](#appendix-a).

<a id="section-4"></a>

## 4. $(-\tfrac12,-\tfrac12)$

제곱 인수 없는 $n=dm$에서 쌍 가중치는 $\mu(n)/\tau(n)$이다. 따라서 $F=1,u$에 대한 합은 정확히 $M(x)-1$, $(M(x)-1)/2$이며, [주정리 §7](main-theorem.md#section-7) 표의 영인 주계수와 일치한다. 여기서 $M(x)=\sum_{n\le x}\mu(n)$이다.

<a id="section-5"></a>

## 5. 거듭제곱 오차와 영점 없는 영역

**명제 2.** 모든 $\varepsilon>0$에 대해 $\sum_{n\le x}(\mu*\mu)(n)\ll x^{\theta+\varepsilon}$이면 $\sigma>\theta$에서 $\zeta(s)\ne0$이다. 역으로 RH이면 $\sum_{n\le x}(\mu*\mu)(n)\ll x^{1/2+\varepsilon}$이다.

**증명.** 부분합산으로 $1/\zeta(s)^2=\sum(\mu*\mu)(n)n^{-s}$가 $\sigma>\theta$에서 수렴하고 해석적이다. 그 영역의 영점은 극이 되므로 있을 수 없다. 역은 RH 아래 $\sigma\ge\tfrac12+\varepsilon$에서 $1/\zeta(s)\ll|t|^\varepsilon$인 것과 Perron 공식으로 따른다. $\square$

따라서 뫼비우스 경우의 거듭제곱 오차 추정은 영점 없는 영역과 직결된다.

<a id="appendix-a"></a>

## 부록 A. 뫼비우스 가중합의 수치 기록

아래 표와 실험 설명은 초기 원고의 보고값이다. 계산 코드와 원자료가 없어 재현하지 않았으며, [전역 공식의 검산](global-formula.md#checks)과는 검증 범위가 다르다.

**보고된 방법.** 제곱 인수 없는 수와 $\omega$를 체로 구하고 약수 쌍을 직접 합했다($x\le10^7$).

**서로소 해제 $(-1,1)$, (3.3):** $A^\natural L/x$의 값이다($10^7$).

| $F$ | 실측 | 1차 $-\tfrac6{\pi^2}F'(0)$ | 2차 |
|:--|--:|--:|--:|
| $u$ | −0.6514 | −0.6079 | −0.6456 |
| $u^2$ | −0.0501 | 0 | −0.0435 |
| $e^u$ | −0.6772 | −0.6079 | −0.6674 |

남은 오차에 $L^2$을 곱한 값은 $10^6$과 $10^7$에서 각각 −1.52→−1.50, −1.70→−1.71, −2.53→−2.54이다.

<a id="references"></a>

## 참고 문헌

- H. L. Montgomery and R. C. Vaughan, [Multiplicative Number Theory I, §6.1](https://personal.science.psu.edu/rcv4/personal/Publications/MNTI/10.0_pp_168_198_The_Prime_Number_Theorem.pdf), 정리 6.6–6.7: 영점 없는 영역과 복소 로그 추정.
- G. Tenenbaum, *Introduction to Analytic and Probabilistic Number Theory*, II.5: Selberg–Delange 방법.
