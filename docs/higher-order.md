# 고차 베타 전개와 국소 특이점

이차 베타 법칙을 고정된 임의 차수로 확장하고, 전개 계수에서 $s=1/2$ 특이점과 모멘트 상쇄가 어떻게 드러나는지 살핀다. **[국소 계산]**은 고정된 특이점의 기여이며, **[휴리스틱]**은 추가 추정이 필요한 논의이다. RH는 리만 가설이다.

- [§1. 고정 차수의 베타 전개](#section-1)
- [§2. 수렴 반경과 모멘트 상쇄](#section-2)
- [§3. 아크사인의 세 점 핵](#section-3)
- [§4. 국소 연장과 전역 오차](#section-4)
- [부록 A. 전개 계수의 수치 기록](#appendix-a)

<a id="notation"></a>

**기호와 전제.** $z,w\in\mathbb C$를 고정하고 $r=z+w$, $\beta=\Re r$, $L=\log x$로 둔다. $\mu$는 뫼비우스 함수, $\omega$는 서로 다른 소인수의 개수이며,

$$
S_x(\xi)=\sum_{\substack{dm\le x\\\mu^2(dm)=1}}
z^{\omega(d)}w^{\omega(m)}e^{i\xi\log d/L},
\qquad
\mathcal B_{z,w}(\xi)=\langle\mathcal K_{z,w},e^{i\xi u}\rangle.
$$

$\mathcal K_{z,w}(u)=g_z(u)g_w(1-u)$는 $g_z(v)=v_+^{z-1}/\Gamma(z)$를 전체해석적으로 연장한 분포로 정의한다. 산술 인자는

$$
H(s,t)=\prod_p(1+zp^{-s}+wp^{-t})(1-p^{-s})^z(1-p^{-t})^w,
\qquad C=H(1,1)
$$

이다. $\widetilde H$는 이 곱의 이차 소수 인수를 분리한 잔여 인자이고, $e_1=z(z+1)/2$, $e_2=w(w+1)/2$이다([주정리 §3.2](main-theorem.md#section-3-2)). 이동 균등 추정과 Hankel 경로 $\mathcal U_Y$는 [주정리 §5–6](main-theorem.md#section-5)의 정의와 증명을 사용한다.

<a id="section-1"></a>

## 1. 고정 차수의 베타 전개

[주정리 부록 A.3](main-theorem.md#appendix-a-3)의 $Z(s)=(s-1)\zeta(s)$와 국소 인수 $G(s_1,s_2)=Z(s_1)^zZ(s_2)^wH(s_1,s_2)/s_2$를 사용한다. $G$를 전개하여

<a id="equation-1-1"></a>

$$
G(s_1,s_2)=\sum_{a,b\ge0}G_{ab}(s_1-1)^a(s_2-1)^b
\tag{1.1}
$$

로 쓴다. $G_{00}=C$, $G_{10}=\partial_1G(1,1)$, $G_{01}=\partial_2G(1,1)$이다.

**명제 1.** 고정된 정수 $J\ge1$과 $|\xi|\le L^{1/4}$에서

<a id="equation-1-2"></a>

$$
S_x(\xi)=xL^{r-1}\sum_{a+b<J}G_{ab}L^{-a-b}\mathcal B_{z-a,w-b}(\xi)+O\!\left(xL^{\beta-1-J}(1+|\xi|)^{B_J}\right).
\tag{1.2}
$$

여기서 $B_J\ge0$와 오차 상수는 고정된 $z,w,J$에 의존한다. $J=2$가 [주정리 (6.2)](main-theorem.md#equation-6-2)이고, 증명은 [주정리 부록 A.3](main-theorem.md#appendix-a-3)의 테일러 전개를 $J$차까지 늘린 것이다. 베타핵 전개의 계수는 $G$의 테일러 계수에서 나오지만, 이 식 자체는 고정 차수의 점근식이다.

$j$차 핵들의 $k$차 모멘트를 함께 모으면 $\sum_{a+b=j}G_{ab}\langle\mathcal K_{z-a,w-b},u^k\rangle=m_{j,k}/\Gamma(r-j+k)$이고

<a id="equation-1-3"></a>

$$
m_{j,k}=\sum_{a+b=j}G_{ab}\,(z-a)_k.
\tag{1.3}
$$

$m_{j,0}$은 $G(1+t,1+t)$의 테일러 계수이다. $m_{j,1}$, $m_{j,2}$에는 대각선 위의 $\partial_1G$, $\partial_1^2G$가 들어간다.

<a id="section-2"></a>

## 2. 국소 인자 $G$의 수렴 반경과 모멘트 상쇄

[주정리 (3.7)](main-theorem.md#equation-3-7)에서 가장 가까운 실축 쪽 특이점 후보는 $s_1=\tfrac12$, $s_1+s_2=1$, $s_2=\tfrac12$이다. 이 후보가 소거되지 않으면 $(1,1)$ 중심의 두 변수 테일러 전개의 공통 수렴 반경은 $\tfrac12$를 넘지 못한다. 더 왼쪽에는 고차 소수 거듭제곱의 층과 $1/s_2$의 극이 있다.

대각선에서는 [주정리 (3.8)](main-theorem.md#equation-3-8)에 의해 세 특이점이 $\zeta(2s)^{-r(r+1)/2}$로 합쳐진다.

- **$r(r+1)/2$가 0 이상의 정수인 경우:** 합친 결과가 $s=\tfrac12$에서 해석적이므로 $m_{j,0}$의 생성함수의 반경은 다음 특이점까지 늘어난다. $r=0,\pm1$에서는 $s=0$의 극까지 거리 $1$이다. 같은 차수의 핵들이 총질량에서 서로 지워지는 경우이다.
- **비상수 모멘트:** 일반적으로 이 상쇄가 유지되지는 않는다. 예를 들어 아크사인에서는 $\partial_1^2G$가 대각선의 $s=\tfrac12$에서 특이하다(§3).

**[휴리스틱]** $r\notin\mathbb Z$에서 $1/\Gamma(r-j)$의 계승 성장은, 반경 $\tfrac12$인 모멘트 급수의 최소항이 $j\approx L/2$에 있고 크기가 대략 $x^{1/2}$임을 시사한다. 최적 절단 오차를 얻으려면 (1.2)의 상수와 나머지를 증가하는 $J$에 대해 제어해야 한다.

<a id="section-3"></a>

## 3. 아크사인의 세 점 핵

$z=w=\tfrac12$이면 $e_1=e_2=\tfrac38$, $zw=\tfrac14$이고, 지수 합이 $1$이다. 대각선에서 $H(s,s)=1/\zeta(2s)$는 $s=\tfrac12$에서 오히려 영점이므로 $m_{j,0}$에는 $s=\tfrac12$이 보이지 않는다. $m_{j,1}$도 $d\leftrightarrow m$ 대칭 때문에 보이지 않는다. 반면 $m_{j,2}$의 생성함수에는 대각선에서 다음 단순극이 생긴다.

$$
\partial_1^2\big[(2s_1-1)^{3/8}(s_1+s_2-1)^{1/4}\big](2s_2-1)^{3/8}\Big|_{s_1=s_2=s}=-\frac34(2s-1)^{-1}.
$$

따라서

<a id="equation-3-1"></a>

$$
m_{j,2}\sim c\,(-2)^j,\qquad c=-\frac3{16}\cdot2Z(\tfrac12)=-\frac{3}{16}\bigl(-\zeta(\tfrac12)\bigr)=-0.27382\ldots
\tag{3.1}
$$

**[수치]** [부록 A](#appendix-a).

비틀림이 있으면 $s-\tfrac12=q/L$에서 $H\approx\frac2L(q-i\xi)^{3/8}(q-\tfrac{i\xi}2)^{1/4}q^{3/8}$이다($\widetilde H(\tfrac12,\tfrac12)=1$). $s=1$에서 뻗은 두 절단을 위아래로 감싸는 Hankel 고리 위에서는 $\zeta(s_1)^{1/2}\zeta(s_2)^{1/2}\approx\zeta(\tfrac12)$이다. 세 점 핵을

<a id="equation-3-2"></a>

$$
\langle\mathcal T,e^{i\xi u}\rangle=\frac1{2\pi i}\int_{\mathcal U_Y}e^q(q-i\xi)^{3/8}\Big(q-\frac{i\xi}2\Big)^{1/4}q^{3/8}\,dq
\tag{3.2}
$$

로 정의하면 $\mathcal T$는 $[0,1]$에 지지된 분포이고

$$
\langle\mathcal T,1\rangle=0,\qquad\langle\mathcal T,u\rangle=0,\qquad\langle\mathcal T,u^2\rangle=-\frac3{16}.
$$

질량은 $q$의 총 차수가 $1$이라서, 1차 모멘트는 상수항의 Hankel 적분이라서 $0$이다. 2차 모멘트는 $(1-t)^{3/8}(1-t/2)^{1/4}$의 $t^2$ 계수 $-3/32$에서 온다.

**[국소 계산]** 고정된 $\xi$에서 $s=\tfrac12$ 근방의 기여는 $4\zeta(\tfrac12)\,x^{1/2}L^{-2}\langle\mathcal T,e^{i\xi u}\rangle+O_\xi(x^{1/2}L^{-3})$이다. $F=u^2$(S 수준)이면 $-\tfrac34\zeta(\tfrac12)x^{1/2}L^{-2}$이고, [영점 기여 §2](zero-contributions.md#section-2)의 정확한 계산과 일치한다. 다음 실축 특이점 후보는 소수 세제곱 항의 $s=\tfrac13$이다.

$r\in\mathbb Z$이면 $1$과 $1+i\alpha$의 절단을 두 점 사이의 선분으로 묶을 수 있다. 두 점을 함께 도는 모노드로미가 $e^{-2\pi ir}=1$이므로 $s=\tfrac12$ 근방의 가지가 단일값이고, 그 극한은 $\zeta(\tfrac12)^r$이다. 이는 $r=-1$에도 적용된다. 일반 정수 $r$의 세 점 핵은 (3.2)의 지수들을 $e_1,zw,e_2$로 바꾸어 얻으며, 총차수는 $r(r+1)/2$이다. 비정수 $r$에서는 절단의 도약을 따로 처리해야 한다.

<a id="section-4"></a>

## 4. 국소 연장과 전역 오차

$(0,1)$의 실수점에서는 $\zeta$가 영이 아니므로, 실축 부근의 국소 연장에 RH는 필요 없다. 이차 소수 인수를 분리한 잔여 인자 $\widetilde H$는 $\Re s,\Re t>1/3$에서 정칙이다. 다만 선도 세 점 핵을 더하는 것만으로 전역 오차가 $x^{1/3}$ 규모로 줄어드는 것은 아니다. 높은 영점의 기여와 나머지 윤곽은 별도로 제어해야 한다.

거듭제곱 오차와 영점 없는 영역의 관계는 [뫼비우스 경우의 명제 2](mobius-cases.md#section-5)를 참조한다.

<a id="appendix-a"></a>

## 부록 A. 아크사인 전개 계수의 수치 기록

아래 표와 실험 설명은 초기 원고의 보고값이다. 계산 코드와 원자료가 없어 재현하지 않았으며, [전역 공식의 검산](global-formula.md#checks)과는 검증 범위가 다르다.

$G$를 다중원 $|s_i-1|=0.3$ 위에서 계산하고 2차원 FFT로 $G_{ab}$를 구했다. 검산으로 $G_{00}=0.6079271$이다.

| $j$ | 개별 $\max\lvert G_{ab}\rvert$ | $m_{j,0}$ | $m_{j,1}$ | $m_{j,2}$ | $m_{j,2}/m_{j-1,2}$ |
|--:|--:|--:|--:|--:|--:|
| 8 | 5.3 | −0.9947 | 0.8949 | −68.89 | −1.981 |
| 10 | 14.5 | −0.9987 | 0.9116 | −273.26 | −1.995 |
| 12 | 42.8 | −0.9997 | 0.9168 | −1094.77 | −2.003 |
| 14 | 134.1 | −0.9999 | 0.9183 | −4401.24 | −2.005 |

- **$m_{j,0}$:** $s=0$의 극에서 예측되는 $(-1)^{j+1}$로 간다.
- **$m_{j,2}$:** 비가 $-2$로 가므로 $s=\tfrac12$의 단순극이 보인다. $m_{14,2}/(-2)^{14}=-0.2686$이 (3.1)의 $-0.2738$로 다가간다. 느리게 다가가는 것은 $s=\tfrac13$의 다음 층 때문으로 보인다.

<a id="references"></a>

## 참고 문헌

- G. Tenenbaum, *Introduction to Analytic and Probabilistic Number Theory*, II.5: Selberg–Delange 방법.
- NIST DLMF, [§5.9](https://dlmf.nist.gov/5.9): Hankel 적분. [§13.6](https://dlmf.nist.gov/13.6): Bessel 함수. [§25.16](https://dlmf.nist.gov/25.16): $\psi$의 명시적 공식.
