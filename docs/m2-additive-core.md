# $M_2$의 가법 핵: 점별 오차, FTA와 분배법칙, 이동 상관

[이차 모멘트](m2-divisor-moment.md)의 잔차 $\mathcal E_2$를 점별로 나누고, 그 오차 부분이 어떤 산술 정보를 담는지 정리한다. 출발점은 [cp 계열 명제 1.3](m2-cp-families.md#proposition-1-3)의 서로소 이차 가법성이다. 결론은 다음과 같다.

1. $M_2=p_2+e_2$로 나누면 $p_2$는 $\log n$만의 함수이고 $\sum e_2$는 톱니항과 상수를 빼면 $\mathcal E_2$와 같다. $e_2$에는 $b$의 편차만 남는다.
2. $e_2$는 로그 보정항이 붙은 코사이클 항등식을 만족하고, 그 핵 $b-\gamma\log=1*(\Lambda(\log-\gamma))$은 정확히 가법적이다.
3. 점별 $e_2$는 소인수 크기 분할의 통계량 $\sum V_i^2$의 편차이지만, 합에서는 가법성 때문에 소수 거듭제곱의 분포로 접힌다. 여기까지는 FTA만 쓰며 Beurling 체계에서도 성립한다.
4. FTA 밖의 정보는 분배법칙에서 온다. 점별로는 $d\mathbb Z$ 위의 이동이 단위 이동으로 환원되고, 합에서는 CRT가 된다. 완전가법 핵 $h$의 분극 항등식은 이동 상관을 $n^2$과 $n(n+1)$에서의 이차 모멘트로 바꾼다.
5. 이동 상관 $\sum h(n)h(n+1)$은 분포 수준 $1$, $\tfrac12$, $>\tfrac12$에 대응하는 세 영역으로 나뉜다. 오차의 톱니는 모듈러 역원으로 꼬여 있으며, 역원이 등분포하는 범위에서는 기약 Hermite 항등식에 의해 $\mathcal E_2$형 톱니로 정확히 바뀐다.

절·명제·수식 번호는 이 문서 안에서 독립적이다. 영점에 관한 결과는 쓰지 않는다.

**기호.** $M_2$, $a_2=(\log n)^2M_2$, $b$, $Q_2$, $Q_b$, $\mathscr L$, $T_2$, $\mathcal P_2$, $\mathcal E_2$는 [이차 모멘트](m2-divisor-moment.md)를 따른다. 그 문서 (5.2)의 상수 $c=\gamma^2+3\gamma_1$을 여기서는 띠 번호 $c$와 구분해 $\kappa:=\gamma^2+3\gamma_1$로 쓴다. 이는 그 문서의 $\kappa_2$, $\kappa_\Lambda$와 다르다. $L=\log n$, $\widetilde B_1(t):=\{t\}-\tfrac12$, $e(\theta)=e^{2\pi i\theta}$이다. 가법함수는 서로소 곱에서 가법적인 함수, 완전가법함수는 모든 곱에서 가법적인 함수이다.

**범위 표시.** 표시가 없는 명제는 가정 없이 증명한다. **[개요]**는 표준 정리를 적용하는 증명의 개요, **[모형]**은 등분포를 가정한 모형 안의 계산, **[휴리스틱]**은 추가 추정이 필요한 논의, **[수치]**는 [검산 스크립트](../scripts/check_m2_additive_core.py)의 계산, **[미해결]**은 알려진 방법이 닿지 않는 부분이다.

| 결과 | 내용 | 범위 |
|:--|:--|:--|
| [명제 1.1, 정리 1.2](#section-1) | $p_2=\mathcal P_2'$, $e_2=(b-\beta)/(4L^2)$, $\sum e_2=\mathcal E_2+\widetilde B_1(x)p_2(x)+C+o(1)$ | 정확 |
| [명제 2.1–2.3](#section-2) | $e_2$의 코사이클, 가법 핵 $b-\gamma\log$, 띠 위의 $e_2$ | 정확 |
| [명제 3.1, 보조정리 3.2](#section-3) | $e_2$의 Poisson–Dirichlet 표현, 가법함수 $=1*$(소수 거듭제곱 위의 함수) | 정확 |
| [§4](#section-4) | FTA만 쓰는 단계의 Beurling 판정 | 관찰과 문헌 |
| [명제 5.1](#section-5) | 점별 분배법칙: $d\mathbb Z$ 위의 이동 $=$ 단위 이동 | 정확 |
| [명제 6.1–6.2](#section-6) | 완전가법 핵 $h=1*\Lambda(\Lambda-\gamma)$, $n^2$ 공식, 분극 항등식 | 정확 |
| [보조정리 7.1, 명제 7.2](#section-7) | 이동 상관의 정확한 개수 공식, 절단 영역의 점근식 | 정확 |
| [§7.3](#section-7-3) | 등차수열 영역과 큰 소수 쌍 영역 | **[개요]**, **[미해결]** |
| [보조정리 8.1, 명제 8.2](#section-8) | 기약 Hermite 항등식, 꼬인 톱니에서 $\mathcal E_2$형 톱니로 | 정확, **[모형]** |

- [§1. 점별 주항과 오차](#section-1)
- [§2. $e_2$의 항등식](#section-2)
- [§3. 두 층: 분할 통계량과 소수 거듭제곱](#section-3)
- [§4. FTA와 Beurling 판정](#section-4)
- [§5. 분배법칙](#section-5)
- [§6. 완전가법 핵과 분극](#section-6)
- [§7. 이동 상관과 세 영역](#section-7)
- [§8. 꼬인 톱니와 Hermite 항등식](#section-8)
- [§9. 정리와 남은 문제](#section-9)
- [검산](#checks)
- [참고 문헌](#references)

<a id="section-1"></a>

## 1. 점별 주항과 오차

<a id="equation-1-1"></a>

$$
\beta(t):=Q_b'(t)=\tfrac12\log^2t+\gamma\log t-\kappa,\qquad
p_2(n):=\frac38+\frac{\gamma}{4\log n}-\frac{\kappa}{4\log^2n},\qquad
e_2(n):=M_2(n)-p_2(n)
\tag{1.1}
$$

로 둔다 ($n\ge2$).

**명제 1.1.** $p_2=\mathcal P_2'$이고

<a id="equation-1-2"></a>

$$
\boxed{e_2(n)=\frac{b(n)-\beta(n)}{4\log^2n}.}
\tag{1.2}
$$

**증명.** $\mathscr Lf$를 미분하면 적분항의 미분 $2f/(x\log^3x)$와 첫 항 미분의 둘째 조각이 상쇄되어 $(\mathscr Lf)'=f'/\log^2x$이다. 따라서 $\mathcal P_2'=Q_2'/\log^2x$이다. [이차 모멘트 (5.3)](m2-divisor-moment.md#section-5-1)을 미분하면 $Q_2'=\tfrac38L^2+\tfrac\gamma4L-\tfrac\kappa4=\tfrac14(L^2+\beta)$이고, (5.9)를 미분하면 $Q_b'=\beta$이다. $4a_2=L^2+b$에서 빼면 (1.2)이다. $\square$

$(\log n)^2$ 부분은 전부 $p_2$로 가고, $e_2$에는 $b$의 편차만 남는다.

<a id="theorem-1-2"></a>

**정리 1.2 (합의 대응).** 상수 $C$가 있어 $x\ge2$에서

<a id="equation-1-3"></a>

$$
\boxed{\sum_{2\le n\le x}e_2(n)=\mathcal E_2(x)+\widetilde B_1(x)\,p_2(x)+C+O\Bigl(\frac1{x\log^2x}\Bigr).}
\tag{1.3}
$$

**증명.** $f=p_2$는 매끄럽고 $f'(t)\ll1/(t\log^2t)$, $f''(t)\ll1/(t^2\log^2t)$이다. Euler–Maclaurin 공식으로

$$
\sum_{2\le n\le x}f(n)=\int_2^xf-\widetilde B_1(x)f(x)+\tfrac12f(2)+\int_2^x\widetilde B_1f'
$$

이다. $\widetilde B_2$를 평균 $0$인 $\widetilde B_1$의 주기적 원시함수라 하면 $\int_x^\infty\widetilde B_1f'=-\widetilde B_2(x)f'(x)-\int_x^\infty\widetilde B_2f''=O(1/(x\log^2x))$이므로 $\int_2^x\widetilde B_1f'$은 상수에 $O(1/(x\log^2x))$를 더한 것이다. $\int_2^xp_2=\mathcal P_2(x)-\mathcal P_2(2)$를 넣고 $T_2$에서 빼면 된다. $\square$

톱니항은 $T_2$가 계단함수이기 때문에 생기며 반정수 $x$에서 사라진다. [정리 8.1](m2-divisor-moment.md#section-8-1)에서 $\mathcal E_2$의 성장 지수는 $\Theta_{\rm simple}\ge\tfrac12$이므로 $\sum e_2$와 $\mathcal E_2$의 성장 지수는 같다. 수치로는 $C$에 해당하는 $\sum p_2-\int p_2+\widetilde B_1p_2$가 $x=10^3,\dots,10^6$에서 $0.266993$으로 안정적이다 **[수치]**.

<a id="section-2"></a>

## 2. $e_2$의 항등식

$\varepsilon(n):=4\log^2n\cdot e_2(n)=b(n)-\beta(n)$으로 둔다.

**명제 2.1 (코사이클).** $\gcd(m,n)=1$, $s=\log m/\log mn$, $t=1-s$이면

<a id="equation-2-1"></a>

$$
\boxed{
\varepsilon(mn)=\varepsilon(m)+\varepsilon(n)-\log m\log n-\kappa,\qquad
e_2(mn)=s^2e_2(m)+t^2e_2(n)-\frac{st}4-\frac{\kappa}{4\log^2mn}.
}
\tag{2.1}
$$

**증명.** $A=\log m$, $B=\log n$으로 두면 $\beta(mn)-\beta(m)-\beta(n)=AB+\kappa$이다($\gamma$ 항은 상쇄된다). $b$는 서로소 곱에서 가법적이므로 $\varepsilon(mn)=b(m)+b(n)-\beta(mn)$에서 첫 식이 나온다. $4\log^2mn$으로 나누면 둘째 식이다. $\square$

**명제 2.2 (가법 핵).** $\varepsilon^\sharp:=b-\gamma\log=\varepsilon+\tfrac12\log^2-\kappa$는 가법적이고

<a id="equation-2-2"></a>

$$
e_2(n)=\frac{\varepsilon^\sharp(n)}{4\log^2n}-\frac18+\frac{\kappa}{4\log^2n}.
\tag{2.2}
$$

**증명.** (2.1)의 결손 $-AB-\kappa$는 $\kappa-\tfrac12\log^2$의 코사이클이므로 $\varepsilon-(\kappa-\tfrac12\log^2)$은 가법적이다. 직접 보면 $b$와 $\log$가 모두 가법적이다. $\square$

$e_2$의 항등식은 $p_2$가 $\log n$만의 함수이므로 [cp 계열 (1.5)](m2-cp-families.md#equation-1-5)에 결정적인 결손항을 붙인 것이고 새로운 산술 정보는 없다. 정확히 가법적인 대상은 $\varepsilon^\sharp$이다.

**명제 2.3 (띠 위의 $e_2$).** $p\nmid c$, $u=\log c$, $v=\log p$, $L=u+v$, $t=u/L$이면

<a id="equation-2-3"></a>

$$
e_2(cp)=t^2e_2(c)+\frac{(1-t)(1-3t)}8-\frac{\gamma(1-t)}{4L},\qquad
4L^2e_2(cp)=\tfrac12(v-u)^2-u^2+b(c)-\gamma L+\kappa.
\tag{2.3}
$$

**증명.** $e_2(p)=\tfrac18-\tfrac\gamma{4v}+\tfrac\kappa{4v^2}$를 (2.1)에 넣으면 $(1-t)^2e_2(p)$의 $\kappa$ 항이 $-\kappa/(4L^2)$와 상쇄된다. 둘째 식은 $b(cp)=b(c)+v^2$와 $v^2-\tfrac12L^2=\tfrac12(v-u)^2-u^2$이다. $\square$

$M_2$의 재귀 $\tfrac12-\tfrac t2+t^2M_2(c)$에서 상수항은 $t$에 대해 단조였지만, $e_2$에서는 상수항 $(1-t)(1-3t)/8$이 $t=\tfrac13$, 즉 $\log p=2\log c$ 근처에서 부호를 바꾼다. $p\to\infty$이면 $e_2(cp)\to\tfrac18$이고 소수에서도 $e_2(p)\to\tfrac18$이다. 따라서 $e_2$는 점별로 작지 않으며, $\mathcal E_2$의 크기는 합성수 쪽 음의 값과의 상쇄로 정해진다.

<a id="section-3"></a>

## 3. 두 층: 분할 통계량과 소수 거듭제곱

**명제 3.1 (점별: 분할 통계량).** $n=p_1p_2\cdots$ (중복 포함), $V_i=\log p_i/\log n$이고 $s(n)$을 [cp 계열 §6.2](m2-cp-families.md#section-6-2)의 제곱인수 보정이라 하면

<a id="equation-3-1"></a>

$$
e_2(n)=\frac14\Bigl(\sum_iV_i^2-\frac12\Bigr)+\frac14s(n)-\frac{\gamma}{4\log n}+\frac{\kappa}{4\log^2n}.
\tag{3.1}
$$

**증명.** (1.2)에서 $e_2=\tfrac14(b/L^2-\tfrac12)-\gamma/(4L)+\kappa/(4L^2)$이고 $b/L^2=\sum V_i^2+s(n)$이다. $\square$

점별로 $e_2$는 소인수 크기 분할의 이차 통계량 $\sum V_i^2$가 Poisson–Dirichlet(1) 평균 $\tfrac12$에서 벗어난 양이다. 이는 인수분해의 한 사영이며, 그 평균은 로그의 모든 차수까지 $p_2$에 흡수된다([이차 모멘트 (9.1)](m2-divisor-moment.md#section-9)).

**보조정리 3.2 (합: 소수 거듭제곱으로 접힘).** $f$가 가법적이고 $f(1)=0$이면 $f=1*g$이다. 여기서 $g$는 소수 거듭제곱에만 값을 가지고 $g(p^k)=f(p^k)-f(p^{k-1})$이다. 따라서

<a id="equation-3-2"></a>

$$
\sum_{n\le x}f(n)=\sum_{p^k\le x}g(p^k)\Bigl\lfloor\frac x{p^k}\Bigr\rfloor.
\tag{3.2}
$$

**증명.** 각 $q^a\parallel n$에서 $f(q^a)=\sum_{k\le a}g(q^k)$로 망원합을 쓰면 $f(n)=\sum_{p^k\mid n}g(p^k)$이다. $\square$

**따름정리 3.3.** $\varepsilon^\sharp=1*\bigl(\Lambda(\log-\gamma)\bigr)$이고, $\Re s>1$에서

<a id="equation-3-3"></a>

$$
\sum_{n\ge1}\frac{b(n)-\gamma\log n}{n^s}=\zeta(s)\Bigl[\Bigl(\frac{\zeta'}{\zeta}\Bigr)'(s)+\gamma\frac{\zeta'}{\zeta}(s)\Bigr].
\tag{3.3}
$$

**증명.** $g(p^k)=k\log^2p-\gamma\log p=\Lambda(p^k)(\log p^k-\gamma)$이다. $\sum\Lambda(n)\log n\,n^{-s}=(\zeta'/\zeta)'$, $\sum\Lambda(n)n^{-s}=-\zeta'/\zeta$이다. $\square$

두 층은 다음과 같다. 점별 요동, 즉 $\sum V_i^2$의 분포는 평균만 남기고 합에서 사라진다. $\sum_{n\le x}\varepsilon^\sharp(n)$에 남는 진동은 가중치 $\Lambda(\log-\gamma)$를 단 소수 거듭제곱의 분포와 격자항 $\lfloor x/d\rfloor=x/d-\widetilde B_1(x/d)-\tfrac12$뿐이다. (3.3)에서 $\zeta$ 인수는 "모든 $n$에 걸쳐 약수 $d$를 센다"는 격자 구조이고, 괄호는 소수 쪽이다. [정리 8.1](m2-divisor-moment.md#section-8-1)이 단순영점만 보는 것은 이 곱에서 $\zeta$ 인수가 $(\zeta'/\zeta)'$의 극 차수를 하나 깎기 때문이다.

<a id="section-4"></a>

## 4. FTA와 Beurling 판정

FTA를 가법적으로 쓰면 $\log=1*\Lambda$이고, 이는 오일러 곱과 동치이다. 보조정리 3.2의 전단사 $f=1*g$, $g=\mu*f$는 "가법함수 ↔ 소수 거듭제곱 위의 함수"이며 이것이 FTA 그 자체이다. 따라서 점별로 소수 거듭제곱 정보와 인수분해 정보는 뫼비우스 반전으로 손실 없이 오간다.

**관찰 4.1.** §1–3의 항등식, 즉 (1.1)–(3.3)과 [cp 계열 (1.5)](m2-cp-families.md#equation-1-5), [이차 모멘트 (5.1), (5.5)](m2-divisor-moment.md#section-5)는 오일러 곱과 $s=1$에서의 로랑 전개만 쓴다. 따라서 정수 계수가 $N(x)=\varrho x+O(x^\theta)$ ($\theta<1$)인 Beurling 일반화 수체계에서도 그 체계의 상수로 그대로 성립한다.

그런데 점근적 규칙성은 FTA를 통해 어느 방향으로도 옮겨지지 않는다.

- Diamond–Montgomery–Vorhauer [1]는 $1/2<\theta<1$에서 $N(x)=\varrho x+O(x^\theta)$이면서 소수 계수가 de la Vallée Poussin 오차 $x\exp(-c\sqrt{\log x})$를 정확히 실현하는 Beurling 체계를 구성했다.
- Broucke–Debruyne–Vindas [2]는 소수 계수가 $\pi(x)=\operatorname{Li}(x)+O(\sqrt x)$인데 정수 계수가 $\Omega_\pm\bigl(x\exp(-c\sqrt{\log x\log\log x})\bigr)$로 진동하는 체계를 구성했다.

**판정 기준.** FTA와 정수 계수의 거듭제곱 절약만 쓰는 논증은 고전적 영점 없는 영역을 넘을 수 없다. 넘는다면 [1]의 체계에서도 같은 결론이 나와 모순이다. 각 단계마다 "이 단계는 [1]의 체계에서도 성립하는가"를 물으면, 그 단계가 RH 쪽 정보를 담을 수 있는지 판정할 수 있다. [1]의 반례를 깨는 정보는 정수의 덧셈 구조에서만 나온다. 후보는 정확한 격자 $N(x)=\lfloor x\rfloor$, 그 결과인 Poisson 합과 함수방정식, 이동 $n\mapsto n+1$과 곱셈의 상호작용이다.

<a id="section-5"></a>

## 5. 분배법칙

**명제 5.1 (점별 분배법칙).** $\gcd(d,m(m+1))=1$이면

<a id="equation-5-1"></a>

$$
\boxed{
\varepsilon^\sharp(dm+d)-\varepsilon^\sharp(dm)=\varepsilon^\sharp(m+1)-\varepsilon^\sharp(m),
}
\qquad
\varepsilon(dm+d)-\varepsilon(dm)=\varepsilon(m+1)-\varepsilon(m)-\log d\,\log\Bigl(1+\frac1m\Bigr).
\tag{5.1}
$$

**증명.** $dm+d=d(m+1)$이고, 양쪽에 명제 2.2의 가법성(또는 (2.1)의 코사이클)을 쓰면 $\varepsilon^\sharp(d)$와 $\kappa$가 소거된다. $\square$

덧셈과 곱셈의 상호작용 중 새로운 것은 $\Delta\varepsilon^\sharp(m)=\varepsilon^\sharp(m+1)-\varepsilon^\sharp(m)$ 하나뿐이다. Beurling 체계에는 $+1$이 없으므로 이 양은 §4의 판정을 통과한다.

이 연산자에는 강한 경직성 정리가 있다. Erdős [3]는 가법함수 $f$가 $f(n+1)-f(n)\to0$을 만족하면 $f=C\log$임을 보였다. 평균판 $\frac1x\sum_{n\le x}|\Delta f(n)|\to0$은 Kátai와 Wirsing이 독립적으로, 밀도 $1$인 집합 위의 판은 Hildebrand가 증명했다([4]의 서술). 따라서 $\log$의 상수배는 $+$와 양립하는 유일한 가법함수이고, $\varepsilon^\sharp$에서 $\log$ 성분을 뺀 나머지가 분배법칙이 감지하는 정보이다.

$\sum_{n\le x}\Delta\varepsilon^\sharp(n)$은 망원합이라 정보가 없다. $\mathcal E_2$ 같은 선형합에 들어가는 분배법칙은 "$d$의 배수가 등간격이다", 즉 $\#\{n\le x:d\mid n\}=\lfloor x/d\rfloor$뿐이다. Beurling 체계에서는 이 자리에 $N(x/d)$가 온다. 첫 번째로 비자명한 대상은 이차 상관이며 §6–8에서 다룬다.

<a id="section-6"></a>

## 6. 완전가법 핵과 분극

$b$는 완전가법이 아니다. $p^a\parallel n$이면 $n^2$에서 그 소수의 몫이 $a(2a+1)$이기 때문이다. 대신

<a id="equation-6-1"></a>

$$
w:=1*\Lambda^2=\sum_{p^a\parallel n}a\log^2p,\qquad
h:=w-\gamma\log=1*G,\qquad G:=\Lambda(\Lambda-\gamma),\qquad
r:=\sum_{p^a\parallel n}\binom a2\log^2p
\tag{6.1}
$$

로 둔다. $w$와 $h$는 완전가법이다.

**명제 6.1.** $b=w+r$이고 $r=(\log n)^2s(n)$이다. 또

<a id="equation-6-2"></a>

$$
\varepsilon=h+r-\tfrac12\log^2n+\kappa,\qquad
a_2(n^2)=4a_2(n)-\tfrac14w(n),\qquad
\varepsilon(n^2)=4\varepsilon(n)-w(n)+2\gamma\log n-3\kappa.
\tag{6.2}
$$

**증명.** $\tfrac12a(a+1)=a+\binom a2$에서 $b=w+r$이다. 첫 식은 $\varepsilon=b-\beta$에 넣은 것이다. $\tfrac12\cdot2a(2a+1)=4\cdot\tfrac12a(a+1)-a$이므로 $b(n^2)=4b(n)-w(n)$이고, $4a_2(n^2)=4L^2+b(n^2)$와 $\beta(n^2)=2L^2+2\gamma L-\kappa$를 쓰면 나머지 두 식이다. $\square$

$r$은 제곱인수 부분에만 있으며 $\sum_{n\le x}s(n)\ll x/\log x$이다([cp 계열 명제 6.1 (ii)](m2-cp-families.md#section-6-2)). 따라서 $\varepsilon$의 핵으로 $h$를 써도 된다.

**명제 6.2 (분극 항등식).**

<a id="equation-6-3"></a>

$$
\boxed{2h(n)h(n+1)=h\bigl(n(n+1)\bigr)^2-\tfrac14h(n^2)^2-\tfrac14h\bigl((n+1)^2\bigr)^2.}
\tag{6.3}
$$

**증명.** $\gcd(n,n+1)=1$이므로 $h(n(n+1))=h(n)+h(n+1)$이고, 완전가법성에서 $h(n^2)=2h(n)$이다. 전개하면 된다. $\square$

선형 수준에서는 $\sum h(n(n+1))=2\sum h(n)+(\text{경계항})$, $h(n^2)=2h(n)$이므로 새 정보가 없고, 둘 다 Beurling 체계에서도 성립한다. 제곱하면 이동 상관이 다항식 값에서의 이차 모멘트로 바뀐다. $\sum h(n^2)^2=4\sum h(n)^2$은 곱셈만으로 정해지는 대조군이고, $d\mid n(n+1)$의 개수에는 분배법칙(CRT)이 들어간다. 국소 밀도는 $\rho_{n^2}(p^k)=p^{\lfloor k/2\rfloor}$, $\rho_{n(n+1)}(p^k)=2$이다.

<a id="section-7"></a>

## 7. 이동 상관과 세 영역

<a id="equation-7-1"></a>

$$
C(x):=\sum_{n\le x}h(n)h(n+1)=\sum_{(d,e)=1}G(d)G(e)\,N(x;d,e),\qquad
N(x;d,e):=\#\{n\le x:\ d\mid n,\ e\mid n+1\}
\tag{7.1}
$$

로 둔다. $d\mid n$, $e\mid n+1$이면 $(d,e)=1$이고, $G$는 소수 거듭제곱에만 값을 가지므로 $d,e\ge2$이다.

<a id="section-7-1"></a>

### 7.1. 정확한 개수

**보조정리 7.1.** $(d,e)=1$, $e\ge2$, $\bar d\equiv d^{-1}\pmod e$이면

<a id="equation-7-2"></a>

$$
\boxed{N(x;d,e)=\frac x{de}-\Bigl[\widetilde B_1\Bigl(\frac x{de}+\frac{\bar d}e\Bigr)-\widetilde B_1\Bigl(\frac{\bar d}e\Bigr)\Bigr].}
\tag{7.2}
$$

**증명.** $n=dm$으로 쓰면 조건은 $m\equiv-\bar d\pmod e$, $1\le m\le x/d$이다. 대표원을 $c\in[1,e-1]$로 잡으면 개수는 $\lfloor(x/d-c)/e\rfloor+1=\frac x{de}+\tfrac12-\tfrac ce-\widetilde B_1\bigl(\tfrac{x/d-c}e\bigr)$이다. $-c/e\equiv\bar d/e\pmod1$이고 $\widetilde B_1$은 정수가 아닌 점에서 홀함수이므로 $\tfrac12-\tfrac ce=\widetilde B_1(\bar d/e)$이다. $\square$

<a id="section-7-2"></a>

### 7.2. 절단 영역

**명제 7.2.** $2\le D\le x$이면

<a id="equation-7-3"></a>

$$
\sum_{\substack{de\le D\\(d,e)=1}}G(d)G(e)N(x;d,e)
=x\Bigl[\frac{(\log D)^4}{24}+O(\log^3D)\Bigr]+O(D\log^3D).
\tag{7.3}
$$

**증명.** (7.2)에서 각 쌍의 오차는 $1$ 이하이다. $|G(d)|\le2\Lambda(d)\log D$이고 Chebyshev 상계로 $\sum_{m\le D}(\Lambda*\Lambda)(m)\ll D\log D$이므로 오차의 합은 $O(D\log^3D)$이다. 주부는 $x\sum G(d)G(e)/(de)$이다. Mertens 공식과 부분합에서 $A(y):=\sum_{d\le y}G(d)/d=\tfrac12\log^2y+O(\log y)$이다. $\Lambda^2$와 $\Lambda\log$의 차는 $k\ge2$인 소수 거듭제곱에서만 생기므로 $O(1)$이고, $\gamma\sum\Lambda(d)/d=O(\log y)$이다. 따라서 $\Lambda=\log D$로 두면 $\sum_{d\le D}\frac{G(d)}dA(D/d)=\int_0^{\Lambda}\tfrac12(\Lambda-u)^2\,u\,du+O(\Lambda^3)=\Lambda^4/24+O(\Lambda^3)$이다. 같은 소수의 쌍은 $\sum_p\log^4p/p^2<\infty$만큼만 기여한다. $\square$

$D=x^{1-\delta}$이면 절단 영역은 $\tfrac{(1-\delta)^4}{24}x\log^4x\,(1+o(1))$이다. 한편 $\sum_{n\le x}w(n)=\sum_{p^k\le x}\log^2p\lfloor x/p^k\rfloor=\tfrac12x\log^2x+O(x\log x)$이므로 $h$의 평균은 $\tfrac12\log^2x$에 가깝다. $n$과 $n+1$의 소인수 구조가 독립이라면 $C(x)\sim\tfrac14x\log^4x$이고 **[휴리스틱]**, 절단 영역은 주항의 많아야 $1/6$만 담는다. 나머지는 $de>x$에 있으므로 이 영역은 오차항이 아니라 주항 영역이다.

수치로 $x=10^6$에서 $C(x)/(x\log^4x)=0.1799$이고 $(\text{$h$의 평균})^2/\log^4x=0.1807$이어서 독립성 예측과 맞는다. 같은 $x$에서 (7.3)의 주합은 $D=10^6$일 때 전체의 약 $14\%$이다. $(S(D)-(\log D)^4/24)/(\log D)^3$은 $D=10^3,\dots,10^6$에서 $-0.19,\,-0.21,\,-0.21,\,-0.22$로 유계이다 **[수치]**.

<a id="section-7-3"></a>

### 7.3. 등차수열 영역과 큰 소수 쌍 영역

**영역 (II): $\min(d,e)\le x^{1/2}(\log x)^{-B}$, $de>D$ [개요].** 예를 들어 $e$가 작으면 $n=dm$에서 $m$이 작고, 남는 조건은 "$d\le x/m$, $d\equiv-\bar m\pmod e$"인 소수 거듭제곱의 개수이다. 가중치 $\log^2$는 부분합으로 처리한다. Bombieri–Vinogradov 정리 [5]에서 이 영역은 국소 밀도 $1/\varphi(e)$의 주항과 $O_A(x(\log x)^{-A})$ 오차로 정리된다. 오차는 로그만큼만 절약되며, 거듭제곱 절약에는 법 $e$의 디리클레 $L$-함수에 대한 정보가 필요하다.

**영역 (III): $d,e>x^{1/2}(\log x)^{-B}$ [미해결].** $n=dm$, $n+1=ek$에서 $m,k$가 $x^{1/2}$ 근처까지 커지므로 분포 수준이 $\tfrac12$보다 커야 한다(Elliott–Halberstam 수준). 이는 "$n$과 $n+1$이 모두 $\sqrt x$ 근처 이상의 소인수 거듭제곱을 가진다"는 사건의 가중 개수이며, 독립성 예측대로라면 주항 크기의 양을 갖는다. 각 쌍 $(d,e)$는 $n$을 많아야 하나 주므로 개수 공식 (7.2)의 주항은 의미가 없다. 이 영역은 $pa+1=qb$ ($p,q$는 큰 소수)를 세는 이항 가법 문제와 같고, 분배법칙의 정보는 여기에 몰려 있다.

<a id="section-8"></a>

## 8. 꼬인 톱니와 Hermite 항등식

보조정리 3.2에서 $\mathcal E_2$ 쪽 톱니는 꼬임이 없다.

<a id="equation-8-1"></a>

$$
\sum_{n\le x}h(n)=x\sum_{d\le x}\frac{G(d)}d-\sum_{d\le x}G(d)\Bigl(\widetilde B_1\Bigl(\frac xd\Bigr)+\frac12\Bigr).
\tag{8.1}
$$

반면 (7.2)의 톱니 $\widetilde B_1\bigl(\frac x{de}+\frac{\bar d}e\bigr)$는 모듈러 역원으로 꼬여 있다. 이 꼬임이 CRT, 즉 분배법칙의 흔적이다.

**보조정리 8.1 (기약 Hermite 항등식).** $e=q^j$이면

<a id="equation-8-2"></a>

$$
\boxed{\sum_{\substack{a\bmod e\\q\nmid a}}\widetilde B_1\Bigl(t+\frac ae\Bigr)=\widetilde B_1(et)-\widetilde B_1(et/q),}
\qquad
\sum_{\substack{a\bmod e\\q\nmid a}}\widetilde B_1\Bigl(\frac ae\Bigr)=0.
\tag{8.2}
$$

**증명.** Hermite 항등식 $\sum_{a=0}^{e-1}\widetilde B_1(t+a/e)=\widetilde B_1(et)$에서 $a=qa'$인 항들의 합 $\sum_{a'\bmod e/q}\widetilde B_1\bigl(t+a'/(e/q)\bigr)=\widetilde B_1(et/q)$를 뺀다. 둘째 식은 $a\leftrightarrow e-a$ 대칭이다. $\square$

**명제 8.2 [모형].** (7.1)에서 $(d,e)$를 $x/(de)$가 거의 일정한 블록으로 묶고, 각 블록에서 $\bar d$가 법 $e$의 기약 잉여류 위에 등분포한다고 하자. 그러면 (7.2)의 톱니 부분은

<a id="equation-8-3"></a>

$$
-\sum_dG(d)\sum_{e=q^j}\frac{G(e)}{\varphi(e)}\Bigl[\widetilde B_1\Bigl(\frac xd\Bigr)-\widetilde B_1\Bigl(\frac x{dq}\Bigr)\Bigr]
=\sum_NW(N)\,\widetilde B_1\Bigl(\frac xN\Bigr)
\tag{8.3}
$$

가 된다. 이는 (8.1)과 같은 꼴이고 가중치만 $W$로 바뀐 것이다. $e$에 대한 합은 §7의 영역 경계에서 자른다.

**계산.** 블록 평균을 기약 잉여류 평균으로 바꾸고 (8.2)에 $t=x/(de)$를 넣으면 $et=x/d$, $et/q=x/(dq)$이다. $\widetilde B_1(\bar d/e)$ 항은 (8.2)의 둘째 식에서 평균이 $0$이다. $\square$

모형과 실제의 차이는 푸리에 전개 $\widetilde B_1(t)=-\sum_{m\ne0}e(mt)/(2\pi im)$로 측정된다. 그 결함은 소수 위의 Kloosterman 분수 합

<a id="equation-8-4"></a>

$$
\sum_{d\sim X}G(d)\,e\Bigl(\frac{m\bar d}e+\frac{mx}{de}\Bigr)
\tag{8.4}
$$

이다. $e$가 작으면 등차수열 문제, 즉 영역 (II)이고, $e$가 크면 영역 (III)의 핵심이다.

<a id="section-9"></a>

## 9. 정리와 남은 문제

- **FTA 층.** §1–3과 §6의 항등식은 오일러 곱만 쓰며 Beurling 체계에서도 성립한다. 띠 위의 Poisson–Dirichlet 요동도 이 층에 속한다.
- **분배법칙 층.** 점별로는 단위 이동 $\Delta\varepsilon^\sharp$ 하나로 줄고(명제 5.1), 합에서는 CRT로 들어간다(보조정리 7.1). $\mathcal E_2$가 실제로 쓰는 것은 가장 약한 형태인 톱니 $\widetilde B_1(x/d)$이다.
- **세 영역.** 이동 상관은 분포 수준 $1$(명제 7.2), $\tfrac12$(영역 II), $>\tfrac12$(영역 III)로 갈린다. 꼬인 톱니가 $\mathcal E_2$형 톱니로 바뀌는 범위(명제 8.2)도 같은 경계를 공유한다.

남은 문제는 다음과 같다.

1. 명제 8.2의 블록 평균을 엄밀하게 하고, 영역 (II)에서 (8.4)의 상쇄를 정량화한다.
2. 영역 (III)에서 소수 위 Kloosterman 분수 합 (8.4)가 상쇄되는가.
3. 이동 상관의 오차와 $\mathcal E_2$를 잇는 정확한 항등식이 있는가. 지금은 같은 꼴의 톱니로 환원되는 모형 수준의 연결만 있다.

<a id="checks"></a>

## 검산

[검산 스크립트](../scripts/check_m2_additive_core.py)는 numpy와 scipy를 사용하며 다음을 확인한다.

- 명제 1.1 ($n<5000$)과 정리 1.2의 상수가 $x=10^3,\dots,10^6$에서 안정적인 것.
- 명제 2.1–2.3, 명제 3.1, 따름정리 3.3, 명제 5.1, 명제 6.1–6.2의 항등식.
- 보조정리 7.1 (무작위 20,000개), 명제 7.2의 정규화된 차, $x=10^6$에서의 $C(x)$.
- 보조정리 8.1과 (8.1).

수치 확인은 증명을 대체하지 않는다.

    python -B scripts/check_m2_additive_core.py

<a id="references"></a>

## 참고 문헌

1. H. G. Diamond, H. L. Montgomery, U. M. A. Vorhauer, *Beurling primes with large oscillation*, Math. Ann. 334 (2006), 1–36.
2. F. Broucke, G. Debruyne, J. Vindas, *Beurling integers with RH and large oscillation*, arXiv:2004.11501.
3. P. Erdős, *On the distribution function of additive functions*, Ann. of Math. 47 (1946), 1–20.
4. A. P. Mangerel, *Additive functions in short intervals, gaps and a conjecture of Erdős*, arXiv:2108.12351 (§1.2에 Kátai, Wirsing, Hildebrand의 결과 서술).
5. E. Bombieri, *On the large sieve*, Mathematika 12 (1965), 201–225.
