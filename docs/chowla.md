# 강한 Chowla 추측: 제곱류 층과 4점 Wick 규칙

이동 상관 $C_h(x)=\sum_{n\le x}\mu(n)\mu(n+h)$의 강한 형태 $C_h(x)\ll x^{1/2+\varepsilon}$을 이산핵 틀에서 다시 본다. 극·영점 쌍으로 분해한 아래에 있는 구조를 찾고, 추측을 그 구조를 뺀 성분에 대한 진술로 줄인다. 결론은 다음과 같다.

1. $\sqrt x$의 진폭과 이동 사이의 공분산은 **제곱류 층**, 즉 $\mu^2$의 산술이 정한다. $C_h(x)=\sum_qA_h(q;x)\mu(q)$로 쓰면 부호는 제곱류 $q$에만 남고, 무작위 곱셈 모형의 공분산은 이 표현의 대각항이다. 제곱류가 일치하는 쌍은 팽창족과 시작점 $h^2/4$ 아래의 유한 개뿐이다. 순수 Chowla 합에는 극극·극영점 성분이 없다.
2. 서로소 부분(원시 부분)으로 나누면 이동 사이에 무상관인 성분이 되고, 수치상 $x$ 방향의 브라운 운동과 구별되지 않는다. 모형에서는 직교성과 고정된 $h$의 중심극한정리가 증명되어 있지만, 성분들이 독립인 것은 아니다.
3. 짧은 구간에서 이동에 대해 평균한 진술 $\sum_{k\ne0}C_k(I)^2\sim Q_I^2$은 짧은 Möbius 디리클레 다항식의 4차 모멘트가 가우스형이라는 진술과 **정확히** 동치이다. 고정된 $h$의 $C_h(x)\ll x^{1/2+\varepsilon}$와 동치인 것은 아니다. 영점 쪽에서는 $1/\zeta'$ 가중 영점합의 4점 Wick 규칙이 되고, 뫼비우스 쪽에서는 "4점 수준에서 μ는 μ² 위의 독립 부호"라는 한 줄이 된다. 소수별 무작위 곱셈 모형에서는 이 규칙이 짧은 구간에서 오차 없이 성립한다.
4. 이산핵·코사인 언어로 쓰면 $\xi\ne0$에서 극이 돌아온다. 그 출처는 끝점 $u\approx1$, 즉 큰 소수이다.
5. 4차를 2k차로 올리면 높이 $T\asymp X,X^2,\dots,X^k$의 평균이 덧셈 에너지 → μ 가중 Vinogradov 평균값 → Steinhaus 모형의 **사다리**를 이룬다. 꼭대기는 $Y\asymp\sqrt X$에서 정확히 대각적이지는 않지만, 모든 $Y\le X^{1-\delta}$에서 점근적으로 대각적이다. M4의 2k차판은 이동한 모멘트의 **순열(permanent) 규칙**이다. 수치상 μ는 6·8차와 μ² 메아리의 3-순환 위상까지 이 규칙을 따르고, 사다리의 모든 단에서 같은 값을 준다.

절·명제·수식 번호는 이 문서 안에서 독립적이다. 기호는 [주정리 §8](main-theorem.md#section-8)의 $d_z$, 이산핵 $\mathscr K_n^{z,w}$, 지수합 $\widetilde S_x^{z,w}(\xi)$와 [M과 ψ](mertens-psi.md)의 소수별 동전을 따른다. $e(\theta)=e^{2\pi i\theta}$이다.

범위 표시는 다음과 같다. **[모형]**은 무작위 곱셈 함수 모형 안의 정리, **[휴리스틱]**은 추가 추정이 필요한 논의, **[추측]**은 제안하는 추측, **[수치]**는 [검산 스크립트](../scripts/check_chowla.py)의 계산이다. 표시가 없는 명제는 가정 없이 증명한다.

| 결과 | 내용 | 범위 |
|:--|:--|:--|
| [명제 1.1–1.2](#section-1-1) | 팽창 항등식과 gcd 분해 | 정확 |
| [명제 1.3](#section-1-2) | 원시 부분의 직교성, 비선형 일치의 유한성 | **[모형]** |
| [명제 1.4–1.6](#section-1-4) | 제곱류 표현, 층별 밀도, 제곱류 일치의 완전한 분류 | 정확 |
| [명제 2.1–2.2](#section-2) | $F\in\mathrm{span}\{1,u\}$에서 극 성분 소멸, $\mathcal K_{-1/2,-1/2}=-\tfrac14\mathcal K_{3/2,3/2}''$ | 정확 |
| [추측 3.1](#section-3) | 원시 부분의 OU 법칙과 LIL형 상계 | **[추측]** |
| [정리 4.1](#section-4-1) | 짧은 구간의 덧셈 에너지 = 디리클레 다항식의 4차 모멘트 | 정확 |
| [추측 4.3 (M4)](#section-4-2) | μ의 4점 덧셈 통계에 대한 Wick 규칙, 정확한 나머지 $R_I(\theta)$ | **[추측]** |
| [명제 4.4](#proposition-4-4) | 독립 부호와 소수별 무작위 곱셈 모형에서 M4가 오차 없이 성립 | **[모형]** |
| [§4.3](#section-4-3) | 영점 쪽 번역: 중간 규모 창의 4점 Wick 규칙과 높이 메아리 | **[휴리스틱]**, RH |
| [정리 4.5](#theorem-4-5) | 2k차 모멘트 사다리: $T\asymp X^r$ 단은 μ 가중 $r$차 Vinogradov 평균값, 꼭대기는 Steinhaus 모형 | 정확 |
| [명제 4.6](#proposition-4-6) | 독립 부호 모형의 정확한 2k차 값 $D_k(Q)$와 이동한 6차의 순열식 | **[모형]** |
| [추측 4.7 (M2k)](#conjecture-4-7) | 순열 규칙과 평평한 사다리 | **[추측]** |
| [명제 4.8](#proposition-4-8) | 꼭대기의 점근적 대각성, $Y=2\sqrt X$의 반례 | 정확 |
| [명제 5.1](#section-5-1) | 이산핵의 부호·모양 분해와 $p^2$의 $\sin^2$ 동전 | 정확 |
| [§5.3](#section-5-3) | $\xi\ne0$에서 극의 귀환과 주호 모형 | **[수치]** |

- [§1. 제곱류 층](#section-1)
- [§2. 극 정보의 부재](#section-2)
- [§3. 원시 성분과 OU 법칙](#section-3)
- [§4. 짧은 구간 사전과 4점 Wick 규칙](#section-4)
- [§5. 이산핵과 코사인 언어](#section-5)
- [§6. 정리와 남은 문제](#section-6)
- [검산](#checks)
- [참고 문헌](#references)

<a id="section-1"></a>

## 1. 제곱류 층

<a id="equation-1-1"></a>

$$
C_h(x)=\sum_{n\le x}\mu(n)\mu(n+h),\qquad
S_2(h)=\prod_p\Big(1-\frac2{p^2}\Big)\prod_{p^2\mid h}\frac{p^2-1}{p^2-2}
\tag{1.1}
$$

로 둔다. $S_2(h)$는 $n$과 $n+h$가 함께 제곱 인수 없을 밀도(Mirsky)이고, 소수 $p$에서 $h$의 $p^2$ 나눔에만 의존한다. 국소 평균을 취하면 $\mathbb E_h S_2(h)=\prod_p(1-p^{-2})^2=(6/\pi^2)^2$이다.

<a id="section-1-1"></a>

### 1.1. 정확한 분해

**명제 1.1 (팽창 항등식).** 소수 $p$와 $h\ge1$에 대해

<a id="equation-1-2"></a>

$$
\sum_{\substack{n\le x\\ p\mid n}}\mu(n)\mu(n+ph)=\sum_{\substack{m\le x/p\\ p\nmid m(m+h)}}\mu(m)\mu(m+h).
\tag{1.2}
$$

**증명.** $n=pm$이면 $\mu(pm)\mu(p(m+h))=\mu(p)^2\mu(m)\mu(m+h)\mathbf 1[p\nmid m(m+h)]$이다. $\square$

즉 $C_{ph}(x)$는 $C_h(x/p)$의 일부를 문자 그대로 포함한다.

**명제 1.2 (gcd 분해).** 다음이 정확히 성립한다.

<a id="equation-1-3"></a>

$$
C_h(x)=\sum_{g\mid h}\mu^2(g)\,P^{[g]}_{h/g}(x/g),\qquad
P^{[g]}_k(y)=\sum_{\substack{m\le y,\ (m,k)=1\\ (g,\,m(m+k))=1}}\mu(m)\mu(m+k).
\tag{1.3}
$$

**증명.** $n$과 $n+h$가 제곱 인수 없으면 $g=\gcd(n,n+h)$는 $h$의 제곱 인수 없는 약수이다. $n=gm$, $k=h/g$로 쓰면 $(m,m+k)=1$이고, $n$과 $n+h$가 제곱 인수 없으므로 $(g,m(m+k))=1$이다. 이때 $\mu(n)\mu(n+h)=\mu(g)^2\mu(m)\mu(m+k)$이고, 대응은 일대일이다. $\square$

$g=1$ 성분 $C^\flat_h(x)=\sum_{n\le x,\,(n,h)=1}\mu(n)\mu(n+h)$을 **원시 부분**이라 한다. $(n,h)=1$이면 $\gcd(n,n+h)=1$이다.

<a id="section-1-2"></a>

### 1.2. 무작위 곱셈 모형에서의 직교성

$f$를 소수마다 독립인 ±1 부호를 곱한 무작위 곱셈 함수라 하자. 제곱 인수 없는 정수에만 지지되며, $\mathbb E[f(n_1)\cdots f(n_4)]=\mathbf 1[n_1n_2n_3n_4=\square]$이다.

**명제 1.3 [모형].** $C^\flat_h(f;x)=\sum_{n\le x,(n,h)=1}f(n)f(n+h)$에 대해

<a id="equation-1-4"></a>

$$
\mathbb E\big[C^\flat_h(f;x)^2\big]=\#\{n\le x:(n,h)=1,\ \mu^2(n)\mu^2(n+h)=1\},\qquad
0\le\mathbb E\big[C^\flat_h C^\flat_{h'}\big]\le\tfrac12\tau(|h^2-h'^2|)\quad(h\ne h').
\tag{1.4}
$$

**증명.** 원시 쌍이면 $n(n+h)$가 제곱 인수 없으므로, 네 수의 곱이 제곱수일 조건은 $n(n+h)=n'(n'+h')$이다. 이는 $(2n+h)^2-(2n'+h')^2=h^2-h'^2$와 같다. $h=h'$이면 $n=n'$이고, $h\ne h'$이면 해는 $|h^2-h'^2|$의 약수쌍으로 정해지므로 많아야 $\tfrac12\tau(|h^2-h'^2|)$개이다. $\square$

전체 $C_h$에서는 $n(n+h)=g^2\,m(m+k)$의 제곱 인수 없는 핵 $m(m+k)$끼리 일치해야 한다. 같은 이유로 그 일치는 선형족 $n'=(b/a)n$을 빼면 $x$와 무관한 유한 개뿐이다. 따라서 모형에서의 공분산은 **팽창 격자**로 주어진다. 서로소이고 제곱 인수 없는 $a,b$에 대해

<a id="equation-1-5"></a>

$$
\mathbb E\big[C_{ak}C_{bk}\big]=\#\Big\{m\le\frac x{\max(a,b)}:\ \mu^2(m)\mu^2(m+k)=1,\ (ab,\,m(m+k))=1\Big\}+O_{a,b,k}(1)
\tag{1.5}
$$

이고, $h'/h$가 제곱 인수 없는 수의 비가 아니면(예: $h'=4h,9h$) 주항이 없다. 오른쪽의 개수는 $\Delta_{ab}(k)\,x/\max(a,b)+o(x)$이고

$$
\Delta_A(k)=\prod_{p\mid A}\Big(1-\frac{2-\mathbf 1_{p\mid k}}p\Big)\prod_{p\nmid A}\Big(1-\frac{2-\mathbf 1_{p^2\mid k}}{p^2}\Big)
$$

이다. 비선형(Pell형) 일치가 별도의 층을 이룰 가능성은 이 명제로 배제되고, §1.4의 명제 1.6은 그 유한 개가 놓일 수 있는 범위까지 정한다.

<a id="section-1-3"></a>

### 1.3. 실측

**[수치]** $N=2^{25}$, $h\le2\cdot10^5$, $x=N-h_{\max}$에서 $z_h=C_h(x)/\sqrt{S_2(h)x}$의 분산은 $1.000$, 첨도는 $+0.015$이다.

| 류 | $\mathbb E[C_h^2/x]$ | $\mathbb E[S_2]$ | 비 |
|:--|--:|--:|--:|
| $h$ 홀수 | 0.3283 | 0.3285 | 0.999 |
| $4\mid h$ | 0.4944 | 0.4928 | 1.003 |
| $9\mid h$ | 0.4145 | 0.4158 | 0.997 |
| $3\parallel h$ | 0.3639 | 0.3638 | 1.000 |

$4\mid h$이면 분산이 1.5배, $9\mid h$이면 8/7배이고, $3\parallel h$는 영향이 없다. 분산은 $p\mid h$가 아니라 $p^2\mid h$에만 반응한다. (1.5)를 정규화한 예측과 실측 $\mathbb E[z_{ak}z_{bk}]$는 다음과 같다.

| $(a,b)$ | (1,2) | (1,3) | (2,3) | (1,5) | (1,4) | (2,9) |
|:--|--:|--:|--:|--:|--:|--:|
| 실측 | 0.1865 | 0.1816 | 0.0712 | 0.1279 | −0.0002 | −0.0020 |
| 예측 | 0.1854 | 0.1824 | 0.0676 | 0.1383 | 0 | 0 |

(1,5)는 약 2 표준오차 낮다. 기본 실행($N=2^{22}$)에서는 (1,2)가 0.156으로 약 16% 낮게 나오는데, 이 크기 의존성은 아직 설명하지 못했다. 같은 크기에서 §3의 로그 자기상관도 큰 지연에서 OU보다 빨리 줄어들어, $x$가 작을 때 증분 사이에 약한 반상관이 있는 것으로 보인다.

<a id="section-1-4"></a>

### 1.4. 제곱류 표현과 일치의 분류

$n$과 $n+h$가 제곱 인수 없을 때 $g=(n,h)$, $q_h(n)=n(n+h)/g^2$로 둔다. $n=gm$, $h=gk$로 쓰면 $(m,m+k)=1$이고 $q_h(n)=m(m+k)$이므로, $q_h(n)$도 제곱 인수가 없다. 이 값을 쌍 $(n,n+h)$의 **제곱류**라 한다.

<a id="proposition-1-4"></a>

**명제 1.4 (제곱류 표현).** $A_h(q;x)=\#\{n\le x:\mu^2(n)\mu^2(n+h)=1,\ q_h(n)=q\}$라 하자. 합을 제곱 인수 없는 $q,r>1$에 대해 취하고 $q\triangle r=qr/(q,r)^2$로 쓰면

<a id="equation-1-6"></a>

$$
C_h(x)=\sum_qA_h(q;x)\,\mu(q),\qquad
C_h(x)C_{h'}(y)=\sum_qA_h(q;x)A_{h'}(q;y)+\sum_{q\ne r}A_h(q;x)A_{h'}(r;y)\,\mu(q\triangle r)
\tag{1.6}
$$

이다. 무작위 곱셈 함수 $f$에서도 $f(n)f(n+h)=f(q_h(n))$이므로, $\mathbb E_f[C_h(f;x)C_{h'}(f;y)]$는 (1.6)의 첫째 합, 즉 대각항과 같다.

**증명.** $\mu(n)\mu(n+h)=\mu(g)^2\mu(m)\mu(m+k)=\mu(q_h(n))$이고, 제곱 인수 없는 $q,r$에 대해 $\mu(q)\mu(r)=\mu(q\triangle r)$이다. 모형에서는 $\mathbb E_f[f(q)f(r)]=\mathbf 1[q=r]$이다. $\square$

$A_h$에는 부호가 없다. 제곱 인수 없는 수의 배치와 gcd 구조만 담고, 부호는 모두 $\mu(q)$에 남는다. 실제 μ에서는 (1.6)의 둘째 합이 남는다. "제곱류 구조가 공분산을 정한다"는 말은 모형에서는 정리이고, 실제 μ에서는 이 둘째 합의 상쇄에 달려 있다.

<a id="proposition-1-5"></a>

**명제 1.5 (층별 밀도).** 명제 1.2의 층 $P^{[g]}_k$가 지지된 $m$의 밀도를 $\delta^{[g]}(k)$라 하면 다음이 성립한다.

<a id="equation-1-7"></a>

$$
\delta^{[g]}(k)=\prod_{p\mid k}\Big(1-\frac1p\Big)\prod_{\substack{p\mid g\\ p\nmid k}}\Big(1-\frac2p\Big)\prod_{p\nmid gk}\Big(1-\frac2{p^2}\Big),\qquad
S_2(h)=\sum_{g\mid h}\frac{\mu^2(g)}g\,\delta^{[g]}(h/g).
\tag{1.7}
$$

**증명.** 첫 식은 소수별 조건 $p\nmid m(m+k)$($p\mid g$) 또는 $(m,k)=1$($p\mid k$)에서 바로 나온다. 둘째 식에서 $g$에 대한 합은 소수별 곱으로 갈라진다.
- $p\nmid h$이면 양변의 국소 인자가 모두 $1-2/p^2$이다.
- $p\parallel h$이면 $p\nmid g$인 층과 $p\mid g$인 층을 더해 $(1-\tfrac1p)+\tfrac1p(1-\tfrac2p)=1-\tfrac2{p^2}$이다.
- $p^2\mid h$이면 $(1-\tfrac1p)+\tfrac1p(1-\tfrac1p)=1-\tfrac1{p^2}$이다. $\square$

원시 층 $\delta^{[1]}=\delta^\flat$는 $p\mid h$에 반응하지만, 층을 모두 더한 $S_2(h)$는 $p^2\mid h$에만 반응한다. §1.3과 §3에서 본 차이의 정확한 이유이다.

<a id="proposition-1-6"></a>

**명제 1.6 (제곱류 일치의 분류).** $n=gu$, $h=gk$와 $n'=g'v$, $h'=g'\ell$를 원시화라 하고 $q_h(n)=q_{h'}(n')$, 즉 $u(u+k)=v(v+\ell)$라 하자.
- $k=\ell$이면 $u=v$이고 $n'/n=h'/h$이다(팽창족).
- $k>\ell$이면 $v=u+r$로 쓸 때 다음이 성립한다.

<a id="equation-1-8"></a>

$$
u(k-\ell-2r)=r(r+\ell),\qquad 1\le r<\frac{k-\ell}2,\qquad u<\frac{k^2-\ell^2}4,\qquad n<\frac{h^2}{4g}\le\frac{h^2}4.
\tag{1.8}
$$

$k<\ell$이면 $n'$ 쪽에 같은 제한이 생긴다. 따라서 두 시작점이 각각 $h^2/4$, $h'^2/4$보다 크면 팽창족만 남는다.

**증명.** $u(u+k)$는 $u$에 대해 증가하므로 $k=\ell$이면 $u=v$이다. $k>\ell$이면 $v>u$이고, $v=u+r$를 넣으면 첫 식이 나온다. 오른쪽이 양수이므로 $k-\ell-2r\ge1$이고, $u\le r(r+\ell)<\tfrac{k-\ell}2\cdot\tfrac{k+\ell}2$이다. 끝으로 $n=gu<gk^2/4=h^2/(4g)$이다. $\square$

명제 1.3의 $\tfrac12\tau(|h^2-h'^2|)$ 상계와 (1.5)의 $O_{a,b,k}(1)$은 이 분류로 구체화된다. 비선형 일치는 $n<h^2/4$에서만 생긴다.

**[수치]** 다음을 확인했다.
- (1.6)의 점별 등식 $\mu(n)\mu(n+h)=\mu(q_h(n))$: $n\le2000$, $h\le30$에서 성립한다.
- (1.7): $h\le2000$에서 $S_2(h)$와 상대오차 $6\cdot10^{-16}$ 안에서 일치한다.
- (1.8): $\ell<k<60$, $u,v<2000$의 해 4194개가 모두 $u<(k^2-\ell^2)/4$를 만족한다.

<a id="section-2"></a>

## 2. 극 정보의 부재

순수 Chowla 합을 각 점의 극·영점 성분의 곱으로 나누면, 극극·극영점 항은 항등적으로 0이다.

**명제 2.1.** $\mathfrak M_x(\xi)=\widetilde S_x^{-1/2,-1/2}(\xi)$의 극 쪽 폐곡선 적분([전역 공식 (1.13)](global-formula.md#equation-1-13))에서 $\xi=0$의 $k$계 도함수의 피적분함수는, $\alpha=\xi/L$에서 오는 $L^{-k}$배를 빼면

<a id="equation-2-1"></a>

$$
(-i)^k\big(\zeta^{-1/2}\big)^{(k)}\zeta^{-1/2}\ \sim\ (s-1)^{1-k}\qquad(s\to1)
\tag{2.1}
$$

이고 일가 유리형이다. $k=0,1$이면 각각 $1/\zeta$, $-\tfrac i2(1/\zeta)'$로 $s=1$에서 해석적이므로 폐곡선 적분은 정확히 0이다. 따라서 $F\in\mathrm{span}\{1,u\}$인 관측량은 $1/L$의 모든 차수에서 극 성분이 없다.

**증명.** $\zeta^{-1/2}=(s-1)^{1/2}Z(s)^{-1/2}$, $Z=(s-1)\zeta$로 쓰면 $k$계 도함수는 $(s-1)^{1/2-k+j}$ 꼴 항의 합이고, $\zeta^{-1/2}$을 곱하면 지수가 정수가 된다. 가장 낮은 차수가 $1-k$이다. $\square$

$k\ge2$부터 극이 생기므로 $\langle\mathcal K_{-1/2,-1/2},u^2\rangle=-\tfrac14\ne0$이다. 모듈러 $q$ 층으로 확장하면 $L(1,\chi)\ne0$이므로 $1/L(s,\chi)$도 $s=1$에서 해석적이다. 극이 생길 수 있는 유일한 경로는 Siegel 영점이다.

**명제 2.2 (극 범함수).** 다음이 성립한다.

<a id="equation-2-2"></a>

$$
\mathcal K_{-1/2,-1/2}=-\tfrac14\,\mathcal K_{3/2,3/2}'',\qquad
\langle\mathcal K_{-1/2,-1/2},F\rangle=-\frac1\pi\int_0^1\sqrt{u(1-u)}\,F''(u)\,du.
\tag{2.2}
$$

**증명.** [주정리 (4.2)](main-theorem.md#equation-4-2)로 양변의 $u^k$ 모멘트를 비교한다. $k\ge2$에서 $(-\tfrac12)_k=-\tfrac14(\tfrac32)_{k-2}$이고 $\Gamma(k-1)=k!/(k(k-1))$이므로 $(-\tfrac12)_k/\Gamma(k-1)=-\tfrac14k(k-1)(\tfrac32)_{k-2}/k!$이다. $k=0,1$에서는 양변이 0이다. 콤팩트 지지 분포는 모멘트로 정해진다. 둘째 식은 $\tfrac14\mathcal K_{3/2,3/2}=\pi^{-1}\sqrt{u(1-u)}$([주정리 (4.5)](main-theorem.md#equation-4-5) 아래)에서 따른다. $\square$

즉 $M$의 극 핵은 반원 핵의 2계 도함수이고, 관측량의 극 정보는 $F''$의 반원 평균이다. 극 정보가 0이 되도록 하려면 이 함수의 값이 0이 되도록 $F$를 사영하면 된다.

**[수치]** $\log x$ 방향 스펙트럼에서 처음 여섯 영점 높이 $\gamma_k$의 파워를 주변 배경과 비교했다($N=2^{25}$). $M(x)=\sum\mu(n)\cdot1$은 상수 1(극뿐인 함수)과의 상관이므로 극영점 상관의 전형이다.

| 대상 | $\gamma_1,\dots,\gamma_6$에서 파워/배경 |
|:--|:--|
| $M(x)/\sqrt x$ | 4603, 1641, 752, 198, 299, 160 |
| $C_h(x)/\sqrt{S_2(h)x}$, $h\le1000$ 평균 | 0.93, 0.98, 0.99, 0.81, 1.16, 1.02 |

<a id="section-3"></a>

## 3. 원시 성분과 OU 법칙

원시 부분의 밀도는 $\delta^\flat(h)=\prod_{p\nmid h}(1-2/p^2)\prod_{p\mid h}(1-1/p)$이다. 원시 부분은 $p^2\mid h$가 아니라 $p\mid h$에 반응한다.

**[수치]** $N=2^{25}$에서 결과는 다음과 같다.
- **분산과 상관 ($h\le10^5$):** $C^\flat_h/\sqrt{\delta^\flat(h)x}$의 분산은 0.998, 첨도는 +0.038이다. 이동 사이 상관 $\mathbb E[z_kz_{ak}]$는 $a=2,3,4,5$에서 $+0.0075,+0.0032,+0.0026,-0.0054$로 모두 표준오차의 2배 안이다. 전체 $C_h$의 $+0.186$, $+0.182$가 사라진다.
- **$x$ 방향 ($h\le1000$의 궤적):** 이웃한 2진 블록 증분의 상관은 $+0.0003$(표준오차 0.010)이고, 같은 통계가 $M(x)$에서는 $-0.845$이다.
- **로그 자기상관:** $\Delta y=0.5,1,2,3$에서 0.764, 0.589, 0.340, 0.204이다. OU 값 $e^{-\Delta y/2}$는 0.779, 0.607, 0.368, 0.223이고, 표준오차는 약 0.02이다.
- **최댓값:** $\sup_{x\le X}|C^\flat_h|/\sqrt{\delta^\flat X}$의 분포가 같은 격자의 브라운 운동과 일치한다(KS p = 0.45).

**[모형] 알려진 것.** 무작위 곱셈 모형에서는 $C^\flat_h(f;x)=\sum_{n\le x}f(n(n+h))$이다. $f$가 제곱 인수 없는 수에만 지지되므로 $(n,h)>1$인 항은 저절로 사라진다. $n(n+h)$가 $n$에 대해 증가하므로

$$
\mathbb E_f\big[C^\flat_h(f;x)\,C^\flat_h(f;y)\big]=\sum_{n\le\min(x,y)}\mu^2(n(n+h))
$$

이고, 따라서 $\log x$에 대한 OU 공분산 $e^{-s/2}$는 모형에서 정확히 성립한다.

그러나 항들이 독립인 것은 아니다. $h=1$, $n=2,5,6,14$의 제곱류는 $6,30,42,210$이고 그 곱이 $1260^2$이므로 $f(6)f(30)f(42)f(210)=1$이 항상 성립한다. 두 개씩은 직교하지만 4점 관계가 남는다. 이런 관계의 총기여가 작다는 것은 별도의 정리이다. Chinis–Shala는 $P(n)=n(n+h)$를 포함하는 다항식에 대해 다음을 증명했다.
- 4차 모멘트: $\mathbb E_f[C^\flat_h(f;x)^4]=3\delta^\flat(h)^2x^2+O_h(x^{2-\eta_h})$ ($\eta_h>0$).
- 중심극한정리: $C^\flat_h(f;x)/\sqrt{\delta^\flat(h)x}\Rightarrow N(0,1)$.

따라서 추측 3.1의 모형판 가운데 고정된 $h$의 일변량 법칙은 정리이다. 남는 것은 실제 μ로의 이전, 여러 $h$와 $x$에 걸친 결합분포와 과정 수준의 수렴, 그리고 (3.1)의 극값이다. 4차 모멘트만으로는 OU 법칙도 LIL도 나오지 않는다.

**[휴리스틱, RH]** 영점쌍 $(\rho,\rho')$가 $C_h$에 주는 항은 $x^{\rho+\rho'-1}$ 규모이고, RH 아래 크기가 $O(1)$이다. 따라서 유한 개의 영점으로는 $\sqrt x$를 만들 수 없다. $\sqrt x$는 높이가 대략 $x/h$인 영점쌍들의 집단적 중심극한 효과이다. 정류 위상으로 보면 쌍 $(\rho,\rho')$는 $(\gamma+\gamma')/\gamma'\approx h/x$인 $x$에서만 기여하므로, 서로 다른 2진 블록은 서로 다른 쌍의 창에 대응한다. 창이 비율 $h/x$에만 의존한다는 점이 명제 1.1의 팽창 구조와 맞는다.

<a id="conjecture-3-1"></a>

**추측 3.1 [추측].** $C^\flat_h(x)/\sqrt{\delta^\flat(h)x}$는 $\log x$에 대해 $h$마다 독립인 Ornstein–Uhlenbeck 과정처럼 행동한다. 특히 다음이 성립한다.

<a id="equation-3-1"></a>

$$
\limsup_{x\to\infty}\frac{|C_h(x)|}{\sqrt{2S_2(h)\,x\log\log x}}=1.
\tag{3.1}
$$

이 추측은 강한 Chowla 추측의 $x^\varepsilon$을 $\sqrt{\log\log x}$로 바꾼다. 낮은 영점이 지배하는 $M(x)$에 대해 Ng가 제안한 $\sqrt x(\log\log\log x)^{5/4}$ 규모와는 다른 법칙이다. 이 범위에서는 $\log\log x\approx2.85$에 불과해서 (3.1) 자체는 수치로 판별할 수 없다.

<a id="section-4"></a>

## 4. 짧은 구간 사전과 4점 Wick 규칙

<a id="section-4-1"></a>

### 4.1. 정확한 다리

$I=(X,X+Y]$, $Y\le\sqrt X/2$로 두고 다음을 정의한다.

<a id="equation-4-1"></a>

$$
P_I(t)=\sum_{n\in I}\mu(n)n^{-it},\quad
S_I(\alpha)=\sum_{n\in I}\mu(n)e(n\alpha),\quad
C_k(I)=\!\!\sum_{n,\,n+k\in I}\!\!\mu(n)\mu(n+k),\quad
Q_I=\sum_{n\in I}\mu^2(n).
\tag{4.1}
$$

평균 $\langle\cdot\rangle_T$의 가중함수 $w$는 푸리에 변환 $\hat w$가 $[-\tfrac12,\tfrac12]$에서 1이고 $[-1,1]$ 밖에서 0인 것으로 잡는다. 이런 $w$는 반드시 음수 값을 가진다. 음이 아닌 적분 가능한 $w$에서는 $\xi\ne0$일 때 $|\hat w(\xi)|<\hat w(0)$이기 때문이다. 따라서 $\langle\cdot\rangle_T$는 부호가 있는 평균이다. 아래의 등식들은 정확한 항등식이지만, 확률 평균의 모멘트로 읽으려면 따로 근거가 필요하다.

**정리 4.1.** $\tfrac{8\pi}3X(1+Y/X)^2<T\le4\pi X$이면 다음이 정확히 성립한다.

<a id="equation-4-2"></a>

$$
\langle|P_I|^2\rangle_T=Q_I,\qquad
\langle|P_I|^4\rangle_T=\sum_{|k|<Y}C_k(I)^2=\int_0^1|S_I(\alpha)|^4\,d\alpha.
\tag{4.2}
$$

**증명.** $\langle(n_1n_2/n_3n_4)^{-it}\rangle_T=\hat w\big(\tfrac T{2\pi}\log\tfrac{n_1n_2}{n_3n_4}\big)$이고 $n_i=X+m_i$로 쓴다.
- $m_1+m_2=m_3+m_4$이면 $|n_1n_2-n_3n_4|=|m_1m_2-m_3m_4|\le Y^2$이므로 $|\log|\le Y^2/X^2$이다. $Y\le\sqrt X/2$이므로 $T\le4\pi X\le\pi X^2/Y^2$이고, 따라서 $\hat w=1$이다.
- 그렇지 않으면 $|n_1n_2-n_3n_4|\ge X-Y^2\ge3X/4$이고 $\max(n_1n_2,n_3n_4)\le(X+Y)^2$이므로 $|\log|\ge\tfrac{3X}{4(X+Y)^2}$이다. 가정한 $T$의 하한에서 $\tfrac T{2\pi}|\log|>1$이고 $\hat w=0$이다.

남는 것은 덧셈 4쌍 전체이고, 그 합은 $\sum_kC_k(I)^2$이다. 마지막 등식은 Parseval이다. $|P_I|^2$도 같은 방식으로 $Q_I$가 된다. $\square$

더 일반적으로 다음이 성립한다. 첫째 식은 $\log((X+m)/(X+m'))\approx(m-m')/X$에 따른 근사이다.

<a id="equation-4-3"></a>

$$
\big\langle|P_I(t)|^2|P_I(t+\tau)|^2\big\rangle_T\approx\sum_kC_k(I)^2e^{-i\tau k/X},\qquad
\int_0^1|S_I(\alpha)|^2|S_J(\alpha+\theta)|^2d\alpha=\sum_kC_k(I)C_k(J)e(k\theta)
\tag{4.3}
$$

둘째 식은 정확하다. 즉 $C_k(I)$는 $|P_I(t)|^2$의 진동수 $k/X$ 성분이다.

**따름정리 4.2.** $\sum_{k\ne0}C_k(I)^2\sim Q_I^2$은 $\langle|P_I|^4\rangle\sim2\langle|P_I|^2\rangle^2$(복소 가우스의 4차 모멘트 값)와 동치이다. $\mathbb E_kS_2(k)=(6/\pi^2)^2$이므로, 왼쪽은 각 $C_k(I)$가 Mirsky 진폭 $\sqrt{S_2(k)(Y-|k|)}$를 가질 때의 크기와 같다. 이는 짧은 구간 안에서 이동에 대해 평균한 진술이다. 고정된 $h$의 $C_h(x)\ll x^{1/2+\varepsilon}$와 동치는 아니며, 이동 제곱평균만으로는 개별 이동의 제곱근 상계가 나오지 않는다.

<a id="section-4-2"></a>

### 4.2. 뫼비우스 언어: 가설 M4

$\widehat{\mu^2}_I(\theta)=\sum_{n\in I}\mu^2(n)e(n\theta)$로 둔다.

**추측 4.3 (M4) [추측].** $Y\to\infty$이고 $Y\le\sqrt X/2$인 모든 구간 $I$, 그리고 $I$와 겹치지 않는 같은 꼴의 구간 $J$에 대해, $\theta$에 균등하게

<a id="equation-4-4"></a>

$$
\int_0^1|S_I(\alpha)|^2|S_I(\alpha+\theta)|^2d\alpha=Q_I^2+|\widehat{\mu^2}_I(\theta)|^2-Q_I+o(Q_I^2),
\tag{4.4}
$$

<a id="equation-4-5"></a>

$$
\int_0^1|S_I(\alpha)|^2|S_J(\alpha+\theta)|^2d\alpha=Q_IQ_J+o(Q_IQ_J).
\tag{4.5}
$$

<a id="proposition-4-4"></a>

**명제 4.4 [모형].** μ를 다음 두 모형의 함수로 바꾸면 (4.4), (4.5)의 좌변의 기댓값은 $o(\cdot)$ 없이 우변과 같다.
1. **독립 부호:** $\varepsilon_n\mu^2(n)$. 여기서 $\varepsilon_n$은 $n$마다 독립인 ±1 부호이다.
2. **소수별 무작위 곱셈 부호:** §1.2의 $f$. 단 $I=(X,X+Y]$와 $J=(X',X'+Y']$가 $Y^2\le X/4$, $Y'^2\le X'/4$를 만족한다.

**증명.** (1) 네 첨자의 짝짓기를 센다. 같은 구간 안의 짝짓기 $\{12\}\{34\}$는 $Q_I^2$, 엇갈린 짝짓기 $\{14\}\{23\}$는 $|\widehat{\mu^2}_I(\theta)|^2$를 준다. 짝짓기 $\{13\}\{24\}$는 $\alpha$ 적분에서 네 항이 모두 같은 경우만 남기고, 네 항이 모두 같은 경우의 중복을 빼면 $-Q_I$가 된다. $I\cap J=\emptyset$이면 엇갈린 짝짓기가 없다.

(2) 좌변을 전개하면 $a-b=d-c=k$인 네 수 $a,b,c,d$에 대한 $\mathbb E_f[f(a)f(b)f(c)f(d)]$의 합이다. 여기서 $a,b$는 $I$, $c,d$는 $I$ 또는 $J$에 있다. $k>0$이면 명제 1.4에 의해 이 기댓값은 $\mathbf 1[q_k(b)=q_k(c)]$이다. $k<0$이면 두 쌍에서 작은 수를 기준으로 같은 논증을 쓴다.
- $k=0$이면 (1)과 같은 짝 $a=b$, $c=d$이다.
- $k\ne0$이고 $q_k(b)=q_k(c)$일 때, 명제 1.6의 팽창족이면 이동이 같으므로 $b=c$이다. 비선형 일치이면 $b$나 $c$가 $k^2/4$보다 작아야 한다. 그런데 $|k|$는 두 구간의 길이보다 작으므로 $k^2/4<X/16$이고, 이런 수는 구간 안에 있을 수 없다.

따라서 같은 구간에서는 엇갈린 짝 $b=c$, $a=d$만 남고, 서로 다른 구간에서는 $k\ne0$인 항이 없다. 남은 항을 세면 (1)과 같은 값이 나온다. $\square$

길이 조건은 필요하다. $(0,15]$에서 $6+10=1+15$이고 $6\cdot10\cdot1\cdot15=30^2$이므로, $f(6)f(10)f(1)f(15)=1$인 짝 아닌 항이 생긴다.

**정확한 나머지.** 확률 모형 없이 짝짓기 항을 분리하면, $D_k(I)=\sum_{n,n+k\in I}\mu^2(n)\mu^2(n+k)$에 대해

$$
\int_0^1|S_I(\alpha)|^2|S_I(\alpha+\theta)|^2d\alpha=Q_I^2+|\widehat{\mu^2}_I(\theta)|^2-Q_I+R_I(\theta),\qquad
R_I(\theta)=\sum_{0<|k|<Y}\big(C_k(I)^2-D_k(I)\big)e(k\theta)
$$

이다. $R_I(\theta)$는 $a+c=b+d$, $\{a,c\}\ne\{b,d\}$인 짝 아닌 4쌍의 부호합 $\sum\mu(a)\mu(b)\mu(c)\mu(d)e((c-d)\theta)$와 같다. 따라서 (4.4)의 내용은 정확히 $\sup_\theta|R_I(\theta)|=o(Q_I^2)$이고, (4.5)의 내용은 $\sup_\theta|\sum_{k\ne0}C_k(I)C_k(J)e(k\theta)|=o(Q_IQ_J)$이다. $Y^2\le X/4$이면 명제 4.4 (2)의 증명에 의해 이 4쌍들의 곱 $abcd$는 제곱수가 아니다. 무작위 모형에서는 각 항의 기댓값이 0이지만, 실제 μ에서는 각 항이 여전히 ±1이다. $Q_I\asymp Y$이면 짝 아닌 4쌍은 $O(Y^3)$개이고, M4는 그 부호합이 $o(Y^2)$라는 요구이다.

따라서 M4는 **4점 덧셈 통계 수준에서 μ는 μ² 위의 독립 부호**라는 진술이다. (4.3)에 따라 각 부분의 뜻은 다음과 같다.
- **(4.4), $\theta=0$:** 따름정리 4.2와 같다.
- **(4.4), $\theta=a/d^2$:** 제곱류 메아리이다. μ²의 국소 푸리에 계수가 $-1/(p^2-1)$이므로 $|\widehat{\mu^2}_I(a/d^2)|^2/Q_I^2\approx\prod_{p\mid d}(p^2-1)^{-2}$이고, 이것이 §1의 $S_2(k)$ 구조와 같다.
- **(4.5):** 서로 다른 구간끼리 Chowla 합이 무상관이라는 것, 즉 블록 독립(이동 평균)이다.
- **팽창:** $S_{pI}(\alpha)$는 $n=pm$ 항들을 통해 $-S^{(p)}_I(p\alpha)$를 포함한다. 진동수 쪽에서 팽창은 $\alpha\mapsto p\alpha$이다.

**[수치]** $N=2^{25}$, $Y=1581$, 구간 14898개에서 결과는 다음과 같다. 왼쪽 값은 $\Sigma(\theta)=\sum_{k\ne0}C_k^2e(k\theta)/Q^2$의 평균, 오른쪽 값은 $(|\widehat{\mu^2}_I|^2-Q)/Q^2$의 평균이다. 둘의 차이가 $R_I(\theta)/Q_I^2$이다.

| θ | 0 | 1/4 | 1/2 | 1/9 | 1/3 | 1/25 | 1/5 | 1/8 | 0.318 |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| $\Sigma(\theta)$ | 0.99857 | 0.10976 | 0.11006 | 0.01433 | 0.01433 | 0.00070 | 0.00108 | −0.00087 | −0.00107 |
| μ² 공식 | 0.99896 | 0.11009 | 0.11009 | 0.01462 | 0.01460 | 0.00075 | 0.00073 | −0.00101 | −0.00098 |

- 499개 θ의 조밀한 격자에서 차이의 rms z는 0.93으로, 순수한 잡음 수준이다.
- 무관한 θ에서 보이는 −0.001은 중복을 빼는 항 $-1/Q_I$이다.
- $1/5=5/25$와 $1/6$도 각각 $d=5$, $d=6$의 메아리이다.
- (4.5)에 해당하는 값은 구간 간격 1, 2, 5에서 $+0.0001$, $+0.0009$, $+0.0002$이다(표준오차 0.0004).
- 디리클레 다항식 쪽: $\langle|P|^2\rangle/Q=1.000$, $\langle|P|^4\rangle/Q^2=2.024$이고, 높이 차 $2\pi X/4$의 메아리는 $+0.112$이다(예측 약 0.11).

현재 기술과의 격차도 적어 둔다. 짧은 구간의 푸리에 균등성(Matomäki–Radziwiłł–Tao)으로 얻는 것은 거의 모든 구간에서 $\int|S_I|^4=o(Y^3)$ 정도이다. M4는 이를 $2Q_I^2\asymp Y^2$, 즉 제곱근 상쇄 수준으로 요구한다.

<a id="section-4-3"></a>

### 4.3. 영점 쪽 번역 [휴리스틱, RH]

RH와 단순 영점 아래 명시공식은 다음을 준다.

<a id="equation-4-6"></a>

$$
P_I(t)\approx\frac Y{\sqrt X}\sum_\rho\frac{X^{i(\gamma-t)}K\big((\gamma-t)Y/X\big)}{\zeta'(\rho)},\qquad K(u)=\frac{e^{iu}-1}{iu}.
\tag{4.6}
$$

즉 $P_I(t)$는 높이 $t\asymp X$ 근처, 폭 $X/Y\approx2\sqrt X$의 중간 규모 창에 있는 영점 약 $\sqrt X\log X/\pi$개의 합이다. 각 영점은 $1/\zeta'(\rho)$로 가중되고 위상 $X^{i\gamma}$를 가진다. 따라서 Montgomery 매개변수는 $\alpha=\log X/\log t\to1$이다.
- **2차 모멘트의 대각항:** $\int|K|^2=2\pi$와 Gonek 상수 $J_{-1}(T)\sim3T/\pi^3$([Ng 쌍대 (5.6)](ng-duality.md#equation-5-6))을 쓰면, 대각항 $\rho=\rho'$만으로 $\frac{Y^2}{X}\cdot\frac3{\pi^3}\cdot\frac XY\cdot2\pi=\frac6{\pi^2}Y$가 나온다. 이는 (4.2)의 $Q_I$와 맞는다.
- **가설 W4:** 이 창에서 $1/\zeta'$ 가중 영점합의 4점 상관이 Wick 규칙을 따른다. 즉 산술 성분을 빼면 연결된 4점 성분이 사라진다. (4.6)과 정리 4.1을 합치면 W4는 따름정리 4.2의 영점 쪽 표현이다.
- **블록 독립:** 위치가 다른 구간끼리는 위상 $X_j^{i(\gamma-\gamma')}$가 쌍 차이 $\gamma-\gamma'$에 대한 적분에서 상쇄되어 무상관이 된다.
- **높이 메아리:** (4.4)의 $\theta=a/d^2$ 성분은, 높이가 $2\pi X\cdot a/d^2$만큼 떨어진 영점 창들이 진폭 수준에서 $\prod_{p\mid d}(p^2-1)^{-2}$의 세기로 상관된다는 예측이 된다. 실제 영점 자료로 반증할 수 있는 형태이다.
- **경계의 위치:** $\alpha\to1$은 소수에서 Goldston–Montgomery가 짧은 구간의 분산과 쌍상관을 잇는 경계와 같다. Farmer–Gonek–Lee–Lester는 이를 여러 $\zeta'/\zeta$의 곱과 고차 상관으로 확장했다. 여기의 W4는 그 Möbius판($1/\zeta'$ 가중)의 4점 경우이다. 가중 2점 상관($J_{-1}(T)\asymp T$)조차 증명되지 않았으므로 W4는 현재 기술 밖에 있다.

<a id="section-4-4"></a>

### 4.4. 2k차 사다리와 순열 규칙

정리 4.1은 $T\asymp X$에서 4차 모멘트를 덧셈 에너지로 바꾸었다. 2k차로 올리면 $T$의 규모마다 다른 산술 통계가 나타난다. $n_i=X+m_i$, $n_i'=X+m_i'$로 쓰고 $\Delta_j=\sum_{i\le k}(m_i^j-m_i'^j)$라 하면

<a id="equation-4-7"></a>

$$
\lambda:=\log\frac{n_1\cdots n_k}{n_1'\cdots n_k'}=\sum_{j\ge1}\frac{(-1)^{j-1}}{jX^j}\,\Delta_j
\tag{4.7}
$$

이다. $T\asymp X^r$의 평균은 $\Delta_1,\dots,\Delta_r$까지만 본다. $1\le r\le k$에 대해

<a id="equation-4-8"></a>

$$
V_{k,r}(I)=\sum_{\substack{n,n'\in I^k\\ \sum_in_i^j=\sum_in_i'^j\ (j\le r)}}\prod_{i=1}^k\mu(n_i)\mu(n_i')
=\int_{[0,1]^r}\Big|\sum_{n\in I}\mu(n)\,e(\alpha_1n+\cdots+\alpha_rn^r)\Big|^{2k}d\alpha
\tag{4.8}
$$

로 둔다. 거듭제곱합 조건은 평행이동에 불변이므로 $n$으로 써도 $m$으로 써도 같다. $V_{k,1}=\int_0^1|S_I|^{2k}d\alpha$는 $k$겹 덧셈 에너지이고, $r\ge2$이면 μ로 가중한 Vinogradov 평균값이다. $r=k$이면 Newton 항등식에 의해 두 다중집합이 같아지므로, $V_{k,k}$는 $I$ 안의 제곱 인수 없는 수로 이루어진 재배열 쌍 $(n,n')$의 개수

<a id="equation-4-9"></a>

$$
D_k(Q)=(k!)^2\,[z^k]\Big(\sum_{m=0}^k\frac{z^m}{(m!)^2}\Big)^{Q},\qquad Q=Q_I
\tag{4.9}
$$

이다. 다중집합마다 배열 수의 제곱을 더한 값이며, $D_2=2Q^2-Q$, $D_3=6Q^3-9Q^2+4Q$, $D_4=24Q^4-72Q^3+82Q^2-33Q$, 일반적으로 $D_k=k!\,Q^k(1+O_k(Q^{-1}))$이다.

<a id="theorem-4-5"></a>

**정리 4.5 (모멘트 사다리).** $k\ge2$, $1\le Y<X/2$로 두고, $\eta_r=\dfrac{krY^{r+1}}{(r+1)(X-Y)}$이라 하자. 평균 $\langle\cdot\rangle_T$는 §4.1과 같고, $A$는 $I$ 안의 제곱 인수 없는 수 전체이다.

1. $1\le r<k$, $\eta_r<\tfrac13$이고 $\dfrac{2\pi rX^r}{1-\eta_r}<T\le\dfrac{\pi rX^r}{\eta_r}$이면 $\langle|P_I|^{2k}\rangle_T=V_{k,r}(I)$이다.
2. $T>2\pi(X+Y)^k$이면 $\langle|P_I|^{2k}\rangle_T=\#\{(n,n')\in A^{2k}:n_1\cdots n_k=n_1'\cdots n_k'\}$이다. 이는 Steinhaus 무작위 곱셈 함수 $f$에 대한 $\mathbb E\,|\sum_{n\in I}\mu(n)f(n)|^{2k}$와 같다. $Y^k\le X$이면 정확히 $D_k(Q_I)$이고, 일반적인 $Y$에서는 명제 4.8의 오차 안에서 $D_k(Q_I)$이다.
3. $a\ne b$, $X\ge e\,(3/2)^{\max(a,b)}$, $T>2\pi$이면 $\langle P_I^{\,a}\,\overline{P_I}^{\,b}\rangle_T=0$이다.

**증명.** (1) $|\Delta_j|<kY^j$이다. $\Delta_1=\cdots=\Delta_r=0$이면 (4.7)에서

$$
|\lambda|<\sum_{j>r}\frac{kY^j}{jX^j}\le\frac{kY^{r+1}}{(r+1)X^r(X-Y)}=\frac{\eta_r}{rX^r}
$$

이므로 $T|\lambda|/2\pi\le\tfrac12$이고 $\hat w=1$이다. 그렇지 않으면 $\Delta_s\ne0$인 가장 작은 $s\le r$에서 같은 꼬리 추정으로 $|\lambda|\ge(1-\eta_s)/(sX^s)\ge(1-\eta_r)/(rX^r)$이다. $\eta_s$와 $sX^s$가 $s$에 대해 증가하기 때문이다. 따라서 $T|\lambda|/2\pi>1$이고 $\hat w=0$이다.

(2) 곱이 다르면 $|\lambda|\ge(X+Y)^{-k}$이므로 $\hat w=0$이다. 곱이 같은 제곱 인수 없는 수들에서는 $\prod\mu(n_i)=(-1)^{\Omega(\prod n_i)}=\prod\mu(n_i')$이므로 가중치가 $1$이다. Steinhaus 함수에서는 $\mathbb E[\prod f(n_i)\overline{f(n_i')}]=\mathbf 1[\prod n_i=\prod n_i']$이다. 곱이 같고 $n_1$이 어느 $n_j'$와도 다르다고 하자. 그러면 $\gcd(n_1,n_j')$는 $0<|n_1-n_j'|<Y$를 나누고, $n_1$은 제곱 인수가 없으므로 $\prod_j\gcd(n_1,n_j')$를 나눈다. 따라서 $n_1<Y^k\le X$가 되어 모순이다. $n_1$을 약분하고 귀납하면 재배열만 남는다.

(3) $a>b$이면 $\lambda\ge(a-b)\log X-b\log(1+Y/X)\ge\log X-b\log\tfrac32\ge1$이고, $a<b$이면 부호만 바뀐다. $\square$

$k=2$, $r=1$에서 $\eta_1=Y^2/(X-Y)$이므로 (1)은 정리 4.1과 같은 범위($T$가 대략 $8\pi X/3$에서 $4\pi X$까지)를 준다. 단 $r$에 오르려면 $Y^{r+1}\ll_kX$가 필요하다. $Y\asymp\sqrt X$이면 바닥($r=1$)과 꼭대기만 쓸 수 있고, $Y^k\ll_kX$이면 사다리 전체가 정확하다.

사다리의 바닥은 덧셈 통계, 꼭대기는 곱셈 통계이다. $T\to\infty$에서 $(p^{-it})_p$가 독립 균등 위상처럼 분포한다는 Bohr의 관찰이 꼭대기의 Steinhaus 모형이다. 짧은 구간에서는 이 모형의 모멘트가 대각값 $D_k$이다. 바닥의 $V_{k,1}$은 따름정리 4.2의 2k차판이다.

<a id="proposition-4-6"></a>

**명제 4.6 [모형].** $\mu(n)$을 $\varepsilon_n\mu^2(n)$으로 바꾸자. $\varepsilon_n$은 $n$마다 독립인 ±1 부호이다.

1. $k=2,3$이면 모든 $r$에서 $\mathbb E\,V_{k,r}=D_k(Q)$이다. $k=4$이면 $\mathbb E\,V_{4,1}=D_4(Q)+9E^*(A)+6\,\mathrm{AP}_3(A)$이고, $r\ge2$에서는 $\mathbb E\,V_{4,r}=D_4(Q)$이다. 여기서 $E^*(A)$는 서로 다른 $a,b,c,d\in A$로 $a+b=c+d$인 순서쌍의 수, $\mathrm{AP}_3(A)$는 $b\ne c$, $b+c=2a$인 $(a,b,c)\in A^3$의 수이다. 둘 다 $O(Q^3)$이다.
2. $\theta_1,\theta_2,\theta_3$와 $G_{jl}=\widehat{\mu^2}_I(\theta_j-\theta_l)$에 대해

<a id="equation-4-10"></a>

$$
\mathbb E\int_0^1\prod_{j=1}^3|S_I(\alpha+\theta_j)|^2d\alpha=\operatorname{perm}G-\sum_{j,l}|G_{jl}|^2+4Q_I
\tag{4.10}
$$

이다. $k=2$에서 같은 계산을 하면 $\operatorname{perm}G-Q_I$, 즉 (4.4)의 오른쪽이 나온다.

**증명.** $\mathbb E\prod\varepsilon=1$인 경우는 각 값이 짝수 번 나올 때뿐이다. 값 $v$가 모두 $c_v$번($c_v$ 짝수), 그중 왼쪽에 $x_v$번 나오면 합 조건은 $\sum_v(2x_v-c_v)v=0$이다. (1) $k\le3$이면 경우를 나누어 보면 이 조건의 해는 $2x_v=c_v$뿐이다. 따라서 양쪽이 재배열이고, 그 수가 $D_k$이다. $k=4$에서는 두 꼴의 해가 새로 생긴다. $\{a,a,b,b\}$ 대 $\{c,c,d,d\}$($a+b=c+d$)와 $\{b,b,c,c\}$ 대 $\{a,a,a,a\}$($b+c=2a$)이고, 순서쌍으로 세면 각각 $9E^*$개, $6\,\mathrm{AP}_3$개이다. 두 꼴 모두 제곱합 조건을 어기므로 $r\ge2$에서는 사라진다. (2) 재배열에 대해 합할 때 안정자의 크기로 나누는 가중치는 $1-\tfrac12\sum_{i<j}\mathbf 1[n_i=n_j]+\tfrac23\mathbf 1[n_1=n_2=n_3]$이다. 각 순열 $\sigma$에서 $\sum_j(\theta_j-\theta_{\sigma(j)})=0$을 쓰면 (4.10)이 나온다. $\square$

<a id="conjecture-4-7"></a>

**추측 4.7 (M2k) [추측].** $k$를 고정하고 $Y\to\infty$, $Y\le\sqrt X$인 짧은 구간을 생각하자.

- **(a) 순열 규칙.** 모든 $\theta_1,\dots,\theta_k$에 대해 다음이 성립한다.

<a id="equation-4-11"></a>

$$
\int_0^1\prod_{j=1}^k|S_I(\alpha+\theta_j)|^2d\alpha=\operatorname{perm}\Big[\widehat{\mu^2}_I(\theta_j-\theta_l)\Big]_{j,l\le k}+o(Q_I^k).
\tag{4.11}
$$

- **(b) 평평한 사다리.** 모든 $1\le r\le k$에서 $V_{k,r}(I)=D_k(Q_I)+o(Q_I^k)$이다.

$k=2$의 (a)는 M4의 (4.4)이다. 서로 다른 구간 $I_j$에 대해 $G_{jl}=\sum_{n\in I_j\cap I_l}\mu^2(n)e(n(\theta_j-\theta_l))$로 쓰면 (4.5)도 포함한다. 복소 가우스 벡터 $Z$에서는 $\mathbb E\prod_j|Z_j|^2=\operatorname{perm}(\operatorname{Cov}Z)$이다. 따라서 (a)는 $\alpha\mapsto S_I(\alpha)$가 공분산 함수 $\widehat{\mu^2}_I$를 갖는 정상 복소 가우스 과정처럼 행동한다는 진술이고, 4점에서 멈추었던 "μ² 위의 독립 부호"를 모든 짝수 차수로 올린 것이다.

정리 4.5의 범위에서 (b)는 $\langle|P_I|^{2k}\rangle_T$가 $T\asymp X$부터 꼭대기까지 $T$에 거의 무관하다는 뜻이다. 꼭대기는 Steinhaus 모형이므로, 높이 $t\asymp X$에서 이미 $P_I(t)$가 $t\to\infty$의 모멘트를 갖는다는 말이 된다. $Y=o(\sqrt X)$에서 모든 $k$에 대해 (b)의 $r=1$이 성립한다고 하자. 그러면 정리 4.5 (1), (3)에 의해 $T\asymp X$에서 부호 있는 평균 $\langle P_I^{\,a}\overline{P_I}^{\,b}\rangle_T/Q_I^{(a+b)/2}$가 표준 복소 가우스의 모멘트 값($a=b$이면 $a!$, 아니면 $0$)으로 수렴한다. $\langle\cdot\rangle_T$는 확률 평균이 아니므로(§4.1), 이를 $P_I(t)$의 분포 수렴으로 읽으려면 양의 가중함수로 사다리를 다시 적어야 한다.

<a id="proposition-4-8"></a>

**명제 4.8 (꼭대기의 점근적 대각성).** $k\ge2$, $1\le Y\le X$, $\varepsilon>0$이고 $A$가 $I$ 안의 제곱 인수 없는 수 전체라 하자. 그러면

<a id="equation-4-12"></a>

$$
0\le\#\{(n,n')\in A^{2k}:n_1\cdots n_k=n_1'\cdots n_k'\}-D_k(Q_I)\ll_{k,\varepsilon}X^\varepsilon\,\frac{Y^{k+1}}X
\tag{4.12}
$$

이다. 따라서 정리 4.5 (2)의 꼭대기는 $Q_I\asymp Y$이면 상대오차 $O(X^\varepsilon Y/X)$로 $D_k(Q_I)$이고, $Y\le X^{1-\delta}$에서 점근적으로 대각적이다.

**증명.** 재배열이 아닌 쌍에서 공통 원소를 다중집합으로 약분하면, 크기 $k'\ge2$이고 공통 원소가 없으며 곱이 같은 쌍 $(u,u')$가 남는다. 약분한 부분을 고르는 방법은 $O_k(Q_I^{k-k'})$가지이다.

$u$를 고정하고, $\prod u_l$의 약수 중 $I$ 안에 있고 $u$에 없는 제곱 인수 없는 수의 개수를 $N(u)$라 하자. $u'$의 원소는 모두 이런 약수이다. 그런 $d$ 하나를 잡고 소수를 나누면 $d=h_1\cdots h_{k'}$, $h_l\mid u_l$로 쓸 수 있다. $h_l$은 0이 아닌 $u_l-d$를 나누므로 $h_l<Y$이다. 따라서

$$
\#\{u\in I^{k'}:d\mid\textstyle\prod u_l,\ d\notin u\}\le\sum_{h_1\cdots h_{k'}=d}\prod_l\Big(\frac Y{h_l}+1\Big)\le(2Y)^{k'}\frac{\tau_{k'}(d)}d
$$

이고, $d\in I$에 대해 더하면 $\sum_uN(u)\ll X^\varepsilon Y^{k'+1}/X$이다. $u'$은 $N(u)$개의 후보에서 $k'-1$개를 고르면 정해지고, $N(u)\le\tau(\prod u_l)\ll X^\varepsilon$이다. 따라서 크기 $k'$의 쌍은 $\ll X^{2\varepsilon}Y^{k'+1}/X$개이고, $k'$에 대해 더하면 (4.12)가 나온다. $\square$

정리 4.5 (2)의 gcd 논증은 점마다 보면 $Y^k\le X$에서 끊기지만, $u$에 대해 평균하면 모든 $Y$에서 통한다. 평균적으로 $k$겹 곱이 $I$ 안에 '남의 약수'를 갖는 경우가 $X^\varepsilon Y/X$개뿐이기 때문이다. 정확한 대각성은 $Y\asymp\sqrt X$에서 깨진다. $X=10^5$, $Y=632=\lfloor2\sqrt X\rfloor$에서

$$
100061\cdot100317\cdot100558=100147\cdot100191\cdot100598
$$

이고, 여섯 수는 모두 제곱 인수가 없으며 공통 원소가 없다. gcd 행렬 $[\gcd(n_i,n_j')]=\begin{pmatrix}43&13&179\\17&21&281\\137&367&2\end{pmatrix}$의 행의 곱이 왼쪽 수, 열의 곱이 오른쪽 수이고, 모든 성분이 $Y$보다 작다.

**[수치]** 정리 4.5와 명제 4.6은 다음과 같이 확인했다.
- **사다리의 이분법:** $(k,r,X,Y)=(2,1,4000,31)$, $(3,1,4000,29)$, $(4,1,2\cdot10^4,45)$, $(3,2,10^6,50)$, $(4,2,10^7,60)$, $(4,3,10^9,25)$에서 모든 다중집합 쌍을 정확한 정수로 계산했다. 남는 쌍의 $|\lambda|X^r$은 최대 $0.065$이고, 지워지는 쌍은 최소 $0.930$이다.
- **꼭대기:** $X=10^6$, $Y=500$에서는 $Y^3>X$이지만 곱이 같은 비자명한 세 쌍이 없다. 제곱 인수 없는 세 쌍 전체를 비교하면 다음과 같다.
  - $Y=\lfloor2\sqrt X\rfloor$: $X=10^5$에서 2개, $X=2.5\cdot10^5$에서 3개의 비자명한 일치가 있다. 모두 공통 원소가 없다.
  - $Y=\sqrt X$: $X=10^6$에서 일치가 없다.
  - 일치의 수는 $D_3(Q)=3.5\cdot10^8$, $1.3\cdot10^9$에 비해 무시할 만하고, (4.12)의 규모 $Y^4/X$보다도 훨씬 작다.
- **명제 4.6:** 작은 집합의 전수 열거로 확인했다.

$N=2^{25}$에서 §4.2와 같은 구간 14898개($Y=1581$)를 쓴다. $\int|S_I|^{2k}$를 명제 4.6의 정확한 모형값으로 나눈 비의 평균(표준오차)과 구간별 표준편차는 다음과 같다.

| $k$ | μ | 독립 부호 | 무작위 곱셈 | 구간별 표준편차 (μ / 독립 / 곱셈) |
|:--|--:|--:|--:|:--|
| 2 | 0.99981 (0.00024) | 1.00023 | 1.00002 | 0.030 / 0.030 / 0.029 |
| 3 | 0.99942 (0.00076) | 1.00071 | 1.00030 | 0.093 / 0.093 / 0.092 |
| 4 | 0.99884 (0.00170) | 1.00118 | 1.00112 | 0.208 / 0.206 / 0.208 |

평균뿐 아니라 구간 사이의 흩어짐도 두 모형과 같다. 이동한 6차 모멘트는 (4.10)과 비교했다. 값은 $Q_I^3$ 단위이다.

| $(\theta_1,\theta_2,\theta_3)$ | $(0,0,0)$ | $(0,0,\tfrac14)$ | $(0,\tfrac14,\tfrac12)$ | $(0,\tfrac19,\tfrac29)$ | $(0,\tfrac14,\tfrac19)$ | 일반 |
|:--|--:|--:|--:|--:|--:|--:|
| 모형값 | 5.9906 | 2.4388 | 1.2555 | 1.0398 | 1.1286 | 0.9974 |
| μ − 모형값 | −0.0035 | −0.0022 | −0.0009 | −0.0004 | −0.0009 | −0.0001 |
| 표준오차 | 0.0046 | 0.0017 | 0.0009 | 0.0006 | 0.0006 | 0.0004 |

- 같은 자료에서 계산했으므로 열들은 서로 상관되어 있다. 독립 부호와 무작위 곱셈 모형도 모든 열에서 표준오차 2배 안에 있다.
- $(0,\tfrac14,\tfrac12)$의 모형값 $1.2555$는 항등·전치 항 $1.3334\approx1+3\cdot\tfrac19$, 두 3-순환 항 $2\Re\,\widehat{\mu^2}_I(-\tfrac14)^2\widehat{\mu^2}_I(\tfrac12)/Q_I^3=-0.0741\approx-\tfrac2{27}$, (4.10)의 낮은 차수 항 $-0.0038$의 합이다. 3-순환 항의 부호는 μ²의 mod 4 구조가 정하는 위상에서 나오며, μ는 이 위상까지 따른다.
- "일반"은 $(0,\,0.3183\ldots,\,0.1234\ldots)$이다.

$k=3$의 사다리 전체는 $Y=80$에서 본다($\eta_2=0.10$, $Y^3\le X$). 1단은 $X\ge10^7$의 모든 구간 294430개, 2단은 그중 21031개에서 계산했다. 모든 값은 $V_{3,r}/D_3$이고, 2단은 $(V_{3,2}-D_3)/D_3$이다.

| 가중치 | 1단 ($T\asymp X$) | 2단 ($T\asymp X^2$) | 3단 ($T\asymp X^3$) |
|:--|--:|--:|--:|
| μ | 0.9989 (0.0007) | $-0.8\times10^{-5}$ ($3.6\times10^{-5}$) | 1 |
| 독립 부호 | 1.0015 (0.0007) | $-2.3\times10^{-5}$ ($3.6\times10^{-5}$) | 1 |
| μ² | 137.15 | $+0.262$ | 1 |

- μ² 가중치에서는 사다리를 내려갈수록 비대각 해가 늘어난다. 2단에서는 이차 Vinogradov 계의 비자명 해가 $26\%$를 더하고, 1단에서는 137배가 된다.
- μ 가중치에서는 모든 단이 대각값이다.
- 1단을 다른 길이의 모든 구간에서 다시 계산하면, $Y=20$에서 μ $1.0002$, 독립 부호 $1.0003$(표준오차 0.0006)이고, $Y=200$에서 $0.9997$, $1.0006$(0.0007)이다.

<a id="section-5"></a>

## 5. 이산핵과 코사인 언어

<a id="section-5-1"></a>

### 5.1. 부호 동전과 모양 동전

**명제 5.1.** 제곱 인수 없는 $n$에서는 모든 $d\mid n$에 대해 $d_{-1/2}(d)d_{-1/2}(n/d)=(-\tfrac12)^{\omega(n)}$이다. 따라서

<a id="equation-5-1"></a>

$$
\mathscr K_n^{-1/2,-1/2}=\mu(n)\,\nu_n,\qquad
\langle\mathscr K_n,e^{i\xi u}\rangle=\mu(n)\,\hat\nu_n(\xi),\quad
\hat\nu_n(\xi)=e^{i\xi/2}\prod_{p\mid n}\cos\frac{\xi u_p}2,\quad u_p=\frac{\log p}{\log n}.
\tag{5.1}
$$

여기서 $\nu_n$은 DDT의 균등 약수측도이다. $p^2\parallel n$이면 $p$의 동전은

<a id="equation-5-2"></a>

$$
\sum_{j=0}^2d_{-1/2}(p^j)d_{-1/2}(p^{2-j})e^{ij\varphi}=-\tfrac18\big(1-e^{i\varphi}\big)^2=\tfrac12e^{i\varphi}\sin^2\frac\varphi2,\qquad\varphi=\frac{\xi\log p}{\log n}
\tag{5.2}
$$

이다.

[M과 ψ (1.4)](mertens-psi.md#equation-1-4)의 동전 $-\cos\theta_p$는 부호 동전 $-1$과 모양 동전 $\cos\theta_p$의 곱으로 갈라진다. 제곱 인수가 있는 소수의 동전은 $\sin^2$ 동전이다. 이 동전은 $\xi=0$에서 2차로 사라지므로 μ 자체에는 보이지 않고, $d_z$ 모형의 추가 원자로만 나타난다.

<a id="section-5-2"></a>

### 5.2. 모양 동전의 제곱과 디리클레 법칙

$|\hat\nu_n(\xi)|^2=\prod_{p\mid n}\frac{1+\cos\xi u_p}2$는 각 소수를 확률 $(\tfrac14,\tfrac12,\tfrac14)$로 세 칸에 나누는 동전이다.

**명제 5.2 [증명 스케치].** 제곱 인수 없는 $n$에 대한 평균은 $\mathbb E\,e^{i\xi(D_1-D_3)}$으로 수렴한다. 여기서 $(D_1,D_2,D_3)\sim\mathrm{Dir}(\tfrac14,\tfrac12,\tfrac14)$이다. 디리클레 분포의 응집 성질로 $D_1-D_3=S(2B-1)$이고, $S\sim\mathrm{Beta}(\tfrac12,\tfrac12)$는 아크사인 분포, $B\sim\mathrm{Beta}(\tfrac14,\tfrac14)$는 $S$와 독립이다.

근거는 [주정리 정리 1](main-theorem.md#section-2)의 3변수판이다. 디리클레 급수 $\prod_p(1+\tfrac14p^{-s_1}+\tfrac12p^{-s_2}+\tfrac14p^{-s_3})$에 같은 Selberg–Delange 논증을 쓴다. 3변수판 자체는 이 문서에서 끝까지 쓰지 않았다.

**[수치]** 수렴은 느리다. $[10^k,2\cdot10^k)$에서의 평균은 다음과 같다.
- **$\xi=3$:** $L=7.3,\dots,16.5$에서 $0.1724\to0.2296$이다. $a+b/L+c/L^2$로 맞춘 외삽값은 0.294이고, 극한은 0.298이다.
- **$\xi=6$:** $0.4081\to0.3120$이고, 외삽값 0.218, 극한 0.233이다. 차이는 더 높은 차수의 오차항으로 보인다.

<a id="section-5-3"></a>

### 5.3. M4의 이산핵판과 극의 귀환

$S_I(\alpha;\xi)=\sum_{n\in I}\mu(n)\hat\nu_n(\xi)\,e(n\alpha)$, $Q_I(\xi)=\sum_{n\in I}\mu^2(n)|\hat\nu_n(\xi)|^2$로 둔다. 명제 5.1에 의해 제곱 인수 없는 $n$에서는 $\mu(n)\hat\nu_n(\xi)=\langle\mathscr K_n,e^{i\xi u}\rangle$이다. 제곱 인수가 있는 $n$의 원자, 예컨대 (5.2)의 $\sin^2$ 동전은 여기서 뺀다. 이 원자까지 넣으려면 Wick 기준항을 실제 계수의 제곱합으로 다시 세워야 한다. 이렇게 두면 (4.4)의 이산핵판은 다음과 같다.

<a id="equation-5-3"></a>

$$
\int_0^1|S_I(\alpha;\xi)|^2|S_I(\alpha+\theta;\eta)|^2d\alpha
=Q_I(\xi)Q_I(\eta)+\Big|\sum_{n\in I}\mu^2(n)\hat\nu_n(\xi)\overline{\hat\nu_n(\eta)}e(n\theta)\Big|^2-\sum_{n\in I}\mu^2(n)|\hat\nu_n(\xi)\hat\nu_n(\eta)|^2+o.
\tag{5.3}
$$

코사인 언어로 읽으면, 4점 수준에서 부호 동전은 $n$마다 독립이고 모양 동전 $\cos$는 그대로 고정된다. $\xi=\eta=0$이면 (4.4)로 돌아간다.

**[수치]** $N=2^{25}$, $Y=1581$에서 결과는 다음과 같다.
- **메아리:** 모양까지 포함한 (5.3)과 대체로 맞는다. $\theta=1/4$에서 $(3,3)$은 0.1270(공식 0.1269), $(6,6)$은 0.0978(공식 0.0959)이다. 한쪽만 $\xi=0$인 $(0,3)$은 $\theta=0,1/4,1/9$ 모두에서 표준오차 1.3배 안이다.
- **$\theta=0$의 초과분:** $\xi\ne0$에서는 Wick 공식을 넘는 초과분이 있다. $\xi=6$에서 $+0.0138$(표준오차 0.0005), $\xi=3$에서 $+0.0029$로 대략 $\xi^2$에 비례한다.
- **진짜 μ만의 성질:** 같은 모양에 독립 부호를 쓰면 $\xi=6$의 초과분이 $-0.0007$, 무작위 곱셈 부호를 쓰면 $-0.0010$이다(표준오차 0.0005). 초과분은 재현되지 않는다.
- **두 진동수를 섞은 경우:** $(\xi,\eta)=(3,6)$에서 $\theta=0$의 값은 $0.0919$이고 공식은 $0.0864$이다.

초과분의 정체는 극 정보이다. (2.2)에 따라 $F=e^{i\xi u}$의 극 정보는 $F''=-\xi^2F$에 비례하므로 $\xi^2$ 차수로 켜진다. 다만 그 극은 $q=1$의 평균에 있지 않다. 평균은 $\xi=6$에서 $+0.0040$으로, [주정리 정리 2](main-theorem.md#section-8-3)의 주항 $(\xi/2)J_1(\xi/2)/L^2\approx0.0036$과 맞을 만큼 작다. 극은 잉여류 쪽에 있다.
- **잉여류 편향:** $a_n=\mu(n)\prod\cos(\xi u_p/2)$의 mod 30 잉여류 평균은 $\xi=6$에서 30과 서로소인 잉여류가 $+0.075$, 나머지가 $-0.074$에서 $+0.056$ 사이이다.
- **출처는 끝점 $u\approx1$:** $n=mP$($P$는 큰 소수)이면 모양 동전이 $-\cos(\xi/2)$이고, $\xi=6$에서 +0.99이다. 잉여류 1 mod 30에서는 소수만으로 평균에 $+0.220$이 들어온다.
- **Hardy–Littlewood형 상관:** 같은 이유로 $\mathbb E[a_na_{n+k}]$가 작은 소수를 많이 가진 $k$에서 커진다. $k=30,60,210$에서 $1.5,2.5,3.8\times10^{-3}$이고, $k=1,29,31$에서는 $10^{-4}$ 이하이다.
- **결정적 확인:** mod 30 잉여류 평균 모형(로그 구간 4개)을 빼면 초과분이 $\xi=6$에서 $0.0138\to0.0006$, $\xi=3$에서 $0.0029\to-0.0001$로 사라진다.

**정리된 진술.** $\xi\ne0$에서 (5.3)은 $\mu(n)\hat\nu_n(\xi)-k^\sharp_n(\xi)$에 적용해야 한다. 여기서 $k^\sharp$는 주호(major arc) 모형, 즉 잉여류 탑의 극이고, 내용은 끝점의 ψ 성분이다. $\xi=0$에서는 $\mu^\sharp=0$이므로 M4가 그대로 남는다. 한쪽만 $\xi=0$인 $(0,3)$에서 (5.3)이 잘 맞는 것도 같은 이유이다. 극극 항이 생기려면 양쪽 모두 극이 있어야 한다. 이는 [M과 ψ §3](mertens-psi.md#section-3)의 "ψ 성분은 끝점에 몰린다"는 관찰이 4점 수준에서 다시 나타난 것이다.

<a id="section-6"></a>

## 6. 정리와 남은 문제

**증명된 것.**
- 팽창 항등식과 gcd 분해(명제 1.1–1.2)
- 무작위 곱셈 모형에서 원시 부분의 직교성과 비선형 일치의 유한성(명제 1.3)
- 제곱류 표현, 층별 밀도, 제곱류 일치의 완전한 분류(명제 1.4–1.6)
- 순수 Chowla 합에서 극 성분의 소멸과 극 범함수(명제 2.1–2.2)
- 짧은 구간의 정확한 다리(정리 4.1)
- 독립 부호 모형과 짧은 구간의 소수별 무작위 곱셈 모형에서 M4가 오차 없이 성립한다는 것(명제 4.4)
- 2k차 모멘트 사다리(정리 4.5)와 독립 부호 모형의 정확한 2k차 값(명제 4.6)
- 사다리 꼭대기의 점근적 대각성과 $Y\asymp\sqrt X$에서 정확한 대각성의 반례(명제 4.8)
- 이산핵의 부호·모양 분해(명제 5.1)

**추측으로 남은 것.**
- M4, 추측 4.3: 정확한 나머지 $R_I(\theta)$의 상쇄 $\sup_\theta|R_I(\theta)|=o(Q_I^2)$. 짧은 구간에서 이동 평균한 Chowla(Mirsky 상수 포함), 제곱류 메아리, 블록 독립이 여기서 따른다. 고정된 $h$의 강한 Chowla는 따로 남는다.
- M2k, 추측 4.7: 순열 규칙과 평평한 사다리. $k=2$가 M4이고, 모든 $k$에서 성립하면 $T\asymp X$에서 부호 있는 평균 $\langle P_I^{\,a}\overline{P_I}^{\,b}\rangle_T/Q_I^{(a+b)/2}$가 복소 가우스의 모멘트 값으로 수렴한다.
- OU 법칙과 (3.1), 추측 3.1: $x^\varepsilon$ 대신 $\sqrt{\log\log x}$를 쓰는 형태. 모형에서는 고정된 $h$의 공분산과 중심극한정리가 알려져 있다(Chinis–Shala). 실제 μ, 여러 이동의 결합분포, 과정 수준의 수렴은 남는다.
- 영점 쪽 W4: 현재 기술로는 가중 2점조차 닿지 않는다.

**다음 단계.**
- 정확한 등식 $\langle|P_I|^2\rangle=Q_I$를 이용해, RH와 $J_{-1}\ll T$ 아래 $\alpha$가 1 가까이까지 성립하는 가중 쌍상관 정리(Montgomery 정리의 Möbius판)를 적는다.
- 실제 영점 자료로 높이 메아리를 측정한다.
- 명제 4.4 (2)를 2k차로 올린다. 짧은 구간에서 $\sum n_i=\sum n_i'$이고 $\prod n_i\prod n_i'$가 제곱수이면 짝이 강제되는지, 아니면 명제 4.8처럼 점근적으로만 그런지 묻는 문제이다. §4.4의 수치에서는 무작위 곱셈 모형이 독립 부호와 같다.
- 양의 가중함수로 사다리를 다시 적어, M2k에서 $P_I(t)$의 분포 수렴(중심극한정리)까지 끌어낸다.
- $N=2^{22}$에서 보인 공분산의 크기 의존성(§1.3)을 설명한다.

<a id="checks"></a>

## 검산

[검산 스크립트](../scripts/check_chowla.py)는 numpy와 scipy를 쓴다.

- **정확한 항등식:** (1.2)와 (1.3)은 $x=2\cdot10^5$에서 정수로 일치한다. (2.2)의 모멘트는 $k<25$에서 기계 정밀도로 일치한다. 정리 4.1의 이분법은 전수 검사($X=4000$, $Y=31$)와 무작위 4쌍($X=10^7,10^{12}$)으로 확인한다. 명제 5.1의 동전 분해도 확인한다.
- **명제 1.3:** $h<h'\le60$, $n\le5000$에서 해의 개수가 $\tfrac12\tau(|h^2-h'^2|)$를 넘지 않는다.
- **§1.4:** (1.6)은 점마다, (1.7)은 $h\le2000$의 모든 $h$에서, (1.8)은 해의 전수 검사로 확인한다. 4점 관계 $6\cdot30\cdot42\cdot210=1260^2$도 확인한다.
- **명제 4.4 (2):** $Y^2\le X/4$인 세 구간에서 각 $k$의 $q_k$가 단사이고 이웃 구간과 겹치는 제곱류가 없음을, $(0,15]$에서는 단사가 아님을 확인한다.
- **§4.4:** 정리 4.5의 이분법은 여섯 경우의 모든 다중집합 쌍에서 정확한 정수로 확인한다. 명제 4.6은 작은 집합의 전수 열거와 부호 전체에 대한 평균으로 확인한다. $D_k$의 다항식도 생성함수 (4.9)와 대조한다.
- **[수치]:** §1.3, §1.4, §2, §3, §4.2, §4.4, §5.2–5.3의 모든 수치를 재현한다. 기본값은 $N=2^{22}$(1분 이내)이고, `--full`은 이 문서의 $N=2^{25}$이다(약 7분). 디리클레 다항식 쪽의 표본 추출, 브라운 운동 비교, 무작위 부호 모형은 정해진 시드를 쓴다.
- **기본 실행에서 달라지는 점:** $N=2^{22}$에서는 §1.3의 공분산이 약 16% 낮게 나오고, §5.3의 초과분은 $Y=500$으로 짧아 약 3 표준오차($\xi=6$에서 $+0.0041$)로 줄어든다. §4.4는 $Y=500$(2k차)과 $Y=50$(사다리)에서 계산하므로 표준오차가 2–3배 크고, 결론은 같다.

이 수치 확인은 증명을 대체하지 않는다.

    python -B scripts/check_chowla.py
    python -B scripts/check_chowla.py --full

<a id="references"></a>

## 참고 문헌

- S. Chowla, *The Riemann Hypothesis and Hilbert's Tenth Problem*, Gordon and Breach (1965): μ의 상관에 관한 추측.
- L. Mirsky, r-free 정수의 패턴에 관한 논문들(1947–1949): 제곱 인수 없는 수 쌍의 밀도 $S_2(h)$.
- N. Ng, The distribution of the summatory function of the Möbius function, Proc. London Math. Soc. (3) 89 (2004), 361–389: $J_{-1}$ 가정, $M(x)$의 극한분포와 $(\log\log\log x)^{5/4}$ 추측.
- N. Ng, The Möbius function in short intervals, in *Anatomy of Integers*, CRM Proc. Lecture Notes 46 (2008), 247–257: 강한 Chowla 가정 아래 짧은 구간의 모멘트.
- J. P. Keating, Z. Rudnick, [Squarefree polynomials and Möbius values in short intervals and arithmetic progressions](https://msp.org/ant/2016/10-2/ant-v10-n2-p04-p.pdf), Algebra Number Theory 10 (2016), 375–420: 분산 $H/\zeta(2)$.
- D. A. Goldston, H. L. Montgomery, Pair correlation of zeros and primes in short intervals, in *Analytic Number Theory and Diophantine Problems*, Progr. Math. 70 (1987), 183–203.
- D. W. Farmer, S. M. Gonek, Y. Lee, S. J. Lester, [Mean values of ζ′/ζ(s), correlations of zeros and the distribution of almost primes](https://nms.kcl.ac.uk/steve.lester/FGLL.pdf) (2013).
- K. Matomäki, M. Radziwiłł, T. Tao, [Fourier uniformity of bounded multiplicative functions in short intervals on average](https://arxiv.org/abs/1812.01224), Invent. Math. 220 (2020), 1–58.
- J. Chinis, B. Shala, [Random Chowla's conjecture for Rademacher multiplicative functions](https://arxiv.org/abs/2409.05952), Trans. Amer. Math. Soc., doi:10.1090/tran/9457: $f(P(n))$의 4차 모멘트(명제 3.1)와 중심극한정리. $P(n)=n(n+h)$를 포함한다.
- J. Bourgain, C. Demeter, L. Guth, Proof of the main conjecture in Vinogradov's mean value theorem for degrees higher than three, Ann. of Math. (2) 184 (2016), 633–682: (4.8)의 거듭제곱합 계와 그 평균값.
