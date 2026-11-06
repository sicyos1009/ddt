# 일반화된 베타핵과 가중 약수합

**DDT의 평균 아크사인 법칙 → 일반화된 베타핵 → 약수별 이산핵 → 영점에 의한 전역 표현**을 다루는 연구 원고다. FTA는 산술의 기본정리, 즉 소인수분해의 유일성을 뜻한다.

## 핵심 문서

| 문서 | 내용 | 자체 번호 |
|:--|:--|:--|
| [주정리](docs/main-theorem.md) | DDT, 이차 베타 법칙, 이산핵과 정확한 변환, 임의 차수 전개 | §1–9, 정리 1–3, 부록 A–B |
| [전역 공식](docs/global-formula.md) | 소수별 합성곱, 영점 곱·가법 공식, 가우스 합산과 이산핵 복원 | §1–3, 정리 1–5 |

주정리를 먼저 읽고 전역 공식으로 이어가면 핵심 흐름을 따라갈 수 있다. 절·정리·명제·수식 번호는 **문서마다 독립적**이며, 다른 문서를 참조할 때는 문서명과 링크를 함께 표시한다. 각 문서 앞에는 필요한 기호와 전제를 정리했다.

$A_x$는 제곱 인수 없는 모형, $\widetilde A_x$는 모든 정수의 $d_z$ 모형이다. 전역 공식은 RH와 단순영점 가정을 쓰지 않는다. 산술합의 평활화는 제거하지만 영점 쪽 가우스 합산 규칙은 유지한다. $O(x^{-1})$ 나머지는 극점·영점 성분을 정확한 윤곽 적분으로 유지한 공식의 추정이다.

## 주제별 후속 연구

| 문서 | 주제 |
|:--|:--|
| [뫼비우스 가중치와 끝점 항](docs/mobius-cases.md) | 정수·반정수 지수, 소멸과 끝점 항, 영점 없는 영역 |
| [고차 베타 전개와 국소 특이점](docs/higher-order.md) | 고차 계수, 수렴 반경, 모멘트 상쇄, 세 점 핵 |
| [영점의 국소 기여와 모멘트·복원](docs/zero-contributions.md) | 영점 중복도, 보조 통계, 영점 높이와 미분값 |
| [RH 아래 보조 통계의 극한분포](docs/limiting-distribution.md) | 극한분포와 평균제곱 오차의 증명 |
| [Ng 공식과 소수·약수 쪽 쌍대](docs/ng-duality.md) | Ng 분포의 가정 약화, 가정 없는 쌍대 극한분포, Mertens 함수와의 연결 |
| [M과 ψ: 코사인 핵과 사인 핵](docs/mertens-psi.md) | $\mu=d_{-1/2}*d_{-1/2}$와 $\varepsilon=d_{-1/2}*d_{1/2}$, ψ의 국소 핵, 네 핵의 아크사인 표현, 소수핵으로 쓴 영점 측도 |
| [약수별 누적합으로 얻는 M(x)](docs/divisor-cumsum.md) | $\mu(n)/\tau(n)$의 약수별 몫과 누적합이 $M(x)$가 되는 항등식, 매끄러운 시험함수의 1차 법칙, 누적 곡선의 $\mathcal K_{-1/2,-1/2}$ 예측 |
| [짝수 중복도 영점과 $d_{1/2}$ 오차항](docs/even-zeros.md) | $\zeta^{1/2}$ 계수합에서 분기절단의 전체 주항을 뺀 잔차 $E$, 짝수 중복도 조건과의 동치, 성장 지수 $=\Theta_{\rm odd}$와 $\Omega_\pm$ |
| [양의 약수분포의 이차 모멘트](docs/m2-divisor-moment.md) | $(1/2,1/2)$ 약수분포의 $M_2(n)$, 전체 주항, 단순영점 기여의 정확식과 로그 전개, $\mathcal E_2$의 성장 지수 $=\Theta_{\rm simple}$ |
| [강한 Chowla 추측: 제곱류 층과 4점 Wick 규칙](docs/chowla.md) | $\sum\mu(n)\mu(n+h)$의 제곱류 층·제곱류 표현과 팽창 격자, 극 성분의 소멸, 원시 성분의 OU 법칙, 짧은 구간의 4차 모멘트 사전, 2k차 사다리와 순열 규칙, 이산핵·코사인 언어와 극의 귀환 |
| [사인과 코사인으로 쓴 RH 판정법](docs/sine-cosine-rh.md) | 로그 동전의 켤레성, Littlewood 상쇄와 임계선 위 켤레 변환으로 쓴 두 판정법, 선 밖 영점의 결손 |

[후속 연구 안내](docs/further-results.md)에 문서별 범위와 수치 기록의 위치를 모았다. **[국소 계산]**, **[휴리스틱]**, **RH** 등의 표시는 각 주장의 범위와 가정을 구분한다.

## 주요 결과 찾기

| 결과 | 위치 |
|:--|:--|
| 이차 베타 법칙과 이산핵 평균 | [주정리 정리 1](docs/main-theorem.md#section-2), [정리 2](docs/main-theorem.md#section-8-3) |
| 임의 차수 베타 법칙, 오차 $O(xL^{\beta-1-J})$ | [주정리 정리 3](docs/main-theorem.md#section-9-1) |
| 영점 곱과 전역 재구성 | [전역 공식 정리 1–2](docs/global-formula.md#section-1-2) |
| 평활화된 가법 공식 | [전역 공식 정리 3](docs/global-formula.md#section-2-3) |
| 원래 정수합의 유한 공식과 정확한 가우스 합산 | [전역 공식 정리 4–5](docs/global-formula.md#section-3-4) |
| RH 아래 극한분포 | [극한분포 명제 1](docs/limiting-distribution.md#section-2) |
| 로그 가중 조건 아래 Mertens 함수의 극한분포 | [Ng 쌍대 명제 1.1](docs/ng-duality.md#section-1) |
| 소수별 아크사인 곱과 $S_1$의 극한분포 (가정 없음) | [Ng 쌍대 정리 1–2](docs/ng-duality.md#section-2-2) |
| 영점 봉우리의 국소 극한과 Ng의 분산 | [Ng 쌍대 명제 5.1–5.2](docs/ng-duality.md#section-5-2) |
| ψ의 극·영점 국소 핵과 M과의 비교 | [M과 ψ 명제 2.1, 표 (2.4)](docs/mertens-psi.md#section-2) |
| 소수핵으로 쓴 영점 측도 (가정 없음) | [M과 ψ 명제 4.1](docs/mertens-psi.md#section-4) |
| 약수별 몫의 누적합 = $M(x)$, 누적 곡선과 베타핵 | [약수별 누적합 명제 1.1, 2.1](docs/divisor-cumsum.md#section-1) |
| Chowla 합의 팽창 항등식과 gcd 분해 | [Chowla 명제 1.1–1.2](docs/chowla.md#section-1-1) |
| 제곱류 표현, 층별 밀도, 제곱류 일치의 분류 | [Chowla 명제 1.4–1.6](docs/chowla.md#section-1-4) |
| 짧은 구간의 덧셈 에너지 = 디리클레 다항식의 4차 모멘트 | [Chowla 정리 4.1](docs/chowla.md#section-4-1) |
| μ의 4점 Wick 규칙(M4)과 제곱류 메아리 | [Chowla 추측 4.3](docs/chowla.md#section-4-2) |
| 짧은 구간의 무작위 곱셈 모형에서 M4가 정확히 성립 | [Chowla 명제 4.4](docs/chowla.md#proposition-4-4) |
| 2k차 모멘트 사다리: 덧셈 에너지에서 Steinhaus 모형까지 | [Chowla 정리 4.5](docs/chowla.md#theorem-4-5) |
| μ의 2k차 순열 규칙(M2k)과 평평한 사다리 | [Chowla 추측 4.7](docs/chowla.md#conjecture-4-7) |
| 사다리 꼭대기의 점근적 대각성과 $Y=2\sqrt X$의 반례 | [Chowla 명제 4.8](docs/chowla.md#proposition-4-8) |
| 짝수 중복도 조건 아래 $d_{1/2}$ 잔차 $E\ll x^\theta\log x$와 그 역 | [짝수 중복도 정리 1, 따름정리 3](docs/even-zeros.md#section-1) |
| $E$의 성장 지수 $=\Theta_{\rm odd}$, 가정 없는 $\Omega_\pm(x^{1/2-\varepsilon})$ | [짝수 중복도 정리 2](docs/even-zeros.md#section-9) |
| 이차 모멘트 잔차 $\mathcal E_2\ll x^{\theta+\varepsilon}$ $\iff$ $\Re s>\theta$에 단순영점 없음 | [이차 모멘트 정리 8.1](docs/m2-divisor-moment.md#section-8-1) |
| 코사인·사인 쪽의 Littlewood 상쇄와 RH의 동치 | [사인·코사인 판정법 1](docs/sine-cosine-rh.md#section-2) |
| 임계선 위 사인 쪽이 코사인 쪽의 켤레 변환인 것과 RH의 동치 | [사인·코사인 판정법 2](docs/sine-cosine-rh.md#section-3) |

## 검산

Python을 사용한다. 핵심 문서의 검산은 mpmath를([주정리 §9](docs/main-theorem.md#section-9)는 numpy와 scipy), [Ng 쌍대](docs/ng-duality.md#checks), [M과 ψ](docs/mertens-psi.md#checks), [Chowla](docs/chowla.md#checks), [약수별 누적합](docs/divisor-cumsum.md#checks), [사인·코사인 판정법](docs/sine-cosine-rh.md#checks)의 검산은 numpy와 scipy를 쓴다.

| 스크립트 | 확인 범위 |
|:--|:--|
| [check_global_formula.py](scripts/check_global_formula.py) | 소수별 합성곱, 유한 영점 곱, 비틀림 없는 윤곽 적분 |
| [check_smoothed_formula.py](scripts/check_smoothed_formula.py) | 비틀린 약수합, 가우스 꼬리 상계, 극점·영점 성분의 분해 |
| [check_desmoothed_formula.py](scripts/check_desmoothed_formula.py) | 정수 격자의 평활화 오차, 경계가 영점 쌍을 가를 때의 윤곽 항등식 |
| [check_beta_expansion.py](scripts/check_beta_expansion.py) | 주정리 §9의 2차 계수 닫힌 꼴, $J\le6$차 전개와 실측 합 비교, 부록 B.1 표의 재현 (`--x 1e7`) |
| [check_ng_tail.py](scripts/check_ng_tail.py) | Ng 쌍대 §1의 Gallagher형 평균제곱 부등식과 로그 가중 부분합 |
| [check_dual_product.py](scripts/check_dual_product.py) | Ng 쌍대 §2의 대칭 비틀기 항등식, 중심 원자 합, 소수별 아크사인 곱 |
| [check_s1_distribution.py](scripts/check_s1_distribution.py) | Ng 쌍대 §3의 Littlewood 공식, 근사 오차, $S_1$의 분포와 왜도 |
| [check_mertens_bridge.py](scripts/check_mertens_bridge.py) | Ng 쌍대 §5의 국소 극한, Plancherel 등식, $\lvert\zeta\rvert^{-2}$ 평균, Gonek 비율 |
| [check_mertens_psi.py](scripts/check_mertens_psi.py) | M과 ψ의 합성곱 항등식, 소수별 동전, 네 국소 핵, $\arg\zeta$의 특성함수, 대역 제한 명시공식과 소수만으로 찾은 영점 |
| [check_chowla.py](scripts/check_chowla.py) | Chowla 합의 팽창·gcd 항등식, 제곱류 표현과 일치 분류, Mirsky 분산과 공분산 격자, 원시 성분, 짧은 구간 이분법, M4와 메아리, 곱셈 모형의 M4, 2k차 사다리와 순열 규칙, 꼭대기의 반례, 이산핵의 동전과 극의 귀환 (`--full`은 $N=2^{25}$) |
| [check_divisor_cumsum.py](scripts/check_divisor_cumsum.py) | 약수별 몫의 합 = $M(x)$, 몫의 분해, $f_x(1)$의 점근식, 핵 누적함수의 닫힌 꼴, 모멘트와 누적 곡선의 수치 (`--figure`는 그림 생성, matplotlib 필요) |
| [check_sine_cosine.py](scripts/check_sine_cosine.py) | 로그 동전, 코사인·사인 쪽의 Littlewood 상쇄, 켤레 변환 공식, 임계선 위 톱니의 재현, 높이와 구간 통계, 그림 1–3 생성(`--figures`) |

    python -B scripts/check_global_formula.py
    python -B scripts/check_smoothed_formula.py
    python -B scripts/check_desmoothed_formula.py
    python -B scripts/check_beta_expansion.py
    python -B scripts/check_ng_tail.py
    python -B scripts/check_dual_product.py
    python -B scripts/check_s1_distribution.py
    python -B scripts/check_mertens_bridge.py
    python -B scripts/check_mertens_psi.py
    python -B scripts/check_chowla.py
    python -B scripts/check_divisor_cumsum.py
    python -B scripts/check_sine_cosine.py

수치 검산은 증명을 대체하지 않는다. 주제별 부록의 기존 수치표는 계산 코드와 원자료가 없어 재현하지 않았으며, 위 스크립트의 검증 범위에 포함되지 않는다. 예외로 주정리 부록 B.1의 표는 `check_beta_expansion.py --x 1e7`이 재현한다.
