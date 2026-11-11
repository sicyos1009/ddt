# 후속 연구 안내

후속 연구는 아래 열세 문서로 나뉜다. 각 문서는 자체 절·명제·수식 번호와 기호 설명을 갖는다. 핵심 흐름은 [주정리](main-theorem.md)에서 [전역 공식](global-formula.md)으로 이어지며, 아래 문서는 필요한 주제별로 읽을 수 있다.

| 문서 | 다루는 문제 | 범위 |
|:--|:--|:--|
| [뫼비우스 가중치와 끝점 항](mobius-cases.md) | 음의 지수의 소멸, 끝점 항, 거듭제곱 오차와 영점 없는 영역 | §1–5, 부록 A |
| [고차 베타 전개와 국소 특이점](higher-order.md) | 고정 차수 전개, 계수의 수렴 반경, 모멘트 상쇄, 세 점 핵 | §1–4, 부록 A |
| [영점의 국소 기여와 모멘트·복원](zero-contributions.md) | 영점 중복도와 모멘트, 보조 통계, 높이·미분값의 복원 | §1–4, 부록 A |
| [RH 아래 보조 통계의 극한분포](limiting-distribution.md) | 극한분포 명제와 이를 뒷받침하는 평균제곱 증명 | §1–3 |
| [Ng 공식과 소수·약수 쪽 쌍대](ng-duality.md) | Ng 분포의 가정 약화, 소수별 아크사인 곱, 누적량 $S_1$, Mertens 함수와의 연결 | §1–5 |
| [M과 ψ: 코사인 핵과 사인 핵](mertens-psi.md) | $M$과 $\psi$의 두 분해, ψ의 국소 핵, 네 핵의 아크사인 표현, $\arg\zeta$의 극한분포, 소수핵으로 쓴 영점 측도 | §1–4 |
| [강한 Chowla 추측: 제곱류 층과 4점 Wick 규칙](chowla.md) | 이동 상관의 제곱류 층, 극 성분의 소멸, 원시 성분과 OU 법칙, 4차 모멘트 사전과 M4, 2k차 사다리와 순열 규칙(M2k), 이산핵·코사인 언어 | §1–6 |
| [약수별 누적합으로 얻는 M(x)](divisor-cumsum.md) | $\mu(n)/\tau(n)$의 약수별 몫, 누적합으로 복원되는 $M(x)$, 누적 곡선과 $\mathcal K_{-1/2,-1/2}$ | §1–3 |
| [짝수 중복도 영점과 $d_{1/2}$ 오차항](even-zeros.md) | 분기절단 주항을 뺀 잔차 $E$, 짝수 중복도 조건 아래 상계와 그 역, 성장 지수 $\Theta_{\rm odd}$와 $\Omega_\pm$, 홀수 영점의 국소 기여 | §1–10 |
| [양의 약수분포의 이차 모멘트](m2-divisor-moment.md) | $(1/2,1/2)$ 약수분포의 $M_2(n)$과 전체 주항, 단순영점 기여의 정확식·상쇄·로그 전개, 단순극을 보존하는 평균 추정, $\mathcal E_2$의 성장 지수 $=\Theta_{\rm simple}$, $b(n)$의 점별 영점 표현, $\psi-x$를 통한 이산–연속 분해와 중복 영점의 소멸, 첫 번째 차로 쓴 RH 판정 | §1–10 |
| [$cp$ 계열의 이차 모멘트](m2-cp-families.md) | 띠 $n=cp$의 정확식과 서로소 이차 가법성, U자 모양, 두 띠의 안정성 부등식과 $M$만의 차, 띠 평균의 1·2차 점근식과 균등판, 소수정리 적분 예측기의 균등한 오차 상계, 띠 분할과 평균 $3/8$의 Poisson–Dirichlet 해석, $\mathcal E_2$로 옮길 때의 장애 | §1–7, 부록 A |
| [$M_2$의 가법 핵](m2-additive-core.md) | 점별 주항·오차 $p_2,e_2$와 $\mathcal E_2$의 대응, $e_2$의 코사이클과 가법 핵, 분할 통계량 층과 소수 거듭제곱 층, FTA와 Beurling 판정, 분배법칙과 분극 항등식, 이동 상관의 세 영역, 꼬인 톱니와 기약 Hermite 항등식 | §1–9 |
| [사인과 코사인으로 쓴 RH 판정법](sine-cosine-rh.md) | 로그 동전의 켤레성, Littlewood 상쇄(판정법 1), 임계선 위 켤레 변환(판정법 2), 선 밖 영점의 결손 | §1–4 |

**[국소 계산]**은 고정된 특이점의 기여로, 무한 영점합의 수렴을 포함하지 않는다. **[휴리스틱]**은 추가 추정이 필요한 논의이다. RH와 LI 등 가정은 사용하는 명제에 명시한다.

<a id="numerics"></a>

## 수치 기록

기존 수치표는 해당 주제의 부록으로 옮겼다.

- [주정리의 이차 법칙과 J차 전개](main-theorem.md#appendix-b)
- [뫼비우스 가중합](mobius-cases.md#appendix-a)
- [고차 전개의 계수](higher-order.md#appendix-a)
- [영점 기여와 복원 실험](zero-contributions.md#appendix-a)

이 기록들은 계산 코드와 원자료가 없어 재현하지 못한 초기 원고의 보고값이다. 예외로 주정리의 표는 [check_beta_expansion.py](../scripts/check_beta_expansion.py)로 재현했다. 현재 실행 가능한 [전역 공식의 검산](global-formula.md#checks), [Ng 쌍대의 검산](ng-duality.md#checks), [M과 ψ의 검산](mertens-psi.md#checks), [Chowla의 검산](chowla.md#checks), [약수별 누적합의 검산](divisor-cumsum.md#checks), [사인·코사인 판정법의 검산](sine-cosine-rh.md#checks), [cp 계열의 검산](m2-cp-families.md#checks), [가법 핵의 검산](m2-additive-core.md#checks)과 구분한다. cp 계열 문서의 부록 A는 그 스크립트가 모두 재현한다.
