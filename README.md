# 일반화된 베타핵과 가중 약수합

**DDT의 평균 아크사인 법칙 → 일반화된 베타핵 → 약수별 이산핵 → 영점에 의한 전역 표현**을 다루는 연구 원고다. FTA는 산술의 기본정리, 즉 소인수분해의 유일성을 뜻한다.

## 핵심 문서

| 문서 | 내용 | 자체 번호 |
|:--|:--|:--|
| [주정리](docs/main-theorem.md) | DDT, 이차 베타 법칙, 이산핵과 정확한 변환 | §1–8, 정리 1–2, 부록 A–B |
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

[후속 연구 안내](docs/further-results.md)에 문서별 범위와 수치 기록의 위치를 모았다. **[국소 계산]**, **[휴리스틱]**, **RH** 등의 표시는 각 주장의 범위와 가정을 구분한다.

## 주요 결과 찾기

| 결과 | 위치 |
|:--|:--|
| 이차 베타 법칙과 이산핵 평균 | [주정리 정리 1](docs/main-theorem.md#section-2), [정리 2](docs/main-theorem.md#section-8-3) |
| 영점 곱과 전역 재구성 | [전역 공식 정리 1–2](docs/global-formula.md#section-1-2) |
| 평활화된 가법 공식 | [전역 공식 정리 3](docs/global-formula.md#section-2-3) |
| 원래 정수합의 유한 공식과 정확한 가우스 합산 | [전역 공식 정리 4–5](docs/global-formula.md#section-3-4) |
| RH 아래 극한분포 | [극한분포 명제 1](docs/limiting-distribution.md#section-2) |

## 검산

Python과 mpmath를 사용한다.

| 스크립트 | 확인 범위 |
|:--|:--|
| [check_global_formula.py](scripts/check_global_formula.py) | 소수별 합성곱, 유한 영점 곱, 비틀림 없는 윤곽 적분 |
| [check_smoothed_formula.py](scripts/check_smoothed_formula.py) | 비틀린 약수합, 가우스 꼬리 상계, 극점·영점 성분의 분해 |
| [check_desmoothed_formula.py](scripts/check_desmoothed_formula.py) | 정수 격자의 평활화 오차, 경계가 영점 쌍을 가를 때의 윤곽 항등식 |

    python -B scripts/check_global_formula.py
    python -B scripts/check_smoothed_formula.py
    python -B scripts/check_desmoothed_formula.py

수치 검산은 증명을 대체하지 않는다. 주제별 부록의 기존 수치표는 계산 코드와 원자료가 없어 재현하지 않았으며, 위 스크립트의 검증 범위에 포함되지 않는다.
