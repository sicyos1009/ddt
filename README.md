# 일반화된 베타 핵과 가중 약수합

DDT의 평균 아크사인 법칙에서 일반화된 베타핵을 거쳐 약수별 이산핵으로 이어지는 연구 원고다. 2026-11-01 개정판의 이차 베타 법칙을 바탕으로, 2026-11-03에 모든 정수의 $d_z$ 가중치와 정확한 Mellin–Fourier 변환을 추가했다.

| 문서 | 내용 | 범위 |
|:--|:--|:--|
| [주정리](docs/main-theorem.md) | DDT 도입, 기존 이차 베타 법칙, 이산핵의 정확한 항등식과 평균의 이차 베타 법칙 | §0–7, 부록 A, B.1 |
| [후속 결과](docs/further-results.md) | 뫼비우스 경우, 고차 베타 전개, 제타함수 영점과 복원, RH 아래 극한분포 | §8–10, 부록 B.2–B.4, C |
| [전역 공식](docs/global-formula.md) | FTA의 소수별 합성곱, 영점 곱 재구성, 유한 가법 공식, 가우스 평활화와 작은 나머지 | §11–12 |

주정리 문서의 §0에서 DDT와의 연결을 읽고, §1–6에서 베타핵과 증명을, §7에서 이산핵과 정확한 변환을 따라갈 수 있다. 기존 제곱 인수 없는 모형의 $A_x$와 새 $d_z$ 모형의 $\widetilde A_x$를 구분한다. §7 추가에 따라 후속 결과의 기존 §7–9는 §8–10으로 옮겼으며, 참조와 링크도 함께 갱신했다. 참고 문헌은 각 문서 끝에 있다.

전역 공식 문서는 이산핵의 정확한 Perron 표현과 Hadamard 곱을 연결한다. 모든 영점을 포함하는 곱 재구성과 유한 절단오차를 제시한다. §12에서는 가우스 평활화와 영점 높이 선택으로, 정확한 극점 기여와 유한 영점 성분합 밖의 나머지를 RH 없이 $O(x^{-1})$로 줄인다. 평활화하지 않은 합으로 돌아갈 때의 오차는 별도로 구분한다.

## 정리와 명제 찾기

| 결과 | 위치 |
|:--|:--|
| 정리 1: 이차 베타 법칙 | [주정리 §1](docs/main-theorem.md#section-1) |
| 보조정리 1: 이동 균등 추정 | [주정리 §5.1](docs/main-theorem.md#section-5-1) |
| 명제 7.1: 이산핵의 총질량과 반사 | [주정리 §7.1](docs/main-theorem.md#section-7-1) |
| 명제 7.2: 정확한 두 변수 변환 | [주정리 §7.2](docs/main-theorem.md#section-7-2) |
| 정리 2: 이산핵 평균의 이차 베타 법칙 | [주정리 §7.3](docs/main-theorem.md#section-7-3) |
| 명제 7.3: 정확한 Perron 표현 | [주정리 §7.5](docs/main-theorem.md#section-7-5) |
| 명제 8.1: 비양의 정수 지수의 소멸 추정 | [후속 결과 §8.1](docs/further-results.md#section-8-1) |
| 명제 9.1: 고차 베타 전개 | [후속 결과 §9.1](docs/further-results.md#section-9-1) |
| 명제 9.2: 뫼비우스 합의 오차와 영점 없는 영역 | [후속 결과 §9.4](docs/further-results.md#section-9-4) |
| 명제 10.1: 모멘트 차수와 영점 중복도 | [후속 결과 §10.3](docs/further-results.md#section-10-3) |
| 명제 10.2: RH 아래 극한분포 | [후속 결과 §10.5](docs/further-results.md#section-10-5) |
| 정리 3: 생성함수의 전역 영점 곱 | [전역 공식 §11.2](docs/global-formula.md#section-11-2) |
| 정리 4: 유한 영점 자료에 의한 재구성 | [전역 공식 §11.3](docs/global-formula.md#section-11-3) |
| 따름정리 4.1: 정확한 전역 재구성 | [전역 공식 (11.12)](docs/global-formula.md#equation-11-12) |
| 명제 11.1: 유한 높이의 가법 공식 | [전역 공식 §11.4](docs/global-formula.md#section-11-4) |
| 명제 11.2: 영점의 국소 반사핵 | [전역 공식 §11.5](docs/global-formula.md#section-11-5) |
| 보조정리 12.1: 영점 쌍을 가르지 않는 절단 높이 | [전역 공식 §12.2](docs/global-formula.md#section-12-2) |
| 정리 5: 평활화된 유한 가법 공식 | [전역 공식 §12.3](docs/global-formula.md#section-12-3) |
| 따름정리 5.1: 로그폭 평활화와 $O(x^{-1})$ 나머지 | [전역 공식 §12.4](docs/global-formula.md#section-12-4) |

**[개요]**는 증명의 뼈대만 제시한 부분이며, **RH**는 리만 가설을 뜻한다. 결과별 가정과 수치 확인은 원문대로 보존했다.

전역 공식의 [검산 스크립트](scripts/check_global_formula.py)는 Python과 mpmath를 사용한다. 소수별 합성곱은 유리수로, 영점 곱과 윤곽 적분은 수치적으로 확인한다.

    python scripts/check_global_formula.py

[평활화 검산](scripts/check_smoothed_formula.py)은 0이 아닌 주파수에서 직접 약수합, 가우스 꼬리 상계, 극점·영점 성분의 정확한 윤곽 분해를 확인한다.

    python scripts/check_smoothed_formula.py
