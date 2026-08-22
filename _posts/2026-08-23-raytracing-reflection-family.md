---
title: 기하광학 3편 — 반사는 왜 하나의 공식으로 충분한가
date: 2026-08-23 19:00:00 +0900
categories:
  - 광학
  - 기하광학
tags:
  - ray-tracing
  - geometric-optics
  - reflection
  - python
  - matlab
description: 반사 벡터 공식이 코드 전체에서 왜 하나로 충분한지 보고, 평면거울·구면거울·빔스플리터·임의형상거울 네 부품에서 실제로 어떻게 쓰이는지 본다. 임의형상거울의 꼭짓점 법선 계산에서 세 점이 서로 달라야 한다는 조건이 깨지면 결과가 얼마나 달라지는지도 확인한다.
math: true
---
2편 마지막에 반사 벡터 공식이 코드 전체에서 여섯 번 나온다고 적었다. `FlatMirror`, `SphericalMirror`, `ArbitraryMirror`, `BeamSplitter`, 그리고 `Prism`과 `SphericalLens`의 TIR 분기. 이번 편에서는 이 공식이 왜 하나로 충분한지 보고, 평면거울·구면거울·빔스플리터·임의형상거울 네 부품에서 실제로 어떻게 쓰이는지 본다.

## 1. 반사 공식이 하나인 이유

굴절은 매질이 바뀔 때마다 다른 각도 관계식이 필요했다. Snell's law 자체가 $n_1 \sin\theta_1 = n_2 \sin\theta_2$로 양쪽 매질의 굴절률을 입력으로 받는다. 렌즈, 프리즘, TIR 분기까지 각도 공식이 매번 조금씩 달랐던 이유다.

반사는 다르다. 입사각과 반사각이 같다는 것뿐이고, 여기엔 굴절률이 아예 등장하지 않는다. 입사 벡터 $\vec v$를 법선 $\hat n$ 방향 성분과 접선 방향 성분으로 나누면

$$ \vec v = (\vec v \cdot \hat n)\hat n + \big[\vec v - (\vec v \cdot \hat n)\hat n\big] $$

반사는 접선 성분은 그대로 두고 법선 성분만 부호를 뒤집는다.

$$ \vec v' = -(\vec v \cdot \hat n)\hat n + \big[\vec v - (\vec v \cdot \hat n)\hat n\big] = \vec v - 2(\vec v \cdot \hat n)\hat n $$

매질 종류가 어디에도 안 들어간다. 그래서 거울이든 프리즘 내부 전반사든 렌즈 출사면 TIR이든 이 식 하나로 끝난다.

```python
def reflect_vector(v, n):
    """법선 n 기준 완전 반사. v, n은 (vx, vy) 튜플/배열, n은 단위벡터.

    부호가 뒤집힌 n을 넣어도 결과는 같다 -- (v.n)이 부호와 함께 뒤집혀서
    상쇄되기 때문이다.
    """
    v = np.asarray(v, dtype=float)
    n = np.asarray(n, dtype=float)
    return v - 2 * np.dot(v, n) * n
```

법선의 부호도 상관없다. 실제로 코드를 보면 함수마다 이 부분을 다르게 다룬다. `FlatMirror`와 `BeamSplitter`는 법선을 그냥 직선의 계수에서 뽑아 쓰고 부호를 신경 쓰지 않는다. `ArbitraryMirror`와 `Prism`은 법선이 입사광 쪽을 향하도록 부호를 맞추는 절차가 따로 있다. 결과는 어느 쪽이든 같다 — 공식 자체가 부호 무관이기 때문이다.

## 2. 법선을 찾는 방법 네 가지

공식이 하나면, 함수마다 다른 부분은 법선을 어떻게 구하느냐뿐이다.

- `FlatMirror`: 거울이 직선이므로 법선도 상수 하나. `(a_e, b_e)`(직선의 계수)를 정규화하면 끝난다.
- `SphericalMirror`: 거울이 원호이므로 법선은 곡률중심 `(cx, cy)`에서 충돌점으로 향하는 반지름 방향.
- `BeamSplitter`: `FlatMirror`와 법선 계산이 같다. 부품을 만드는 함수부터 `make_beam_splitter`가 `make_flat_mirror`를 그대로 호출해서 만든다 — 빔스플리터는 물리적으로 반투명한 평면거울이라 모양을 정의하는 코드까지 같은 게 자연스럽다. 다만 반사광 하나만 만들지 않고, 입사 방향 그대로인 투과광까지 같이 반환한다.
- `ArbitraryMirror`: 경계가 다각형(선분들의 모음)이라, 변에 부딪히면 그 변에 수직인 방향이 법선이다. 꼭짓점에 정확히 부딪히면 얘기가 달라진다(3절에서 다룬다).

네 경우 모두 법선 기준 입사각과 반사각이 실제로 같은지 확인했다.

```python
def ang_err(vin, vout, n):
    ang_in = np.rad2deg(np.arccos(abs(np.dot(vin, n))))
    ang_out = np.rad2deg(np.arccos(abs(np.dot(vout, n))))
    return ang_in, ang_out, ang_in - ang_out
```

```text
FlatMirror        입사각=20.000000°  반사각=20.000000°  차이=0.0e+00°
SphericalMirror    입사각=7.662256°  반사각=7.662256°  차이=-1.9e-13°
BeamSplitter(반사)  입사각=45.000000° 반사각=45.000000°  차이=3.6e-14°
ArbitraryMirror(변) 입사각=26.565051° 반사각=26.565051°  차이=-1.4e-14°
```

네 함수 다 부동소수점 오차 수준에서 반사법칙을 만족한다.

<img src="/assets/img/posts/raytracing-reflection-family/fig1-reflection-family.png" alt="네 반사 부품의 광선 경로" width="700">
_그림1. 네 반사 부품의 광선 경로 — 반사 공식은 하나, 법선 찾는 법만 다르다_

## 3. 다각형 꼭짓점의 국소 곡률 근사와 그 오차

`ArbitraryMirror`는 매끈한 곡면 거울도 다각형으로 잘게 쪼개서 표현한다. 문제는 꼭짓점이다. 변에 수직인 법선은 변마다 값이 다르므로, 꼭짓점에서는 어느 쪽 변의 법선을 써야 할지가 애매하다. 그래서 이 함수는 꼭짓점에서 인접한 3개 점(이전 점, 꼭짓점, 다음 점)을 지나는 원을 하나 잡고, 그 원의 중심에서 꼭짓점으로 향하는 방향을 법선으로 쓴다. 점 3개는 원 하나를 유일하게 정하므로, 추가 파라미터 없이 국소 곡률을 근사할 수 있다.

```python
def _vertex_normal(x0, y0, x1, y1, x2, y2, xc, yc, dup_point):
    """꼭짓점(xc,yc)에서의 국소 곡률 법선 -- 인접 3점을 지나는 원의 중심에서 구한다.

    dup_point=True는 P에서 세 번째 점(x2,y2) 대신 두 번째 점(x1,y1)을 다시 쓴다.
    세 점이 서로 달라야 원이 유일하게 정해지므로, 이 경우 그 조건이 깨진다.
    dup_point=False가 세 번째 점을 올바르게 쓰는 버전이다.
    """
    M = np.array([[x0, y0, 1], [x1, y1, 1], [x2, y2, 1]], dtype=float)
    if dup_point:
        P = -np.array([x0**2 + y0**2, x1**2 + y1**2, x1**2 + y1**2], dtype=float)
    else:
        P = -np.array([x0**2 + y0**2, x1**2 + y1**2, x2**2 + y2**2], dtype=float)
    J = np.linalg.solve(M, P)
    acx, acy = -J[0] / 2, -J[1] / 2
    n = np.array([acx - xc, acy - yc])
    return n / np.linalg.norm(n)
```

이 함수는 `P`에서 세 번째 점 `(x2, y2)` 대신 `(x1, y1)`을 한 번 더 쓸 수 있게 `dup_point`로 분기해뒀다. 원이 유일하게 정해지려면 세 점이 서로 달라야 하는데, `dup_point=True`는 이 조건을 깬다. 두 변의 길이·기울기가 얼마나 다른지(대칭 → 비대칭 → 더 비대칭)를 바꿔가며, 두 경우가 법선 방향을 실제로 얼마나 다르게 만드는지 비교했다(`verify_vertex_normal.py`).

```text
대칭 꼭짓점      : n_dup=[-0.210, -0.978]  n_uniq=[0.000, -1.000]  차이=12.095°
비대칭 꼭짓점     : n_dup=[-0.057, -0.998]  n_uniq=[0.310, -0.951]  차이=21.319°
더 비대칭한 꼭짓점: n_dup=[0.161, -0.987]  n_uniq=[0.685, -0.728]  차이=33.983°
```

예상과 반대다. 대칭 꼭짓점이 가장 크게 어긋날 거라 생각했는데, 오히려 가장 작게 어긋나고 비대칭이 심해질수록 차이가 12° → 21° → 34°로 커진다. 이유를 따져보면 자연스럽다 — `dup_point=True`는 원을 정할 때 실제 세 번째 점 대신 꼭짓점 자신을 재사용하는데, 이 대체가 만드는 왜곡의 크기는 결국 "실제 세 번째 점이 꼭짓점에서 얼마나 멀리, 얼마나 다른 방향으로 떨어져 있는가"에 달려 있다. 두 변의 길이·기울기 차이가 커질수록(비대칭할수록) 대체된 점(꼭짓점 자신)과 실제 세 번째 점의 차이도 커지므로, 법선도 그만큼 더 크게 틀어진다. 대칭 꼭짓점에 수직으로 쏜 광선은 왔던 길을 그대로 되짚어 나가야 하는데(그림2), 세 번째 점을 재사용한 법선으로 계산하면 12° 어긋난 방향으로 튄다 — 비대칭 꼭짓점에서는 이 왜곡이 더 크다는 뜻이다.

<img src="/assets/img/posts/raytracing-reflection-family/fig2-vertex-normal-typo.png" alt="꼭짓점 법선 계산의 영향" width="600">
_그림2. 세 번째 점을 재사용하면(빨간 실선) 반사 방향이 12.1° 달라진다 — 세 점을 모두 쓰면(파란 점선) 입사 경로를 그대로 되짚는다_

## 4. 한계

반사 공식 자체는 매질과 무관해서 하나로 충분하지만, 정확성은 전부 법선을 찾는 로직에 달려 있다. `ArbitraryMirror`에서 봤듯이 변에 부딪히는 경우와 꼭짓점에 부딪히는 경우는 법선을 구하는 계산 방식 자체가 다르고, 그만큼 결과도 계산 방식에 민감하다.

반사율도 다루지 않는다. `BeamSplitter`조차 반사광과 투과광을 세기 배분 없이 각각 완전한 광선으로 만든다 — 실제 빔스플리터라면 반사:투과 비율(예: 50:50)만큼 두 광선의 세기가 나뉘어야 한다.

시리즈에서 아직 다루지 않은 부분도 있다. `Aperture`/`Stops`/`FlatSensor`는 반사도 굴절도 없이 광선을 막거나 세기만 기록하는 종류라 이번 반사 계열과는 성격이 다르다. `CompoundSphericalLens`(N면 렌즈)는 `SphericalLens`의 일반화라 이미 다룬 것과 원리가 겹친다.
