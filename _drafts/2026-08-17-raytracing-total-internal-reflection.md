---
title: 기하광학 2편 — 전반사(TIR) 판정과 렌즈에 적용하기
date: 2026-08-17 19:00:00 +0900
categories:
  - 광학
  - 기하광학
tags:
  - ray-tracing
  - geometric-optics
  - snells-law
  - python
  - matlab
description: 프리즘의 전반사(TIR) 판정 로직을 검증하고, 같은 로직을 구면 렌즈 굴절 함수의 출사면에도 적용했다. 반사 공식이 코드 전체에서 몇 번이나 반복되는지까지 정리한다.
math: true
---
1편에선 렌즈 굴절함수를 구현한 내용을 기록했다. 그리고 구현 한계점으로 이렇게 적었다.
"이번에 검증한 렌즈 굴절 함수는 전반사(total internal reflection, TIR) 여부를 판정하지 않는다."
Snell's law에서 $\sin\theta_2 = (n_1/n_2)\sin\theta_1$이 1을 넘으면 `arcsin`이 정의역을 벗어나는데, 렌즈 함수는 그 경우를 걸러내지 않고 `NaN`을 그냥 흘려보낸다. 반면 프리즘 함수(`CalculateRayPath_Prism`)는 애초에 이 경우를 처리하도록 짜여 있었다. 이번 편에서는 그 프리즘의 TIR 판정 로직을 먼저 검증하고, 같은 로직을 렌즈 함수의 출사면에도 옮겨 붙인다.

## 1. 프리즘의 전반사 판정 — 다각형 경계 위의 반사법칙

렌즈 함수와 프리즘 함수는 광선-면 교차를 찾는 방식부터 다르다. 렌즈는 면이 원호라서 `arccos`/`arcsin`으로 교차각을 직접 풀었지만(1편 참고), 프리즘의 경계는 꼭짓점을 잇는 선분들의 모음이라 직선-직선 교차와 벡터 연산만으로 충분하다. 그리고 렌즈 함수는 입사면 1번, 출사면 1번으로 끝나지만, 프리즘 함수는 더 이상 어떤 경계와도 만나지 않을 때까지 `while` 루프를 돈다 — 그 안에서 전반사가 몇 번이든 일어날 수 있다는 뜻이다.

```python
def trace_ray_through_prism(ray, prism, n_lens, n_air, max_bounces=20, verbose=False):
    """프리즘 경계와 더 이상 만나지 않을 때까지 광선을 추적한다."""
    xs, ys, vx, vy = map(float, ray)
    boundary = prism["BOUNDARY"]
    air2glass = True
    path = [[xs, ys, vx, vy]]

    for _ in range(max_bounces):
        # ... (가장 가까운 경계 교차점 (xc, yc)를 찾는 부분, 생략) ...

        # 면의 법선 -- 선분에 수직인 방향, 광선이 들어오는 쪽을 향하도록 부호를 맞춘다
        nx, ny = -(y2 - y1), (x2 - x1)
        nnorm = np.hypot(nx, ny)
        nx, ny = nx / nnorm, ny / nnorm
        if vx * nx + vy * ny < 0:
            nx, ny = -nx, -ny

        theta_in = np.arccos(np.clip((vx * nx + vy * ny) / np.hypot(vx, vy), -1, 1))
        n1, n2 = (n_air, n_lens) if air2glass else (n_lens, n_air)
        sin_out = n1 * np.sin(theta_in) / n2

        if abs(sin_out) <= 1:
            theta_out = np.arcsin(sin_out)
            # ... Snell's law로 굴절 방향 계산 ...
            air2glass = not air2glass
        else:
            vout = np.array([vx, vy]) - 2 * (vx * nx + vy * ny) * np.array([nx, ny])

        vx, vy = vout / np.linalg.norm(vout)
        xs, ys = xc + 0.5 * vx, yc + 0.5 * vy  # 자기 자신과 재충돌 방지
        path.append([xc, yc, vx, vy])
    return np.array(path)
```

`abs(sin_out) <= 1`로 정의역을 먼저 확인하고, 벗어나면 굴절 대신 그 면의 법선 $\hat n$을 기준으로 반사시킨다.

$$ \vec v' = \vec v - 2(\vec v \cdot \hat n)\hat n $$

이 코드로 정삼각형 N-BK7 프리즘(한 변 600mm)에 수직 입사광을 쏴봤다.

<img src="/assets/img/posts/raytracing-total-internal-reflection/fig1-prism-tir-path.png" alt="정삼각형 프리즘에서의 전반사 경로" width="600">
_그림1. 수직 입사 광선이 정삼각형 프리즘 안에서 한 번 전반사한다_

처음엔 당연히 프리즘을 그냥 통과할 거라고 예상했는데, 실제로 돌려보면 안쪽 면에서 전반사가 걸리고 나서야(내부입사각 60°) 다른 면으로 무굴절 출사한다(내부입사각 0°, 마침 그 면에 수직으로 부딪힌다).

```text
face=(0.0,-300.0)-(0.0,300.0)  theta_in=0.000deg  투과  theta_out=0.000deg
face=(0.0,300.0)-(519.6,0.0)  theta_in=60.000deg  전반사(TIR)
face=(519.6,0.0)-(0.0,-300.0)  theta_in=0.000deg  투과  theta_out=0.000deg
```

이유는 간단하다. N-BK7의 임계각은 $\theta_c = \arcsin(1/n) = 41.41°$인데, 정삼각형의 꼭지각은 60°다. 입사면에 수직으로 들어온 광선은 굴절 없이 그대로 직진하므로, 다음 면에는 정확히 $60° - 0° = 60°$의 내부입사각으로 부딪힌다. $60° > 41.41°$이니 전반사할 수밖에 없다. 오히려 프리즘을 그냥 투과시키려면 입사각을 충분히 크게 줘야 한다는 뜻이다. 입사각 $\theta_1$, 굴절각 $\theta_1'$, 꼭지각 $A$ 사이에는 $\theta_1' + \theta_2' = A$ 관계가 성립하므로, 두 번째 면에서 임계각을 넘지 않으려면

$$ \theta_1' > A - \theta_c, \qquad \theta_1 > \arcsin\big(n \sin(A - \theta_c)\big) $$

이 조건을 만족해야 한다. 여기 대입하면 $\theta_1 > 28.8131°$가 나온다. 코드로 입사각을 0°부터 이분법으로 스윕해서 실제 코드가 관측하는 전반사↔투과 경계를 찾아보면(`verify_tir.py`) 28.813100°로, 이론값과 오차 7.6e-10° 이내로 일치한다.

<img src="/assets/img/posts/raytracing-total-internal-reflection/fig2-prism-angle-sweep.png" alt="입사각 스윕에 따른 전반사/투과" width="600">
_그림2. 입사각이 임계값(28.8°)을 넘으면 전반사 대신 투과한다_

이건 실은 낯선 현상이 아니다. 거울 코팅 없이 프리즘 내부의 전반사만으로 빛의 경로를 꺾는 porro 프리즘(쌍안경에 쓰이는 그 프리즘)이 정확히 이 원리로 동작한다 — 45°로 꺾인 두 면에 빛을 수직으로 넣으면 45° > 41.4°(BK7 임계각)라 반드시 전반사하고, 코팅 없이도 반사면 역할을 한다.

## 2. 렌즈에도 같은 로직을 넣으면? — 출사면 전반사 처리

렌즈 함수의 출사면 굴절 계산은 이렇게 생겼다(1편에서 그대로 가져온 부분).

```python
phi2 = np.arctan2(vy1, vx1)
theta_t2 = np.arcsin(n_lens * np.sin(w2 - phi2) / n_air)  # 도메인 체크가 없다
vx2, vy2 = -np.cos(theta_t2 + w2), -np.sin(theta_t2 + w2)
```

`n_lens * sin(w2 - phi2) / n_air`가 1을 넘으면, 즉 내부입사각이 임계각을 넘으면 이 `arcsin`은 애초에 실수해가 없다. 원본 코드는 이 경우를 걸러내지 않고 그대로 넘기기 때문에, 임계각을 넘는 광선에서는 조용히 `NaN`이 나온다.

이 함수는 광선이 렌즈의 어느 면으로 먼저 들어오는지 스스로 판단하지 않는다 — 호출하는 쪽이 `entry_surface`를 지정해야 한다. R1/R2는 렌즈를 정의할 때의 매개변수 순서일 뿐 광원에서 봤을 때 어느 면이 먼저인지와는 무관하므로(1편 참고), 항상 `find_entry_surface`로 먼저 판정한 값을 넘겨야 한다 — 그렇지 않으면 전반사와 무관한 자리(출사면 교차점 계산)에서 `NaN`이 날 수도 있다.

평철렌즈(R1≈평면, R2=-50mm, T=20mm, D=78mm)의 y0=36mm 광선으로 확인해보면:

```text
find_entry_surface: 1번 면

수정 전:
[[-9.99935200e+00  3.60000000e+01  1.00000000e+00 -1.21878769e-05]
 [-5.30123745e+00  3.59999427e+01             nan             nan]]
-> RuntimeWarning: invalid value encountered in arcsin
```

프리즘과 같은 패턴을 적용해 고쳤다 — 정의역을 벗어나면 굴절 대신, 출사면의 법선 $(\cos w_2, \sin w_2)$(원의 매개변수 각을 그대로 쓴다, 접선의 방향이 곧 반지름 방향이니까) 기준으로 반사시킨다.

```python
sin_t2 = n_lens * np.sin(w2 - phi2) / n_air
if abs(sin_t2) <= 1:
    theta_t2 = np.arcsin(sin_t2)
    vx2, vy2 = -np.cos(theta_t2 + w2), -np.sin(theta_t2 + w2)
else:
    n2x, n2y = np.cos(w2), np.sin(w2)
    dot2 = vx1 * n2x + vy1 * n2y
    vx2, vy2 = vx1 - 2 * dot2 * n2x, vy1 - 2 * dot2 * n2y
```

```text
수정 후:
[[-9.99935200e+00  3.60000000e+01  1.00000000e+00 -1.21878769e-05]
 [-5.30123745e+00  3.59999427e+01  3.68088814e-02 -9.99322324e-01]]
```

이 결과가 진짜 반사인지 두 가지로 검증했다(`verify_tir.py`). 먼저 법선 기준 반사법칙(입사각=반사각)이 성립하는지: 46.05508420° vs 46.05508420°, 오차 7.1e-15°. 다음으로, TIR이 시작되는 경계(y0)를 코드로 이분법 스윕해서 찾은 뒤, 그 경계에서의 내부입사각을 이론 임계각과 비교했다: 코드가 찾은 경계는 y0=33.072057mm이고, 그 지점의 내부입사각은 41.410389°로 이론값 $\theta_c=41.410389°$와 오차 2.7e-10° 이내로 일치한다.

<img src="/assets/img/posts/raytracing-total-internal-reflection/fig3-lens-marginal-tir.png" alt="평철렌즈 마진광선의 전반사" width="600">
_그림3. 평철렌즈(R2=-50mm) 마진광선 39개 중 6개는 출사면에서 전반사한다_

전체 조리개(D=78mm)의 절반(39mm) 중, 실제로 전반사가 걸리는 구간은 가장자리 33mm 지점부터다 — 조리개의 바깥쪽 15% 정도에서만 나타나는 현상이라는 뜻이다. 곡률이 완만한 렌즈라면 이 구간 자체가 존재하지 않거나 조리개 훨씬 바깥에 있어서 실무에서는 잘 안 보이지만, 이번 렌즈처럼 R2가 D/2에 가까운(근반구형에 가까운) 강한 곡률에서는 무시할 수 없다.

## 3. 한계

프리즘 함수와 다르게, 이번에 고친 렌즈 함수는 여전히 "진입 1회 + 퇴장 1회"만 모델링한다. 출사면에서 전반사한 광선의 방향은 이제 정확히 계산되지만, 그 광선이 렌즈 안에서 다음에 어디로 가는지는 이 함수에서 구현되지 않았다. — 프리즘처럼 경계와 더 이상 만나지 않을 때까지 도는 루프가 없다. 실제로 이 반사된 광선을 계속 추적하려면 렌즈 함수 자체를 프리즘처럼 다시 짜야 한다.

마지막으로 하나 더 확인한 게 있다. 이번 편에서 반사에 두 번(프리즘의 내부 반사, 렌즈의 TIR) 쓴 공식 $\vec v' = \vec v - 2(\vec v\cdot\hat n)\hat n$은 사실 이 라이브러리 전체에서 정확히 같은 형태로 여섯 번 나온다 — `FlatMirror`, `SphericalMirror`, `ArbitraryMirror`, `BeamSplitter`(반사 성분), 그리고 이번에 고친 `Prism`/`SphericalLens`의 TIR 분기까지. 함수마다 다른 건 오직 충돌점과 법선을 어떻게 찾느냐(평면 대수 / 원의 후보해 / 다각형 위 국소 곡률)와, 한 번만 반사하고 끝나는지 여러 번 반사할 때까지 도는지 뿐이다. 굴절은 매번 다른 각도 공식이 필요했지만, 반사는 거울이든 프리즘이든 렌즈든 이 벡터 공식 하나로 끝난다. 3편에서는 이 반사 공식이 코드 전체에서 정말로 일관되게 구현돼 있는지, `FlatMirror`/`SphericalMirror`/`ArbitraryMirror`/`BeamSplitter`를 한 번에 검증한다.
