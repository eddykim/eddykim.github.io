---
title: 기하광학 1편 — 광선의 구면 굴절 계산과 검증
date: 2026-08-15 19:00:00 +0900
categories:
  - 광학
  - 기하광학
tags:
  - ray-tracing
  - geometric-optics
  - snells-law
  - python
  - matlab
description: MATLAB 2D 광선추적기의 구면 렌즈 굴절 계산을 Python으로 재구성했다. 광선과 구면의 교차점을 구하는 방법부터, 광학 가역성과 thick-lens 공식으로 결과를 검증하는 과정까지 정리한다.
math: true
---
학교나 연구실에서 레이트레이싱 혹은 광학 설계 툴을 지원하지 않는 경우가 종종 있다. 내가 속했던 연구실이 그러했는데 필요성이 그리 크진 않았다. 설계하는 광학계가 Infinite corrected(설명필요) 렌즈를 사용하기도 하고 특정 파장 영역, 특히 사용하는 가시광 영역을 보장하는 광 부품을 사용하기도 했기에 그러했다. 하지만 점점 더 정밀한 측정 데이터가 필요해지고 고배율 렌즈의 후초평면(back focal plane) 을 관찰할 일이 많아짐에 따라 연구실에서는 자체적으로 사용할, 지그 제작 및 초점거리를 특정할 레이트레이싱 툴을 구현하고 사용했었다.

기본 골자는 MATLAB 오픈소스 코드였고, 이를 Python으로 옮기고 다듬는 것이 내가 했던 일이다. 함수 수가 60개가 넘고 렌즈·미러·프리즘·빔스플리터·조리개까지 부품 종류도 많아서, 다 옮기고 나서 든 생각은 "그래서 이게 제대로 옮겨진 걸 어떻게 확인하지"였다. 겉보기엔 그럴듯한 광선 그림이 나와도, 굴절 각도가 미묘하게 틀렸거나 렌즈 가장자리 근처에서 엉뚱한 점을 골랐다면 눈으로는 못 잡아낸다. 이번 포스트들에서는 기하광학에 관한 간략한 소개와 함께, 부품별 구현 원리와 그 검증 방법을 하나씩 정리한다.

그래서 가장 기초가 되는 개념— 광선이 구면(spherical surface)을 지나며 굴절 — 을 가지고 원리부터 다시 짚어가며 검증한 과정을 정리한다. 이 글의 코드는 검증에 필요한 부분만(`make_spherical_lens`, `refract_through_lens` 등) 뽑아 재구성 했다(`raytracing-spherical-lens-refraction`).

## 1. 광선과 구면의 교차 — 두 개의 후보해

렌즈 하나를 통과하는 광선을 계산하려면 제일 먼저 "광선이 렌즈의 어느 면을 어디서 만나는가"부터 알아야 한다. 코드에서 광선은 `[x, y, vx, vy]` 네 값으로 표현되는데, $(x_0, y_0)=(x, y)$는 지금 계산하는 면을 향해 출발하는 위치, $(v_x, v_y)$는 그 방향벡터다. 렌즈의 첫 면(광원에서 먼저 만나는 면)을 계산할 때는 $(x_0, y_0)$가 광원의 위치이고, 렌즈를 통과해 두 번째 면을 계산할 때는 첫 면에서 방금 굴절되어 나온 점이 새로운 $(x_0, y_0)$가 된다 — 즉 광원 자체가 아니라 "이 면을 향해 지금 광선이 출발하는 자리"다. 이 시작점과 방향벡터를 직선의 일반형으로 쓰면 다루기 편하다.

$$ ax + by = c, \qquad a = -v_y,\ \ b = v_x,\ \ c = a x_0 + b y_0 $$

방향벡터에 수직인 $(a,b)=(-v_y, v_x)$를 계수로 쓰면 광선이 수직이든 수평이든 예외 없이 이 형태로 표현된다. 렌즈의 한쪽 면은 중심 $(c_x, c_y)$, 곡률반경 $R$인 원의 일부이므로, 원 위의 점을 $(R\cos w + c_x,\ R\sin w + c_y)$로 매개변수화해서 직선의 방정식에 대입하면 각도 $w$에 대한 삼각방정식이 나온다. 직선과 원은 최대 두 점에서 만나므로, 이 방정식도 보통 두 개의 해를 갖는다.

문제는 이 두 해 중 실제로 광선이 렌즈 면에 부딪히는 지점은 하나뿐이라는 것이다. 나머지 하나는 광선을 반대로 연장했을 때 만나는 점이거나, 렌즈의 물리적 지름 $D$ 바깥(구면의 연장선 위, 렌즈가 실제로 존재하지 않는 자리)에 있는 점이다. 그래서 이 둘을 정답이 아니라 "후보해"라고 부르고, 뒤에서 렌즈 경계 안에 있는지를 기준으로 진짜 교차점을 가려낸다. 코드에서는 이 두 후보해를 `arccos`와 `arcsin` 두 가지 삼각함수 항등식으로 각각 구한다.

```python
def find_entry_surface(ray, lens):
    """광선이 렌즈의 두 면 중 어느 쪽에 먼저 닿는지 찾는다.

    광선을 ax+by=c 직선으로, 각 면의 경계를 선분들의 모음으로 보고 교차를
    검사한다. 반환값은 (면 번호(1 또는 2), 교차까지 거리) 또는 못 맞았으면 None.
    """
    x, y, vx, vy = ray
    a_r, b_r = -vy, vx
    c_r = a_r * x + b_r * y

    best = None
    for surf_idx, boundary in enumerate(lens["BOUNDARY"], start=1):
        for i in range(len(boundary) - 1):
            x1, y1 = boundary[i]
            x2, y2 = boundary[i + 1]
            a_e, b_e = y2 - y1, x1 - x2
            c_e = a_e * x1 + b_e * y1

            M = np.array([[a_r, b_r], [a_e, b_e]])
            if abs(np.linalg.det(M)) < 1e-8:
                continue
            J = np.linalg.solve(M, np.array([c_r, c_e]))
            ...
```

여기서는 렌즈 경계를 직선(원의 근사 다각형)들의 모음으로 보고 교차를 검사해서 "어느 면에 먼저 닿는가"만 결정한다. 실제 굴절 계산(`refract_through_lens`)에서는 원 방정식을 정확히 풀어 두 후보각 $w_c$(arccos), $w_s$(arcsin)를 구하고, 그중 렌즈의 물리적 지름 $D$ 안에 들어오는 쪽을 고른다.

```python
def _pick_candidate(boundary, offset, wc, xc, yc, ws, xs, ys, ref_x, ref_y):
    """arccos/arcsin 두 후보해 중 렌즈 물리 경계 안에 있는 쪽을 고른다.

    두 후보 다 경계 안이면(보통의 경우) 기준점에 더 가까운 쪽을 쓴다.
    """
    def inside(px, py):
        return (
            px <= np.max(boundary[:, 0]) + offset and px >= np.min(boundary[:, 0]) - offset
            and py <= np.max(boundary[:, 1]) + offset and py >= np.min(boundary[:, 1]) - offset
        )

    in_c, in_s = inside(xc, yc), inside(xs, ys)
    if in_c and in_s:
        if np.hypot(xc - ref_x, yc - ref_y) < np.hypot(xs - ref_x, ys - ref_y):
            return wc, xc, yc
        return ws, xs, ys
    elif in_c:
        return wc, xc, yc
    else:
        return ws, xs, ys
```

구면은 평면이 아니므로, 렌즈 정점(vertex)에서의 접평면을 기준으로 재면 광축에서 멀어질수록 그 평면으로부터 점점 더 파고 들어간다. 이 파고 들어간 깊이를 sag(처짐량)라고 부르고, 렌즈 가장자리($D/2$)에서는 $\mathrm{sag} = R - \sqrt{R^2 - (D/2)^2}$로 계산된다. `offset`은 이 sag에 작은 배율(0.01)을 곱한 여유값이다. 부동소수점 계산 오차 때문에 렌즈 경계에 딱 걸쳐 있는 점이 "경계 밖"으로 잘못 판정되는 걸 막아준다.

교차점을 찾고 나면 그 지점에서 Snell's law를 적용한다.

$$ n_1 \sin\theta_1 = n_2 \sin\theta_2 $$

여기서 $\theta_1, \theta_2$는 각각 입사·굴절 광선이 그 지점의 표면 법선과 이루는 각이다. 코드는 원의 매개변수 각 $w$를 이용해 법선 방향을 얻고, 입사 방향벡터와의 상대각을 구해 $\theta_1$을 만든 다음, Snell's law로 $\theta_2$를 구해 새 방향벡터로 되돌린다.

```python
theta_t1 = np.arcsin(n_air * np.sin(np.pi - w1 + np.arctan2(vy, vx)) / n_lens)
vx1, vy1 = np.cos(theta_t1 - (np.pi - w1)), np.sin(theta_t1 - (np.pi - w1))
```

이 각도항의 세부 유도(원의 어느 매개변수화를 썼는지에 따라 달라지는 부분)까지 여기서 다시 풀지는 않는다. 대신 이 식이 실제로 Snell's law를 만족하는 결과를 내는지는 2절에서 광학 가역성과 이론값 대조로 수치 검증한다.

이 렌즈를 실제로 통과시켜 보면 이렇게 나온다.

<img src="/assets/img/posts/raytracing-spherical-lens-refraction/fig1-lens-focusing.png" alt="평행광선이 렌즈를 지나 초점에 모이는 모습" width="600">
_그림1. 평행광선 11개가 양볼록 렌즈(R1=1000mm, R2=-1000mm, t=100mm, D=500mm, N-BK7)를 지나 한 점으로 모인다_

렌즈는 원점에 놓았고, 광원은 $x=1000$mm에서 $-x$ 방향으로 진행하는 평행광 11가닥(지름 200mm)이다. 여기서 처음엔 헷갈렸던 점 하나 — 광원이 $+x$에서 $-x$로 진행하니 당연히 "R1 면"(첫 번째로 정의한 면)을 먼저 만날 거라 생각했는데, 실제로 코드를 돌려보면 `find_entry_surface`가 반환하는 진입면은 2번(R2)이다. R1/R2는 렌즈를 정의할 때 매개변수 순서일 뿐 "광원에서 봤을 때 먼저"라는 물리적 좌우와 직접 대응하지 않는다 — 이 렌즈는 R1=1000(볼록), R2=-1000(볼록)인데, 두 면의 정점(vertex) 좌표를 직접 계산해보면 R1 면이 $x=-t/2$, R2 면이 $x=+t/2$에 있다. 즉 $+x$ 쪽에서 오는 광선은 R2 면을 먼저 만난다. 이 부호 규약은 2절의 가역성 검증에서 다시 나온다.

## 2. 검증 — 가역성과 이론값 대조

### 광학 가역성

광학 가역성(optical reversibility)이란 광선이 지나온 경로를 반대 방향으로 그대로 되짚어 갈 수 있다는 원리다. Snell's law가 시간 역전에 대해 대칭이기 때문에 성립한다 — 굴절된 광선의 진행 방향을 반대로 뒤집어 쏘면, 원래 들어왔던 경로를 정확히 되짚어 나가야 한다. 이건 `refract_through_lens`의 두 분기(`entry_surface==1`일 때와 그 반대일 때)가 서로 모순 없이 같은 물리를 구현했는지 보는 좋은 교차검증이 된다 — 정방향 진입은 한쪽 분기를, 역추적은 반대쪽 분기를 타기 때문이다.

```python
ray_fwd = np.array([1000.0, 100.0, -1.0, 0.0])
entry_fwd = find_entry_surface(ray_fwd, lens)[0]
path_fwd = refract_through_lens(ray_fwd, entry_fwd, lens, n_lens, N_AIR)

# 출사 광선을 렌즈 밖으로 500mm 더 진행시킨 지점에서, 방향만 반대로 뒤집는다
x2, y2, vx2, vy2 = path_fwd[-1]
far_point = np.array([x2 + 500 * vx2, y2 + 500 * vy2])
ray_retrace = np.array([far_point[0], far_point[1], -vx2, -vy2])
entry_retrace = find_entry_surface(ray_retrace, lens)[0]
path_retrace = refract_through_lens(ray_retrace, entry_retrace, lens, n_lens, N_AIR)
```

이렇게 정방향과 역추적을 나란히 돌려보면 다음과 같이 나온다.

```text
정방향 진입면=2, 역추적 진입면=1 (반대쪽 면이어야 함)
원래 입사 광선 시작점: (1000.000, 100.000), 방향 (-1.000, 0.000)
역추적 도착점(렌즈 첫 굴절점): (44.987, 100.000), 방향 (1.000, -0.000)
위치 오차: 0.00e+00 mm
```

진입면이 정방향과 정확히 반대(2 → 1)로 나오고, 역추적한 광선이 원래 입사 광선의 시작 $y$좌표(100mm)와 방향(부호만 반대)을 부동소수점 오차 없이 그대로 복원한다.

<img src="/assets/img/posts/raytracing-spherical-lens-refraction/fig2-reversibility.png" alt="가역성 검증" width="600">
_그림2. 출사 광선을 반대로 쏘면 입사 경로를 정확히 되짚는다 (파란 실선과 빨간 점선이 렌즈~x=1000 구간에서 완전히 겹친다)_

### thick-lens 공식과 대조 — 근축은 거의 정확히, 가장자리는 어긋난다

가역성 검증은 `refract_through_lens`의 두 분기가 서로 일관되게 동작하는지는 보여주지만, 그 결과가 실제 렌즈의 물리량과 맞는지는 확인해주지 않는다 — 두 분기가 똑같이 틀렸다면 가역성은 그래도 성립하기 때문이다. 그래서 이 광선추적 코드와 완전히 무관하게 유도된 독립적인 이론식과 대조해본다. 렌즈의 두께까지 고려한 thick-lens 공식으로 초점거리와 후방초점거리(back focal distance)를 구할 수 있다.

$$ \frac{1}{f} = (n-1)\left[\frac{1}{R_1} - \frac{1}{R_2} + \frac{(n-1)t}{n R_1 R_2}\right], \qquad \mathrm{BFD} = f\left[1 - \frac{(n-1)t}{n R_1}\right] $$

이 공식은 근축광선(paraxial ray, 광축에 아주 가까운 광선) 근사 위에서 유도된 것이다. 광선추적 시뮬레이션은 근사 없이 실제 기하로 계산하므로, 광축에서 조금이라도 떨어진(marginal) 광선은 이 공식과 정확히 일치할 이유가 없다 — 오히려 정확히 일치한다면 그게 이상한 것이다(구면수차, spherical aberration).

```python
n = n_lens  # N-BK7, 750nm
inv_f = (n - 1) * (1 / R1 - 1 / R2 + (n - 1) * T / (n * R1 * R2))
f = 1 / inv_f
bfd = f * (1 - (n - 1) * T / (n * R1))
expected_focus_x = -T / 2 - bfd
```

이 렌즈에 계산해보면 이론값과 시뮬레이션이 이렇게 갈린다.

```text
n(N-BK7, 750nm) = 1.511835
thick-lens 이론 초점(근축): f=993.698mm, 후방초점거리=960.056mm, 초점 x=-1010.056
시뮬레이션 근축 초점(|y0|<=20mm) x=-1010.229, 이론 대비 -0.172mm (-0.017%)
시뮬레이션 최외곽 초점(|y0|=100mm) x=-995.010, 근축 대비 +15.219mm (구면수차)
```

광축에서 20mm 이내(렌즈 반지름 250mm의 8%)로 들어오는 근축 광선은 이론값과 0.017% 차이로 사실상 정확히 일치한다. 반면 조리개 가장자리(광축에서 100mm)로 들어오는 광선은 근축 초점보다 15.2mm 렌즈 쪽으로 당겨져 초점을 맺는다 — 단순 양볼록 렌즈에서 흔히 보는 undercorrected 구면수차 방향(가장자리 광선일수록 초점이 짧아짐)과 일치한다.

<img src="/assets/img/posts/raytracing-spherical-lens-refraction/fig3-focus-vs-theory.png" alt="초점 위치 비교" width="600">
_그림3. 근축 이론값(검은 점선)과 시뮬레이션의 근축 초점(파란 점선)은 거의 겹치고, 최외곽 광선 초점(빨간 점선)만 구면수차만큼 어긋난다_

## 3. 한계

이번 포스트에서는 광선이 구면을 지나며 굴절하는 계산 — 직선-원 교차로 렌즈 면과의 진입점을 찾고, 그 지점에서 Snell's law로 굴절 방향을 구하는 부분(`refract_through_lens`) — 을 Python으로 재구성하고, 가역성 검증과 thick-lens 이론값 대조라는 두 가지 서로 다른 방식으로 그 결과가 맞는지 확인했다. 가역성 검증에서는 정방향과 역추적 경로가 부동소수점 오차 없이 일치했고, thick-lens 이론값과의 대조에서는 근축 초점이 0.017% 오차로 사실상 일치하면서 가장자리 광선에서는 예상된 방향의 구면수차만큼만 벗어났다. 두 검증이 서로 독립적인 방법으로 같은 결론(구현이 맞다)에 도달했다는 점에서, `refract_through_lens`의 구면 굴절 계산은 신뢰할 만하다고 본다.

다만 이번에 검증한 `refract_through_lens`는 전반사(total internal reflection) 여부를 판정하지 않는다. Snell's law에서 $\sin\theta_2 = (n_1/n_2)\sin\theta_1$이 1을 넘어서면 `arcsin`이 정의역을 벗어나는데, 이 함수는 그 경우를 걸러내지 않고 그냥 NumPy가 계산한 `NaN`을 다음 계산으로 흘려보낸다. MATLAB 원본도 이 구면렌즈 함수에서는 마찬가지였다 — 전반사를 명시적으로 판정하는 건 프리즘 계산 함수 쪽뿐이었다. 렌즈는 대개 공기보다 굴절률이 높은 매질로 들어가는 각이 완만해서 전반사 조건에 잘 안 걸리지만, 곡률이 크고 조리개가 넓은 조합에서는 걸릴 수 있다. 다음 편에서 프리즘의 전반사 판정 로직과, 그걸 프리즘이 아닌 렌즈에도 적용하면 어떻게 되는지를 다룬다.
