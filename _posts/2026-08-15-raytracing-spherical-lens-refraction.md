---
title: 기하광학 1편 — 구면 렌즈 굴절, 광선추적 검증하기
date: 2026-08-15 19:00:00 +0900
categories: [광학, 기하광학]
tags: [ray-tracing, geometric-optics, snells-law, python, matlab]
description: 20년 된 MATLAB 2D 광선추적기를 Python으로 옮기면서, 구면 렌즈 굴절 하나를 붙잡고 "제대로 옮겼다"를 어떻게 확신할 수 있는지 정리한다. 실행 중 만난 버그, 소스 대조로 찾은 버그, 광학 가역성과 thick-lens 공식으로 한 교차검증까지 담았다.
math: true
---

예전에 쓰던 2D 광선추적(ray tracing) MATLAB 코드를 Python으로 옮길 일이 있었다. 함수 수가 60개가 넘고, 렌즈·미러·프리즘·빔스플리터·조리개까지 부품 종류도 많다. 다 옮기고 나서 든 생각은 "그래서 이게 제대로 옮겨진 걸 어떻게 확인하지"였다. 겉보기엔 그럴듯한 광선 그림이 나와도, 굴절 각도가 미묘하게 틀렸거나 렌즈 가장자리 근처에서 엉뚱한 점을 골랐다면 눈으로는 못 잡아낸다.

그래서 가장 기초가 되는 조각 하나 — 광선이 구면(spherical surface) 하나를 지나며 굴절하는 계산 — 을 붙잡고, 원리부터 다시 짚어가며 검증한 과정을 정리한다. 이 글의 코드는 원본 MATLAB 프로젝트 전체가 아니라 이 검증에 필요한 부분만(`make_spherical_lens`, `refract_through_lens` 등) 뽑아 재구성한 것이다(`_code/raytracing-spherical-lens-refraction/`).

## 1. 광선과 구면의 교차 — 두 개의 후보해

렌즈 하나를 통과하는 광선을 계산하려면 제일 먼저 "광선이 렌즈의 어느 면을 어디서 만나는가"부터 알아야 한다. 광선은 위치 $(x_0, y_0)$와 방향벡터 $(v_x, v_y)$로 주어지는데, 이걸 직선의 일반형으로 쓰면 다루기 편하다.

$$ ax + by = c, \qquad a = -v_y,\ \ b = v_x,\ \ c = a x_0 + b y_0 $$

방향벡터에 수직인 $(a,b)=(-v_y, v_x)$를 계수로 쓰면 광선이 수직이든 수평이든 예외 없이 이 형태로 표현된다. 렌즈의 한쪽 면은 중심 $(c_x, c_y)$, 곡률반경 $R$인 원의 일부이므로, 원 위의 점을 $(R\cos w + c_x,\ R\sin w + c_y)$로 매개변수화해서 직선의 방정식에 대입하면 각도 $w$에 대한 삼각방정식이 나온다. 직선과 원은 최대 두 점에서 만나므로, 이 방정식도 보통 두 개의 해를 갖는다.

코드에서는 이 두 해를 `arccos`와 `arcsin` 두 가지 삼각함수 항등식으로 각각 구한다.

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

`offset`은 렌즈 가장자리($D/2$)에서 구면과 그 접평면 사이의 sag(처짐량)에 비례하는 여유값이다. 부동소수점 계산 오차 때문에 렌즈 경계에 딱 걸쳐 있는 점이 "경계 밖"으로 잘못 판정되는 걸 막아준다.

교차점을 찾고 나면 그 지점에서 Snell's law를 적용한다.

$$ n_1 \sin\theta_1 = n_2 \sin\theta_2 $$

여기서 $\theta_1, \theta_2$는 각각 입사·굴절 광선이 그 지점의 표면 법선과 이루는 각이다. 코드는 원의 매개변수 각 $w$를 이용해 법선 방향을 얻고, 입사 방향벡터와의 상대각을 구해 $\theta_1$을 만든 다음, Snell's law로 $\theta_2$를 구해 새 방향벡터로 되돌린다.

```python
theta_t1 = np.arcsin(n_air * np.sin(np.pi - w1 + np.arctan2(vy, vx)) / n_lens)
vx1, vy1 = np.cos(theta_t1 - (np.pi - w1)), np.sin(theta_t1 - (np.pi - w1))
```

이 각도항의 세부 유도(원의 어느 매개변수화를 썼는지에 따라 달라지는 부분)까지 여기서 다시 풀지는 않는다. 대신 이 식이 실제로 Snell's law를 만족하는 결과를 내는지는 3절에서 광학 가역성과 이론값 대조로 수치 검증한다.

이 렌즈를 실제로 통과시켜 보면 이렇게 나온다.

![평행광선이 렌즈를 지나 초점에 모이는 모습](/assets/img/posts/raytracing-spherical-lens-refraction/fig1-lens-focusing.png){: width="600" }
_그림1. 평행광선 11개가 양볼록 렌즈(R1=1000mm, R2=-1000mm, t=100mm, D=500mm, N-BK7)를 지나 한 점으로 모인다_

렌즈는 원점에 놓았고, 광원은 $x=1000$mm에서 $-x$ 방향으로 진행하는 평행광 11가닥(지름 200mm)이다. 여기서 처음엔 헷갈렸던 점 하나 — 광원이 $+x$에서 $-x$로 진행하니 당연히 "R1 면"(첫 번째로 정의한 면)을 먼저 만날 거라 생각했는데, 실제로 코드를 돌려보면 `find_entry_surface`가 반환하는 진입면은 2번(R2)이다. R1/R2는 렌즈를 정의할 때 매개변수 순서일 뿐 "광원에서 봤을 때 먼저"라는 물리적 좌우와 직접 대응하지 않는다 — 이 렌즈는 R1=1000(볼록), R2=-1000(볼록)인데, 두 면의 정점(vertex) 좌표를 직접 계산해보면 R1 면이 $x=-t/2$, R2 면이 $x=+t/2$에 있다. 즉 $+x$ 쪽에서 오는 광선은 R2 면을 먼저 만난다. 이 부호 규약은 3절의 가역성 검증에서 다시 나온다.

## 2. 처음 돌렸을 때 만난 버그 둘

### IndexError — BD_OFFSET을 2차원으로 착각

`make_spherical_lens`가 만드는 `BD_OFFSET`은 두 면의 sag 여유값을 담은 1차원 배열이다.

```python
bd_offset = 0.01 * np.array(
    [np.abs(R1) - np.sqrt(R1**2 - (D / 2) ** 2), np.abs(R2) - np.sqrt(R2**2 - (D / 2) ** 2)]
)
```

그런데 MATLAB 원본은 `BD_OFFSET(1,1)`, `BD_OFFSET(1,2)`처럼 1×2 행벡터로 인덱싱한다. 처음 포팅할 때 이 표기를 그대로 옮기면서 `BD_OFFSET[0, 0]`(2차원 인덱싱)을 썼다. 최소 재현 코드로 실행하면 이렇게 난다.

```python
BD_OFFSET = 0.01 * np.array([1.234, 2.345])  # make_spherical_lens가 만드는 형태와 동일
offset = BD_OFFSET[0, 0]
```

```text
Traceback (most recent call last):
  File "bug_demo.py", line 21, in <module>
    offset = BD_OFFSET[0, 0]
             ~~~~~~~~~^^^^^^
IndexError: too many indices for array: array is 1-dimensional, but 2 were indexed
```

MATLAB은 스칼라도 기본적으로 1×1 행렬로 다루기 때문에 `(1,1)` 인덱싱이 자연스럽지만, NumPy 배열은 실제로 선언된 차원만큼만 인덱싱을 받는다. 고치는 건 인덱스를 하나만 쓰면 된다.

```python
offset = BD_OFFSET[0]  # 1차원 배열이므로 인덱스 하나만 필요
```

### 소스 대조로 찾은 버그 — 출구면 경계조건 누락

두 번째는 실행하다 걸린 게 아니라, MATLAB 원본과 한 줄씩 대조하다 찾았다. 렌즈에 처음 들어가는 면(R1 또는 R2)에서는 `_pick_candidate`가 두 후보 중 렌즈 경계 안에 있는 쪽을 정확히 골랐는데, 렌즈를 빠져나가는 반대쪽 면에서는 이 경계조건 검사 없이 그냥 "기준점에 더 가까운 후보"만으로 정했다.

두 후보 $w_c$, $w_s$는 직선-원 교차의 두 근(가까운 해/먼 해에 해당)인데, 보통은 "더 가까운 쪽"이 실제 렌즈 위의 점과 일치한다. 하지만 이건 어디까지나 경험적 휴리스틱이고, 렌즈 곡률과 광선 각도에 따라 렌즈 바깥(구면의 연장선 위, 물리적으로 렌즈가 끝난 자리)에 있는 후보가 오히려 더 가까운 경우가 있을 수 있다 — 경계조건 검사는 바로 이런 경우를 걸러내려고 있는 안전장치다. 처음 포팅한 코드는 광선이 들어가는 면(entry surface)에서는 이 검사를 했지만, 렌즈를 빠져나가는 반대쪽 면에서는 검사 없이 그냥 거리 비교만 했다.

```python
# 처음 포팅한 코드 -- 출구면 쪽은 inside() 검사가 아예 빠져 있었다
ddc2 = np.linalg.norm([x2_c - x1, y2_c - y1])
dds2 = np.linalg.norm([x2_s - x1, y2_s - y1])
if ddc2 < dds2:
    w2, x2, y2 = w2_c, x2_c, y2_c
else:
    w2, x2, y2 = w2_s, x2_s, y2_s
```

이 안전장치가 렌즈의 절반에서만 작동하는 셈이 된다. 지금 쓰는 예시 렌즈·광선 조합에서는 이 차이가 눈에 띄는 정도로 나타나지는 않았지만(그래서 실행 중에는 못 잡았다), 원본과 동작이 다르다는 것 자체가 문제라 고쳤다. 지금 버전(위 `_pick_candidate` 전체)은 두 면 모두에서 같은 `inside()` 검사를 거치도록 통일한 것이다.

## 3. 검증 — 가역성과 이론값 대조

### 광학 가역성

Snell's law는 시간 역전에 대해 대칭이다. 즉 굴절된 광선의 진행 방향을 반대로 뒤집어 쏘면, 원래 들어왔던 경로를 정확히 되짚어 나가야 한다. 이건 `refract_through_lens`의 두 분기(`entry_surface==1`일 때와 그 반대일 때)가 서로 모순 없이 같은 물리를 구현했는지 보는 좋은 교차검증이 된다 — 정방향 진입은 한쪽 분기를, 역추적은 반대쪽 분기를 타기 때문이다.

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

실행 결과:

```text
정방향 진입면=2, 역추적 진입면=1 (반대쪽 면이어야 함)
원래 입사 광선 시작점: (1000.000, 100.000), 방향 (-1.000, 0.000)
역추적 도착점(렌즈 첫 굴절점): (44.987, 100.000), 방향 (1.000, -0.000)
위치 오차: 0.00e+00 mm
```

진입면이 정방향과 정확히 반대(2 → 1)로 나오고, 역추적한 광선이 원래 입사 광선의 시작 $y$좌표(100mm)와 방향(부호만 반대)을 부동소수점 오차 없이 그대로 복원한다.

![가역성 검증](/assets/img/posts/raytracing-spherical-lens-refraction/fig2-reversibility.png){: width="600" }
_그림2. 출사 광선을 반대로 쏘면 입사 경로를 정확히 되짚는다 (파란 실선과 빨간 점선이 렌즈~x=1000 구간에서 완전히 겹친다)_

### thick-lens 공식과 대조 — 근축은 거의 정확히, 가장자리는 어긋난다

렌즈의 두께까지 고려한 thick-lens 공식으로 초점거리와 후방초점거리(back focal distance)를 구할 수 있다.

$$ \frac{1}{f} = (n-1)\left[\frac{1}{R_1} - \frac{1}{R_2} + \frac{(n-1)t}{n R_1 R_2}\right], \qquad \mathrm{BFD} = f\left[1 - \frac{(n-1)t}{n R_1}\right] $$

이 공식은 근축광선(paraxial ray, 광축에 아주 가까운 광선) 근사 위에서 유도된 것이다. 광선추적 시뮬레이션은 근사 없이 실제 기하로 계산하므로, 광축에서 조금이라도 떨어진(marginal) 광선은 이 공식과 정확히 일치할 이유가 없다 — 오히려 정확히 일치한다면 그게 이상한 것이다(구면수차, spherical aberration).

```python
n = n_lens  # N-BK7, 750nm
inv_f = (n - 1) * (1 / R1 - 1 / R2 + (n - 1) * T / (n * R1 * R2))
f = 1 / inv_f
bfd = f * (1 - (n - 1) * T / (n * R1))
expected_focus_x = -T / 2 - bfd
```

```text
n(N-BK7, 750nm) = 1.511835
thick-lens 이론 초점(근축): f=993.698mm, 후방초점거리=960.056mm, 초점 x=-1010.056
시뮬레이션 근축 초점(|y0|<=20mm) x=-1010.229, 이론 대비 -0.172mm (-0.017%)
시뮬레이션 최외곽 초점(|y0|=100mm) x=-995.010, 근축 대비 +15.219mm (구면수차)
```

광축에서 20mm 이내(렌즈 반지름 250mm의 8%)로 들어오는 근축 광선은 이론값과 0.017% 차이로 사실상 정확히 일치한다. 반면 조리개 가장자리(광축에서 100mm)로 들어오는 광선은 근축 초점보다 15.2mm 렌즈 쪽으로 당겨져 초점을 맺는다 — 단순 양볼록 렌즈에서 흔히 보는 undercorrected 구면수차 방향(가장자리 광선일수록 초점이 짧아짐)과 일치한다.

![초점 위치 비교](/assets/img/posts/raytracing-spherical-lens-refraction/fig3-focus-vs-theory.png){: width="600" }
_그림3. 근축 이론값(검은 점선)과 시뮬레이션의 근축 초점(파란 점선)은 거의 겹치고, 최외곽 광선 초점(빨간 점선)만 구면수차만큼 어긋난다_

## 4. 한계

이번에 검증한 `refract_through_lens`는 전반사(total internal reflection) 여부를 판정하지 않는다. Snell's law에서 $\sin\theta_2 = (n_1/n_2)\sin\theta_1$이 1을 넘어서면 `arcsin`이 정의역을 벗어나는데, 이 함수는 그 경우를 걸러내지 않고 그냥 NumPy가 계산한 `NaN`을 다음 계산으로 흘려보낸다. MATLAB 원본도 이 구면렌즈 함수에서는 마찬가지였다 — 전반사를 명시적으로 판정하는 건 프리즘 계산 함수 쪽뿐이었다. 렌즈는 대개 공기보다 굴절률이 높은 매질로 들어가는 각이 완만해서 전반사 조건에 잘 안 걸리지만, 곡률이 크고 조리개가 넓은 조합에서는 걸릴 수 있다. 다음 편에서 프리즘의 전반사 판정 로직과, 그걸 프리즘이 아닌 렌즈에도 적용하면 어떻게 되는지를 다룬다.
