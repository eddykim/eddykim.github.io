---
title: "광선추적 2편 — 렌즈 가장자리로 들어간 빛은 왜 밖으로 나오지 못하는가"
date: 2026-09-08 20:00:00 +0900
categories: [광학, 기하광학]
page_id: raytracing-total-internal-reflection
tags: [ray-tracing, geometric-optics, snells-law, python, matlab]
description: 프리즘의 전반사 판정을 검증하고, 같은 판정을 구면 렌즈의 출사면에 넣었다. 전반사한 빛이 어디로 가는지와 반사 공식이 몇 번 반복되는지도 정리한다.
math: true
---

[1편](/posts/raytracing-spherical-lens-refraction/)에서는 렌즈 굴절 함수를 구현한 내용을 기록했고, 구현의 한계로 이렇게 적었다. "이번에 검증한 `refract_through_lens`는 전반사(total internal reflection) 여부를 판정하지 않는다." 스넬의 법칙(Snell's law)에서 $\sin\theta_2 = (n_1/n_2)\sin\theta_1$이 1을 넘으면 `arcsin`이 정의역을 벗어난다. 렌즈 함수는 그 경우를 걸러내지 않고 `NaN`을 그대로 다음 계산에 넘긴다. 반면 프리즘 함수(`CalculateRayPath_Prism`)는 애초에 이 경우를 처리하도록 짜여 있었다. 이번 편에서는 프리즘 함수의 전반사(total internal reflection, TIR) 판정을 먼저 검증하고, 같은 판정을 렌즈 함수의 출사면에도 넣는다.

## 1. 프리즘 함수는 전반사를 어떻게 판정하는가

프리즘 함수는 `arcsin`이 정의역을 벗어나는 순간을 어떻게 알아채는가? 굴절각을 구하기 전에 사인 값이 1을 넘는지 먼저 확인하고, 넘으면 굴절 대신 그 면에서 반사시킨다. 다각형 경계의 각 선분 위에서 반사법칙을 적용하는 셈이다.

렌즈 함수와 프리즘 함수는 광선과 면의 교차점을 찾는 방식부터 다르다. 렌즈는 면이 원호라서 `arccos`/`arcsin`으로 교차각을 직접 풀었다([1편](/posts/raytracing-spherical-lens-refraction/) 참고). 프리즘의 경계는 꼭짓점을 잇는 선분들의 모음이라 직선과 직선의 교차와 벡터 연산만으로 충분하다. 또 렌즈 함수는 입사면 1번, 출사면 1번으로 끝나지만, 프리즘 함수는 광선이 더 이상 어떤 경계와도 만나지 않을 때까지 반복문을 돈다. 그래서 그 안에서 전반사가 몇 번이든 일어날 수 있다. 아래 발췌는 반복 횟수에 상한(`max_bounces=20`)을 둔 `for` 반복문이다.

```python
def trace_ray_through_prism(ray, prism, n_lens, n_air, max_bounces=20, verbose=False):
    """프리즘 경계와 더 이상 만나지 않을 때까지 광선을 추적한다."""
    xs, ys, vx, vy = map(float, ray)
    boundary = prism["BOUNDARY"]
    air2glass = True
    path = [[xs, ys, vx, vy]]

    for _ in range(max_bounces):
        # ... (가장 가까운 경계 교차점 (xc, yc)를 찾는 부분, 생략) ...

        # 면의 법선 -- 선분에 수직인 방향, 광선 진행 방향과 같은 쪽(v·n ≥ 0)을 향하도록 부호를 맞춘다
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
            # ... 스넬의 법칙으로 굴절 방향 계산 ...
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

여기서 $\vec v$는 입사 방향, $\vec v'$은 반사 방향, $\hat n$은 면의 단위 법선이다. 즉, 법선 방향 성분만 부호를 뒤집고 면에 나란한 성분은 그대로 둔다. $\hat n$ 대신 $-\hat n$을 넣어도 두 번 곱해지므로 결과가 같다. 반사에는 법선이 어느 쪽을 향하는지가 상관없다는 뜻이다([3편](/posts/raytracing-reflection-family/)에서 다시 다룬다).

이 코드로 한 변이 600mm인 정삼각형 N-BK7 프리즘에 수직 입사광을 넣었다.

<img src="/assets/img/posts/raytracing-total-internal-reflection/fig1-prism-tir-path.png" alt="정삼각형 프리즘 안의 전반사 경로" width="600">
_그림1. 한 변 600mm 정삼각형 N-BK7 프리즘(파장 750nm)의 한 면에 수직으로 들어온 광선. 두 번째 면에 내부입사각 60°로 닿아 전반사하고, 세 번째 면에는 수직으로 닿아 꺾이지 않고 나간다._

수직으로 들어온 광선이니 프리즘을 그대로 통과할 것 같지만, 실행해 보면 광선은 두 번째 면에서 먼저 전반사하고 세 번째 면에 수직으로 닿아 꺾이지 않고 나간다. 여기서 내부입사각은 유리 안쪽에서 면의 법선과 잰 입사각이다(아래 출력의 `theta_in`).

```text
face=(0.0,-300.0)-(0.0,300.0)  theta_in=0.000deg  투과  theta_out=0.000deg
face=(0.0,300.0)-(519.6,0.0)  theta_in=60.000deg  전반사(TIR)
face=(519.6,0.0)-(0.0,-300.0)  theta_in=0.000deg  투과  theta_out=0.000deg
```

이유는 간단하다. N-BK7의 임계각(critical angle)은 $\theta_c = \arcsin(1/n) = 41.41°$인데(750nm에서 $n=1.511835$), 정삼각형의 꼭지각(apex angle)은 60°다. 입사면에 수직으로 들어온 광선은 굴절 없이 그대로 직진하므로, 다음 면에는 정확히 $60° - 0° = 60°$의 내부입사각으로 부딪힌다. $60° > 41.41°$이니 전반사할 수밖에 없다.

임계각이 무엇인지는 그림으로 보면 분명해진다. 굴절률이 큰 쪽에서 작은 쪽으로 빛이 나갈 때 굴절각은 입사각보다 크다. 입사각을 키우면 굴절광은 점점 면 쪽으로 눕고, 굴절각이 90°에 이르면 면을 스치듯 나간다(그림2). 즉, 임계각은 굴절광이 면을 스치는 입사각이다. 이보다 비스듬히 들어오면 나갈 굴절광이 없으니 빛이 전부 되돌아온다. 코드에서 `abs(sin_out) > 1`이 바로 이 상황이다.

<img src="/assets/img/posts/raytracing-total-internal-reflection/ext-critical-angle-diagram.png" alt="입사각이 커질수록 굴절광이 면을 스치다가 전반사로 바뀌는 개념도" width="700">
_그림2. 굴절률이 큰 매질(아래, n1)에서 작은 매질(위, n2)로 나가는 광선의 입사각을 키워 가는 모습. 입사각 θ1에서는 굴절각 θ2로 꺾여 나가고, 임계각 θc에서는 굴절광이 면을 따라 90°로 스치며, 그보다 크면 전부 반사된다(total reflection). 옅은 선은 면에서 일부만 반사되는 빛이다. (출처: Jfmelero (adapted by Gavin R Putland), [ReflexionTotal en](https://commons.wikimedia.org/wiki/File:ReflexionTotal_en.svg), [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0))_

임계각은 굴절률이 클수록 작다. 589nm에서 물(n=1.333)은 48.6°, 용융 석영(1.4584)은 43.3°, N-BK7(d선 1.5168)은 41.2°, 다이아몬드(2.417)는 24.4°다. 굴절률은 Hecht의 표 4.1과 SCHOTT 데이터시트의 값이다.

그렇다면 이 프리즘을 전반사 없이 투과시키려면 어떻게 해야 하는가? 입사각을 충분히 크게 주면 된다. 입사각 $\theta_1$, 첫 면의 굴절각 $\theta_1'$, 두 번째 면의 내부입사각 $\theta_2'$, 꼭지각 $A$ 사이에는 $\theta_1' + \theta_2' = A$ 관계가 성립한다. 따라서 두 번째 면에서 임계각을 넘지 않으려면 다음 조건을 만족해야 한다.

$$ \theta_1' > A - \theta_c, \qquad \theta_1 > \arcsin\big(n \sin(A - \theta_c)\big) $$

즉, 첫 면에서 충분히 꺾여 들어가야 두 번째 면에서 각이 임계각 아래로 내려온다. 이 식은 광선이 첫 면의 법선보다 꼭지각에서 먼 쪽(그림3처럼 아래쪽)에서 비스듬히 들어올 때의 이야기다. 꼭지각 쪽에서 들어오면 $\theta_2' = A + \theta_1'$이 되어 입사각과 상관없이 전반사한다. 값을 대입하면 $\theta_1 > 28.8131°$가 나온다. 입사각을 0°부터 이분법(bisection)으로 훑어(sweep) 실제 코드가 관측하는 전반사↔투과 경계를 찾으면(`verify_tir.py`) 28.813100°로, 이론값과 $7.6\times10^{-10}$° 이내로 일치한다.

<img src="/assets/img/posts/raytracing-total-internal-reflection/fig2-prism-angle-sweep.png" alt="입사각에 따른 전반사와 투과" width="600">
_그림3. 같은 프리즘에 입사각 0°~45°로 넣은 광선 여섯 개. 경계값(28.8131°)의 바로 아래인 28.813°까지는 두 번째 면에서 전반사하고, 35°와 45°는 두 번째 면을 투과한다._

그림3에서 투과한 광선의 출사각을 보면 경계의 뜻이 드러난다. 두 번째 면에서 법선과 잰 출사각을 계산하면 입사각 45°에서 53.5°, 35°에서 67.6°, 30°에서 80.3°이고, 경계 바로 위인 28.82°에서는 89.3°다. 경계에 다가갈수록 출사광이 면을 스치듯 눕는다. 출사각이 90°에 이르는 지점이 곧 전반사의 경계이고, 그림2에서 본 임계각의 정의가 그대로 나타난다.

이는 낯선 현상이 아니다. 쌍안경에 들어가는 Porro 프리즘(Porro prism)은 거울 코팅 없이 프리즘 내부의 전반사만으로 빛의 경로를 꺾는다. 직각 프리즘의 빗변에 빛을 수직으로 넣으면, 빗변과 45°를 이루는 두 짧은 면에 각각 45°의 내부입사각으로 부딪힌다. 45° > 41.4°(앞에서 구한 750nm의 N-BK7 임계각. d선 587.6nm에서는 41.2°)라 반드시 전반사하므로 코팅 없이도 반사면 역할을 한다. 빛은 두 번 반사한 뒤 들어온 방향의 반대로 되돌아 나간다(그림4).

<img src="/assets/img/posts/raytracing-total-internal-reflection/ext-porro-prism-photo.jpg" alt="직각 프리즘 안에서 두 번 전반사해 되돌아 나오는 빛줄기" width="600">
_그림4. 직각 프리즘의 빗변에 수직으로 들어간 빛줄기. 빛은 두 짧은 면에서 차례로 전반사해 들어온 방향의 반대로 나가며, 짧은 면 밖으로 새어 나오는 빛은 보이지 않는다. (출처: Zátonyi Sándor, (ifj) Fizped, [Visszatükrözés prizmán 2](https://commons.wikimedia.org/wiki/File:Visszat%C3%BCkr%C3%B6z%C3%A9s_prizm%C3%A1n_2.jpg), [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0))_

45° 면이 전반사하려면 $\sin 45° > 1/n$, 곧 $n > \sqrt{2} \approx 1.414$면 된다. 앞에서 든 용융 석영(1.4584)도 N-BK7도 이 조건을 넘지만, 물(1.333)은 넘지 못한다. 같은 이유로 반사면에 물이 묻으면 전반사가 깨진다. 유리와 물의 경계에서는 임계각이 $\arcsin(1.333/1.5168) = 61.5°$로 커지므로, 45°로 닿은 빛은 물이 묻은 자리에서 물 쪽으로 빠져나간다. 전반사면은 코팅이 필요 없는 대신 깨끗하게 유지해야 한다.

## 2. 렌즈 출사면에서는 어느 광선부터 전반사하는가

같은 판정이 렌즈에 없으면 어느 광선에서 무엇이 깨지는가? 렌즈 함수의 출사면 굴절 계산은 다음과 같다([1편](/posts/raytracing-spherical-lens-refraction/)에서 그대로 가져온 부분).

```python
phi2 = np.arctan2(vy1, vx1)
theta_t2 = np.arcsin(n_lens * np.sin(w2 - phi2) / n_air)  # 정의역 확인이 없다
vx2, vy2 = -np.cos(theta_t2 + w2), -np.sin(theta_t2 + w2)
```

`n_lens * sin(w2 - phi2) / n_air`가 1을 넘으면, 즉 내부입사각이 임계각을 넘으면 이 `arcsin`은 애초에 실수해가 없다. 원본 코드는 이 경우를 걸러내지 않는다. 그래서 임계각을 넘는 광선에서는 예외 없이 경고(`RuntimeWarning`) 한 줄만 남기고 `NaN`이 나온다.

이 함수는 광선이 어느 면으로 먼저 들어오는지 스스로 판단하지 않으므로, 호출하는 쪽이 `find_entry_surface`로 판정한 진입면(`entry_surface`)을 넘겨야 한다([1편](/posts/raytracing-spherical-lens-refraction/) 참고). 잘못 넘기면 전반사와 무관하게 틀린 경로가 나오거나 `NaN`이 날 수 있다.

평볼록 렌즈(plano-convex lens)로 확인했다. R1은 $10^6$mm로 사실상 평면이고, R2=-50mm, T=20mm, D=78mm다. 평면 쪽으로 높이 y0=36mm인 평행광선을 넣은 결과는 다음과 같다.

```text
find_entry_surface: 1번 면

수정 전:
[[-9.99935200e+00  3.60000000e+01  1.00000000e+00 -1.21878769e-05]
 [-5.30123745e+00  3.59999427e+01             nan             nan]]
-> RuntimeWarning: invalid value encountered in arcsin
```

출력의 마지막 줄이 앞에서 말한 경고다. 프리즘과 같은 방식으로 고쳤다. 정의역을 벗어나면 굴절 대신 출사면의 법선 $(\cos w_2, \sin w_2)$를 기준으로 반사시킨다. 구면 위의 한 점에서 법선은 반지름 방향이므로, 원의 매개변수 각 $w_2$를 그대로 쓰면 된다. 이 렌즈처럼 R2가 음수이면 이 벡터는 그 점에서 중심 쪽을 향하지만, 앞에서 본 대로 반사 공식은 법선의 부호와 무관하다.

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

이 결과가 진짜 반사인지 두 가지로 검증했다(`verify_tir.py`). 먼저 법선을 기준으로 반사법칙(입사각=반사각)이 성립하는지 보았다. 입사각과 반사각이 모두 46.05508420°이고 차이는 $7.1\times10^{-15}$°다. 다음으로, 전반사가 시작되는 광선 높이를 이분법으로 찾은 뒤 그 높이의 내부입사각을 이론 임계각과 비교했다. 코드가 찾은 경계는 y0=33.072057mm이고, 그 지점의 내부입사각은 41.410389°로 이론값 $\theta_c=41.410389°$와 $2.7\times10^{-10}$° 이내로 일치한다.

이 경계는 식으로도 바로 나온다. 평면 쪽으로 들어온 평행광은 평면에서 꺾이지 않으므로, 높이 $h$인 광선이 구면에 닿을 때 법선과 이루는 각은 $\sin\theta = h/\lvert R_2\rvert$를 만족한다. 이 값이 $1/n$을 넘는 순간, 곧 $h = \lvert R_2\rvert/n = 33.07$mm부터 전반사가 시작된다.

<img src="/assets/img/posts/raytracing-total-internal-reflection/fig3-lens-marginal-tir.png" alt="평볼록 렌즈 최외곽 광선의 전반사와 그 뒤 경로" width="600">
_그림5. 평볼록 렌즈(R2=-50mm, T=20mm, N-BK7, 750nm)의 평면 쪽에 2mm 간격으로 넣은 평행광선 39개. 높이 33.07mm 바깥의 6개(빨간선)가 출사면에서 전반사한다. 렌즈 함수는 반사 뒤 경로를 계산하지 않으므로, 빨간선의 반사 뒤 경로는 그림용으로 따로 짠 반복 추적(옆면은 매끄러운 면으로 가정)으로 그렸다._

전체 조리개(D=78mm)의 절반인 39mm 가운데 전반사가 일어나는 구간은 높이 33mm부터 가장자리까지다. 반지름으로는 바깥쪽 15%, 원형 조리개의 면적으로는 28%다. 그림5의 빨간선은 전반사한 빛이 어디로 가는지도 보여 준다. 높이 36mm와 38mm의 네 광선은 반대쪽 곡면에 처음과 같은 각(46.06°, 49.46°)으로 다시 닿아 한 번 더 전반사한 뒤, 평면을 지나 광원 쪽으로 되돌아간다. 원 안에서 현의 양 끝이 반지름과 이루는 각은 같기 때문이다. 높이 34mm의 두 광선은 한 번 전반사한 뒤 렌즈 옆면으로 빠져나간다. 어느 쪽도 초점으로 가지 않는다. 즉, 전반사한 빛은 상을 만드는 데 쓰이지 못한다.

그렇다면 어떤 렌즈에서 이 현상을 걱정해야 하는가? 경계 높이 $\lvert R_2\rvert/n$이 반구경 D/2보다 작으면 전반사 구간이 조리개 안에 들어온다. 평볼록 렌즈의 초점거리는 $f = \lvert R_2\rvert/(n-1)$이므로, 이 조건을 F수로 바꾸면 $f/D < n/\big(2(n-1)\big)$이다. N-BK7(750nm)이면 이 값은 약 1.48이다. 이번 렌즈는 f=97.7mm, F/1.25라 이 조건에 든다. 즉, 평면 쪽으로 평행광을 받는 N-BK7 평볼록 렌즈라면 F수가 약 1.48보다 작을 때 가장자리에서 전반사가 생긴다고 보아도 무난하다. 곡률이 완만한 렌즈에서 이 현상을 보기 어려운 이유다. 볼록면을 광원 쪽으로 돌린 경우는 여기서 계산하지 않았다.

## 3. 고친 렌즈 함수에는 무엇이 남았는가

프리즘 함수와 달리, 이번에 고친 렌즈 함수는 여전히 진입 1회와 출사 1회만 다룬다. 출사면에서 전반사한 광선의 방향은 이제 정확히 계산되지만, 그 광선이 렌즈 안에서 다음에 어디로 가는지는 이 함수가 계산하지 않는다. 프리즘처럼 경계와 더 이상 만나지 않을 때까지 도는 반복문이 없기 때문이다. 그림5의 반사 뒤 경로를 그리려고 그림 코드에 그런 반복 추적을 따로 넣었고, 그 첫 교차점과 방향은 렌즈 함수의 출사 결과와 좌표·방향 성분 모두 $2.5\times10^{-14}$ 이내로 같았다. 반사된 광선을 라이브러리 안에서 계속 추적하려면 렌즈 함수 자체를 프리즘처럼 다시 짜야 한다.

## 정리

전반사 판정은 `arcsin`에 넣을 값이 1을 넘는지 하나로 끝난다. 넘으면 굴절 대신 $\vec v' = \vec v - 2(\vec v\cdot\hat n)\hat n$으로 반사시키면 된다. 정삼각형 N-BK7 프리즘(750nm)은 입사각 28.8° 아래에서 두 번째 면에 반드시 전반사한다. 평면 쪽으로 평행광을 받는 평볼록 렌즈는 높이 $\lvert R_2\rvert/n$ 바깥에서 출사면 전반사가 생기고, 이번 렌즈에서는 그 높이가 33.07mm였다. 전반사한 빛은 초점으로 가지 않고 광원 쪽이나 렌즈 옆면으로 빠져나간다.

이번 편에서 반사에 두 번(프리즘의 내부 반사, 렌즈의 전반사) 쓴 공식은 사실 이 라이브러리 전체에서 정확히 같은 형태로 여섯 번 나온다. `FlatMirror`, `SphericalMirror`, `ArbitraryMirror`, `BeamSplitter`(반사 성분), 그리고 이번에 고친 `Prism`과 `SphericalLens`의 전반사 분기다. 함수마다 다른 것은 충돌점과 법선을 어떻게 찾느냐(평면 대수, 원의 후보해, 다각형 위 국소 곡률)와, 한 번만 반사하고 끝나는지 여러 번 반사할 때까지 도는지뿐이다. 굴절은 매번 다른 각도 공식이 필요했지만, 반사는 거울이든 프리즘이든 렌즈든 이 벡터 공식 하나로 끝난다. [3편](/posts/raytracing-reflection-family/)에서는 이 반사 공식이 코드 전체에서 정말로 일관되게 구현돼 있는지, `FlatMirror`/`SphericalMirror`/`ArbitraryMirror`/`BeamSplitter`를 한 번에 검증한다.

## 참고자료

- E. Hecht, *Optics*, 5th ed. (Pearson, 2017), 4.7절(전반사와 임계각), 표 4.1(589nm 굴절률), 5.5절(프리즘의 꼭지각 관계, Porro 프리즘).
- SCHOTT, N-BK7 광학 유리 데이터시트(d선 굴절률 $n_d = 1.5168$). [SCHOTT 광학 유리](https://www.schott.com/en-us/products/optical-glass-p1000267)
