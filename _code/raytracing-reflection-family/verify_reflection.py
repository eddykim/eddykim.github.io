"""네 반사 함수의 반사법칙(입사각=반사각) 검증.

실행: python verify_reflection.py
"""
import numpy as np

from reflection_optics import (
    make_arbitrary_mirror,
    make_beam_splitter,
    make_flat_mirror,
    make_spherical_mirror,
    reflect_arbitrary_mirror,
    reflect_beam_splitter,
    reflect_flat_mirror,
    reflect_spherical_mirror,
)


def ang_err(vin, vout, n):
    ang_in = np.rad2deg(np.arccos(abs(np.dot(vin, n))))
    ang_out = np.rad2deg(np.arccos(abs(np.dot(vout, n))))
    return ang_in, ang_out, ang_in - ang_out


# FlatMirror
mirror = make_flat_mirror(200, 0, 0, np.deg2rad(20))
ray = np.array([-100.0, 30.0, 1.0, 0.0])
out = reflect_flat_mirror(ray, mirror)
(x1, y1), (x2, y2) = mirror["BOUNDARY"][0]
n = np.array([-(y2 - y1), x2 - x1]); n /= np.linalg.norm(n)
ai, ao, err = ang_err(ray[2:4], out[0, 2:4], n)
print(f"FlatMirror        입사각={ai:.6f}°  반사각={ao:.6f}°  차이={err:.1e}°")

# SphericalMirror
mirror = make_spherical_mirror(-300, 200, 0, 0, 0)
ray = np.array([-100.0, 40.0, 1.0, 0.0])
out = reflect_spherical_mirror(ray, mirror)
cx, cy = mirror["cx"], mirror["cy"]
n = np.array([out[0, 0] - cx, out[0, 1] - cy]); n /= np.linalg.norm(n)
ai, ao, err = ang_err(ray[2:4], out[0, 2:4], n)
print(f"SphericalMirror    입사각={ai:.6f}°  반사각={ao:.6f}°  차이={err:.1e}°")

# BeamSplitter (reflected component)
mirror = make_beam_splitter(60, 0, 0, np.deg2rad(45))
ray = np.array([-100.0, 0.0, 1.0, 0.0])
r, t = reflect_beam_splitter(ray, mirror)
(x1, y1), (x2, y2) = mirror["BOUNDARY"][0]
n = np.array([-(y2 - y1), x2 - x1]); n /= np.linalg.norm(n)
ai, ao, err = ang_err(ray[2:4], r[0, 2:4], n)
print(f"BeamSplitter(반사)  입사각={ai:.6f}° 반사각={ao:.6f}°  차이={err:.1e}°")

# ArbitraryMirror (segment hit, not a vertex)
mirror = make_arbitrary_mirror([-100, 0, 100, 200], [0, 50, 0, 50])
ray = np.array([-50.0, 100.0, 0.0, -1.0])
out = reflect_arbitrary_mirror(ray, mirror)
(x1, y1), (x2, y2) = mirror["BOUNDARY"][0]
n = np.array([-(y2 - y1), x2 - x1]); n /= np.linalg.norm(n)
ai, ao, err = ang_err(ray[2:4], out[1, 2:4], n)
print(f"ArbitraryMirror(변) 입사각={ai:.6f}° 반사각={ao:.6f}°  차이={err:.1e}°")
