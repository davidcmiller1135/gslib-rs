"""Find a kt3d configuration where GSLIB's picksup skips a super-block that
holds a point inside the search ellipsoid of an actual grid node.

Geometry follows vendor/gslib-2.0/src: setsupr (nxsup = min(nx, MAXSBX = 21),
xsizsup = nx*xsiz/nxsup, xmnsup = xmn - xsiz/2 + xsizsup/2), picksup (a block
offset is searched iff the smallest sqdist over the 3x3x3 lattice
(i1-i2)/2*size + offset, i1,i2 in {-1,1}, is <= radsqd) and setrot (2-D:
major axis along the azimuth, minor scaled by radius/radius1).
Usage: python3 -I picksup_kt3d.py
"""
import math, itertools

def rows(az, r, r1, r2):
    a = math.radians(90.0 - az)
    c, s = math.cos(a), math.sin(a)
    return [(c, s, 0.0), (-s * r / r1, c * r / r1, 0.0), (0.0, 0.0, r / r2)]

def sq(m, d):
    return sum((m[i][0] * d[0] + m[i][1] * d[1] + m[i][2] * d[2]) ** 2 for i in range(3))

def lattice_min(m, off, size):
    best = math.inf
    for u in itertools.product((-1, 0, 1), repeat=3):
        d = [u[k] * size[k] + off[k] for k in range(3)]
        best = min(best, sq(m, d))
    return best

def box_min(m, lo, hi, n=200):
    """Smallest q over the box [lo, hi] (2-D, z = 0), and where."""
    best, where = math.inf, None
    for i in range(n + 1):
        for j in range(n + 1):
            p = (lo[0] + (hi[0] - lo[0]) * i / n, lo[1] + (hi[1] - lo[1]) * j / n, 0.0)
            q = sq(m, p)
            if q < best:
                best, where = q, p
    return best, where

# Grid: 42 x 42 nodes of 5, so super-blocks of 10 hold 2 x 2 nodes, each node
# 2.5 off its block centre.
nx, xmn, xsiz = 42, 2.5, 5.0
size = (10.0, 10.0, 10.0)
found = []
for az in range(5, 86, 5):
    for ratio in (0.1, 0.15, 0.2, 0.25):
        for r in (15.0, 20.0, 25.0, 30.0):
            m = rows(az, r, r * ratio, r)
            for bi in range(-4, 5):
                for bj in range(-4, 5):
                    lm = lattice_min(m, (bi * 10.0, bj * 10.0, 0.0), size)
                    if lm <= r * r:
                        continue
                    # Node offsets from its block centre: +-2.5 on each axis.
                    for tx, ty in itertools.product((-2.5, 2.5), repeat=2):
                        lo = (bi * 10.0 - 5.0 - tx, bj * 10.0 - 5.0 - ty)
                        hi = (bi * 10.0 + 5.0 - tx, bj * 10.0 + 5.0 - ty)
                        q, where = box_min(m, lo, hi, 60)
                        if q < 0.97 * r * r:
                            found.append((round(math.sqrt(q) / r, 4), az, ratio, r, bi, bj, tx, ty, where))
found.sort()
print(len(found))
for f in found[:8]:
    print(f)
