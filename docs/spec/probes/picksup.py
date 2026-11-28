# Does GSLIB's picksup (min over a 3x3x3 lattice of the block-difference box)
# miss super-blocks whose true minimum anisotropic distance is within the
# search radius? 2-D check: azimuth rotation with anisotropy, z ignored.
import math, itertools
def metric(az_deg, anis):
    # GSLIB setrot for 2-D: azimuth clockwise from north, major axis along az.
    a = math.radians(90.0 - az_deg)
    c, s = math.cos(a), math.sin(a)
    # rows: major (scaled 1), minor (scaled 1/anis)
    return [(c, s), (-s / anis, c / anis)]
def sq(m, dx, dy):
    return sum((r[0]*dx + r[1]*dy) ** 2 for r in m)
def lattice_min(m, xo, yo, xs, ys):
    return min(sq(m, xo + u*xs, yo + v*ys) for u in (-1,0,1) for v in (-1,0,1))
def true_min(m, xo, yo, xs, ys, n=400):
    best = math.inf
    for i in range(n+1):
        t = -1 + 2*i/n
        for (dx, dy) in ((xo + t*xs, yo - ys), (xo + t*xs, yo + ys), (xo - xs, yo + t*ys), (xo + xs, yo + t*ys)):
            best = min(best, sq(m, dx, dy))
    # interior minimum (origin inside the difference box)
    if abs(xo) <= xs and abs(yo) <= ys: best = 0.0
    return best
found = []
for az in (30, 45, 60, 20):
  for anis in (0.2, 0.3, 0.5):
    m = metric(az, anis)
    xs = ys = 10.0
    for radius in (25.0, 40.0, 55.0):
      for i in range(-12, 13):
        for j in range(-12, 13):
          lm = lattice_min(m, i*xs, j*ys, xs, ys)
          tm = true_min(m, i*xs, j*ys, xs, ys)
          if tm <= radius**2 < lm:
            found.append((az, anis, radius, i, j, round(math.sqrt(tm),3), round(math.sqrt(lm),3)))
print(len(found)); print(found[:6])
