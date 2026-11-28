"""Data for the kt3d/picksup_probe fixture (docs/spec/search.md step 3).

40 points scattered over the 42 x 42 node grid (spacing 5, origin 2.5) by an
additive recurrence, with a smooth value. Deterministic; run from the
repository root to regenerate: python3 -I docs/spec/probes/picksup_data.py
"""
import math
lines = ["picksup probe: 40 scattered points for kt3d", "4", "X", "Y", "Z", "V"]
for i in range(40):
    x = 5.0 + ((i * 0.6180339887) % 1.0) * 200.0
    y = 5.0 + ((i * 0.4142135624) % 1.0) * 200.0
    v = 2.0 + math.sin(x / 40.0) + 0.5 * math.cos(y / 30.0)
    lines.append(f"{x:.4f} {y:.4f} 0.0000 {v:.6f}")
with open("tests/fixtures/gslib/kt3d/picksup_probe/probe.dat", "w") as f:
    f.write("\n".join(lines) + "\n")
