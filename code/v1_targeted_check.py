"""v1-B addendum: exact distance of the SPECIFIC targeted logical after
partial inflation of the Bring code (lifts the coarse operator into the
mixed cellulation, then ILP-minimizes over its coset).

Expected output: targeted logical d = 3 -> 5.
"""
import numpy as np
from cellulation import bring_code, gf2_rank, gf2_nullspace
from v1_partial import refine_partial
from scipy.optimize import milp, LinearConstraint, Bounds

b = bring_code()
kerX = gf2_nullspace(b.HX())
stack = b.HZ().copy()
Lz0 = None
for v in kerX:
    t = np.vstack([stack, v[None, :]])
    if gf2_rank(t) > gf2_rank(stack):
        Lz0 = v
        break
support = np.flatnonzero(Lz0).tolist()
touching = [fi for fi, f in enumerate(b.faces) if set(support) & set(f)]
mixed, info = refine_partial(b, touching)

# lift Lz0: each SUPPORT edge -> its half-edges in the mixed cell
lift = np.zeros(mixed.nE, dtype=np.uint8)
for e in support:
    for h in info["halves"][e]:
        lift[h] ^= 1
assert not ((mixed.HX().astype(int) @ lift.astype(int)) % 2).any(), \
    "lifted operator must commute with all X checks"

HZ = mixed.HZ().astype(float)
n, m = mixed.nE, HZ.shape[0]
nvar = n + m + 1
c = np.zeros(nvar)
c[:n] = 1
A, lb, ub = [], [], []
for i in range(m):
    r = np.zeros(nvar)
    r[:n] = HZ[i]
    r[n + i] = -2
    A.append(r); lb.append(0); ub.append(0)
r = np.zeros(nvar)
r[:n] = lift
r[n + m] = -2
A.append(r); lb.append(1); ub.append(1)
res = milp(c=c, constraints=LinearConstraint(np.array(A), lb, ub),
           bounds=Bounds(np.zeros(nvar),
                         np.concatenate([np.ones(n),
                                         np.full(m, HZ.sum(1).max()), [n]])),
           integrality=np.ones(nvar))
print(f"targeted logical: d = 3 -> {int(round(res.fun))} "
      f"after inflating its {len(touching)} touching faces "
      f"(+{info['new_data'] + info['new_anc']} atoms)")
