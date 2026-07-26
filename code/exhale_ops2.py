"""Addendum: exact min weight of the TARGET operator's own coset
(operator + row-space of HZ^T, i.e. products of Z stabilizers) in the
corridor-deflated code. This is the literal 'weight of the logical
operation' -- physical gates in an applied logical, patch-boundary size in
surgery, qubits read in a transversal measurement.

Also: sanity-locate the weight-3 representative inside the corridor.
"""
import numpy as np
from cellulation import gf2_rank, gf2_nullspace, code_k, torus, bring_code
from v1_partial import refine_partial
from scipy.optimize import milp, LinearConstraint, Bounds


def coset_min_weight(cell, op):
    """min |z| s.t. z = op XOR HZ^T y over GF(2), via ILP.
    Row i: z_i + sum_j HZt[i,j] y_j - 2 s_i = op_i."""
    HZt = cell.HZ().T.astype(float)          # nE x nF
    n, m = HZt.shape
    nvar = n + m + n                          # z, y, slack s
    c = np.zeros(nvar); c[:n] = 1
    A, lb, ub = [], [], []
    for i in range(n):
        r = np.zeros(nvar)
        r[i] = 1
        r[n:n + m] = HZt[i]
        r[n + m + i] = -2
        A.append(r); lb.append(float(op[i])); ub.append(float(op[i]))
    ymax = float(HZt.sum(0).max())
    smax = (1 + HZt.sum(1).max()) / 2 * np.ones(n)
    res = milp(c=c, constraints=LinearConstraint(np.array(A), lb, ub),
               bounds=Bounds(np.zeros(nvar),
                             np.concatenate([np.ones(n), np.full(m, 1.0),
                                             smax])),
               integrality=np.ones(nvar))
    assert res.success, res.message
    z = np.round(res.x[:n]).astype(np.uint8)
    return int(round(res.fun)), z


# ---------------- Bring corridor ----------------
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
corridor = [fi for fi, f in enumerate(b.faces) if set(support) & set(f)]
allb = list(range(b.nF))

# fully fine reference
fine_b, finfo = refine_partial(b, allb)
lift_f = np.zeros(fine_b.nE, dtype=np.uint8)
for e in support:
    for h in finfo["halves"][e]:
        lift_f[h] ^= 1
w_fine, _ = coset_min_weight(fine_b, lift_f)

# corridor-deflated
refined = [fi for fi in allb if fi not in corridor]
mixed, minfo = refine_partial(b, refined)
lift_m = np.zeros(mixed.nE, dtype=np.uint8)
for e in support:
    for h in minfo["halves"][e]:
        lift_m[h] ^= 1
assert not ((mixed.HX().astype(int) @ lift_m.astype(int)) % 2).any()
w_mixed, zmin = coset_min_weight(mixed, lift_m)

print("Bring target Z-bar operator (own coset min weight):")
print(f"  fully fine      : {w_fine}")
print(f"  corridor exhale : {w_mixed}   (lift itself has weight "
      f"{int(lift_m.sum())}; support edges unsubdivided: "
      f"{all(len(minfo['halves'][e]) == 1 for e in support)})")
sup_m = np.flatnonzero(zmin).tolist()
# is the min rep inside the corridor? (each of its edges only in corridor faces)
inside = all(all(fi in corridor for fi, f in enumerate(mixed.faces)
                 if ei in f) for ei in sup_m)
print(f"  min rep support size {len(sup_m)}, entirely inside corridor: {inside}")

# ---------------- torus 1-column annulus ----------------
L = 3
base = torus(L)
allt = list(range(base.nF))
fine_t, tinfo = refine_partial(base, allt)
# target: Z-logical with a representative inside face-column 0
# (extract Z-basis of the coarse torus, take one whose support is vertical col)
kerXt = gf2_nullspace(base.HX())
stack = base.HZ().copy()
Lzs = []
for v in kerXt:
    t = np.vstack([stack, v[None, :]])
    if gf2_rank(t) > gf2_rank(stack):
        stack = t
        Lzs.append(v)
refined_t = [i * L + j for j in range(1, L) for i in range(L)]  # deflate col 0
mixed_t, mtinfo = refine_partial(base, refined_t)
print("\nTorus, deflate face-column 0; per coarse Z-basis element:")
for idx, Lz in enumerate(Lzs):
    sup = np.flatnonzero(Lz).tolist()
    lf = np.zeros(fine_t.nE, dtype=np.uint8)
    lm = np.zeros(mixed_t.nE, dtype=np.uint8)
    for e in sup:
        for h in tinfo["halves"][e]:
            lf[h] ^= 1
        for h in mtinfo["halves"][e]:
            lm[h] ^= 1
    wf, _ = coset_min_weight(fine_t, lf)
    wm, _ = coset_min_weight(mixed_t, lm)
    print(f"  Z-logical {idx} (coarse support {sup}): fine {wf} -> corridor {wm}")
