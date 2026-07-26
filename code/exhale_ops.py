"""Exhale-to-operate: partial DEFLATION as a logical-operator-weight dial.

v1-B showed partial INFLATION raises a targeted logical's distance
(d_target = 3+c on the torus) iff the inflated region is a transversal cut.
This script runs the inverse question: starting from the FINE code, deflate
(exhale) a sub-region and measure how far the min weight of logical
operators drops -- because lattice-surgery / logical-measurement cost
scales with the weight w_log of the operated logical, this is a runtime
dial for the COST OF LOGICAL OPERATIONS.

Duality conjecture (min-vs-max):
  - inflation must block EVERY representative (transversal cut) to help;
  - deflation only needs to expose ONE representative to help.
So targeted deflation should be guaranteed-effective where targeted
inflation could be dodged. Bystanders may suffer "collateral exposure"
(the reverse of v1-B's collateral armor).

All numbers exact: GF(2) rank for k, ILP (HiGHS) for per-logical weights.
"""
import numpy as np
from cellulation import (Cellulation, gf2_rank, gf2_nullspace, code_k,
                         torus, bring_code)
from v1_partial import refine_partial, per_logical_distances
from scipy.optimize import milp, LinearConstraint, Bounds


def min_weight_anticommuting(cell, Lz):
    """ILP: min |x|, HZ x = 0 mod 2, <Lz, x> = 1 mod 2  (weight of the
    cheapest X-logical operation that acts on the class dual to Lz)."""
    HZ = cell.HZ().astype(float)
    n, m = cell.nE, HZ.shape[0]
    nvar = n + m + 1
    c = np.zeros(nvar); c[:n] = 1
    A, lb, ub = [], [], []
    for i in range(m):
        r = np.zeros(nvar); r[:n] = HZ[i]; r[n + i] = -2
        A.append(r); lb.append(0); ub.append(0)
    r = np.zeros(nvar); r[:n] = Lz; r[n + m] = -2
    A.append(r); lb.append(1); ub.append(1)
    res = milp(c=c, constraints=LinearConstraint(np.array(A), lb, ub),
               bounds=Bounds(np.zeros(nvar),
                             np.concatenate([np.ones(n),
                                             np.full(m, HZ.sum(1).max()),
                                             [n]])),
               integrality=np.ones(nvar))
    assert res.success, res.message
    return int(round(res.fun))


def atoms(cell):
    return cell.nE + cell.nF + cell.nV          # data + Z-anc + X-anc


def hops(cell):
    return int(cell.HZ().sum() + cell.HX().sum())


# ================================================================ torus
L = 3
base = torus(L)
allfaces = list(range(base.nF))
fine_t, _ = refine_partial(base, allfaces)
fine_t.validate()
FT_ATOMS, FT_HOPS = atoms(fine_t), hops(fine_t)

print("=" * 74)
print("PART 1  torus L=3: deflate c face-columns of the FINE [[72,2,6]] code")
print("        (refined set = remaining columns; c=3 -> fully coarse)")
print("=" * 74)
print(f"{'deflated':>9s} {'n':>4s} {'k':>2s} {'per-logical w':>14s} "
      f"{'atoms back':>10s} {'hops/rd':>8s}")
for c in range(0, L + 1):
    refined = [i * L + j for j in range(c, L) for i in range(L)]
    if refined:
        cellm, info = refine_partial(base, refined)
        cellm.validate()
    else:
        cellm = base
    ds = per_logical_distances(cellm)
    print(f"{c:>7d}c  {cellm.nE:4d} {code_k(cellm):2d} {str(sorted(ds)):>14s} "
          f"{FT_ATOMS - atoms(cellm):10d} {hops(cellm):8d}")

# one face-row corridor (annulus PARALLEL to the target cycle)
refined = [i * L + j for i in range(1, L) for j in range(L)]
cellr, _ = refine_partial(base, refined)
cellr.validate()
ds = per_logical_distances(cellr)
print(f"{'1 row':>8s}  {cellr.nE:4d} {code_k(cellr):2d} {str(sorted(ds)):>14s} "
      f"{FT_ATOMS - atoms(cellr):10d} {hops(cellr):8d}")

# ================================================================ Bring
print()
print("=" * 74)
print("PART 2  Bring {5,5}: corridor exhale of the FINE [[120,8,6]] code")
print("=" * 74)
b = bring_code()
allb = list(range(b.nF))
fine_b, _ = refine_partial(b, allb)
fine_b.validate()
FB_ATOMS, FB_HOPS = atoms(fine_b), hops(fine_b)

# coarse weight-3 target logical (same extraction as v1_targeted_check)
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
refined = [fi for fi in allb if fi not in corridor]
print(f"target: coarse weight-{len(support)} logical; corridor = its "
      f"{len(corridor)} touching faces (deflated); {len(refined)} faces stay fine")

mixed, info = refine_partial(b, refined)
mixed.validate()
assert code_k(mixed) == 8

# lift the coarse target into the mixed cellulation and certify it
lift = np.zeros(mixed.nE, dtype=np.uint8)
for e in support:
    for h in info["halves"][e]:
        lift[h] ^= 1
assert not ((mixed.HX().astype(int) @ lift.astype(int)) % 2).any(), \
    "lifted operator must commute with all X checks"

print("\nfine per-logical w :", sorted(per_logical_distances(fine_b)))
w_target_fine = None  # target class weight in fully-fine code, via lift there
fine_lift = np.zeros(fine_b.nE, dtype=np.uint8)
_, finfo = refine_partial(b, allb)
for e in support:
    for h in finfo["halves"][e]:
        fine_lift[h] ^= 1
assert not ((fine_b.HX().astype(int) @ fine_lift.astype(int)) % 2).any()
w_target_fine = min_weight_anticommuting(fine_b, fine_lift)
ds_mixed = per_logical_distances(mixed)
w_target = min_weight_anticommuting(mixed, lift)
print("mixed per-logical w:", sorted(ds_mixed))
print(f"TARGET class operator weight: {w_target_fine} -> {w_target}")
print(f"atoms returned to reservoir: {FB_ATOMS - atoms(mixed)} "
      f"(of {FB_ATOMS - atoms(b)} for a full exhale)")
print(f"hops/rd: fine {FB_HOPS} -> corridor {hops(mixed)} -> coarse {hops(b)}")

# ============================================================ economics
print()
print("=" * 74)
print("PART 3  cost of ONE logical operation (KT3 surgery model)")
print("        merge = d_t rounds x c_s*2*(w_log + d_t) hops; c_s = 1")
print("=" * 74)
moved = FB_ATOMS - atoms(mixed)          # corridor exhale, 1 hop/atom/way
for name, w, d_t, morph in [
        ("operate in fine (no morph)", 6, 6, 0),
        ("corridor exhale -> operate -> re-inhale", 3, 3, 2 * moved),
        ("full exhale -> operate -> re-inhale", 3, 3,
         2 * (FB_ATOMS - atoms(b)))]:
    surg = d_t * 1 * 2 * (w + d_t)
    print(f"  {name:42s} surgery {surg:4d} + morph {morph:4d} "
          f"= {surg + morph:4d} hops | transversal readout: {w} qubits")
print("\n(protection price: the operated class runs at d = its own weight "
      "for the window;\n bystander distances = 'mixed per-logical w' above)")
