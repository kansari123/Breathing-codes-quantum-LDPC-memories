"""v1-B: partial breathing. Inflate a sub-region; does the targeted
logical's distance rise, or does homology route around the armor?

Torus test: inflate 0..3 full columns of faces (a column is a transversal
cut for the horizontal logical class -- every horizontal cycle must cross it).
Bring control: inflate only the faces touching one logical representative
(NOT a transversal cut) -- homology should dodge.
"""
import numpy as np
from cellulation import (Cellulation, gf2_rank, gf2_nullspace, code_k,
                         torus, bring_code)
from scipy.optimize import milp, LinearConstraint, Bounds


def refine_partial(cell, faceset):
    """Quad-refine only faces in faceset; neighbors keep their face with
    subdivided edges replaced by their halves. Returns mixed cellulation
    and cost info."""
    faceset = set(faceset)
    sub_edges = set()
    for fi in faceset:
        sub_edges |= set(cell.faces[fi])
    mid = {}
    nv = cell.nV
    for e in sorted(sub_edges):
        mid[e] = nv
        nv += 1
    cen = {}
    for fi in sorted(faceset):
        cen[fi] = nv
        nv += 1
    edges, eidx = [], {}

    def add(u, v):
        key = frozenset({u, v})
        if key not in eidx:
            eidx[key] = len(edges)
            edges.append(key)
        return eidx[key]

    halves = {}
    for e in range(cell.nE):
        u, v = sorted(tuple(cell.edges[e]))
        if e in sub_edges:
            halves[e] = (add(u, mid[e]), add(mid[e], v))
        else:
            halves[e] = (add(u, v),)
    faces = []
    for fi, f in enumerate(cell.faces):
        # order edges into a vertex cycle (as in refine())
        fedges = list(f)
        e0 = cell.edges[fedges[0]]
        remaining = fedges[1:]
        u0, v0 = tuple(e0)
        cyc_v, cyc_e = [u0, v0], [fedges[0]]
        while remaining:
            cur = cyc_v[-1]
            for e in remaining:
                a, b = tuple(cell.edges[e])
                if a == cur or b == cur:
                    remaining.remove(e)
                    cyc_e.append(e)
                    cyc_v.append(b if a == cur else a)
                    break
        cyc_v = cyc_v[:-1]
        p = len(cyc_e)
        if fi not in faceset:
            fe = []
            for e in cyc_e:
                fe += list(halves[e])
            faces.append(tuple(fe))
        else:
            spoke = {e: add(cen[fi], mid[e]) for e in cyc_e}
            for i in range(p):
                e_i, e_j = cyc_e[i], cyc_e[(i + 1) % p]
                w = cyc_v[(i + 1) % p]
                faces.append((add(w, mid[e_i]), add(w, mid[e_j]),
                              spoke[e_i], spoke[e_j]))
    mixed = Cellulation(nv, edges, faces, name=cell.name + f"_part{len(faceset)}")
    new_data = mixed.nE - cell.nE
    new_anc = (mixed.nF + mixed.nV) - (cell.nF + cell.nV)
    hops = int(mixed.HZ().sum() + mixed.HX().sum())
    return mixed, dict(new_data=new_data, new_anc=new_anc, hops=hops, halves=halves)


def per_logical_distances(cell):
    """ILP min weight per Z-logical class (order = extraction order)."""
    HZ = cell.HZ().astype(float)
    kerX = gf2_nullspace(cell.HX())
    stack = cell.HZ().copy()
    Lzs = []
    for v in kerX:
        test = np.vstack([stack, v[None, :]])
        if gf2_rank(test) > gf2_rank(stack):
            stack = test
            Lzs.append(v)
    n, m = cell.nE, HZ.shape[0]
    out = []
    for Lz in Lzs:
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
        r[:n] = Lz
        r[n + m] = -2
        A.append(r); lb.append(1); ub.append(1)
        res = milp(c=c, constraints=LinearConstraint(np.array(A), lb, ub),
                   bounds=Bounds(np.zeros(nvar),
                                 np.concatenate([np.ones(n),
                                                 np.full(m, HZ.sum(1).max()),
                                                 [n]])),
                   integrality=np.ones(nvar))
        out.append(int(round(res.fun)))
    return out


if __name__ == "__main__":
    L = 3
    base = torus(L)
    print("=== torus: inflate 0..3 full columns (transversal cuts) ===")
    print(f"{'cols':>4s} {'n':>4s} {'k':>2s} {'per-logical d':>16s} "
          f"{'+atoms':>7s} {'hops/rd':>8s}")
    base_hops = int(base.HZ().sum() + base.HX().sum())
    for ncols in range(0, L + 1):
        faceset = [i * L + j for j in range(ncols) for i in range(L)]
        if ncols == 0:
            cellm, info = base, dict(new_data=0, new_anc=0, hops=base_hops)
        else:
            cellm, info = refine_partial(base, faceset)
            cellm.validate()
        k = code_k(cellm)
        ds = per_logical_distances(cellm)
        extra = info["new_data"] + info["new_anc"]
        print(f"{ncols:4d} {cellm.nE:4d} {k:2d} {str(sorted(ds)):>16s} "
              f"{extra:7d} {info['hops']:8d}")

    print("\n=== Bring control: inflate faces touching one logical rep ===")
    b = bring_code()
    kerX = gf2_nullspace(b.HX())
    stack = b.HZ().copy()
    Lz0 = None
    for v in kerX:
        t = np.vstack([stack, v[None, :]])
        if gf2_rank(t) > gf2_rank(stack):
            Lz0 = v
            break
    support = set(np.flatnonzero(Lz0).tolist())
    touching = [fi for fi, f in enumerate(b.faces) if support & set(f)]
    print(f"logical rep weight {len(support)}; inflating {len(touching)} faces")
    mixed, info = refine_partial(b, touching)
    mixed.validate()
    ds0 = per_logical_distances(b)
    ds1 = per_logical_distances(mixed)
    print(f"k: {code_k(b)} -> {code_k(mixed)}")
    print(f"per-logical d before: {sorted(ds0)}")
    print(f"per-logical d after : {sorted(ds1)}  "
          f"(+{info['new_data'] + info['new_anc']} atoms)")
