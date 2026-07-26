"""v1-A: deep breaths. Does the shuttling cost concentrate on the seams?

Layout model: after refinement, each ORIGINAL face's interior sub-lattice is
laid out as a local flat patch (interleaved ancillas, Rydberg gates, 0 hops).
A check is 'interior' iff every qubit it touches is interior to the SAME
original face; otherwise it is a 'seam' check and its ancilla tours its w
qubits (w hops).

Iterate the 1->4 quad refinement with provenance labels:
  edge label = ('seam', original_edge) if it lies on an original edge,
               ('int', original_face) if it lies strictly inside one.
"""
import numpy as np
from cellulation import Cellulation, refine, gf2_rank, code_k, torus, bring_code


def refine_tracked(cell, edge_labels, face_orig):
    """One quad refinement, propagating provenance.
    edge_labels: list per edge; face_orig: list per face -> original face id."""
    fine, info = refine(cell)
    new_labels = [None] * fine.nE
    for e in range(cell.nE):
        new_labels[info["halfA"][e]] = edge_labels[e]
        new_labels[info["halfB"][e]] = edge_labels[e]
    for (fi, e), s in info["spokes"].items():
        new_labels[s] = ("int", face_orig[fi])
    new_face_orig = [None] * fine.nF
    for fi, quads in enumerate(info["face_quads"]):
        for q in quads:
            new_face_orig[q] = face_orig[fi]
    assert all(l is not None for l in new_labels)
    return fine, new_labels, new_face_orig


def classify_hops(cell, edge_labels):
    """(interior_checks, seam_checks, seam_hops_per_round) over both sectors."""
    def owner(e):
        kind, idx = edge_labels[e]
        return idx if kind == "int" else None
    n_int = n_seam = hops = 0
    for f in cell.faces:                       # Z checks
        owners = {owner(e) for e in f}
        if None not in owners and len(owners) == 1:
            n_int += 1
        else:
            n_seam += 1
            hops += len(f)
    incid = [[] for _ in range(cell.nV)]       # X checks (vertex stars)
    for i, e in enumerate(cell.edges):
        u, v = tuple(e)
        incid[u].append(i)
        incid[v].append(i)
    for star in incid:
        owners = {owner(e) for e in star}
        if None not in owners and len(owners) == 1:
            n_int += 1
        else:
            n_seam += 1
            hops += len(star)
    return n_int, n_seam, hops


if __name__ == "__main__":
    D_T = {1: 3, 2: 6, 4: 12, 8: 24}   # d = 3*l (exact at l<=2 via ILP;
                                       # upper bound = lifted logical at l>2)
    print(f"{'code':10s} {'l':>2s} {'n':>5s} {'k':>2s} {'d':>3s} "
          f"{'naive hops/rd':>13s} {'seam hops/rd':>12s} {'interior %':>10s}")
    results = {}
    for base in [torus(3), bring_code()]:
        cell = base
        labels = [("seam", e) for e in range(cell.nE)]
        forig = list(range(cell.nF))
        l = 1
        rows = []
        while l <= 8:
            naive = int(cell.HZ().sum() + cell.HX().sum())
            n_int, n_seam, seam_hops = classify_hops(cell, labels)
            k = code_k(cell) if cell.nE <= 2000 else -1
            rows.append((l, cell.nE, k, D_T[l], naive, seam_hops,
                         100 * n_int / (n_int + n_seam)))
            print(f"{base.name:10s} {l:2d} {cell.nE:5d} {k:2d} {D_T[l]:3d} "
                  f"{naive:13d} {seam_hops:12d} {rows[-1][-1]:9.1f}%")
            if l == 8:
                break
            cell, labels, forig = refine_tracked(cell, labels, forig)
            l *= 2
        results[base.name] = rows

    # ---- kill-test 3 rematch at matched d=12, protect all 8 logicals ----
    print("\n--- KT3 rematch: Bring block at l=4 (d=12) vs 8 vaults d=12 ---")
    rows = results["bring_30_8"]
    l4 = rows[2]
    n4, seam4 = l4[1], l4[5]
    anc4 = None
    # recompute atoms: fine data + fine ancillas vs coarse
    cell = bring_code()
    labels = [("seam", e) for e in range(cell.nE)]
    forig = list(range(cell.nF))
    for _ in range(2):
        cell, labels, forig = refine_tracked(cell, labels, forig)
    extra_atoms_A = (cell.nE - 30) + (cell.nF + cell.nV - 24)
    patch = 2 * 12 * 12 - 1
    print(f"breathe l=4: extra atoms={extra_atoms_A}, seam hops/rd={seam4} "
          f"(naive would be {l4[4]})")
    print(f"8 vaults d=12: atoms={8 * patch}, hops/rd=0, "
          f"place+retrieve={16 * patch}, surgeries(2x8, nominal)="
          f"{8 * 2 * 12 * 2 * (3 + 12)}")
    W = np.arange(1, 201)
    A = 2 * extra_atoms_A + W * seam4
    B_od = 16 * patch + 8 * 2 * (12 * 2 * (3 + 12))
    cross = W[np.where(B_od < A)[0]]
    print(f"hop crossover (on-demand vaults win): W >= "
          f"{cross[0] if len(cross) else '>200'}")
    print(f"space: breathe {extra_atoms_A} vs vaults {8 * patch} "
          f"({8 * patch / extra_atoms_A:.1f}x)")
