"""Cellulations of closed surfaces, quad refinement, and code parameters.

A cellulation is (V, E, F) with:
  edges: list of frozenset({u, v}) vertex pairs (qubits live on edges)
  faces: list of tuples of edge indices (each face = cycle of edges)

Surface code convention:
  H_X rows = vertex stars   (X-type checks, detect Z errors)
  H_Z rows = face boundaries (Z-type checks, detect X errors)
We simulate the X-error / Z-check sector only (CSS symmetry).
"""
import numpy as np
from itertools import combinations
from scipy.optimize import milp, LinearConstraint, Bounds


class Cellulation:
    def __init__(self, n_vertices, edges, faces, name=""):
        self.nV = n_vertices
        self.edges = [frozenset(e) for e in edges]
        self.faces = [tuple(f) for f in faces]
        self.name = name
        self.nE = len(self.edges)
        self.nF = len(self.faces)

    def euler(self):
        return self.nV - self.nE + self.nF

    def HZ(self):
        """Face-edge incidence (Z checks)."""
        H = np.zeros((self.nF, self.nE), dtype=np.uint8)
        for i, f in enumerate(self.faces):
            for e in f:
                H[i, e] ^= 1
        return H

    def HX(self):
        """Vertex-edge incidence (X checks)."""
        H = np.zeros((self.nV, self.nE), dtype=np.uint8)
        for i, e in enumerate(self.edges):
            u, v = tuple(e)
            H[u, i] ^= 1
            H[v, i] ^= 1
        return H

    def validate(self):
        # every edge in exactly 2 faces (closed surface)
        cnt = np.array(self.HZ().sum(axis=0)).flatten()
        assert np.all(cnt == 2), f"{self.name}: edge-face incidence != 2: {cnt}"
        # CSS commutation: HX @ HZ.T = 0 mod 2
        M = (self.HX().astype(int) @ self.HZ().T.astype(int)) % 2
        assert not M.any(), f"{self.name}: HX HZ^T != 0"
        return True


def gf2_rank(A):
    A = A.copy().astype(np.uint8)
    r, rows, cols = 0, A.shape[0], A.shape[1]
    for c in range(cols):
        piv = None
        for rr in range(r, rows):
            if A[rr, c]:
                piv = rr
                break
        if piv is None:
            continue
        A[[r, piv]] = A[[piv, r]]
        for rr in range(rows):
            if rr != r and A[rr, c]:
                A[rr] ^= A[r]
        r += 1
        if r == rows:
            break
    return r


def gf2_nullspace(A):
    """Basis of the nullspace of A over GF(2)."""
    A = A.copy().astype(np.uint8)
    rows, cols = A.shape
    pivots = {}
    r = 0
    for c in range(cols):
        piv = None
        for rr in range(r, rows):
            if A[rr, c]:
                piv = rr
                break
        if piv is None:
            continue
        A[[r, piv]] = A[[piv, r]]
        for rr in range(rows):
            if rr != r and A[rr, c]:
                A[rr] ^= A[r]
        pivots[c] = r
        r += 1
    free = [c for c in range(cols) if c not in pivots]
    basis = []
    for fc in free:
        v = np.zeros(cols, dtype=np.uint8)
        v[fc] = 1
        for pc, pr in pivots.items():
            if A[pr, fc]:
                v[pc] = 1
        basis.append(v)
    return np.array(basis, dtype=np.uint8) if basis else np.zeros((0, cols), np.uint8)


def code_k(cell):
    n = cell.nE
    return n - gf2_rank(cell.HX()) - gf2_rank(cell.HZ())


def logical_X_reps(cell):
    """Representatives of X-logicals: nullspace(HZ... ) careful:
    X-error e is undetected iff HZ e = 0; harmless iff e in rowspace(HX^T)
    (products of X stars). Logical classes = ker(HZ)/im(HX^T).
    Returns basis vectors of ker(HZ) that are independent mod im(HX^T)."""
    ker = gf2_nullspace(cell.HZ())
    HX = cell.HX()
    reps = []
    stack = HX.copy()
    base_rank = gf2_rank(stack)
    for v in ker:
        test = np.vstack([stack, v[None, :]])
        if gf2_rank(test) > gf2_rank(stack):
            stack = test
            reps.append(v)
    # first (rank(stack_final) - base_rank) additions are the logical reps
    return np.array(reps, dtype=np.uint8)


def distance_X_ilp(cell, verbose=False):
    """Min weight X-error with HZ e = 0 and e not in im(HX^T), via ILP.
    Solve for each anti-commuting Z-logical constraint: min |e| s.t.
      HZ e = 0 mod 2   and   <Lz, e> = 1 mod 2
    where Lz ranges over a basis of Z-logicals (ker(HX)/im(HZ^T)).
    d_X = min over Lz."""
    HZ = cell.HZ().astype(float)
    # Z logicals: ker(HX)/im(HZ^T)
    kerX = gf2_nullspace(cell.HX())
    HZm = cell.HZ()
    stack = HZm.copy()
    Lzs = []
    for v in kerX:
        test = np.vstack([stack, v[None, :]])
        if gf2_rank(test) > gf2_rank(stack):
            stack = test
            Lzs.append(v)
    n = cell.nE
    best = None
    for Lz in Lzs:
        m = HZ.shape[0]
        # vars: x (n binary), y (m int, HZ x = 2y), z (1 int, Lz x = 2z + 1)
        nvar = n + m + 1
        c = np.zeros(nvar)
        c[:n] = 1.0
        A_rows, lb, ub = [], [], []
        for i in range(m):
            row = np.zeros(nvar)
            row[:n] = HZ[i]
            row[n + i] = -2.0
            A_rows.append(row)
            lb.append(0.0)
            ub.append(0.0)
        row = np.zeros(nvar)
        row[:n] = Lz.astype(float)
        row[n + m] = -2.0
        A_rows.append(row)
        lb.append(1.0)
        ub.append(1.0)
        cons = LinearConstraint(np.array(A_rows), lb, ub)
        bounds = Bounds(
            np.zeros(nvar),
            np.concatenate([np.ones(n), np.full(m, np.ceil(HZ.sum(1).max() / 2)), [n / 2]]),
        )
        res = milp(c=c, constraints=cons, bounds=bounds,
                   integrality=np.ones(nvar))
        assert res.success, res.message
        w = int(round(res.fun))
        if verbose:
            print(f"  logical -> min weight {w}")
        best = w if best is None else min(best, w)
    return best


# ---------------------------------------------------------------- builders

def torus(L):
    """L x L square-tiled torus. Qubits on edges."""
    def vid(i, j):
        return (i % L) * L + (j % L)
    edges, eidx = [], {}
    for i in range(L):
        for j in range(L):
            for (u, v) in [(vid(i, j), vid(i, j + 1)), (vid(i, j), vid(i + 1, j))]:
                e = frozenset({u, v})
                if e not in eidx:
                    eidx[e] = len(edges)
                    edges.append(e)
    # careful: L=2 would collapse edges; require L>=3
    assert L >= 3
    faces = []
    for i in range(L):
        for j in range(L):
            f = [eidx[frozenset({vid(i, j), vid(i, j + 1)})],
                 eidx[frozenset({vid(i + 1, j), vid(i + 1, j + 1)})],
                 eidx[frozenset({vid(i, j), vid(i + 1, j)})],
                 eidx[frozenset({vid(i, j + 1), vid(i + 1, j + 1)})]]
            faces.append(tuple(f))
    return Cellulation(L * L, edges, faces, name=f"torus_L{L}")


def icosahedron_graph():
    """Vertices/edges of the icosahedron via golden-ratio coordinates."""
    phi = (1 + 5 ** 0.5) / 2
    verts = []
    for s1 in (+1, -1):
        for s2 in (+1, -1):
            verts += [(0, s1 * 1, s2 * phi), (s1 * 1, s2 * phi, 0), (s2 * phi, 0, s1 * 1)]
    verts = [np.array(v, float) for v in verts]
    n = len(verts)
    edge_len = 2.0  # nearest-neighbor distance
    edges = []
    for i, j in combinations(range(n), 2):
        if abs(np.linalg.norm(verts[i] - verts[j]) - edge_len) < 1e-6:
            edges.append(frozenset({i, j}))
    assert n == 12 and len(edges) == 30, (n, len(edges))
    return verts, edges


def bring_code():
    """[[30,8,3]] hyperbolic {5,5} code on Bring's surface (genus 4).

    Combinatorics of the great dodecahedron:
      vertices = icosahedron vertices (12)
      edges    = icosahedron edges (30)   <- qubits
      faces    = for each icosahedron vertex v, the pentagon formed by
                 the 5 edges among v's neighbors (12 pentagons)
    """
    verts, edges = icosahedron_graph()
    eidx = {e: i for i, e in enumerate(edges)}
    adj = {i: set() for i in range(12)}
    for e in edges:
        u, v = tuple(e)
        adj[u].add(v)
        adj[v].add(u)
    faces = []
    for v in range(12):
        nb = list(adj[v])
        # order neighbors into the 5-cycle they form
        cyc = [nb[0]]
        while len(cyc) < 5:
            cur = cyc[-1]
            nxt = [w for w in nb if w in adj[cur] and w not in cyc]
            cyc.append(nxt[0])
        f = [eidx[frozenset({cyc[i], cyc[(i + 1) % 5]})] for i in range(5)]
        faces.append(tuple(f))
    return Cellulation(12, edges, faces, name="bring_30_8")


# ------------------------------------------------------------- refinement

def refine(cell):
    """1->p quad refinement ("inhale").

    Each edge e = (u,v) splits at a midpoint m_e into halves; each face f
    gets a center c_f joined by spokes to the midpoints of its edges.
    Each p-gon face becomes p quads.

    Qubit bookkeeping for the physical morph:
      half A of edge e  -> KEEPS the old physical qubit of e   (index map old_qubit)
      half B of edge e  -> NEW qubit, initialized |0>
      spokes            -> NEW qubits, initialized |0>

    Returns (fine_cell, info) where info maps structure between levels.
    """
    # new vertex ids: 0..nV-1 old, then midpoints (per edge), then centers (per face)
    mid = {i: cell.nV + i for i in range(cell.nE)}
    cen = {i: cell.nV + cell.nE + i for i in range(cell.nF)}
    nV = cell.nV + cell.nE + cell.nF

    edges = []
    eidx = {}

    def add(u, v):
        e = frozenset({u, v})
        if e not in eidx:
            eidx[e] = len(edges)
            edges.append(e)
        return eidx[e]

    halfA, halfB = {}, {}  # old edge idx -> fine edge idx
    for i, e in enumerate(cell.edges):
        u, v = sorted(tuple(e))
        halfA[i] = add(u, mid[i])
        halfB[i] = add(mid[i], v)

    spokes = {}  # (face, old_edge) -> fine edge idx
    faces = []
    face_quads = []  # per old face: list of fine face indices
    for fi, f in enumerate(cell.faces):
        # order the face's edges into a cycle of (vertex, edge) steps
        fedges = list(f)
        # build vertex cycle
        e0 = cell.edges[fedges[0]]
        # find orientation: walk
        remaining = fedges[1:]
        u0, v0 = tuple(e0)
        cyc_v = [u0, v0]
        cyc_e = [fedges[0]]
        while remaining:
            cur = cyc_v[-1]
            nxt = None
            for e in remaining:
                a, b = tuple(cell.edges[e])
                if a == cur:
                    nxt, other = e, b
                    break
                if b == cur:
                    nxt, other = e, a
                    break
            assert nxt is not None, "face edges do not form a cycle"
            remaining.remove(nxt)
            cyc_e.append(nxt)
            cyc_v.append(other)
        assert cyc_v[-1] == cyc_v[0]
        cyc_v = cyc_v[:-1]
        p = len(cyc_e)
        for e in cyc_e:
            spokes[(fi, e)] = add(cen[fi], mid[e])
        quads = []
        for i in range(p):
            e_i, e_j = cyc_e[i], cyc_e[(i + 1) % p]
            shared_v = cyc_v[(i + 1) % p]
            q = (add(shared_v, mid[e_i]), add(shared_v, mid[e_j]),
                 spokes[(fi, e_i)], spokes[(fi, e_j)])
            quads.append(len(faces))
            faces.append(q)
        face_quads.append(quads)

    fine = Cellulation(nV, edges, faces, name=cell.name + "_fine")
    info = dict(halfA=halfA, halfB=halfB, spokes=spokes, face_quads=face_quads)
    return fine, info


def lift_cycle(cycle_edges, info):
    """A coarse Z-logical (set of coarse edge indices) -> fine edge set."""
    out = []
    for e in cycle_edges:
        out += [info["halfA"][e], info["halfB"][e]]
    return out


if __name__ == "__main__":
    for cell in [torus(3), bring_code()]:
        cell.validate()
        fine, info = refine(cell)
        fine.validate()
        k, kf = code_k(cell), code_k(fine)
        print(f"{cell.name}: V={cell.nV} E={cell.nE} F={cell.nF} chi={cell.euler()}  k={k}")
        print(f"{fine.name}: V={fine.nV} E={fine.nE} F={fine.nF} chi={fine.euler()}  k={kf}")
        assert cell.euler() == fine.euler() and k == kf
        d = distance_X_ilp(cell)
        df = distance_X_ilp(fine)
        print(f"  d_X coarse = {d},  d_X fine = {df}")
