"""CIRCUIT-LEVEL NOISE ATTACK on the breathing protocol (v2 roadmap item #1).

Replaces phenomenological MPP measurements with explicit ancilla-based
CNOT extraction. Uniform depolarizing circuit noise, strength p:
  - R  (|0>) on any qubit  -> X_ERROR(p) after
  - RX (|+>) on X-ancillas -> Z_ERROR(p) after
  - every CX               -> DEPOLARIZE2(p) after
  - M  (Z basis)           -> X_ERROR(p) before
  - MX                     -> Z_ERROR(p) before
  - every live qubit idle in a tick -> DEPOLARIZE1(p)
Both check sectors are present (X errors AND Z errors live), so ancilla
hook errors exist:
  - X fault on an X-check ancilla (control) mid-extraction propagates
    X onto the remaining data targets  -> correlated X chains (the threat
    to the Z-basis logicals we track).
Gate order within each round = greedy edge-coloring layers of the combined
(check, data) bipartite graph. Arbitrary but fixed and identical across
strategies -> fair relative comparison.

Morph bookkeeping (mirrors circuits.py, plus the X sector):
  inhale:  fold detector per OLD face: quads(f) xor last coarse Z(f)
           X sector: goes GAUGE across the morph. Old-vertex fine stars are
           NOT the coarse operator (halfA sits at the smaller endpoint, so
           fine stars contain fresh B atoms) -> first post-morph X outcomes
           are random at all vertices. Discovered by stim's determinism
           gate; an honest protocol cost (one blind X round per morph).
  exhale:  Z fold: first coarse Z(f) xor last quads(f) xor b_e outcomes;
           X sector gauge again (same reason, reversed).
           B-halves AND spokes destructively Z-measured with noise
           (spoke outcomes unused, as in v0).
stim's DEM construction certifies every fold deterministic = correctness gate.
"""
import time
import numpy as np
import stim
from cellulation import torus, bring_code, refine
from circuits import logical_Z_reps, Rec


# ---------------------------------------------------------------- scheduling

def color_layers(pairs):
    """Greedy edge coloring. pairs = list of (kind, anc_qid, data_qid).
    Returns list of layers; within a layer no qubit appears twice."""
    used = {}
    layers = {}
    for kind, a, d in pairs:
        ua = used.setdefault(a, set())
        ud = used.setdefault(d, set())
        c = 0
        while c in ua or c in ud:
            c += 1
        ua.add(c)
        ud.add(c)
        layers.setdefault(c, []).append((kind, a, d))
    return [layers[c] for c in sorted(layers)]


# ---------------------------------------------------------------- builder

class CL:
    def __init__(self, p, order_seed=None):
        self.c = stim.Circuit()
        self.rec = Rec()
        self.p = p
        self.order_seed = order_seed

    def idle(self, busy, live):
        idle = sorted(set(live) - set(busy))
        if idle:
            self.c.append("DEPOLARIZE1", idle, self.p)

    def reset_data(self, qubits):
        self.c.append("R", qubits)
        self.c.append("X_ERROR", qubits, self.p)

    def _phase(self, checks, kind, live):
        """One extraction phase (all Z checks, or all X checks).
        Sequential phases keep detectors deterministic on arbitrary graphs:
        any cross-sector stabilizer overlap is even (closed surface), so
        backward-propagated sensitivities cancel at phase boundaries. An
        interleaved schedule needs parity-matched gate orders per adjacent
        check pair (the surface-code 'Z-order' trick), which has no known
        canonical solution on curved {5,5} graphs.
        Returns measurement indices parallel to checks."""
        c, p = self.c, self.p
        anc = [a for a, _ in checks]
        if not anc:
            return []
        phase_live = set(live) | set(anc)
        # reset tick
        if kind == "z":
            c.append("R", anc)
            c.append("X_ERROR", anc, p)
        else:
            c.append("RX", anc)
            c.append("Z_ERROR", anc, p)
        self.idle(set(anc), phase_live)
        # CNOT layers
        pairs = [(kind, a, d) for a, ds in checks for d in ds]
        if self.order_seed is not None:
            # FIXED schedule: same permutation every round (fresh rng per call)
            np.random.default_rng(self.order_seed).shuffle(pairs)
        for layer in color_layers(pairs):
            busy = set()
            for _, a, d in layer:
                if kind == "z":
                    c.append("CX", [d, a])   # control data, target ancilla
                else:
                    c.append("CX", [a, d])   # control ancilla, target data
                c.append("DEPOLARIZE2", [a, d], p)
                busy.add(a)
                busy.add(d)
            self.idle(busy, phase_live)
        # measure tick
        if kind == "z":
            c.append("X_ERROR", anc, p)
            c.append("M", anc)
        else:
            c.append("Z_ERROR", anc, p)
            c.append("MX", anc)
        m = self.rec.take(len(anc))
        m = m if isinstance(m, list) else [m]
        self.idle(set(anc), phase_live)
        return m

    def sub_round(self, z_checks, x_checks, live_data):
        """One full syndrome round: Z phase then X phase."""
        zm = self._phase(z_checks, "z", live_data)
        xm = self._phase(x_checks, "x", live_data)
        return zm, xm

    def det(self, rec_indices):
        self.c.append("DETECTOR", [self.rec.rec(i) for i in rec_indices])


# ------------------------------------------------------------ check builders

def cell_checks(cell, qubit_of_edge, z_anc0, x_anc0):
    """(anc_qid, [data_qids]) lists for a cellulation's Z (face) and X
    (vertex-star) checks. qubit_of_edge maps edge idx -> physical qubit."""
    zc = [(z_anc0 + i, [qubit_of_edge[e] for e in f])
          for i, f in enumerate(cell.faces)]
    stars = [[] for _ in range(cell.nV)]
    for ei, e in enumerate(cell.edges):
        u, v = tuple(e)
        stars[u].append(qubit_of_edge[ei])
        stars[v].append(qubit_of_edge[ei])
    xc = [(x_anc0 + i, s) for i, s in enumerate(stars)]
    return zc, xc


# ------------------------------------------------------------------ circuits

def static_cl(cell, rounds, p, order_seed=None):
    n = cell.nE
    b = CL(p, order_seed)
    data = list(range(n))
    zc, xc = cell_checks(cell, list(range(n)), n, n + cell.nF)
    b.reset_data(data)
    prev_z = prev_x = None
    for r in range(rounds):
        zm, xm = b.sub_round(zc, xc, data)
        for i in range(cell.nF):
            b.det([zm[i]] + ([prev_z[i]] if prev_z else []))
        if prev_x:
            for i in range(cell.nV):
                b.det([xm[i], prev_x[i]])
        prev_z, prev_x = zm, xm
    b.c.append("X_ERROR", data, p)
    b.c.append("M", data)
    m = b.rec.take(n)
    m = m if isinstance(m, list) else [m]
    for i, f in enumerate(cell.faces):
        b.det([m[e] for e in f] + [prev_z[i]])
    for li, log in enumerate(logical_Z_reps(cell)):
        b.c.append("OBSERVABLE_INCLUDE", [b.rec.rec(m[e]) for e in log], li)
    return b.c


def breathing_cl(cell, R1, RF, R2, p, order_seed=None):
    assert R1 >= 1 and RF >= 1 and R2 >= 1
    fine, info = refine(cell)
    A, B = info["halfA"], info["halfB"]
    spokes = sorted(set(info["spokes"].values()))
    nq = fine.nE
    # ancilla layout after the data block
    cz0 = nq
    cx0 = cz0 + cell.nF
    fz0 = cx0 + cell.nV
    fx0 = fz0 + fine.nF

    coarse_data = [A[e] for e in range(cell.nE)]
    new_qubits = sorted([B[e] for e in range(cell.nE)] + spokes)
    fine_data = list(range(fine.nE))

    czc, cxc = cell_checks(cell, {e: A[e] for e in range(cell.nE)}, cz0, cx0)
    fzc, fxc = cell_checks(fine, list(range(fine.nE)), fz0, fx0)

    logicals = logical_Z_reps(cell)
    b = CL(p, order_seed)

    b.reset_data(coarse_data)
    prev_cz = prev_cx = None
    for r in range(R1):
        zm, xm = b.sub_round(czc, cxc, coarse_data)
        for i in range(cell.nF):
            b.det([zm[i]] + ([prev_cz[i]] if prev_cz else []))
        if prev_cx:
            for i in range(cell.nV):
                b.det([xm[i], prev_cx[i]])
        prev_cz, prev_cx = zm, xm

    # ---- INHALE ----
    b.reset_data(new_qubits)
    zm, xm = b.sub_round(fzc, fxc, fine_data)
    for i in range(cell.nF):                       # Z fold per old face
        b.det([zm[q] for q in info["face_quads"][i]] + [prev_cz[i]])
    # X sector: NO fold. Old-vertex fine stars contain fresh B atoms
    # (halfA attaches to the smaller endpoint), so first fine X outcomes
    # are gauge at ALL vertices. One round of X-sector blindness per morph.
    prev_fz, prev_fx = zm, xm
    for r in range(RF - 1):
        zm, xm = b.sub_round(fzc, fxc, fine_data)
        for i in range(fine.nF):
            b.det([zm[i], prev_fz[i]])
        for v in range(fine.nV):
            b.det([xm[v], prev_fx[v]])
        prev_fz, prev_fx = zm, xm

    # ---- EXHALE ----
    b_list = [B[e] for e in range(cell.nE)]
    b.c.append("X_ERROR", b_list, p)
    b.c.append("M", b_list)
    mB = b.rec.take(cell.nE)
    mB = mB if isinstance(mB, list) else [mB]
    b.c.append("X_ERROR", spokes, p)               # noisy, outcomes unused
    b.c.append("M", spokes)
    b.rec.take(len(spokes))
    zm, xm = b.sub_round(czc, cxc, coarse_data)
    for i, f in enumerate(cell.faces):             # Z fold back
        b.det([zm[i]] + [prev_fz[q] for q in info["face_quads"][i]] +
              [mB[e] for e in f])
    # X sector gauge again after exhale (same reason, reversed).
    prev_cz, prev_cx = zm, xm
    for r in range(R2 - 1):
        zm, xm = b.sub_round(czc, cxc, coarse_data)
        for i in range(cell.nF):
            b.det([zm[i], prev_cz[i]])
        for i in range(cell.nV):
            b.det([xm[i], prev_cx[i]])
        prev_cz, prev_cx = zm, xm

    for li, log in enumerate(logicals):
        b.c.append("OBSERVABLE_INCLUDE",
                   [b.rec.rec(mB[e]) for e in log], li)

    b.c.append("X_ERROR", coarse_data, p)
    b.c.append("M", coarse_data)
    m = b.rec.take(cell.nE)
    m = m if isinstance(m, list) else [m]
    for i, f in enumerate(cell.faces):
        b.det([m[e] for e in f] + [prev_cz[i]])
    for li, log in enumerate(logicals):
        b.c.append("OBSERVABLE_INCLUDE", [b.rec.rec(m[e]) for e in log], li)
    return b.c


def fine_static_cl(cell, rounds, p, order_seed=None):
    fine, _ = refine(cell)
    return static_cl(fine, rounds, p, order_seed)


# ------------------------------------------------------------------ analysis

def verify(name, circ):
    t0 = time.time()
    dem = circ.detector_error_model(decompose_errors=True)
    print(f"  [ok] {name}: dets={circ.num_detectors} obs={circ.num_observables} "
          f"dem_ok ({time.time()-t0:.1f}s)")
    return dem


def fault_distance(name, circ, autopsy=False):
    t0 = time.time()
    err = circ.shortest_graphlike_error(canonicalize_circuit_errors=autopsy)
    fd = len(err)
    print(f"  fd[{name}] = {fd}  ({time.time()-t0:.1f}s)")
    return fd, err


if __name__ == "__main__":
    P_STRUCT = 1e-3   # any nonzero p; fd is structural
    for cell in [torus(3), bring_code()]:
        print(f"== {cell.name} : determinism gate ==")
        verify("static coarse", static_cl(cell, 4, P_STRUCT))
        verify("static fine  ", fine_static_cl(cell, 4, P_STRUCT))
        verify("breath 2/3/2 ", breathing_cl(cell, 2, 3, 2, P_STRUCT))
        verify("breath 1/1/1 ", breathing_cl(cell, 1, 1, 1, P_STRUCT))
