"""Stim circuits for static and breathing surface-code memories.

Phenomenological noise, X-error / Z-check sector:
  each round: X_ERROR(p) on live data, then noisy MPP of every Z face check.
Initial state |0...0> (Z-basis), final ideal Z readout of live data.

Breathing circuit structure:
  R1 coarse rounds
  INHALE: reset new qubits (B-halves + spokes), fine round with per-OLD-face
          fold detectors (product of the old face's quads vs last coarse msmt)
  RF-1 more fine rounds (quad vs quad detectors)
  EXHALE: noisy Z-measure of B-halves (spokes ideal; unused in bookkeeping),
          coarse round with fold detectors
          coarse[f] (+) last quads(f) (+) B-outcomes of f's edges
  R2-1 more coarse rounds
  final readout in the coarse code
Logical observables (all k of them) accumulate the B-half exhale outcomes on
their cycle edges, then the final data readout.
"""
import stim
import numpy as np
from cellulation import (Cellulation, refine, gf2_rank, gf2_nullspace,
                         lift_cycle)


def logical_Z_reps(cell):
    """Basis of Z-logicals: ker(HX) mod im(HZ). X errors flip these."""
    kerX = gf2_nullspace(cell.HX())
    HZ = cell.HZ()
    stack = HZ.copy()
    reps = []
    for v in kerX:
        test = np.vstack([stack, v[None, :]])
        if gf2_rank(test) > gf2_rank(stack):
            stack = test
            reps.append(v)
    return [np.flatnonzero(r).tolist() for r in reps]


class Rec:
    """Measurement record bookkeeping."""
    def __init__(self):
        self.n = 0

    def take(self, k=1):
        out = list(range(self.n, self.n + k))
        self.n += k
        return out if k > 1 else out[0]

    def rec(self, idx):
        return stim.target_rec(idx - self.n)


def _round(c, rec, data, faces_qubits, p, meas_p=None):
    """One noisy syndrome round. Returns list of msmt indices, per face."""
    if meas_p is None:
        meas_p = p
    c.append("X_ERROR", data, p)
    idxs = []
    for fq in faces_qubits:
        targets = []
        for q in fq:
            targets.append(stim.target_z(q))
            targets.append(stim.target_combiner())
        c.append("MPP", targets[:-1], meas_p)
        idxs.append(rec.take())
    return idxs


def _face_qubits(cell, qubit_of_edge):
    return [[qubit_of_edge[e] for e in f] for f in cell.faces]


def static_circuit(cell, rounds, p):
    c = stim.Circuit()
    rec = Rec()
    n = cell.nE
    data = list(range(n))
    fq = _face_qubits(cell, list(range(n)))
    c.append("R", data)
    prev = None
    for r in range(rounds):
        cur = _round(c, rec, data, fq, p)
        for i in range(cell.nF):
            if prev is None:
                c.append("DETECTOR", [rec.rec(cur[i])])
            else:
                c.append("DETECTOR", [rec.rec(cur[i]), rec.rec(prev[i])])
        prev = cur
    # final ideal readout
    c.append("X_ERROR", data, p)
    c.append("M", data)
    m0 = rec.take(n)
    m0 = m0 if isinstance(m0, list) else [m0]
    for i, f in enumerate(cell.faces):
        targs = [rec.rec(m0[e]) for e in f] + [rec.rec(prev[i])]
        c.append("DETECTOR", targs)
    for li, log in enumerate(logical_Z_reps(cell)):
        c.append("OBSERVABLE_INCLUDE", [rec.rec(m0[e]) for e in log], li)
    return c


def breathing_circuit(cell, R1, RF, R2, p):
    """Coarse R1 -> inhale -> fine RF -> exhale -> coarse R2 -> readout."""
    assert R1 >= 1 and RF >= 1 and R2 >= 1
    fine, info = refine(cell)
    c = stim.Circuit()
    rec = Rec()

    # qubit ids = fine edge ids; coarse edge e lives on atom halfA[e]
    A = info["halfA"]
    B = info["halfB"]
    spokes = sorted(set(info["spokes"].values()))
    coarse_data = [A[e] for e in range(cell.nE)]
    new_qubits = sorted([B[e] for e in range(cell.nE)] + spokes)
    fine_data = list(range(fine.nE))
    coarse_fq = _face_qubits(cell, {e: A[e] for e in range(cell.nE)})
    fine_fq = _face_qubits(fine, list(range(fine.nE)))

    logicals = logical_Z_reps(cell)

    c.append("R", coarse_data)
    prev_c = None
    for r in range(R1):
        cur = _round(c, rec, coarse_data, coarse_fq, p)
        for i in range(cell.nF):
            t = [rec.rec(cur[i])] + ([rec.rec(prev_c[i])] if prev_c else [])
            c.append("DETECTOR", t)
        prev_c = cur

    # ---- INHALE ----
    c.append("R", new_qubits)
    cur_f = _round(c, rec, fine_data, fine_fq, p)
    for i in range(cell.nF):  # fold detectors per OLD face
        t = [rec.rec(cur_f[q]) for q in info["face_quads"][i]]
        t.append(rec.rec(prev_c[i]))
        c.append("DETECTOR", t)
    prev_f = cur_f
    for r in range(RF - 1):
        cur_f = _round(c, rec, fine_data, fine_fq, p)
        for i in range(fine.nF):
            c.append("DETECTOR", [rec.rec(cur_f[i]), rec.rec(prev_f[i])])
        prev_f = cur_f

    # ---- EXHALE ----
    b_list = [B[e] for e in range(cell.nE)]
    c.append("X_ERROR", b_list, p)      # error just before destructive msmt
    c.append("M", b_list)
    mB_first = rec.take(cell.nE)
    mB = {e: mB_first[e] for e in range(cell.nE)} if isinstance(mB_first, list) \
        else {0: mB_first}
    c.append("M", spokes)               # ideal; outcomes unused
    rec.take(len(spokes))
    cur_c = _round(c, rec, coarse_data, coarse_fq, p)
    for i, f in enumerate(cell.faces):
        t = [rec.rec(cur_c[i])]
        t += [rec.rec(prev_f[q]) for q in info["face_quads"][i]]
        t += [rec.rec(mB[e]) for e in f]
        c.append("DETECTOR", t)
    prev_c = cur_c
    for r in range(R2 - 1):
        cur_c = _round(c, rec, coarse_data, coarse_fq, p)
        for i in range(cell.nF):
            c.append("DETECTOR", [rec.rec(cur_c[i]), rec.rec(prev_c[i])])
        prev_c = cur_c

    # observables absorb B outcomes on their cycles
    for li, log in enumerate(logicals):
        c.append("OBSERVABLE_INCLUDE", [rec.rec(mB[e]) for e in log], li)

    # final ideal readout in coarse code
    c.append("X_ERROR", coarse_data, p)
    c.append("M", coarse_data)
    mfin = rec.take(cell.nE)
    m = {e: mfin[e] for e in range(cell.nE)}
    for i, f in enumerate(cell.faces):
        targs = [rec.rec(m[e]) for e in f] + [rec.rec(prev_c[i])]
        c.append("DETECTOR", targs)
    for li, log in enumerate(logicals):
        c.append("OBSERVABLE_INCLUDE", [rec.rec(m[e]) for e in log], li)
    return c


def fine_static_circuit(cell, rounds, p):
    fine, _ = refine(cell)
    return static_circuit(fine, rounds, p)


if __name__ == "__main__":
    from cellulation import torus, bring_code
    for cell in [torus(3), bring_code()]:
        for name, circ in [
            ("static coarse", static_circuit(cell, 5, 0.01)),
            ("static fine  ", fine_static_circuit(cell, 5, 0.01)),
            ("breathing    ", breathing_circuit(cell, 2, 3, 2, 0.01)),
        ]:
            # this raises if any detector is non-deterministic (bookkeeping bug)
            dem = circ.detector_error_model(decompose_errors=True)
            print(f"{cell.name} | {name} | detectors={circ.num_detectors} "
                  f"obs={circ.num_observables} | DEM ok")
