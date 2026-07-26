"""Kill-test 1: does the morph create a fault path shorter than min(d)?"""
import stim
from cellulation import torus, bring_code, refine
from circuits import static_circuit, fine_static_circuit, breathing_circuit
from circuits import logical_Z_reps, Rec, _round, _face_qubits


def inhale_only_circuit(cell, R1, RF, p):
    """Coarse R1 -> inhale -> fine RF -> readout in FINE code.
    Isolates the inhale transition (no exhale)."""
    fine, info = refine(cell)
    c = stim.Circuit()
    rec = Rec()
    A = info["halfA"]
    spokes = sorted(set(info["spokes"].values()))
    coarse_data = [A[e] for e in range(cell.nE)]
    new_qubits = sorted([info["halfB"][e] for e in range(cell.nE)] + spokes)
    fine_data = list(range(fine.nE))
    coarse_fq = _face_qubits(cell, {e: A[e] for e in range(cell.nE)})
    fine_fq = _face_qubits(fine, list(range(fine.nE)))

    c.append("R", coarse_data)
    prev_c = None
    for r in range(R1):
        cur = _round(c, rec, coarse_data, coarse_fq, p)
        for i in range(cell.nF):
            t = [rec.rec(cur[i])] + ([rec.rec(prev_c[i])] if prev_c else [])
            c.append("DETECTOR", t)
        prev_c = cur
    c.append("R", new_qubits)
    cur_f = _round(c, rec, fine_data, fine_fq, p)
    for i in range(cell.nF):
        t = [rec.rec(cur_f[q]) for q in info["face_quads"][i]]
        t.append(rec.rec(prev_c[i]))
        c.append("DETECTOR", t)
    prev_f = cur_f
    for r in range(RF - 1):
        cur_f = _round(c, rec, fine_data, fine_fq, p)
        for i in range(fine.nF):
            c.append("DETECTOR", [rec.rec(cur_f[i]), rec.rec(prev_f[i])])
        prev_f = cur_f
    c.append("X_ERROR", fine_data, p)
    c.append("M", fine_data)
    mfin = rec.take(fine.nE)
    m = {e: mfin[e] for e in range(fine.nE)}
    for i, f in enumerate(fine.faces):
        c.append("DETECTOR", [rec.rec(m[e]) for e in f] + [rec.rec(prev_f[i])])
    coarse_logs = logical_Z_reps(cell)
    from cellulation import lift_cycle
    for li, log in enumerate(coarse_logs):
        lifted = lift_cycle(log, info)
        c.append("OBSERVABLE_INCLUDE", [rec.rec(m[e]) for e in lifted], li)
    return c


def fault_distance(circ):
    err = circ.shortest_graphlike_error()
    return len(err)


if __name__ == "__main__":
    p = 0.01
    for cell in [torus(3), bring_code()]:
        print(f"\n=== {cell.name} (d_coarse=3, d_fine=6) ===")
        rows = [
            ("static coarse (7 rounds)", static_circuit(cell, 7, p)),
            ("static fine   (7 rounds)", fine_static_circuit(cell, 7, p)),
            ("inhale-only  (2c + 5f, ends fine)", inhale_only_circuit(cell, 2, 5, p)),
            ("full breath  (2c + 3f + 2c)", breathing_circuit(cell, 2, 3, 2, p)),
            ("full breath  (1c + 5f + 1c)", breathing_circuit(cell, 1, 5, 1, p)),
        ]
        for name, circ in rows:
            print(f"  {name:36s} fault distance = {fault_distance(circ)}")
