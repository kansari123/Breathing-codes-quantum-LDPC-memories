"""Kill-test 2: is a breath worth its cost, and after what hold time?

Strategies over T total noisy rounds:
  coarse : static coarse code, T rounds          (cheap, weak)
  breathe: 3 coarse + (T-6) fine + 3 coarse      (morph in, hold, morph out)
  fine   : static fine code, T rounds            (reference ceiling, 4x qubits)
"""
import numpy as np
import stim
import pymatching
import json
from cellulation import torus, bring_code, refine
from circuits import static_circuit, fine_static_circuit, breathing_circuit

SHOTS = 50_000
RNG = np.random.default_rng(7)


def logical_block_error(circ, shots=SHOTS):
    dem = circ.detector_error_model(decompose_errors=True)
    matching = pymatching.Matching.from_detector_error_model(dem)
    sampler = circ.compile_detector_sampler(seed=int(RNG.integers(2**31)))
    dets, obs = sampler.sample(shots, separate_observables=True)
    pred = matching.decode_batch(dets)
    fails = np.any(pred != obs, axis=1)
    P = fails.mean()
    se = np.sqrt(P * (1 - P) / shots)
    return P, se


def run(cell, ps, Ts):
    out = {}
    nf = refine(cell)[0].nE
    for p in ps:
        for T in Ts:
            row = {}
            row["coarse"] = logical_block_error(static_circuit(cell, T, p))
            row["fine"] = logical_block_error(fine_static_circuit(cell, T, p))
            if T >= 8:
                row["breathe"] = logical_block_error(
                    breathing_circuit(cell, 3, T - 6, 3, p))
            out[(p, T)] = row
            msg = f"{cell.name} p={p} T={T}: " + "  ".join(
                f"{k}={v[0]:.4f}±{v[1]:.4f}" for k, v in row.items())
            print(msg, flush=True)
    return out


if __name__ == "__main__":
    ps = [0.01, 0.02]
    Ts = [8, 12, 16, 24, 32, 48]
    results = {}
    for cell in [torus(3), bring_code()]:
        results[cell.name] = run(cell, ps, Ts)
    ser = {name: {f"{p}|{T}": {k: list(v) for k, v in row.items()}
                  for (p, T), row in res.items()}
           for name, res in results.items()}
    with open("results.json", "w") as f:
        json.dump(ser, f, indent=1)
    print("saved results.json")
