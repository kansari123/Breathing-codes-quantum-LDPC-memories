"""Reproduce the ENTIRE circuit-level attack (CIRCUIT_LEVEL.md) in one run.

  pip install -r requirements.txt && python run_circuit_level.py

Sections (total ~3-5 min laptop):
  1. determinism gate      - stim certifies every fold detector under real gates
  2. fault distances       - default schedule, both codes (KT1-CL headline)
  3. schedule scan + hernia table (Bring) + torus control table
  4. autopsies             - what the weight-2 mechanisms physically are
  5. corroboration         - independent (non-graphlike) search
  6. economics Monte Carlo - KT2-CL, writes ../data/results_circuit_level.json
                             and ../figures/breathing_cl.png
Seeds fixed; MC numbers vary at statistical level only.
"""
import json
import os
import time

import numpy as np
import stim
import pymatching

from cellulation import torus, bring_code, refine
from circuit_level import (static_cl, fine_static_cl, breathing_cl,
                           fault_distance, verify)

_BASE = os.path.dirname(os.path.abspath(__file__))  # outputs resolve relative to this file, not the CWD

P = 1e-3
SHOTS = 20000
ECON_TS = [8, 12, 16, 24, 48]
BRING_SAFE_SEED = 8          # breathing-safe schedule found in section 3
GOOD_COARSE_SEEDS = [3, 4, 8, 12, 15, 18, 19, 24]


def banner(s):
    print(f"\n{'=' * 70}\n{s}\n{'=' * 70}", flush=True)


# ------------------------------------------------------------ 1. determinism
banner("1. DETERMINISM GATE (raises if any fold detector is wrong)")
for cell in [torus(3), bring_code()]:
    print(f"-- {cell.name}")
    verify("static coarse", static_cl(cell, 4, P))
    verify("static fine  ", fine_static_cl(cell, 4, P))
    verify("breath 2/3/2 ", breathing_cl(cell, 2, 3, 2, P))
    verify("breath 1/1/1 ", breathing_cl(cell, 1, 1, 1, P))

# ------------------------------------------------- 2. fd, default schedule
banner("2. FAULT DISTANCES, default schedule (KT1-CL headline)")
for cell in [torus(3), bring_code()]:
    print(f"-- {cell.name}")
    for T in (3, 7):
        fault_distance(f"static coarse T={T}", static_cl(cell, T, P))
        fault_distance(f"static fine   T={T}", fine_static_cl(cell, T, P))
    fault_distance("breath 1/1/1      ", breathing_cl(cell, 1, 1, 1, P))
    fault_distance("breath 2/3/2      ", breathing_cl(cell, 2, 3, 2, P))

# ------------------------------------------------- 3. schedules + hernia
banner("3. SCHEDULE SCAN (fixed per-round orders)")
cell = bring_code()
res = {}
for seed in range(40):
    fd = len(static_cl(cell, 3, P, order_seed=seed).shortest_graphlike_error())
    res.setdefault(fd, []).append(seed)
print("Bring static coarse T=3, 40 seeds, fd->count:",
      {k: len(v) for k, v in sorted(res.items())})
print("fd=3 seeds:", res.get(3, []))

print("\nBring hernia table (T=7 statics vs breath 2/3/2):")
print(" seed | static_c | static_f | breath | hernia?")
hern = 0
for s in GOOD_COARSE_SEEDS:
    fc = len(static_cl(cell, 7, P, order_seed=s).shortest_graphlike_error())
    ff = len(fine_static_cl(cell, 7, P, order_seed=s).shortest_graphlike_error())
    fb = len(breathing_cl(cell, 2, 3, 2, P, order_seed=s).shortest_graphlike_error())
    h = fb < min(fc, ff)
    hern += h
    print(f"  {s:3d}  |    {fc}     |    {ff}     |   {fb}    |"
          f" {'YES' if h else 'no'}")
print(f"-> hernia in {hern}/{len(GOOD_COARSE_SEEDS)} good-coarse schedules; "
      f"safe schedules exist (e.g. seed {BRING_SAFE_SEED})")

print("\nTorus control (10 seeds):")
tor = torus(3)
for s in range(10):
    fc = len(static_cl(tor, 7, P, order_seed=s).shortest_graphlike_error())
    ff = len(fine_static_cl(tor, 7, P, order_seed=s).shortest_graphlike_error())
    fb = len(breathing_cl(tor, 2, 3, 2, P, order_seed=s).shortest_graphlike_error())
    print(f"  seed {s}: c={fc} f={ff} breath={fb}"
          f" {'HERNIA' if fb < min(fc, ff) else ''}")

# --------------------------------------------------------------- 4. autopsy
banner("4. AUTOPSIES (what the weight-2 mechanisms are)")
fine, info = refine(cell)
B_set = {info["halfB"][e] for e in range(cell.nE)}
A_set = {info["halfA"][e] for e in range(cell.nE)}


def qname(q):
    if q < 120:
        return f"data{q}({'A' if q in A_set else 'B' if q in B_set else 'spoke'})"
    if q < 132:
        return f"cZanc{q - 120}"
    if q < 144:
        return f"cXanc{q - 132}"
    if q < 204:
        return f"fZanc{q - 144}"
    return f"fXanc{q - 204}"


def autopsy(name, circ):
    err = circ.shortest_graphlike_error(canonicalize_circuit_errors=True)
    print(f"-- {name}: weight {len(err)}")
    for i, e in enumerate(err):
        loc = e.circuit_error_locations[0]
        qs = [t.gate_target.qubit_value for t in loc.flipped_pauli_product]
        inst = loc.instruction_targets
        tr = [t.gate_target.qubit_value for t in inst.targets_in_range]
        print(f"   mech {i}: {inst.gate} fault, flipped "
              f"{[qname(q) for q in qs if q is not None]} "
              f"(gate on {[qname(q) for q in tr if q is not None]})")


autopsy("Bring static coarse, default sched (coarse-phase hooks)",
        static_cl(cell, 3, P))
autopsy("Bring breath 2/3/2, default sched (same coarse hooks)",
        breathing_cl(cell, 2, 3, 2, P))
autopsy("Bring breath 2/3/2, seed 3 (THE HERNIA: fine-phase hooks)",
        breathing_cl(cell, 2, 3, 2, P, order_seed=3))

# ---------------------------------------------------------- 5. corroboration
banner("5. INDEPENDENT SEARCH CORROBORATION")
for nm, c in [("Bring breath seed 3 (hernia)",
               breathing_cl(cell, 2, 3, 2, P, order_seed=3)),
              ("Bring static coarse default",
               static_cl(cell, 3, P))]:
    e2 = c.search_for_undetectable_logical_errors(
        dont_explore_detection_event_sets_with_size_above=5,
        dont_explore_edges_with_degree_above=5,
        dont_explore_edges_increasing_symptom_degree=True)
    print(f"  {nm}: search finds weight {len(e2)}")

# -------------------------------------------------------------- 6. economics
banner("6. ECONOMICS MONTE CARLO (KT2-CL)")


def run_mc(circ, seed):
    dem = circ.detector_error_model(decompose_errors=True)
    m = pymatching.Matching.from_detector_error_model(dem)
    dets, obs = circ.compile_detector_sampler(seed=seed).sample(
        SHOTS, separate_observables=True)
    pred = m.decode_batch(dets)
    return float(np.any(pred != obs, axis=1).mean())


results = {}
t0 = time.time()
for cc, sched in [(torus(3), None), (bring_code(), BRING_SAFE_SEED)]:
    tag = cc.name
    results[tag] = {"schedule_seed": sched, "p": P, "shots": SHOTS, "data": {}}
    for T in ECON_TS:
        row = {
            "coarse": run_mc(static_cl(cc, T, P, order_seed=sched), 11 + T),
            "breathe": run_mc(breathing_cl(cc, 3, T - 6, 3, P,
                                           order_seed=sched), 12 + T),
            "fine": run_mc(fine_static_cl(cc, T, P, order_seed=sched), 13 + T),
        }
        results[tag]["data"][T] = row
        print(f"  {tag} T={T:2d}: coarse={row['coarse']:.4f} "
              f"breathe={row['breathe']:.4f} fine={row['fine']:.4f} "
              f"[{time.time() - t0:.0f}s]", flush=True)

results["hernia_price_T48"] = {
    "breathe_seed3_hernia": run_mc(
        breathing_cl(bring_code(), 3, 42, 3, P, order_seed=3), 99),
    "breathe_seed8_safe": results["bring_30_8"]["data"][48]["breathe"],
}
print("  hernia price @T=48:", results["hernia_price_T48"])

json.dump(results, open(os.path.join(_BASE, "..", "data", "results_circuit_level.json"), "w"), indent=1)
print("  wrote ../data/results_circuit_level.json")

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.6))
    for ax, tag, title in [
            (axes[0], "torus_L3", "torus [[18,2,3]] (flat control)"),
            (axes[1], "bring_30_8",
             "Bring {5,5} [[30,8,3]] (curved, safe schedule)")]:
        D = results[tag]["data"]
        Ts = sorted(D)
        for strat, mk, lbl in [("coarse", "o-", "coarse (never breathe)"),
                               ("breathe", "s-", "breathe 3c+(T-6)f+3c"),
                               ("fine", "^-", "fine (always inflated)")]:
            ax.semilogy(Ts, [D[t][strat] for t in Ts], mk, label=lbl, ms=4)
        ax.set_title(title, fontsize=9.5)
        ax.set_xlabel("rounds T")
        ax.grid(alpha=.3, which="both")
    axes[0].set_ylabel("block logical error (any of k wrong)")
    axes[0].legend(fontsize=7.5)
    fig.suptitle("Circuit-level noise, p=1e-3", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(_BASE, "..", "figures", "breathing_cl.png"), dpi=160)
    print("  wrote ../figures/breathing_cl.png")
except ImportError:
    print("  (matplotlib not installed; skipped figure)")

banner("DONE — compare against CIRCUIT_LEVEL.md")
