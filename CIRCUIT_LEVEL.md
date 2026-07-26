# CIRCUIT-LEVEL ATTACK — v2 item #1 executed (July 2026)

**Question attacked:** does the breathing protocol survive real gates —
explicit CNOT extraction with depolarizing noise on every gate, idle, reset,
and measurement, including during the morph rounds?
**Code:** `code/circuit_level.py` · **Data:** `data/results_circuit_level.json`
· **Figure:** `figures/breathing_cl.png`

## Verdict in one line
The v0 "no hernia, proven exactly" claim is **falsified as stated** — generic
gate schedules on the curved code produce a real, morph-located distance dip
(3→2) — but the protocol **survives with one new mandatory ingredient**:
breathing-safe schedules exist, are found by cheap search, and are certifiable
exactly at compile time. Economics improve at circuit level.

## Model
Sequential two-phase rounds (all Z-extraction, then all X-extraction);
CX data→ancZ / ancX→data; DEPOLARIZE2(p) per CX, X/Z_ERROR(p) at
reset/measure, DEPOLARIZE1(p) per idle tick; gate order = greedy
edge-coloring of a fixed per-round permutation (`order_seed`). Both check
sectors live, so ancilla hook errors exist. stim's DEM determinism check
certifies every fold detector under real gates (the correctness gate; it
caught both findings below).

## Two protocol facts discovered by the determinism gate
1. **Interleaved X/Z schedules are invalid on arbitrary curved graphs**
   unless CNOT orders are parity-matched per adjacent check pair (the
   surface-code "Z-order" trick has no known {5,5} analog). Sequential
   phases restore determinism generically.
2. **The X sector must run gauge (one blind round) across every morph, both
   directions.** Old-vertex X-stars are NOT preserved operators: each split
   edge's surviving atom sits at one endpoint, so fine stars contain fresh
   B atoms. Invisible phenomenologically; unavoidable at circuit level.

## KT1-CL — fault distances (exact, graphlike; independent search corroborates)
Default schedule: torus c/f/breath = 3/6/3 (clean). Bring = 2/6/2 — the
**static** coarse code loses a unit to hooks on its weight-5 X-checks
(autopsy: two X-ancilla faults; in the breath they sit in the plain coarse
rounds). Schedule scan (fixed per-round order, 40 seeds): coarse fd=3 in
10/40; fine wobbles 4–6.

**The hernia:** under schedules where BOTH statics hold (coarse 3, fine 5–6),
`breath 2/3/2` drops to **fd = 2 in 5/8 cases** — mechanisms are hooks on
**fine X-ancillas during the breath** (one touches a fresh B atom), harmless
in the static fine circuit, lethal inside the morph structure. This is the
predicted Tier-1 failure mode, real at circuit level.

**The rescue:** 3/8 good-static schedules (seeds 8, 15, 19) give
**breath fd = 3 = min(statics)** — hernia-free compilations exist. The flat
torus never herniates (10/10 schedules). New design rule: **schedules must be
verified jointly (coarse + fine + morph)**; per-phase-good ≠ breath-good.
The verification is exact and takes seconds.

## KT2-CL — economics at p = 1e-3, 20k shots (torus default; Bring seed 8)
| T | torus c / b / f | Bring c / b / f |
|---|---|---|
| 8  | .016 / .019 / .001 | .095 / .111 / .002 |
| 12 | .024 / .022 / .001 | .140 / .119 / .003 |
| 16 | .031 / .021 / .001 | .184 / .114 / .005 |
| 24 | .046 / .021 / .002 | .263 / .116 / .006 |
| 48 | .089 / .024 / .006 | .456 / .121 / .013 |

τ* ≈ **10 rounds** both codes (was 14–20 phenomenologically); T=48 win
**3.7× (torus) / 3.8× (Bring)** (was 2.3–2.5×); breathe curve stays flat
(fixed toll, near-free holding — signature intact). **Hernia price at
p=1e-3, T=48: unmeasurable** (0.1216 vs 0.1213 at 20k shots) — the fd-2
mechanism is a two-fault coincidence; it caps low-p scaling, not near-term
rates. It still must be compiled away for the distance story to hold.

## Limitations
Graphlike fd is primary (DEM decomposes cleanly; independent
`search_for_undetectable_logical_errors` corroborates key circuits at
weight 2 / ≥3). 40 random schedules ≠ schedule optimization. Uniform
depolarizing only; no leakage, loss, or morph-transport noise. Smallest
codes, one live inhale level. Single-session code; the exact fd layer is the
most trustworthy.

## What this changes in the paper
Replace "no hernia (exact)" with: hernia exists generically at circuit level
on curved codes, is localized to fine-phase hooks, and is removed by
certified breathing-safe schedules (exact check). Add the two protocol facts
above. Economics section strengthens.
