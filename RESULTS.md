# Breathing codes — v0 falsification run

**Question:** can a quantum memory change its fine-graining level ("curvature")
at runtime — coarse for dense storage, fine for protection — without the morph
itself killing the benefit?

**Setup.** Phenomenological noise (X errors + measurement flips, rate p),
Z-check sector, stim + pymatching (MWPM), 50k shots/point. One full breath:
3 coarse rounds → inhale (reset new qubits, fold detectors per old face) →
fine hold → exhale (destructive Z-readout of new qubits, fold detectors) →
3 coarse rounds → readout.

**Codes** (both verified: commutation, Euler characteristic, k via GF(2) rank,
d via ILP):

| code | coarse | after 1 inhale (1→4 quad refinement) |
|---|---|---|
| torus (flat control) | [[18, 2, 3]] | [[72, 2, 6]] |
| Bring {5,5} (curved, genus 4) | [[30, 8, 3]] | [[120, 8, 6]] |

k preserved by topology in both; d exactly doubles. The curved code stores
8/30 = 0.27 logicals per physical qubit when exhaled — 2.4× the flat control.

## Kill-test 1: mid-morph distance dip ("hernia") — SURVIVED (phenomenological; REVISED at circuit level, see below)

Exact circuit fault distances (stim `shortest_graphlike_error`):

| circuit | torus | Bring |
|---|---|---|
| static coarse | 3 | 3 |
| static fine | 6 | 6 |
| inhale-only (ends fine) | 3 | 3 |
| full breath (2c+3f+2c) | 3 | 3 |
| full breath (1c+5f+1c) | 3 | 3 |

Every breathing variant sits exactly at min(d_coarse, d_fine). The morph
never opens a fault path shorter than the weakest static phase.

## Kill-test 2: breath economics — SURVIVED, with a measured price

Logical block error P_L (any of the k logicals wrong) vs total duration T.
Break-even hold time τ* (breathe beats stay-coarse):

| code, p | τ* (rounds) | P_L at T=48: coarse vs breathe |
|---|---|---|
| torus, p=0.01 | ≈ 16 | 0.105 vs 0.046 (2.3×) |
| torus, p=0.02 | ≈ 20 | 0.393 vs 0.292 (1.3×) |
| Bring, p=0.01 | ≈ 14 | 0.460 vs 0.181 (2.5×) |
| Bring, p=0.02 | ≈ 14 | 0.918 vs 0.651 (1.4×) |

Signature confirmed: the breathe curve is nearly **flat in T** — a fixed
"toll" (the two morphs + coarse bookends at d=3) followed by near-free
holding at d=6. Below τ*, don't breathe; above it, breathing strictly wins,
while occupying the 4× qubit count only during the hold window.

## Registered verdicts

1. Hernia: **not observed** (exact, not statistical). Fold detectors suffice.
2. Economics: **breath pays for itself** beyond τ* ≈ 14–20 rounds at
   p = 1–2%. The toll is dominated by time spent at d=3, not by the morph
   steps themselves.
3. Not yet tested: atom-move budget vs teleport-to-processor baseline
   (kill-criterion 3) — needs a movement-cost model; circuit-level noise;
   larger ℓ; spatially partial breathing (inhale only a subregion).

## Files

- `cellulation.py` — surfaces, quad refinement, k/d verification (ILP)
- `circuits.py` — static + breathing stim circuits (fold-detector morphs)
- `hernia.py` — kill-test 1
- `economics.py` — kill-test 2 (Monte Carlo)
- `results_torus.json`, `results_bring.json` — raw numbers
- `breathing_v0.png` — 4-panel summary figure

## Kill-test 3: moves vs the teleport-to-processor baseline — PARTIAL KILL

Exact hop counts from the code structures (coarse block 120 ancilla-hops/round,
fine 480; inhale/exhale each move 180 atoms; a d=6 planar patch = 71 atoms,
~0 transport/round). Surgery cost swept 0.5×–4× nominal.

| scenario | movement (atom-hops) | space (peak extra atoms) |
|---|---|---|
| protect 1 of 8 logicals | teleport wins from W ≥ 1–2 | teleport wins (71 vs 180) |
| protect all 8 (burst) | teleport wins from W ≥ 5–22 | **breathe wins 3.2×** (180 vs 568) |

Verdict: as a general single-qubit processor substitute, breathing loses the
movement budget — criterion 3 fires. It survives as (a) whole-block burst
armor under atom scarcity and (b) short windows. Move-induced error is
negligible at 1e-4/hop, significant at 1e-3/hop.

Rescue hypothesis (untested): the hop penalty exists because the fine code is
non-planar *everywhere at ℓ=2*. At larger ℓ the interior of every face is a
flat patch needing no transport — curvature (and hop cost) concentrates on the
seams. Partial breathing + large-ℓ layouts attack exactly this term.

**Scoreboard: 3 registered kill-criteria → 2 clean survivals, 1 partial kill.**
Honest shape of the result: breathing codes = dense storage with on-demand
burst armor, not a processor replacement.

## v1: the rescue hypotheses — BOTH CONFIRMED

**v1-A, deep breaths (seam locality).** Provenance-tracked refinement at
ℓ = 1, 2, 4, 8 with flat-interior layouts (interior checks = 0 transport):

| ℓ | d | naive hops/rd (Bring) | seam-only hops/rd | interior checks |
|---|---|---|---|---|
| 1 | 3 | 120 | 120 | 0% |
| 2 | 6 | 480 | 420 | 10.5% |
| 4 | 12 | 1920 | 1140 | 40.5% |
| 8 | 24 | 7680 | 2580 | 66.5% |

Shuttling grows like the perimeter (×2.3 per doubling), not the area (×4).
Hops per unit distance saturate instead of growing — deep breaths get
cheaper per unit of armor. KT3 rematch at matched d=12, protect-all-8:
vaults still hop-cheaper beyond W≈8 rounds, but breathing keeps a 2.6× space
win (900 vs 2296 atoms) and locality roughly doubled its favorable window.

**v1-B, partial breathing.** Mixed cellulations (inflate a face subset)
validate, preserve k, and give exact ILP per-logical distances:

- Torus, inflating 1/2/3 transversal columns: targeted logical d = 4/5/6
  while the other stays 3 (until full breath) — a graded **distance dial**,
  +1 per column at 42/78/108 extra atoms.
- Bring (curved), inflating only the 6 faces touching one weight-3 logical:
  **target 3→5, six bystanders 3→4 for free, one class dodges (stays 3).**
  Homology can route around local armor, but hyperbolic density makes
  dodging hard — local inflation leaks protection broadly.

Design rule extracted: guarantees require inflating a *transversal cut* of
the target class; on curved codes, local inflation buys wide collateral
protection but not universal guarantees.


## Circuit-level attack (v2 item #1) — executed; changes one headline

Full record: CIRCUIT_LEVEL.md. Compressed: under explicit CNOT extraction
with full depolarizing noise, a genuine morph-located hernia (fd 3→2, hook
faults on fine X-check ancillas, one touching a fresh B atom) fires in 5/8
otherwise-good gate schedules on the curved code and 0/10 on the flat one;
breathing-safe schedules exist (3/8, Bring seeds 8/15/19), are certified
exactly in seconds, and become a compile-time rule (verify schedules jointly
coarse+fine+morph). Static coarse Bring is itself schedule-fragile (fd=3 in
10/40). Two forced protocol facts: sequential X/Z extraction phases on
curved graphs; one gauge X-round per morph, both directions. Economics
improve: τ*≈10 rounds at p=1e-3, 3.7–3.8× block-error win at T=48; the
hernia's near-term rate cost is unmeasurable (low-p scaling threat only).

## Hardware (two-run record) — morph primitive validated on silicon
kingston: machinery ✓, economics ✗ (device above threshold — statics
invert). marrakesh: 4/4 ✓ incl. hardware τ*≤8 and monotone dial; run-1
anomaly resolved to transpiler scheduling (fixed by barriers+DD). Devices
bracket the theory's below-threshold precondition in both directions.
