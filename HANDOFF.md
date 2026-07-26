# BREATHING CODES — Project Handoff

**Status:** v0+v1+v1-C complete. Concept survived pre-registered falsification
(2 clean, 1 partial kill). Manuscript at v5 (9 pp) with the
exhale-to-operate result folded in.
**Owner:** Kamran (Stanford). **This document:** everything a collaborator
(or future you) needs to pick the project up cold.

---

## 0. Quick verdict

A quantum memory whose curvature/fine-graining level is a **runtime
variable** — dense curved phase for storage, flatter high-distance phase for
compute windows or noise bursts — works, with one mandatory ingredient found
at circuit level. Phenomenologically the morph is exactly safe; under real
gates a genuine morph-located hernia (fd 3→2, hook faults on fine X-ancillas)
fires in 5/8 otherwise-good schedules on the curved code (never on the flat
control) and is eliminated by **certified breathing-safe schedules** (3/8
exist; exact check, seconds — a compile-time rule: verify schedules jointly
across coarse+fine+morph). The breath has a measured price (τ* ≈ 14–20 rounds
phenomenological, ≈10 at circuit level p=1e-3; wins grow to 3.7–3.8× at
T=48), loses long single-qubit holds to the teleport baseline on transport,
but wins whole-block protection at ~1/3 the atom footprint — and two v1
mechanisms (seam-local deep breaths; partial breathing with a per-logical
distance dial) widen the niche. Novel phenomena found: the breath toll τ*,
the morph hernia + its certification rule, collateral armor on curved codes,
the homology dodger, the transversal-cut placement rule, two forced
protocol facts (sequential X/Z phases on curved graphs; one gauge X-round
per morph), and (v1-C) the exhale-to-operate dial: corridor deflation
provably lowers a targeted logical's operator weight ell*d -> d,
guaranteed by a min-max duality that the inflation direction lacks, with
collateral exposure as its price on curved codes.

---

## 1. Origin and intellectual context

- Grew out of the "computation as curvature" research program (curvature as
  a controllable computational resource), transplanted from neural-network
  substrates onto stabilizer codes.
- Trigger question: *has anyone used non-Euclidean geometry / curvature to
  dynamically change code distance?* Literature answer: curvature is a huge
  **static** design tool (hyperbolic codes, semi-hyperbolic dial), and
  distance is dynamically changed only in **flat** codes (Q3DE patch growth,
  Surf-Deformer, Floquet schedules). Dimension even has a runtime
  gauge-fixing switch. **Curvature never did.** That gap is this project.
- Methodology: registered-falsification style. Three kill-criteria stated
  before any simulation; experiments built to fire them.

## 2. The idea, precisely

One homological code block on a closed surface; qubits on edges; the
fine-graining level ℓ of the tiling is time-dependent.

- **Exhaled (coarse, curved):** high k/n, low d — storage mode.
- **Inhaled (fine ℓ, flatter):** n → Θ(ℓ²)n, k unchanged (topology),
  d → ≥ ℓ·d — armor mode.
- **Inhale protocol:** fetch new atoms (edge B-halves + face spokes) from a
  reservoir, init |0⟩, begin fine checks. First fine round is gauge data;
  one deterministic *fold detector* per old face:
  D_in(f) = m_coarse_last(f) ⊕ ⨁_{q∈quads(f)} m_fine_1(q).
- **Exhale protocol:** destructively Z-measure B-halves (+spokes; spokes
  cancel in all folds), resume coarse checks; fold detector:
  D_out(f) = m_coarse(f) ⊕ ⨁ m_fine_last(quads f) ⊕ ⨁_{e∈∂f} b_e.
- Logical Z̄ cycles lift to their subdivided images at inhale
  (deterministic; B qubits start |0⟩) and absorb {b_e} at exhale.
- Both directions are code deformation = gauge fixing in a parent subsystem
  code (Vuillot et al. 2019). The natural hardware host: reconfigurable
  neutral-atom arrays (mid-run rearrangement is routine there).

**Glossary:** breath = inhale+hold+exhale · hernia = hypothetical mid-morph
distance dip · τ* = break-even hold time · seam = check touching an
original edge/vertex (needs transport) vs interior (flat, transport-free) ·
dodger = logical class that routes around a partial inflation ·
transversal cut = inflated region every representative of a class must cross ·
corridor = both faces adjacent to each support edge of one representative
(deflate it to expose that representative at coarse weight) ·
collateral exposure = bystander weight drops caused by a corridor exhale
(mirror of collateral armor).

## 3. Registered kill-criteria and verdicts

| # | Criterion | Test | Verdict |
|---|---|---|---|
| 1 | Mid-morph distance dip ("hernia") | exact fault distance of full breathing circuits, phenomenological AND circuit level | **REVISED** — phenomenological: always = min(d) = 3 (exact). Circuit level: real hernia (fd 3→2, fine-X hooks) in 5/8 good-static schedules on Bring, 0/10 on torus; retired by certified breathing-safe schedules (3/8; exact, cheap). Rule: verify schedules jointly coarse+fine+morph. See CIRCUIT_LEVEL.md |
| 2 | Morph cost > hold benefit | Monte Carlo, 3 strategies vs duration T | **SURVIVED** — τ* ≈ 14–20 rounds at p=1–2%; 2.3–2.5× win at T=48; breathe curve flat in T (fixed toll, near-free holding) |
| 3 | More atom transport than teleport-to-patch | exact structural hop counts, assumptions swept | **PARTIAL KILL** — teleport wins single-logical long holds (W≥1–2); breathing wins whole-block footprint 3.2× (180 vs 568 atoms) and short bursts (W ≲ 5–22) |

## 4. Test systems (all properties machine-verified)

| code | coarse | fine (ℓ=2) | construction |
|---|---|---|---|
| torus L=3 (flat control) | [[18,2,3]] | [[72,2,6]] | standard square tiling |
| Bring {5,5} (curved, genus 4) | [[30,8,3]] | [[120,8,6]] | great-dodecahedron combinatorics: icosahedron vertices+edges; one pentagon per icosahedral vertex = the 5-cycle among its neighbors |

Verified per code: every edge in exactly 2 faces; H_X·H_Zᵀ = 0; χ and k by
GF(2) rank; ALL distances by ILP (HiGHS via scipy.milp). Fine-graining law
n′=ℓ²-ish·n, k′=k, d′≥ℓd met with equality at ℓ=2 (d exactly doubles).
Rate story: Bring stores 8/30 = 0.27 logicals/qubit exhaled — 2.4× the flat
control.

## 5. Full numerical results

> Circuit-level results (KT1-CL hernia + schedule scan, KT2-CL economics at
> p=1e-3: τ*≈10, 3.7–3.8× at T=48, two protocol facts) live in
> **CIRCUIT_LEVEL.md** — the authoritative record for v2 item #1.

### 5.1 KT1 — fault distances (exact)
static coarse 3 | static fine 6 | inhale-only 3 | breath(2c+3f+2c) 3 |
breath(1c+5f+1c) 3 — identical on torus and Bring. No hernia.

### 5.2 KT2 — economics (phenomenological p, 50k shots/point, MWPM)
Strategies over T rounds: coarse×T; breathe = 3c+(T−6)f+3c; fine×T.
Block error = any of k logicals wrong.

| code, p | τ* | P_L at T=48 (coarse → breathe → fine) |
|---|---|---|
| torus, 0.01 | ≈16 | 0.105 → 0.046 → 0.013 |
| torus, 0.02 | ≈20 | 0.393 → 0.292 → 0.182 |
| Bring, 0.01 | ≈14 | 0.460 → 0.181 → 0.019 |
| Bring, 0.02 | ≈14 | 0.918 → 0.651 → 0.297 |

Raw numbers: data/results_torus.json, data/results_bring.json.

### 5.3 KT3 — movement (exact structure, model assumptions swept)
Model: non-planar checks = ancilla tours w qubits (w hops)/round both
sectors; planar patch = 0 transport; reservoir place/retrieve = 1 hop/atom;
surgery merge = d rounds × c_s·2(w_log+d) hops, c_s ∈ {0.5,1,2,4}.
Numbers: coarse 120 hops/rd; fine 480; breath moves 180 atoms each way;
d=6 rotated patch = 71 atoms, 0 hops/rd.
Crossovers (teleport becomes hop-cheaper): protect-1: W≥1–2; protect-all-8
on-demand: W≥5/7/12/22 for c_s=0.5/1/2/4; parked: W≥2/4/9/19.
Footprint: 180 (any n) vs 71·n. Move-error: negligible at 1e-4/hop,
significant at 1e-3/hop.

### 5.4 v1-A — seam locality (provenance-tracked refinement)
Interior check = all its qubits interior to one original face → 0 transport.

| ℓ | d | naive hops/rd (Bring/torus) | seam-only | interior % |
|---|---|---|---|---|
| 1 | 3 | 120/72 | 120/72 | 0 |
| 2 | 6 | 480/288 | 420/252 | 10.5 |
| 4 | 12 | 1920/1152 | 1140/684 | 40.5 |
| 8 | 24 | 7680/4608 | 2580/1548 | 66.5 |

Seam hops ×≈2.3 per ℓ-doubling (perimeter), vs ×4 naive (area).
KT3 rematch at d=12 protect-all-8: vaults hop-cheaper from W≈8, breathing
keeps 2.6× footprint win (900 vs 2296 atoms).
(d at ℓ≥4 from construction bound d′≥ℓd + lifted-logical upper bound;
ILP-exact only at ℓ≤2 — flagged in paper limitations.)

### 5.5 v1-B — partial breathing (mixed cellulations, exact ILP distances)
Torus, inflate c transversal columns: **d_target = 3+c** (3→4→5→6),
other class stays 3 until full breath; cost 42/78/108 atoms,
156/228/288 hops/rd. k=2 preserved; all mixed cells validate.
Bring, inflate the 6 faces touching one weight-3 logical (+108 atoms,
half the surface): **target 3→5**, six bystanders 3→4 (collateral armor),
**one class dodges** (stays 3). Design rule: guaranteed per-logical gains
require inflating a transversal cut of the target class; on curved codes
local inflation buys wide but unguaranteed collateral protection.


### 5.6 v1-C — exhale-to-operate (partial deflation as an operator-weight dial)
Question inverted from v1-B: instead of inflating to RAISE a target's
distance, deflate (exhale) a sub-region of the FINE code to LOWER a
target's logical operator weight — because logical-measurement/surgery
cost scales with representative weight w (Theta(w) ancillas in
state-of-the-art qLDPC surgery). All numbers ILP-exact; lift certified by
commutation assertion before every optimization (v1_targeted_check
pattern). No stim required (structural layer only).

**Corridor rule:** deflate both faces adjacent to each support edge of one
coarse representative (6 faces on Bring, 2L on torus). The representative
survives unsubdivided => coset-min weight returns to coarse value.

| system | move | target op weight | bystanders / protection |
|---|---|---|---|
| torus fine [[72,2,6]] | 1 face-annulus | transverse X-bar 6->3 | parallel 6->5; c annuli: 6-c |
| torus | Z0's 2-sided corridor (6/9 faces) | Z0 6->3 exact | Z1 6->4; X-sector d [3,4]; 66/108 atoms back |
| Bring fine [[120,8,6]] | 6-face corridor of weight-3 rep | target Z 6->3 exact | X-sector [3,3,3,3,4,4,4,4]; conjugate X-bar 6->4; 72/180 atoms back; hold 480->336 hops/rd |

**Min-max duality (the theory point):** distance/weight is a MIN over
representatives. Inflation helps only by blocking ALL representatives
(transversal cut; dodgers exist — v1-B). Deflation helps by exposing ONE.
=> targeted weight reduction is GUARANTEED; no dodger is possible in this
direction. Cleanest theorem-shaped statement in the project.

**Collateral exposure:** the mirror of collateral armor. Bring corridor
(half the surface) pulls 4 of 8 classes to weight 3, rest to 4. A corridor
opened for one cheap op weakens+cheapens many bystanders — batching bonus
AND protection hazard; corridor placement is a scheduling problem.

**Economics (KT3 model, c_s=1):** fine merge = 144 hops, w=6 boundary, 6
rounds. Corridor-exhale->operate->re-inhale = 36 + 2*72 = 180 hops, w=3,
3 rounds, transversal readout 3 vs 6 qubits. Amortizes at M>=2 ops per
corridor window (108 hops/op at M=2, 72 at M=4); collateral exposure means
one corridor cheapens ops for half the block. Protection price: operated
class runs at d = its own weight during the window (same bargain as the
breath itself — keep windows short).

**Novelty check (2026-07-25, web, 4 targeted searches):** distinct from
(a) weight-reduction / qLDPC surgery (Hastings; Cohen et al. Sci Adv 2022;
gauge-fixed QLDPC surgery arXiv:2407.18393) — those reduce CHECK/ancilla
weight around a static code and take representative weight w as given
(their Theta(w) cost is exactly what our dial feeds); (b) dynamic
deformation (Q3DE, Surf-Deformer, defect-adaptive codes) — deform to
RESTORE distance, never deliberately lower operator weight for cheap ops;
(c) code switching / dimensional jump — changes the available GATE SET,
not operator weight (closest conceptual rhyme; cited); (d) Fu-Gottesman
dynamical distance — analysis framework only; (e) logical spectroscopy
(arXiv:2607.05386, the uploaded tarball) — STATIC addressable low-weight
bases; complementary (tells you where operators are, we change how heavy
they are at runtime). Surface-code patch shrink/grow is spiritually
adjacent but planar, single-logical, and not framed as a per-logical
weight tool; our claim is scoped to targeted per-logical reduction on
shared-surface multi-logical blocks + corridor rule + duality + collateral
exposure. Caveat: single-session review; independent pass recommended.

**Files:** code/exhale_ops.py (dial + economics), code/exhale_ops2.py
(coset ILP for the operator's own weight), figures/breathing_exhale.png.
**Untested:** everything circuit-level. Corridor morphs inherit the hernia
risk; joint schedule certification required before trusting fault
distance (flagged in paper limitations).

## 6. Methods notes (for whoever re-runs or extends)

- Noise: phenomenological, X-error/Z-check sector only (CSS symmetry);
  X_ERROR(p) on live data each round; check outcomes flipped w.p. p
  (MPP prob arg); B-half exhale readout flipped w.p. p; final data
  readout ideal (identical across strategies → fair).
- Simulation: stim 1.16 circuits with manual measurement-record
  bookkeeping (circuits.py:Rec). stim's DEM construction **certifies every
  fold detector deterministic** — this was the main correctness gate for
  the morph bookkeeping. Every error mechanism triggers ≤2 detectors
  (graphlike) → PyMatching MWPM on decompose_errors=True DEM;
  fault distances via Circuit.shortest_graphlike_error().
- ILP distances: min-weight coset member per logical; mod-2 handled via
  integer slack (Hx = 2y; L·x = 2z+1). Same method as the
  hyperbolic-codes literature.
- Reproducibility: `pip install -r requirements.txt && ./run_all.sh`
  (~10–15 min laptop). economics.py is the long step. MC seeds fixed
  (rng seed 7) but stim sampler seeds drawn from it — expect
  statistical-level variation only.
- One dev-time bug worth knowing: an inline targeted-lift script once
  lifted ALL edges instead of the logical's support; caught by the
  commutation assertion. The saved v1_targeted_check.py is correct and
  asserts the lift commutes before optimizing.

## 7. File map

```
code/
  cellulation.py      surfaces (torus, Bring), quad refinement, GF2 rank,
                      ILP distances, logical extraction   [foundation]
  circuits.py         static + breathing stim circuits, fold detectors
  hernia.py           KT1 + inhale-only variant (exact fault distances)
  economics.py        KT2 Monte Carlo (3 strategies × T × p)
  kt3_moves.py        KT3 hop/footprint accounting + sensitivity sweep
  v1_deep.py          provenance-tracked ℓ=1..8 refinement, seam
                      classification, KT3 rematch at d=12
  v1_partial.py       partial refinement, per-logical ILP distances,
                      torus column dial + Bring control
  v1_targeted_check.py  lifts the specific targeted logical → d=3→5
  exhale_ops.py       v1-C: deflation dials, per-logical ILP, economics
  exhale_ops2.py      v1-C: coset ILP (operator's own min weight)
  requirements.txt / run_all.sh
data/    results_torus.json, results_bring.json  (raw KT2 numbers)
figures/ breathing_v0.png (KT2), kt3_moves.png (KT3), breathing_v1.png (v1)
paper/   breathing_codes.pdf  (compiled, 6 pp)
RESULTS.md            running lab-notebook version of Secs. 3–5
CIRCUIT_LEVEL.md      circuit-level attack: full science record (v2 item #1)
  code/circuit_level.py, code/run_circuit_level.py,
  data/results_circuit_level.json, figures/breathing_cl.png
hardware/             IBM notebook (clean + executed) + breathing_hw.py
overleaf_package.zip  upload-ready Overleaf project (see §9)
```

## 8. Limitations / threats to validity (be upfront with reviewers)

1. ~~Phenomenological noise only~~ — DONE (CIRCUIT_LEVEL.md): explicit CNOT
   extraction with full depolarizing noise incl. morph rounds. Remaining:
   uniform depolarizing only — no leakage, no atom loss (now the TOP open
   threat: fold detectors assume every outcome exists), no morph-transport
   noise; 40 random schedules = search, not optimization.
2. Smallest possible codes (d=3 base); one inhale level simulated live.
3. ℓ≥4 distances by construction bound, not ILP.
4. KT3 rests on a stated transport model (perfect patch locality; surgery
   cost swept 8× but not derived from a real AOD schedule).
5. All code written and run in a single session; independent review
   recommended before submission. Exact results (KT1 table, all ILP
   distances, hop counts) are properties of the circuits/codes and are the
   most trustworthy layer.

## 9. Paper + Overleaf

- paper/breathing_codes.pdf: 6 pages, two-column, 3 figures, 22 refs
  (all anchored to sources verified during the literature phase).
- Red placeholders to fill: author surname; acknowledgments
  (advisor/funding).
- overleaf_package.zip: upload directly at overleaf.com → New Project →
  Upload Project. Main file breathing_codes.tex, compiler pdfLaTeX
  (default). Figures included at project root. Two compile passes resolve
  references (Overleaf does this automatically).
- Suggested pre-submission steps: (a) supervised pass on Sec. 6 movement
  assumptions, (b) run the v2 items below that a referee will request.

## 10. v2 roadmap (prioritized)

1. ~~Circuit-level noise~~ — **DONE** (see CIRCUIT_LEVEL.md; hernia found,
   autopsied, retired by schedule certification; economics improved).
1b. **Atom loss during morphs** — new top threat: fold detectors assume
   every outcome exists; a lost atom returns none. No decoder story yet.
1c. **Constructive schedule optimization** — find breathing-safe orders by
   design, not 3-in-8 luck; a hard constraint for the curvature scheduler.
2. **Live partial-breath dynamics**: v0 stim machinery + v1_partial
   cellulations → simulate inflating one stripe and watch the targeted
   logical's error drop. (Machinery exists; ~1 day.) Now also the
   deflation direction: live corridor-exhale → surgery → re-inhale, with
   joint schedule certification (v1-C is structural only; hernia risk
   untested there).
3. **Transversal-cut theory**: which face subsets are simultaneous
   transversal cuts for the most classes per atom? Systolic optimization;
   the collateral-armor effect suggests hyperbolic surfaces are unusually
   favorable.
4. **Breathing bicycle codes**: fine-graining-like weight-preserving morphs
   of bivariate bicycle codes → connects to IBM-roadmap hardware.
5. **Real AOD movement model** for KT3 (replace the swept assumptions).
6. **Curvature scheduler**: compiler pass allocating ℓ(region, t) under an
   atom budget given an algorithm timeline + noise map; τ*, the paired
   placement rules (inflate a transversal cut to protect / exhale a
   corridor to operate), and collateral-exposure management are its first
   constraints. v1-C supplies the gate-cost term.
7. **Fuse corridor morph with surgery**: both touch the same boundary;
   can the exhale and the merge share rounds/moves?

## 11. Key literature anchors

Static curvature: Freedman–Meyer–Luo (Z2-systolic freedom); Breuckmann &
Terhal QST 2, 035007 (2017) [semi-hyperbolic dial]; Delfosse ISIT'13
kd²≤C(log k)²n; Guth–Lubotzky JMP 55, 082202 (2014); golden codes
arXiv:1712.08578. Good qLDPC: Panteleev–Kalachev STOC'22;
Leverrier–Zémor FOCS'22; Bravyi et al. Nature 627, 778 (2024) [gross code].
Expanders are negatively curved: Salez arXiv:2101.08242. Dynamic-but-flat:
Q3DE (MICRO'22); Surf-Deformer arXiv:2405.06941; Fu–Gottesman
arXiv:2403.04163 [dynamical distance]. Switching machinery: Bombín NJP 17,
083002 (2015) & NJP 18, 043038 (2016) [dimensional jump]; Vuillot et al.
NJP 21, 033028 (2019) [deformation = gauge fixing]; Vasmer–Kubica PRXQ 3,
030319 (2022); Higgott–Breuckmann PRX 11, 031039 (2021) [subsystem
semi-hyperbolic]. Baseline architecture: Xu et al. Nat. Phys. 20, 1084
(2024) [zoned qLDPC atom arrays]; Bluvstein et al. Nature 626, 58 (2024).
Tools: Gidney (Stim) Quantum 5, 497 (2021); Higgott–Gidney (PyMatching)
arXiv:2303.15933. Fine-graining formulas cross-checked against
arXiv:2602.10423.

## 12. IBM hardware notebook (hardware/)

`breathing_codes_ibm.ipynb` runs the **1D morph primitive** on real devices
— a repetition code that inhales d=3→6 (CNOT extension + fold detectors),
holds, and exhales (tail measured into the decoding frame). Heavy-hex and
the campaign's measured decoherence walls rule out the 2D curved codes;
this validates the protocol layer and says so in its receipt.

Four pre-registered checks, all PASS in the device-like Aer validation
(executed copy ships alongside with outputs embedded):
E1 fold-detector rate ≤2.5× bulk (measured 0.98) · E2 breathe beats coarse
at longest hold, z≥3 (6.9; τ*≤8 rounds) · E3 a GHZ logical superposition
survives a full breath, ⟨X̄⟩ 0.481 vs static-6's 0.503 · E4 partial dial
monotone, z≥3 (6.3).

Operational: single top control panel (`RUN_ON_HARDWARE`, token/CRN,
`IBM_BACKEND=None` → least-busy device with ≥11 qubits, pinning supported),
three-rung pip ladder, pass-manager self-heal, Open-Plan job mode (no
Session), QPU-window guard, pre-submit routed-2Q assert (<300 vs the
measured walls), noiseless correctness gate before anything submits.
~15 circuits × 8192 shots ≈ 1–2 min QPU. Workflow: flip the switch, run
top-to-bottom, paste `out_breathing/receipt.json` back for the bar-check
before the paper's hardware subsection is written.

### §12 status update — two-run hardware record (COMPLETE)
Run 1 ibm_kingston: 2/4 (E1, E3 pass; E2/E4 fail, localized to
above-pseudo-threshold operation — statics invert z=+3.8); exposed the
scheduling anomaly (breathe detectors 1.7–5.5× static-fine) + two owned
instrumentation flaws (relative E1 bar; no raw capture). Run 2
ibm_marrakesh (barriers+DD+twirl, raw npz, W1 control, absolute E1):
**4/4** — fold 1.28, toll z=5.9 with hardware τ*≤8, dial z=4.0, coherent
transport 0.429(10); anomaly resolved to scheduling (1.1–1.2×); one
pathological ancilla (~0.45 rate) dominates absolute rates on all arms
equally — future nicety: add readout error to best_chain cost. E3 nuance
reported honestly: breath free on kingston, ~33% coherence cost on
marrakesh (double exposure to the bad ancilla; unresolved). Cross-device
punchline: economics pass exactly where statics permit, fail exactly
where they invert. Hardware section folded into paper (Sec. hw, 8 pp).

### Paper v4 (final form)
Restructured to standard layout (contributions + Results subsections; the
pre-specified-analysis history lives in this handoff, not the paper),
field-standard terminology (no informal coinages beyond the named
protocol), em dashes removed, author line: Kamran Ansari. Figure
kt3_moves.png regenerated with a neutral title. 8 pp, clean-room
verified; overleaf_package.zip current.
Title finalized: 'Breathing codes: quantum LDPC memories that morph their curvature at runtime' (LDPC added to title + first two abstract mentions).

### Paper v5 (exhale-to-operate integrated)
New Sec. 'Exhale to operate' (results, after partial refinement) + figure
breathing_exhale.png; abstract, contributions (now vii items), discussion
(4th regime; scheduler gate-cost term; paired placement rules), and
limitations (corridor morphs structural-only, hernia caveat) updated.
Three refs added: Cohen et al. Sci Adv 2022, gauge-fixed QLDPC surgery
arXiv:2407.18393, logical spectroscopy arXiv:2607.05386. 9 pp, clean
compile, overleaf_package updated. Novelty check recorded in Sec 5.6.
