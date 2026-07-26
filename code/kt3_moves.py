"""Kill-test 3: is a breath cheaper than shipping qubits to a processor?

Task: give logical qubit(s) stored in the coarse Bring block distance-6
protection for a window of W rounds, then return to dense storage.

Protocol A (breathe): inflate the whole block coarse->fine, hold W, exhale.
Protocol B (teleport): keep block coarse; lattice-surgery-teleport the
  target logical(s) into dedicated d=6 planar surface patch(es), hold W
  there, teleport back (the zoned-architecture baseline).

Cost model (stated, swept):
- Non-planar codes need transport every syndrome round: each check ancilla
  physically tours its w data qubits => "atom-hops" per round = sum of all
  check weights (both X and Z sectors). Parallel AOD steps per round ~ 2*w_max.
- Planar surface patches: ancillas interleaved in the lattice, local Rydberg
  gates, 0 transport per round.
- Placement/retrieval from reservoir: 1 hop per atom each way.
- Surgery (one merge: d rounds of joint checks): bridge checks touch the
  patch boundary (local) and the LDPC logical support (non-local).
  hops per merge = d_t * m_bridge, with m_bridge = c_s * 2*(w_log + d_t).
  Teleport = 2 merges (in + out counted per direction below).
- Common costs (coarse block's own 120 hops/round) excluded: identical in
  both protocols.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from cellulation import bring_code, refine

cell = bring_code()
fine, info = refine(cell)

# ---- structural counts (exact, from the objects) ----
inc = lambda c: int(c.HZ().sum() + c.HX().sum())      # ancilla hops / round
HOPS_COARSE = inc(cell)          # common baseline (excluded from marginals)
HOPS_FINE = inc(fine)
ANC_COARSE = cell.nF + cell.nV
ANC_FINE = fine.nF + fine.nV
NEW_DATA = fine.nE - cell.nE
NEW_ANC = ANC_FINE - ANC_COARSE
PLACED = NEW_DATA + NEW_ANC      # atoms fetched at inhale, returned at exhale

D_T = 6                          # target distance
W_LOG = 3                        # coarse logical weight (=d_coarse)
PATCH_ATOMS = D_T * D_T + (D_T * D_T - 1)   # rotated surface code, data+anc

print(f"Bring block: coarse hops/round={HOPS_COARSE}, fine hops/round={HOPS_FINE}")
print(f"Inhale places {NEW_DATA} data + {NEW_ANC} ancilla = {PLACED} atoms")
print(f"d=6 patch: {PATCH_ATOMS} atoms parked, 0 hops/round")


def breathe_hops(W):
    """Marginal hops vs staying coarse."""
    return 2 * PLACED + W * (HOPS_FINE - HOPS_COARSE)


def teleport_hops(W, n_log, c_s=1.0, park=False):
    """n_log logicals shipped to n_log patches. park=True: patches
    pre-exist (no placement); else placed/retrieved on demand."""
    merge = D_T * c_s * 2 * (W_LOG + D_T)
    place = 0 if park else 2 * n_log * PATCH_ATOMS
    return place + n_log * 2 * merge  # hold = 0 hops


def peak_extra_atoms_A():
    return PLACED


def peak_extra_atoms_B(n_log):
    return n_log * PATCH_ATOMS


if __name__ == "__main__":
    Ws = np.arange(1, 121)
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))

    # -- panel 1: hops, protect ONE logical
    ax = axes[0]
    ax.plot(Ws, [breathe_hops(w) for w in Ws], "tab:blue", lw=2,
            label="breathe whole block")
    for c_s, ls in [(0.5, ":"), (1.0, "-"), (2.0, "--"), (4.0, "-.")]:
        ax.plot(Ws, [teleport_hops(w, 1, c_s) for w in Ws], "tab:red",
                ls=ls, lw=1.4, label=f"teleport x1 (surgery x{c_s})")
    ax.set_title("Protect 1 of 8 logicals — atom-hops")
    ax.set_xlabel("hold window W (rounds)")
    ax.set_ylabel("marginal atom-hops")
    ax.legend(fontsize=7.5)
    ax.grid(alpha=0.3)

    # -- panel 2: hops, protect ALL 8 (noise burst / whole-block compute)
    ax = axes[1]
    ax.plot(Ws, [breathe_hops(w) for w in Ws], "tab:blue", lw=2,
            label="breathe whole block")
    for park, ls, lab in [(False, "-", "teleport x8 (on-demand patches)"),
                          (True, "--", "teleport x8 (patches parked)")]:
        ax.plot(Ws, [teleport_hops(w, 8, 1.0, park) for w in Ws], "tab:red",
                ls=ls, lw=1.6, label=lab)
    ax.set_title("Protect all 8 logicals — atom-hops")
    ax.set_xlabel("hold window W (rounds)")
    ax.legend(fontsize=7.5)
    ax.grid(alpha=0.3)

    # -- panel 3: peak extra atom footprint
    ax = axes[2]
    ns = np.arange(1, 9)
    ax.bar(ns - 0.18, [peak_extra_atoms_A()] * 8, 0.36, color="tab:blue",
           label="breathe (whole block, any n)")
    ax.bar(ns + 0.18, [peak_extra_atoms_B(n) for n in ns], 0.36,
           color="tab:red", label="teleport (n patches)")
    ax.set_title("Peak extra atoms needed")
    ax.set_xlabel("logicals protected simultaneously")
    ax.set_ylabel("extra atoms")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y")

    fig.suptitle("Kill-test 3: breathe-in-place vs teleport-to-processor "
                 "(Bring [[30,8,3]] block, d=6 protection)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig("kt3_moves.png", dpi=160)

    # crossovers
    print("\nCrossovers (teleport becomes hop-cheaper than breathing):")
    for n_log, park in [(1, False), (8, False), (8, True)]:
        for c_s in [0.5, 1.0, 2.0, 4.0]:
            b = np.array([breathe_hops(w) for w in Ws])
            t = np.array([teleport_hops(w, n_log, c_s, park) for w in Ws])
            idx = np.where(t < b)[0]
            w0 = Ws[idx[0]] if len(idx) else None
            print(f"  n={n_log} park={park} surgery x{c_s}: "
                  f"teleport wins from W >= {w0}")
    print(f"\nPeak extra atoms: breathe={peak_extra_atoms_A()} (any n); "
          f"teleport: {[peak_extra_atoms_B(n) for n in [1,2,4,8]]} for n=1,2,4,8")
    # move-induced error sanity: extra effective error per fine round
    for eps in [1e-4, 1e-3]:
        per_anc = 4 * eps  # ~4 hops per ancilla per round
        print(f"move error eps={eps}: adds ~{per_anc:.1e} per check per round "
              f"(vs p=1e-2 measurement error -> {'negligible' if per_anc < 2e-3 else 'significant'})")
