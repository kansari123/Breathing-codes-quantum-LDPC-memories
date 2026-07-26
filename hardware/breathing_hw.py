"""breathing_hw.py -- the breathing-code morph primitive on IBM hardware (1D).

A repetition code (bit-flip sector) that BREATHES: start at d=3, inhale to
d=6 mid-run (CNOT extension + new checks with fold detectors), hold, exhale
(Z-measure the tail, fold outcomes into the decoding frame), finish at d=3.

Four pre-registered experiments (PASS bars in the receipt):
  E1 FOLD  -- morph-event detector fire rates consistent with bulk rounds
              (<= 2.5x bulk median): the gauge-fixing bookkeeping survives
              real correlated noise.
  E2 TOLL  -- breathe vs stay-coarse vs stay-fine memory curves; a breath
              has a measurable toll and a break-even hold time tau*.
              PASS if P_L(breathe) < P_L(coarse) at the longest hold, 3 sigma.
  E3 COHER -- a logical superposition (GHZ = |+bar>) survives a full breath:
              prepare GHZ-3, inhale to GHZ-6, one check round, exhale with
              X-basis tail readout folded into the frame, read Xbar.
              PASS if <Xbar>_breathed > 0 at 3 sigma.
  E4 DIAL  -- partial breathing: extend by c in {0,1,2,3} qubits mid-run;
              P_L strictly decreases with c. PASS if P_L(3) < P_L(0), 3 sigma.

Scope (honest): 1D has no curvature -- this validates the PROTOCOL layer of
breathing codes (morph safety, fold detectors, toll, coherent transport,
graded dial) on real devices; the 2D curved-code claims stay in simulation.
Feasibility is set by measured device constants from this hardware campaign:
retention ~0.994/routed-2Q, decoherence walls ~348 gates @10q (marrakesh);
every circuit here stays under ~130 routed 2Q on <= 11 qubits in a line.
"""
import json, os, time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister

# ----------------------------------------------------------------- schedule
# Layout (fine): d0 a0 d1 a1 d2 a2 d3 a3 d4 a4 d5   (11 qubits in a line)
# Coarse phase uses d0..d2 with a0,a1. Inhale extends by EXT qubits.

class Sched:
    """Builds one breathing (or static) circuit AND its decoding structure
    from the same event list, so circuit and decoder cannot disagree."""

    def __init__(self, n_core=3, ext=3, tag=""):
        self.n_core, self.ext, self.tag = n_core, ext, tag
        self.n_data = n_core + ext
        self.qd = QuantumRegister(self.n_data, "d")
        self.qa = QuantumRegister(max(1, self.n_data - 1), "a")
        self.events = []          # for receipts
        self.meas = []            # (kind, idx, round) per classical bit
        self.detectors = []       # list of (list of meas indices, label)
        self.edges = []           # (u, v_or_None, meas_flip?, qubit, label)
        self.obs_meas = []        # meas indices whose XOR = observable
        self.qc = QuantumCircuit(self.qd, self.qa)
        self.last = {}            # ancilla -> last meas index
        self._round = 0
        self._nm = 0

    # -- helpers
    def _measure(self, q, kind, idx):
        self.qc.measure(q, self._creg_bit())
        self.meas.append((kind, idx, self._round))
        self._nm += 1
        return self._nm - 1

    def _creg_bit(self):
        if not hasattr(self, "cr"):
            self.cr = ClassicalRegister(300, "m")   # trimmed at finalize
            self.qc.add_register(self.cr)
        return self.cr[self._nm]

    # -- physical events
    def prep_zero(self, n):
        pass  # |0...0> is the default

    def prep_ghz(self, n):
        self.qc.h(self.qd[0])
        for i in range(n - 1):
            self.qc.cx(self.qd[i], self.qd[i + 1])

    def extend(self, n_from, n_to):
        """Inhale: CNOT-copy the boundary into fresh |0> qubits."""
        self.qc.barrier()
        for i in range(n_from, n_to):
            self.qc.cx(self.qd[i - 1], self.qd[i])
        self.qc.barrier()
        self.events.append(("extend", n_from, n_to, self._round))

    def synd_round(self, n_active, morph=False):
        """One check round on ancillas 0..n_active-2. Detector for each
        ancilla vs its previous outcome; first-ever outcome is its own
        detector (deterministic 0 from |0..0> / GHZ / post-extension)."""
        A = range(n_active - 1)
        self.qc.barrier()
        for a in A:
            self.qc.cx(self.qd[a], self.qa[a])
        for a in A:
            self.qc.cx(self.qd[a + 1], self.qa[a])
        for a in A:
            m = self._measure(self.qa[a], "synd", a)
            lbl = ("morph" if (morph or (a in self.first_time(a))) and morph
                   else "bulk")
            if a in self.last:
                self.detectors.append(([m, self.last[a]],
                                       "morph" if morph else "bulk"))
            else:
                self.detectors.append(([m], "morph" if morph else "bulk"))
            self.last[a] = m
        for a in A:
            self.qc.reset(self.qa[a])
        self._round += 1

    def first_time(self, a):
        return []  # (labelling handled inline above)

    def retire_tail_z(self, n_keep):
        """Exhale (memory-Z): Z-measure data n_keep..end; their outcomes are
        terminal readouts for the tail -- each retired ancilla's last synd
        outcome folds with the measured parities (fold detector), and the
        boundary ancilla (n_keep-1) folds with the first retired qubit."""
        outs = {}
        self.qc.barrier()
        for i in range(n_keep, self.n_used):
            outs[i] = self._measure(self.qd[i], "dataZ", i)
        for a in range(n_keep - 1, self.n_used - 1):
            group = [self.last.pop(a)]
            if a == n_keep - 1:                # boundary ancilla: core side
                group += [outs[a + 1]]
                self._pending_boundary = (a, group)  # closes at final readout
                continue
            group += [outs[a], outs[a + 1]]
            self.detectors.append((group, "morph"))
        self.n_used = n_keep
        self.events.append(("retire", n_keep, self._round))

    def final_z(self):
        outs = {}
        self.qc.barrier()
        for i in range(self.n_used):
            outs[i] = self._measure(self.qd[i], "dataZ", i)
        for a in range(self.n_used - 1):
            if a in self.last:
                self.detectors.append(([self.last.pop(a), outs[a],
                                        outs[a + 1]], "final"))
        if getattr(self, "_pending_boundary", None):
            a, group = self._pending_boundary
            self.detectors.append((group + [outs[a]], "morph"))
            self._pending_boundary = None
        self.obs_meas.append(outs[0])

    def finalize(self):
        used = self._nm
        qc = QuantumCircuit(self.qd, self.qa, ClassicalRegister(used, "m"))
        for inst in self.qc.data:
            if inst.operation.name == "measure":
                qc.measure(inst.qubits[0], qc.cregs[0][self.qc.cregs[0].index(
                    inst.clbits[0])])
            else:
                qc.append(inst.operation, inst.qubits, inst.clbits)
        qc.metadata = {"tag": self.tag}
        return qc


def breathe_circuit(R1, W, R2, ext=3, tag=""):
    """|0bar> memory: R1 coarse rounds -> inhale(ext) -> W fine rounds
    (first is the morph round) -> exhale -> R2 coarse rounds -> readout."""
    s = Sched(ext=ext, tag=tag)
    s.n_used = s.n_core
    for _ in range(R1):
        s.synd_round(s.n_used)
    if ext > 0 and W > 0:
        s.extend(s.n_core, s.n_core + ext)
        s.n_used = s.n_core + ext
        s.synd_round(s.n_used, morph=True)
        for _ in range(W - 1):
            s.synd_round(s.n_used)
        s.retire_tail_z(s.n_core)
        s.synd_round(s.n_used, morph=True)
    for _ in range(R2 - 1):
        s.synd_round(s.n_used)
    s.final_z()
    return s


def static_circuit(n, T, tag=""):
    s = Sched(ext=n - 3 if n > 3 else 0, tag=tag)
    s.n_core = n
    s.n_used = n
    for _ in range(T):
        s.synd_round(n)
    s.final_z()
    return s


def ghz_breath(ext=3, tag="ghz_breathe"):
    """|+bar> through a full breath; returns (circuit, parse_fn)."""
    s = Sched(ext=ext, tag=tag)
    s.prep_ghz(3)
    s.n_used = 3
    s.synd_round(3)                       # ZZ checks stabilize GHZ
    s.extend(3, 3 + ext)
    s.n_used = 3 + ext
    s.synd_round(s.n_used, morph=True)
    x_tail = []
    for i in range(3, s.n_used):          # exhale in X basis
        s.qc.h(s.qd[i])
        x_tail.append(s._measure(s.qd[i], "dataX", i))
    x_core = []
    for i in range(3):                    # Xbar readout
        s.qc.h(s.qd[i])
        x_core.append(s._measure(s.qd[i], "dataX", i))
    s.x_bits = x_tail + x_core
    return s


def ghz_static(n, tag):
    s = Sched(ext=max(0, n - 3), tag=tag)
    s.prep_ghz(n)
    s.n_used = n
    s.synd_round(n)
    s.x_bits = []
    for i in range(n):
        s.qc.h(s.qd[i])
        s.x_bits.append(s._measure(s.qd[i], "dataX", i))
    return s


# ------------------------------------------------------------------ decode
def build_matcher(sched, p_data=0.01, p_meas=0.02):
    """Spacetime matching graph from the schedule's detector structure.
    Nodes = detectors; time edges = shared synd meas; space edges = shared
    data qubit at one round (incl. boundaries); observable frame = edges
    touching data qubit 0's final readout / left boundary."""
    import pymatching
    m = pymatching.Matching()
    det_of_meas = {}
    for di, (group, _) in enumerate(sched.detectors):
        for g in group:
            det_of_meas.setdefault(g, []).append(di)
    n_det = len(sched.detectors)
    w_meas = float(np.log((1 - p_meas) / p_meas))
    w_data = float(np.log((1 - p_data) / p_data))
    added = set()

    def link(dets, w, frame):
        dets = sorted(set(dets))
        key = (tuple(dets), frame)
        if key in added:
            return
        added.add(key)
        fid = {0} if frame else set()
        if len(dets) == 1:
            m.add_boundary_edge(dets[0], weight=w, fault_ids=fid,
                                merge_strategy="smallest-weight")
        elif len(dets) == 2:
            m.add_edge(dets[0], dets[1], weight=w, fault_ids=fid,
                       merge_strategy="smallest-weight")

    # measurement errors: every measurement appearing in exactly 2 detectors
    for mi, dets in det_of_meas.items():
        kind = sched.meas[mi][0]
        if len(dets) == 2:
            link(dets, w_meas if kind == "synd" else w_data,
                 frame=(kind == "dataZ" and sched.meas[mi][1] == 0))
        elif len(dets) == 1:
            link(dets, w_meas if kind == "synd" else w_data,
                 frame=(kind == "dataZ" and sched.meas[mi][1] == 0))
    # space errors: data qubit i between ancillas i-1, i at each round:
    # detectors indexed per round via the synd measurements
    synd_by_round = {}
    for mi, (kind, a, r) in enumerate(sched.meas):
        if kind == "synd":
            synd_by_round.setdefault(r, {})[a] = det_of_meas[mi][-1]
    for r, amap in synd_by_round.items():
        anc = sorted(amap)
        nq = max(anc) + 2
        for i in range(nq):
            dets = [amap[a] for a in (i - 1, i) if a in amap]
            link(dets, w_data, frame=(i == 0))
    m.set_boundary_nodes(set())  # boundaries already explicit
    return m, n_det


def shots_to_records(sched, bitstrings):
    arr = np.array([[int(b) for b in s[::-1]] for s in bitstrings],
                   dtype=np.uint8)
    dets = np.zeros((arr.shape[0], len(sched.detectors)), dtype=np.uint8)
    for di, (group, _) in enumerate(sched.detectors):
        for g in group:
            dets[:, di] ^= arr[:, g]
    obs = np.zeros(arr.shape[0], dtype=np.uint8)
    for g in sched.obs_meas:
        obs ^= arr[:, g]
    return arr, dets, obs


def logical_error(sched, bitstrings, p_data=0.01, p_meas=0.02):
    matcher, _ = build_matcher(sched, p_data, p_meas)
    arr, dets, obs = shots_to_records(sched, bitstrings)
    pred = matcher.decode_batch(dets).astype(np.uint8)[:, 0]
    fails = (pred ^ obs).mean()
    se = float(np.sqrt(max(fails * (1 - fails), 1e-12) / len(obs)))
    labels = [lab for _, lab in sched.detectors]
    rates = {lab: float(dets[:, [i for i, l in enumerate(labels) if l == lab]]
                        .mean()) if lab in labels else 0.0
             for lab in ("bulk", "morph", "final")}
    return float(fails), se, rates


def ghz_xbar(sched, bitstrings):
    arr = np.array([[int(b) for b in s[::-1]] for s in bitstrings],
                   dtype=np.uint8)
    par = np.zeros(arr.shape[0], dtype=np.uint8)
    for g in sched.x_bits:
        par ^= arr[:, g]
    v = 1.0 - 2.0 * par.mean()
    return float(v), float(2 * np.std(par) / np.sqrt(len(par)))


# --------------------------------------------------------------- run modes
def _aer_backend(noisy=True):
    from qiskit_aer import AerSimulator
    if not noisy:
        return AerSimulator()
    from qiskit_aer.noise import NoiseModel, depolarizing_error, ReadoutError
    nm = NoiseModel()
    nm.add_all_qubit_quantum_error(depolarizing_error(0.012, 2),
                                   ["cx", "cz", "ecr"])
    nm.add_all_qubit_quantum_error(depolarizing_error(3e-4, 1),
                                   ["x", "sx", "h", "rz", "reset"])
    nm.add_all_qubit_readout_error(ReadoutError([[0.97, 0.03],
                                                 [0.05, 0.95]]))
    return AerSimulator(noise_model=nm)


def _pm(backend, initial_layout=None):
    """Pass-manager with the version-skew self-heal chain. Layout is bound
    HERE, at construction -- StagedPassManager.run() takes no layout kwarg."""
    from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
    kw = dict(optimization_level=1, initial_layout=initial_layout)
    try:
        return generate_preset_pass_manager(backend=backend, **kw)
    except Exception:
        try:
            return generate_preset_pass_manager(
                target=backend.target, translation_method="translator", **kw)
        except Exception:
            return generate_preset_pass_manager(
                basis_gates=backend.configuration().basis_gates,
                coupling_map=backend.configuration().coupling_map, **kw)


def best_chain(backend, L):
    """Greedy best line of L qubits by summed 2Q error (job-mode friendly)."""
    try:
        t = backend.target
        err = {}
        for g in ("cz", "ecr", "cx"):
            if g in t.operation_names:
                for q, props in t[g].items():
                    if props and props.error is not None:
                        err[tuple(sorted(q))] = min(
                            err.get(tuple(sorted(q)), 1.0), props.error)
        adj = {}
        for (a, b), e in err.items():
            adj.setdefault(a, []).append((b, e))
            adj.setdefault(b, []).append((a, e))
        best, best_cost = None, None
        for start in adj:
            path, cost, seen = [start], 0.0, {start}
            while len(path) < L:
                nxt = [(e, q) for q, e in adj.get(path[-1], [])
                       if q not in seen]
                if not nxt:
                    break
                e, q = min(nxt)
                path.append(q); seen.add(q); cost += e
            if len(path) == L and (best_cost is None or cost < best_cost):
                best, best_cost = path, cost
        return best
    except Exception:
        return None



def _interleaved_layout(chain, n_qubits):
    """Physical layout d0 a0 d1 a1 ... along the chain. Circuit register
    order is d0..d_{nd-1}, a0..a_{nd-2} with n_qubits = 2*nd - 1, so data
    sits on even chain slots and ancillas between their two data qubits."""
    nd = (n_qubits + 1) // 2
    if not chain or len(chain) < 2 * nd - 1:
        return None
    return ([chain[2 * i] for i in range(nd)] +
            [chain[2 * j + 1] for j in range(nd - 1)])


def _compile_for(backend, circuits):
    """The exact hardware compile path (shared with the smoke test): best
    chain by calibrated 2Q error, interleaved layout per circuit width,
    per-width pass managers, no-layout fallback on any failure."""
    chain = best_chain(backend, 11)
    pms, isa = {}, []
    for qc in circuits:
        lay = _interleaved_layout(chain, qc.num_qubits)
        key = None if lay is None else tuple(lay)
        if key not in pms:
            try:
                pms[key] = _pm(backend, initial_layout=lay)
            except Exception:
                pms[key] = _pm(backend)
        try:
            isa.append(pms[key].run(qc))
        except Exception:
            if None not in pms:
                pms[None] = _pm(backend)
            isa.append(pms[None].run(qc))
    n2q = [sum(1 for i in c.data if i.operation.num_qubits == 2)
           for c in isa]
    return isa, n2q, chain


def _build_jobs(R1=2, R2=2):
    jobs = []   # (sched, kind, meta)
    for T in (4, 8, 12, 16):
        jobs.append((static_circuit(3, T, f"coarse_T{T}"), "mem",
                     dict(strategy="coarse", T=T)))
    for T in (8, 16):
        jobs.append((static_circuit(6, T, f"fine_T{T}"), "mem",
                     dict(strategy="fine", T=T)))
    for W in (2, 4, 8, 12):
        jobs.append((breathe_circuit(R1, W, R2, tag=f"breathe_W{W}"), "mem",
                     dict(strategy="breathe", T=R1 + W + R2, W=W)))
    for c in (1, 2):
        jobs.append((breathe_circuit(R1, 8, R2, ext=c, tag=f"dial_c{c}"),
                     "mem", dict(strategy="dial", c=c, T=R1 + 8 + R2)))
    jobs.append((breathe_circuit(R1, 1, R2, tag="breathe_W1"), "mem",
                 dict(strategy="nullbreath", T=R1 + 1 + R2, W=1)))
    jobs.append((ghz_static(3, "ghz3"), "ghz", dict(which="static3")))
    jobs.append((ghz_static(6, "ghz6"), "ghz", dict(which="static6")))
    jobs.append((ghz_breath(tag="ghz_breathe"), "ghz", dict(which="breathed")))
    return jobs


def run_all(mode="aer", shots=8192, backend_name=None,
            outdir="out_breathing", noisy=True, seed=7):
    os.makedirs(outdir, exist_ok=True)
    jobs = _build_jobs()

    circuits = [s.finalize() for s, _, _ in jobs]
    n2q = []
    if mode == "ibm":
        from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2
        service = QiskitRuntimeService()
        if backend_name:
            backend = service.backend(backend_name)
        else:
            backend = service.least_busy(operational=True, simulator=False,
                                         min_num_qubits=11)
        print("backend:", backend.name,
              "(pinned)" if backend_name else "(least busy, >=11 qubits)")
        isa, n2q, chain = _compile_for(backend, circuits)
        if chain:
            print("qubit chain (interleaved d/a):", chain)
        est = len(isa) * shots * 90e-6 + 8 * len(isa)
        print(f"pre-submit: {len(isa)} circuits x {shots} shots, "
              f"routed-2Q per circuit {min(n2q)}-{max(n2q)}, "
              f"projected ~{est:.0f} s QPU")
        assert max(n2q) < 300, "circuit exceeds measured decoherence wall"
        sampler = SamplerV2(mode=backend)      # job mode: Open-Plan legal
        try:
            sampler.options.dynamical_decoupling.enable = True
            sampler.options.twirling.enable_measure = True
        except Exception:
            pass
        job = sampler.run(isa, shots=shots)
        print("job id:", job.job_id())
        res = job.result()
        all_bits = [r.data.m.get_bitstrings() for r in res]
    else:
        backend = _aer_backend(noisy=noisy)
        pm = _pm(backend)
        isa = [pm.run(qc) for qc in circuits]
        n2q = [sum(1 for i in c.data if i.operation.num_qubits == 2)
               for c in isa]
        job = backend.run(isa, shots=shots, memory=True, seed_simulator=seed)
        result = job.result()
        all_bits = [result.get_memory(i) for i in range(len(isa))]

    # ---------------- analysis (+ raw detector capture for autopsies)
    mem, ghz, raw = {}, {}, {}
    for (s, kind, meta), bits in zip(jobs, all_bits):
        if kind == "mem":
            P, se, rates = logical_error(s, bits)
            mem[s.tag] = dict(meta, P=P, se=se, rates=rates)
            _, dets, obs = shots_to_records(s, bits)
            metaarr = []
            for group, lab in s.detectors:
                rr, aa = -1, -1
                for g in group:
                    kk, ai, ri = s.meas[g]
                    if kk == "synd":
                        rr, aa = ri, ai
                        break
                metaarr.append([rr, aa, {"bulk": 0, "morph": 1,
                                         "final": 2}[lab]])
            raw[f"{s.tag}__dets"] = dets
            raw[f"{s.tag}__obs"] = obs
            raw[f"{s.tag}__meta"] = np.array(metaarr, dtype=np.int16)
        else:
            v, se = ghz_xbar(s, bits)
            ghz[meta["which"]] = dict(x=v, se=se)
    np.savez_compressed(f"{outdir}/raw_detectors.npz", **raw)

    # E1 fold consistency: morph rate vs STATIC-FINE bulk (absolute bar --
    # run-1 lesson: a breathe-internal baseline can be inflated and make
    # the relative bar vacuous)
    fine_bulk = [m["rates"]["bulk"] for m in mem.values()
                 if m["strategy"] == "fine"]
    morph = [m["rates"]["morph"] for m in mem.values()
             if m["strategy"] in ("breathe", "nullbreath")]
    e1_ratio = float(np.median(morph) / max(np.mean(fine_bulk), 1e-9))
    e1 = e1_ratio <= 2.5
    # E2 toll / tau*
    cW = {m["T"]: (m["P"], m["se"]) for m in mem.values()
          if m["strategy"] == "coarse"}
    bW = {m["T"]: (m["P"], m["se"]) for m in mem.values()
          if m["strategy"] == "breathe"}
    Tmax = max(bW)
    z2 = (cW[Tmax][0] - bW[Tmax][0]) / max(
        np.hypot(cW[Tmax][1], bW[Tmax][1]), 1e-9)
    e2 = z2 >= 3.0
    tau = None
    for T in sorted(bW):
        if T in cW and bW[T][0] < cW[T][0]:
            tau = T
            break
    # E3 coherent breath
    e3 = ghz["breathed"]["x"] > 3 * ghz["breathed"]["se"]
    # E4 dial
    dial = {0: cW[12], 3: bW[12]}
    for m in mem.values():
        if m["strategy"] == "dial":
            dial[m["c"]] = (m["P"], m["se"])
    z4 = (dial[0][0] - dial[3][0]) / max(
        np.hypot(dial[0][1], dial[3][1]), 1e-9)
    e4 = z4 >= 3.0

    receipt = dict(mode=mode, shots=shots,
                   backend=backend.name if mode == "ibm"
                   else f"aer(noisy={noisy})",
                   routed_2q=dict(min=int(min(n2q)), max=int(max(n2q))),
                   mem=mem, ghz=ghz, tau_star=tau,
                   checks=dict(E1_fold=dict(ratio=e1_ratio, passed=bool(e1)),
                               E2_toll=dict(z=float(z2), passed=bool(e2)),
                               E3_coherent=dict(passed=bool(e3)),
                               E4_dial=dict(z=float(z4), passed=bool(e4))))
    with open(f"{outdir}/receipt.json", "w") as f:
        json.dump(receipt, f, indent=1)

    # ---------------- figure
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.8))
    ax = axes[0]
    for strat, col in (("coarse", "tab:red"), ("breathe", "tab:blue"),
                       ("fine", "tab:gray")):
        pts = sorted((m["T"], m["P"], m["se"]) for m in mem.values()
                     if m["strategy"] == strat)
        if pts:
            T, P, S = zip(*pts)
            ax.errorbar(T, P, yerr=S, marker="o", color=col, label=strat)
    if tau:
        ax.axvline(tau, ls="--", c="k", alpha=0.5)
        ax.text(tau, ax.get_ylim()[0], f" tau*<= {tau}", fontsize=8)
    ax.set_xlabel("total rounds T"); ax.set_ylabel("logical error P_L")
    ax.set_title("E2: breath toll and tau*"); ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax = axes[1]
    cs = sorted(dial)
    ax.errorbar(cs, [dial[c][0] for c in cs], yerr=[dial[c][1] for c in cs],
                marker="s", color="tab:blue")
    ax.set_xlabel("inserted qubits c (d = 3+c)"); ax.set_ylabel("P_L")
    ax.set_title("E4: partial-breath dial"); ax.grid(alpha=0.3)
    ax = axes[2]
    names = ["static3", "static6", "breathed"]
    ax.bar(names, [ghz[n]["x"] for n in names],
           yerr=[ghz[n]["se"] for n in names],
           color=["tab:gray", "tab:gray", "tab:blue"])
    ax.set_ylabel(r"$\langle \bar X\rangle$")
    ax.set_title("E3: coherent breath (GHZ)"); ax.grid(alpha=0.3, axis="y")
    fig.suptitle(f"Breathing morph primitive -- {receipt['backend']}, "
                 f"{shots} shots", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fp = f"{outdir}/breathing_hw.png"
    fig.savefig(fp, dpi=150)
    return receipt, fp


def correctness_gate():
    """Noiseless Aer: zero detector rate, zero logical error, GHZ <Xbar>=1."""
    r, _ = run_all(mode="aer", shots=1024, noisy=False,
                   outdir="out_gate")
    for tag, m in r["mem"].items():
        assert m["P"] == 0.0, (tag, m)
        assert all(v == 0.0 for v in m["rates"].values()), (tag, m)
    for k, g in r["ghz"].items():
        assert abs(g["x"] - 1.0) < 1e-9, (k, g)
    print("[OK] correctness gate: all detectors deterministic, P_L = 0, "
          "<Xbar> = 1 noiselessly")


def hardware_compile_smoke():
    """Exercise the EXACT hardware compile branch (chain pick, interleaved
    layout, per-width PMs, routing) on a fake heavy-hex backend -- zero QPU.
    This is the test that would have caught the initial_layout kwarg bug."""
    try:
        from qiskit_ibm_runtime.fake_provider import FakeSherbrooke
        backend = FakeSherbrooke()
    except Exception as e:
        print("compile smoke skipped (no fake_provider):", type(e).__name__)
        return
    circuits = [s.finalize() for s, _, _ in _build_jobs()]
    isa, n2q, chain = _compile_for(backend, circuits)
    assert len(isa) == len(circuits)
    assert max(n2q) < 300, f"routed-2Q {max(n2q)} exceeds the wall budget"
    print(f"[OK] hardware compile smoke (FakeSherbrooke): {len(isa)} circuits"
          f" transpiled, routed-2Q {min(n2q)}-{max(n2q)}, chain head "
          f"{chain[:6] if chain else None}")


def autopsy(outdir="out_breathing"):
    """Zero-QPU post-mortem from raw_detectors.npz: detector rates by
    phase label, by round, and by ancilla; hottest cells; comparison of
    every breathe circuit against the static-fine baseline."""
    z = np.load(f"{outdir}/raw_detectors.npz")
    tags = sorted({k.split("__")[0] for k in z.files})
    ref = None
    if "fine_T16__dets" in z.files:
        ref = float(z["fine_T16__dets"].mean())
        print(f"reference: static fine_T16 mean detector rate = {ref:.4f}")
    fig, ax = plt.subplots(figsize=(8, 4))
    hot = []
    for tag in tags:
        dets, meta = z[f"{tag}__dets"], z[f"{tag}__meta"]
        rate = dets.mean(0)
        lab = {0: "bulk", 1: "morph", 2: "final"}
        by_lab = {l: float(rate[meta[:, 2] == c].mean())
                  for c, l in lab.items() if (meta[:, 2] == c).any()}
        line = f"{tag:14s} " + "  ".join(f"{l}={v:.4f}"
                                         for l, v in by_lab.items())
        if ref:
            line += f"   (x{by_lab.get('bulk', 0)/max(ref,1e-9):.1f} fine-bulk)"
        print(line)
        for di in np.argsort(rate)[::-1][:3]:
            hot.append((float(rate[di]), tag, int(meta[di, 0]),
                        int(meta[di, 1]), lab[int(meta[di, 2])]))
        if tag.startswith(("breathe", "fine_T16")):
            rounds = sorted(set(meta[meta[:, 0] >= 0, 0]))
            rr = [float(rate[meta[:, 0] == r].mean()) for r in rounds]
            ax.plot(rounds, rr, marker="o", ms=3, label=tag)
    print("\nhottest detectors (rate, tag, round, ancilla, label):")
    for h in sorted(hot, reverse=True)[:8]:
        print(f"  {h[0]:.4f}  {h[1]:14s} round={h[2]:3d} anc={h[3]:2d} {h[4]}")
    ax.set_xlabel("round"); ax.set_ylabel("mean detector rate")
    ax.set_title("per-round detector rates (breathe vs static fine)")
    ax.legend(fontsize=7); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(f"{outdir}/autopsy.png", dpi=140)
    print(f"figure: {outdir}/autopsy.png")
