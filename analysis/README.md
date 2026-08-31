# Low-compute extension of the DNA-origami scissor nanoswitch

A zero-new-MD analysis chain that extends the M.S. thesis
*"In silico rational design of DNA origami-based allosteric metamaterials"*
(Campagna & Carnevale, 2025) using only its existing simulation outputs. Every
tool here runs on a laptop in seconds — no LAMMPS, no GPU, no new trajectories.

## The device

The nanoswitch is a flat **"X" (scissors)**: two rigid DNA-origami beams
(~8,700 nucleotides each, ~100 nm long) crossing at a central **crux** bridged
by a flexible **poly-T linker** that acts as the hinge. The thesis varied the
linker length (10/20/30/40 nt) and measured the **standard deviation σ of the
hinge bond angle**, finding the counter-intuitive result that *longer linkers
reduce* angular fluctuation (σ: 6.24° → 2.27°, exponential decay).

This chain reinterprets that result as mechanics, measures the scissor geometry
directly, and turns the hinge angle into a quantitative biosensor readout —
addressing the two gaps flagged at the thesis defense: (1) put a *functional
element* on the hinge to make the cross a real biosensor, and (2) characterize
the *global open/close* behavior, not just the local angle jitter.

## Pipeline

```
  thesis σ table ──▶ stiffness_from_sigma.py ─┐
                                              ├─▶ binding_to_signal.py ──▶ sensor design specs
  snapshot .xyz  ──▶ hinge_geometry.py ──▶ optical_readout.py ─┘
```

| Tool | Input | What it computes | Outputs |
|------|-------|------------------|---------|
| `stiffness_from_sigma.py` | thesis Table 1 (σ) | hinge torsional stiffness `k_θ = k_BT/σ²`, exponential fit, actuation-energy table | `stiffness_table.csv`, `stiffness_from_sigma.png` |
| `hinge_geometry.py` | a snapshot `.xyz` | plane (PCA) → EM two-line **beam split** → hinge angle, crux, beam length, **tip–tip reporter distances** | `hinge_<name>.json`, `hinge_<name>.png` |
| `optical_readout.py` | Step-A JSONs | angle → reporter separation `r(θ)=2ρ·sin(θ/2)` → **FRET** + **plasmon-ruler** signal | `optical_readout_sweep.csv`, `optical_readout.png` |
| `binding_to_signal.py` | k_θ + optical models | **ligand ΔG → Δθ → Δsignal**, mechanical SNR, per-ligand deflection | `binding_to_signal.csv`, `binding_to_signal.png` |

## Key results

**1. The linker stiffens the hinge (thesis result, reframed as mechanics).**
`k_θ = k_BT/σ²` turns the σ table into a stiffness curve: the hinge stiffens
**7.5×** from 10→40 nt. Fit: `σ(L) = 1.96° + 10.35·exp(−L/11.3 nt)` — an ~11 nt
stiffening length scale and a ~2° angular noise floor.

**2. The structure is a clean planar X at ~80° open.**
The geometry tool splits the two beams almost perfectly (8748 vs 8740 atoms),
locates the crux at the center, and measures hinge angles of **81°/78°** (two
frames of the open state) and **32°** (a closed snapshot) — resolving genuine
open-vs-closed conformations. Beam length ~100 nm; tip–tip reporter distance
~65 nm when open.

**3. Two transducers, two placements.**
FRET (R₀≈5.4 nm) only works with dyes **near the crux** (optimal radius
ρ≈5.3 nm), where efficiency swings **E = 0.92 → 0.08** across closed→open. At
the tips FRET is dead (~10⁻⁷). The **plasmon ruler** (Au, D=30 nm) works **at
the tips**, giving a **~52 nm spectral shift** between states — the right
readout for the tens-of-nm tip regime.

**4. The biosensor loop, fully quantitative.** A **10 kcal/mol** (nM-aptamer)
ligand deflects a soft 10-nt hinge **36°** (ΔFRET 0.42) but a stiff 40-nt hinge
only **13°** (ΔFRET 0.24). The emergent design law:

> **Mechanical SNR is linker-independent** — `Δθ/σ = √(2ΔG/k_BT)` (stiffness
> cancels; SNR = 5.8 for every linker at 10 kcal/mol), and reaching a 3σ
> deflection costs a **universal 2.68 kcal/mol (4.5 k_BT)** regardless of length.

So the linker doesn't set detectability per unit binding energy — it sets the
**absolute angular scale**, i.e. *which reporter geometry can transduce it*.
Pick a short soft hinge (large angles) for the tip plasmon ruler; a long stiff
hinge (small precise angles, low noise) for near-crux FRET. The thesis'
counter-intuitive stiffening becomes a tunable design feature.

## Running

```bash
# from the repo root; venv already provisioned (numpy, matplotlib, pypdf)
.venv/bin/python analysis/stiffness_from_sigma.py
.venv/bin/python analysis/hinge_geometry.py [snapshot.xyz]   # default: switch_origami_structure.xyz
.venv/bin/python analysis/optical_readout.py                 # reads hinge_*.json
.venv/bin/python analysis/binding_to_signal.py               # reads optical_readout.py
```

`hinge_geometry.py` auto-detects three snapshot formats: element-letter XYZ
(A/T/G/C), MDAnalysis XYZ, and LAMMPS-dump lines (`id type x y z …`).

## Units & assumptions

- Lengths in oxDNA reduced units; 1 unit ≈ 0.8518 nm.
- Simulation T = 0.1 (oxDNA) ≈ 300 K, since 1 oxDNA energy unit ≈ 4.142×10⁻²⁰ J
  gives k_BT = 0.1 units; RT = 0.596 kcal/mol.
- Harmonic hinge (equipartition) near the free-energy minimum.
- FRET: `E = 1/(1+(r/R₀)⁶)`. Plasmon ruler: `Δλ/λ₀ = 0.18·exp(−(r−D)/(0.23·D))`
  (Jain et al., *Nano Lett.* 2007).

## Honest limitations

- **Single frames, not distributions.** The production trajectory dumps are no
  longer on disk, so hinge angles here are *instantaneous* geometries — signal-
  vs-angle is a design map, not a re-measured σ. Recovering σ needs the dumps
  or a re-run.
- **Chiral-plasmonic CD is out of scope.** A circular-dichroism readout needs an
  out-of-plane twist the current planar X lacks; `optical_readout.py` keeps a
  dormant `chirality_proxy()` to motivate a deliberately non-planar redesign.
- **Coarse-grained.** Inherits oxDNA2's abstractions (no sequence-specific
  stacking beyond seqdep, implicit ions via Debye–Hückel).

## Natural next steps (higher compute)

- Free-energy profile **F(θ)** of the scissors via metadynamics on a *truncated*
  hinge model (feasible overnight on a laptop) — the true global open/close
  landscape and anharmonicity beyond the harmonic k_θ used here.
- Model an **aptamer/strand-displacement clamp** bridging the beams and compute
  the two-state allosteric coupling directly.
