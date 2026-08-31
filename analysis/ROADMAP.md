# Project Roadmap — DNA-origami scissor nanoswitch biosensor

**Goal:** take the thesis extension from "nice analysis" to something that
catches the eye of a top structural-DNA-nanotech lab (e.g., ASU Biodesign —
Hao Yan / Stephanopoulos: origami, chiral plasmonics, allosteric devices) and a
translational diagnostics startup. Two currencies: **a validated, publishable
insight** (academic) and **a demoable product with a number attached**
(startup). Both are within reach at low compute.

Status: `[x]` done · `[~]` partial · `[ ]` to-do

---

## Compute tiers (reiterated)

| Tier | Cost | Definition |
|------|------|------------|
| 0 | free | arithmetic / geometry on existing outputs, no MD |
| 1 | seconds | Python modeling, no MD |
| 2 | overnight on a laptop | one MD run each, **truncated** hinge model (linker + beam stubs, 5–10× smaller) |
| 3 | cluster | full-origami enhanced sampling, lattices, ML datasets |

### Tier 0 — done
- `[x]` σ → hinge stiffness `k_θ=k_BT/σ²` (7.5× stiffening, actuation-energy table) — `stiffness_from_sigma.py`
- `[x]` single-frame scissor geometry (hinge angle, crux, beam length, tip distances) — `hinge_geometry.py`
- `[~]` angle *distribution* / σ recompute — **blocked**: trajectory dumps not retained (needs a re-run)

### Tier 1 — mostly done
- `[x]` optical readout: FRET (crux) + plasmon ruler (tips) — `optical_readout.py`
- `[x]` binding ΔG → Δθ → Δsignal, universal-SNR law — `binding_to_signal.py`
- `[ ]` **dose–response / limit-of-detection** (Langmuir isotherm → signal vs. analyte concentration) — *buildable now*
- `[ ]` **inverse-design helper** (target signal + ligand Kd → recommended linker length + reporter radius) — *buildable now*

### Tier 2 — the rigor upgrade
- `[ ]` **F(θ) free-energy landscape** via metadynamics on a truncated hinge (PLUMED/colvars on the existing oxDNA2 stack) — the real global open/close curve; *validates* the harmonic k_θ and the universal-detectability claim
- `[ ]` salt (`rhos`) and temperature sweeps → stimulus-response curves
- `[ ]` sequence-dependent linkers (mixed / GC / hairpin) beyond poly-T

### Tier 3 — cluster-scale (shelve for a laptop)
- `[ ]` umbrella sampling on the full origami (error-barred F(θ))
- `[ ]` metamaterial lattice → auxetic / negative-Poisson response (**cash in the "metamaterials" title**)
- `[ ]` ML surrogate + inverse-design dataset

---

## New to-dos — for an academic nanotech lab (ASU-style)

1. `[ ]` **F(θ) landscape (Tier 2)** — converts the harmonic estimate and the
   linker-independent-SNR law into a *validated* result. This is the
   publishable core; everything else supports it.
2. `[ ]` **Preprint (bioRxiv / arXiv)** — polish `extended_thesis.pdf` into a
   methods-complete preprint with a reproducibility statement. A citable
   artifact is the fastest way to get a lab's attention.
3. `[ ]` **Open-source, pip-installable analysis package + notebook** — the
   four tools as `origami-hinge` with a Jupyter walkthrough. Tools get cited and
   adopted; it signals engineering competence.
4. `[ ]` **Chiral-plasmonic redesign** — deliberately introduce an out-of-plane
   twist so the two beams form a chiral pair → circular-dichroism readout.
   Directly aligned with ASU's chiral-plasmonics strength; the current planar X
   gives CD ≈ 0 (already flagged by `chirality_proxy()`). Design + Tier-2 sim.
5. `[ ]` **Experimental validation plan** — oxView design export, staple set,
   and a concrete assay (gel/TEM for folding, FRET or DLS for actuation) a
   collaborating wet lab could run. Turns a dry project into a collaboration ask.
6. `[ ]` **Error bars + sequence dependence** — replicas and mixed-sequence
   linkers; standard reviewer-proofing.

## New to-dos — for a translational diagnostics startup

*(Inferring the audience is a translational DNA-nanotech / diagnostics startup;
tell me if "solidSF" is a specific company and I'll tailor to it.)*

1. `[ ]` **Pick a target + state an LOD** — e.g., miRNA-21 (oncology) or cardiac
   troponin; the dose–response model yields a predicted limit of detection. A
   product needs a number.
2. `[ ]` **Interactive web design-calculator (Artifact)** — sliders for linker
   length / reporter placement / ligand Kd → live Δθ, FRET & plasmon signal,
   SNR, LOD. A 10-second shareable demo that sells the concept. *Buildable now.*
3. `[ ]` **One-page product brief** — assay format (homogeneous mix-and-read),
   readout hardware (plate reader, or a smartphone camera for the plasmonic
   color shift), reagent bill-of-materials, cost, time-to-result.
4. `[ ]` **Novelty / IP memo** (technical, not legal) — what's defensible:
   linker-tuned transducer matching, orthogonal dual-channel (crux-FRET +
   tip-plasmon) readout, and the design rules themselves.
5. `[ ]` **Manufacturability note** — scalable origami production (phage
   scaffold, one-pot folding), lyophilized reagent stability for point-of-care.
6. `[ ]` **Positioning scan** — vs. molecular beacons, aptamer switches, and
   existing origami sensors: what this device does that they don't (large
   mechanical gain, dual readout, quantitative design rules).

---

## Recommended next three (highest leverage / lowest compute)

1. **Dose–response + LOD model** (Tier 1, now) — the single most persuasive
   addition for a diagnostics audience; makes the sensor claim quantitative.
2. **Interactive design-calculator Artifact** (now) — the eye-catching demo;
   reuses the exact physics already coded.
3. **F(θ) metadynamics pipeline** (Tier 2 — prep now, run overnight) — the
   rigor that earns an academic's attention and validates the headline law.
