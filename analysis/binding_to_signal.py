#!/usr/bin/env python3
"""
Close the biosensor loop: ligand-binding free energy -> hinge deflection ->
optical signal, using Tier-0 stiffness and the Step-B optical models together.

Physics chain
-------------
    dG_bind  --(mechanical work)-->  d(theta)  --(transducer)-->  d(signal)

* Tier-0 gives the hinge stiffness k_theta(L) = k_B T / sigma^2 (kcal/mol/deg^2).
* Bending the hinge by d(theta) from its rest angle costs the harmonic work
      dG_deform = 1/2 * k_theta * d(theta)^2 .
  A ligand can drive at most d(theta)_max = sqrt(2 * dG_bind / k_theta).
* The thermal angular noise is sigma = sqrt(k_B T / k_theta) -- i.e. exactly the
  thesis sigma (the chain round-trips), so the mechanical (angular) SNR is
      SNR = d(theta)_signal / sigma = sqrt(2 * dG_bind / k_B T),
  which is INDEPENDENT of linker length. Equivalently the energy to reach a
  target SNR = n is universal:  dG = 1/2 n^2 k_B T.  The linker sets the
  absolute angular scale (which reporter/placement transduces it), not the
  per-energy detectability.
* Step B's transducers map the achievable d(theta) to a real observable:
  FRET (dyes near the crux) and the Au plasmon ruler (particles at the tips).

Reuses fret_efficiency / plasmon_shift / tip_distance from optical_readout.py.
"""

import os
import numpy as np
from optical_readout import fret_efficiency, plasmon_shift, tip_distance

# ---- Tier-0 stiffness (mirrors stiffness_from_sigma.py) -------------------- #
L      = np.array([10, 20, 30, 40], float)
SIGMA  = np.array([6.235704832, 3.736856603, 2.675925874, 2.272207295])  # deg
EPS_J, KB_J, T_SIM = 4.142e-20, 1.380649e-23, 0.1
T_K    = EPS_J * T_SIM / KB_J                 # ~300 K
RT     = 1.987204e-3 * T_K                    # kcal/mol (~0.596)
DEG    = np.pi / 180.0

K_THETA = RT / (SIGMA * DEG)**2 * DEG**2      # kcal/mol/deg^2  (== RT/sigma_deg^2)
#   note: k in kcal/mol/deg^2 reduces to RT / sigma_deg^2

# ---- reference ligand affinities (dG = -RT ln Kd) -------------------------- #
REF = {
    "6-bp toehold (~9)":      9.0,
    "small-mol aptamer uM (~8)":  8.2,
    "protein aptamer nM (~12)":  12.3,
    "antibody pM (~16)":         16.5,
}

# ---- transducer operating points (from Step B) ----------------------------- #
RHO_FRET, R0, TH_FRET   = 5.3, 5.4, 61.0     # dyes near crux; E~0.5 here
RHO_TIP, D_AUNP, TH_PLAS = 50.0, 30.0, 47.0   # AuNP at tips; steep plasmon here


def dtheta_max(dG, k):
    return np.sqrt(2.0 * dG / k)


def signal_change(theta0, dtheta, channel):
    """Magnitude of the optical signal change for an opening deflection."""
    if channel == "fret":
        s0 = fret_efficiency(tip_distance(theta0, RHO_FRET), R0)
        s1 = fret_efficiency(tip_distance(theta0 + dtheta, RHO_FRET), R0)
        return abs(s1 - s0)                      # dE (unitless)
    else:
        s0 = plasmon_shift(tip_distance(theta0, RHO_TIP), D_AUNP)
        s1 = plasmon_shift(tip_distance(theta0 + dtheta, RHO_TIP), D_AUNP)
        return abs(s1 - s0) * 530.0              # absolute nm shift (lam0=530)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    print(f"T = {T_K:.0f} K, RT = {RT:.3f} kcal/mol\n")

    # universal detectability law
    print("Universal detectability (linker-independent):")
    for n in (1, 3, 5):
        print(f"  SNR={n}:  dG = 1/2 n^2 kT = {0.5*n*n*RT:.2f} kcal/mol")
    print()

    # main coupling table at a representative ligand budget
    DG = 10.0
    print(f"=== Ligand budget dG_bind = {DG:.0f} kcal/mol (typical nM aptamer) ===")
    hdr = (f"{'L(nt)':>5} {'k_theta':>9} {'sigma':>7} {'dtheta_max':>10} "
           f"{'SNR':>5} {'dFRET@61':>9} {'dlambda@47(nm)':>14}")
    print(hdr); print("-"*len(hdr))
    rows = []
    for i in range(len(L)):
        k = K_THETA[i]; sig = SIGMA[i]
        dth = dtheta_max(DG, k)
        snr = dth / sig
        dE = signal_change(TH_FRET, dth, "fret")
        dl = signal_change(TH_PLAS, dth, "plas")
        rows.append((L[i], k, sig, dth, snr, dE, dl))
        print(f"{L[i]:5.0f} {k:9.4f} {sig:7.2f} {dth:9.1f}  {snr:5.1f} "
              f"{dE:9.2f} {dl:14.1f}")
    print("\n(SNR is identical across L: sqrt(2*dG/kT) = "
          f"{np.sqrt(2*DG/RT):.1f} -- linker sets the angular scale, not SNR.)\n")

    # inverse: energy to reach a detectable 3-sigma deflection, per linker
    print("Energy to drive a 3-sigma (SNR=3) deflection -- what a sensor needs:")
    for i in range(len(L)):
        dth3 = 3*SIGMA[i]
        dG3 = 0.5*K_THETA[i]*dth3**2
        dE = signal_change(TH_FRET, dth3, "fret")
        print(f"  L={L[i]:.0f} nt: dtheta={dth3:4.1f} deg  costs {dG3:.2f} kcal/mol"
              f"  -> dFRET={dE:.2f}")
    print("  (cost is ~2.68 kcal/mol for every L: 4.5 kT.)\n")

    # per-ligand achievable deflection, per linker
    print("Max hinge deflection dtheta (deg) achievable by common ligands:")
    lhdr = f"{'ligand':28s}" + "".join(f"{int(l):>7}nt" for l in L)
    print(lhdr)
    for name, dG in REF.items():
        line = f"{name:28s}" + "".join(f"{dtheta_max(dG,K_THETA[i]):9.1f}"
                                        for i in range(len(L)))
        print(line)

    # CSV
    import csv
    with open(os.path.join(here, "binding_to_signal.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["L_nt","k_theta_kcal_deg2","sigma_deg",
                    "dtheta_max_deg@10kcal","SNR","dFRET","dlambda_nm"])
        for r in rows:
            w.writerow([r[0], f"{r[1]:.5f}", f"{r[2]:.3f}", f"{r[3]:.2f}",
                        f"{r[4]:.2f}", f"{r[5]:.3f}", f"{r[6]:.2f}"])
    print(f"\nWrote {os.path.join(here,'binding_to_signal.csv')}")

    # figure: energy landscape + achievable signal
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        dth = np.linspace(0, 40, 300)
        fig, ax = plt.subplots(1, 2, figsize=(12, 5))
        colors = ["#e63946", "#f4a261", "#2a9d8f", "#1d3557"]
        for i in range(len(L)):
            ax[0].plot(dth, 0.5*K_THETA[i]*dth**2, color=colors[i],
                       label=f"{int(L[i])} nt (k={K_THETA[i]:.3f})")
        for name, dG in REF.items():
            ax[0].axhline(dG, ls=":", color="gray", lw=.8)
            ax[0].text(0.3, dG+0.3, name.split(" (")[0], fontsize=7, color="gray")
        ax[0].set_xlabel(r"hinge deflection $\Delta\theta$ (deg)")
        ax[0].set_ylabel(r"binding energy required $\Delta G$ (kcal/mol)")
        ax[0].set_title("Actuation cost vs deflection (Tier-0 stiffness)")
        ax[0].set_ylim(0, 20); ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)

        # resulting FRET change vs deflection at the FRET operating point
        dE = [signal_change(TH_FRET, d, "fret") for d in dth]
        ax[1].plot(dth, dE, color="#e07a5f", label="FRET $|\\Delta E|$ @61 deg")
        for i in range(len(L)):
            dmax = dtheta_max(10.0, K_THETA[i])
            ax[1].scatter([dmax], [signal_change(TH_FRET, dmax, "fret")],
                          color=colors[i], zorder=5,
                          label=f"{int(L[i])} nt @10 kcal/mol")
        ax[1].set_xlabel(r"hinge deflection $\Delta\theta$ (deg)")
        ax[1].set_ylabel(r"FRET change $|\Delta E|$")
        ax[1].set_title("Transduced signal (dots = reachable by 10 kcal/mol)")
        ax[1].legend(fontsize=7); ax[1].grid(alpha=.3)
        fig.tight_layout()
        png = os.path.join(here, "binding_to_signal.png")
        fig.savefig(png, dpi=140)
        print(f"Wrote {png}")
    except ImportError:
        print("(matplotlib missing -> figure skipped)")


if __name__ == "__main__":
    main()
