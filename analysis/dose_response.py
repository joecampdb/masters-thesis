#!/usr/bin/env python3
"""
Tier-1 (zero MD): predicted analyte dose-response and limit of detection (LOD)
for the scissor nanoswitch used as an allosteric biosensor.

Model
-----
* 1:1 binding (Langmuir): bound fraction f = [A] / ([A] + Kd).
* Two-state transduction: unbound devices rest at theta_open; a bound device
  is actuated (closed) by the binding free energy dG = -RT ln(Kd) into a
  harmonic hinge of stiffness k_theta (Tier-0), giving a deflection
  d(theta) = sqrt(2 dG / k_theta) (capped by geometry).
* Ensemble signal is population-weighted:
      S([A]) = S_unbound + (S_bound - S_unbound) * f          (sigmoid in log[A])
  with S from the Step-B FRET or plasmon transducer.
* Blank noise from the thermal angular fluctuation sigma (= the thesis sigma),
  propagated through the transducer:  sigma_S = |dS/dtheta|(theta_open) * sigma.
* LOD (IUPAC 3-sigma): net signal = 3 sigma_S  =>  f_LOD = 3 sigma_S / DR,
  [A]_LOD = Kd * f_LOD / (1 - f_LOD),  where DR = |S_bound - S_unbound|.

Reuses fret_efficiency / plasmon_shift / tip_distance from optical_readout.py.
"""

import os
import numpy as np
from optical_readout import fret_efficiency, plasmon_shift, tip_distance

# ---- Tier-0 stiffness (mirrors stiffness_from_sigma.py) -------------------- #
L      = np.array([10, 20, 30, 40], float)
SIGMA  = np.array([6.235704832, 3.736856603, 2.675925874, 2.272207295])  # deg
RT     = 1.987204e-3 * 300.0                         # kcal/mol (~0.596)
K_THETA = RT / SIGMA**2                              # kcal/mol/deg^2

# ---- transducer settings (from Step B) ------------------------------------- #
THETA_OPEN = 90.0                # unbound rest angle (deg)
RHO_FRET, R0 = 5.3, 5.4          # dyes near crux
RHO_TIP, D_AUNP, LAM0 = 50.0, 30.0, 530.0   # AuNP at tips

def dG_from_Kd(Kd_M):
    return -RT * np.log(Kd_M)     # kcal/mol (Kd in molar)

def signal(theta, channel):
    if channel == "fret":
        return fret_efficiency(tip_distance(theta, RHO_FRET), R0)
    s = plasmon_shift(tip_distance(theta, RHO_TIP), D_AUNP)
    return np.nan_to_num(s) * LAM0                    # nm shift

def dS_dtheta(theta, channel, h=0.5):
    return (signal(theta + h, channel) - signal(theta - h, channel)) / (2 * h)

def response(Kd_M, k_theta, sigma_deg, channel):
    dG = dG_from_Kd(Kd_M)
    dtheta = min(np.sqrt(2 * dG / k_theta), THETA_OPEN)     # cap at full close
    theta_bound = THETA_OPEN - dtheta
    S_un = signal(THETA_OPEN, channel)
    S_bd = signal(theta_bound, channel)
    DR = abs(S_bd - S_un)
    sigma_S = abs(dS_dtheta(THETA_OPEN, channel)) * sigma_deg
    f_LOD = 3 * sigma_S / DR if DR > 0 else np.inf
    if f_LOD >= 1:
        LOD = np.inf
    else:
        LOD = Kd_M * f_LOD / (1 - f_LOD)
    return dict(dtheta=dtheta, theta_bound=theta_bound, S_un=S_un, S_bd=S_bd,
                DR=DR, sigma_S=sigma_S, f_LOD=f_LOD, LOD=LOD)

def main():
    here = os.path.dirname(os.path.abspath(__file__))
    Kd = 1e-9                       # 1 nM aptamer (protein biomarker regime)
    print(f"Analyte binding Kd = {Kd*1e9:.0f} nM  "
          f"(dG = {dG_from_Kd(Kd):.1f} kcal/mol), rest angle {THETA_OPEN:.0f} deg\n")

    for ch, unit in (("fret", "dE"), ("plas", "nm")):
        print(f"=== {ch.upper()} channel ===")
        hdr = (f"{'L(nt)':>5} {'sigma':>6} {'dtheta':>7} {'DR':>8} "
               f"{'sigma_S':>9} {'LOD':>12}")
        print(hdr); print("-"*len(hdr))
        for i in range(len(L)):
            r = response(Kd, K_THETA[i], SIGMA[i], ch)
            lod = ("inf" if not np.isfinite(r["LOD"])
                   else f"{r['LOD']*1e12:.1f} pM")
            print(f"{L[i]:5.0f} {SIGMA[i]:6.2f} {r['dtheta']:7.1f} "
                  f"{r['DR']:8.3f} {r['sigma_S']:9.4f} {lod:>12}")
        print()

    # dose-response figure for a representative linker (20 nt)
    i = 1
    conc = np.logspace(-12, -5, 300)            # 1 pM .. 10 uM
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(12, 5))
        for ch, c, lab in (("fret", "#e07a5f", "FRET |E|"),
                           ("plas", "#1d3557", "plasmon shift (nm)")):
            r = response(Kd, K_THETA[i], SIGMA[i], ch)
            f = conc / (conc + Kd)
            S = r["S_un"] + (r["S_bd"] - r["S_un"]) * f
            axx = ax[0] if ch == "fret" else ax[0].twinx()
            axx.plot(conc, S, color=c, label=lab)
            axx.set_ylabel(lab, color=c)
            axx.tick_params(axis="y", labelcolor=c)
            if np.isfinite(r["LOD"]):
                axx.axvline(r["LOD"], color=c, ls=":", lw=1)
        ax[0].axvline(Kd, color="gray", ls="--", lw=1)
        ax[0].text(Kd, ax[0].get_ylim()[1]*0.9, " Kd", fontsize=8, color="gray")
        ax[0].set_xscale("log"); ax[0].set_xlabel("analyte concentration (M)")
        ax[0].set_title(f"Dose-response, {int(L[i])} nt linker (dotted = LOD)")
        ax[0].grid(alpha=.3, which="both")

        # LOD vs linker, both channels
        for ch, c in (("fret", "#e07a5f"), ("plas", "#1d3557")):
            lods = [response(Kd, K_THETA[j], SIGMA[j], ch)["LOD"]*1e12
                    for j in range(len(L))]
            ax[1].plot(L, lods, "o-", color=c, label=f"{ch}")
        ax[1].set_xlabel("linker length (nt)"); ax[1].set_ylabel("LOD (pM)")
        ax[1].set_title(f"Predicted LOD vs linker (Kd={Kd*1e9:.0f} nM)")
        ax[1].legend(); ax[1].grid(alpha=.3)
        fig.tight_layout()
        png = os.path.join(here, "dose_response.png")
        fig.savefig(png, dpi=140); print(f"Wrote {png}")
    except ImportError:
        print("(matplotlib missing -> figure skipped)")

    # CSV
    import csv
    with open(os.path.join(here, "dose_response.csv"), "w", newline="") as fcsv:
        w = csv.writer(fcsv)
        w.writerow(["channel","L_nt","sigma_deg","dtheta_deg","dynamic_range",
                    "sigma_S","LOD_M"])
        for ch in ("fret", "plas"):
            for j in range(len(L)):
                r = response(Kd, K_THETA[j], SIGMA[j], ch)
                w.writerow([ch, L[j], SIGMA[j], f"{r['dtheta']:.2f}",
                            f"{r['DR']:.4f}", f"{r['sigma_S']:.5f}",
                            f"{r['LOD']:.3e}"])
    print(f"Wrote {os.path.join(here,'dose_response.csv')}")

if __name__ == "__main__":
    main()
