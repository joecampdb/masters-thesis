#!/usr/bin/env python3
"""
Tier-0 (zero new MD): reinterpret the thesis bond-angle standard deviations
(Table 1) as the mechanical stiffness of the DNA-origami scissor hinge.

Physics
-------
Near the bottom of its free-energy well the hinge behaves as a harmonic
angular spring, U(theta) = 1/2 k_theta (theta - theta0)^2. Equipartition then
ties the measured angular fluctuation directly to the spring constant:

        <(d theta)^2> = k_B T / k_theta      =>   k_theta = k_B T / sigma^2

So each sigma in Table 1 is, with no extra simulation, an effective torsional
stiffness. Longer linker -> smaller sigma -> STIFFER hinge, which reframes the
thesis' counter-intuitive result in the language of mechanics/metamaterials.

Unit note
---------
The production runs used oxDNA reduced units at T = 0.1, and 1 oxDNA energy
unit epsilon ~ 4.142e-20 J, so k_B T = 0.1 epsilon corresponds to ~300 K
(k_B * 300 K = 4.14e-21 J = 0.1 epsilon). We report k_theta in three ways:
  * k_B T / rad^2   (dimensionless; this is just 1/sigma_rad^2)
  * kcal/mol / rad^2 (multiply by RT = 0.593 kcal/mol at 300 K)
  * kcal/mol / deg^2 (the per-degree form, often handier for design)

We also fit sigma(L) = sigma_inf + A*exp(-L/lambda) (the exponential decay the
thesis reports) with a dependency-light method, and report the actuation
energy needed to open the scissors by a target angle -- the number that sets
how much ligand-binding free energy a downstream biosensor must supply.
"""

import numpy as np

# ---- Thesis Table 1: linker length (nt) -> sigma of bond-angle dist (deg) ----
L     = np.array([10, 20, 30, 40], dtype=float)
sigma = np.array([6.235704832, 3.736856603, 2.675925874, 2.272207295])

# ---- constants / unit mapping ----
T_SIM        = 0.1            # oxDNA reduced temperature used in the runs
EPS_J        = 4.142e-20      # 1 oxDNA energy unit in Joules
KB_J         = 1.380649e-23   # Boltzmann constant, J/K
T_KELVIN     = EPS_J * T_SIM / KB_J           # physical T implied by T_sim=0.1
RT_KCAL      = 1.987204e-3 * T_KELVIN         # kcal/mol at that temperature
DEG2RAD      = np.pi / 180.0

def exp_decay_fit(x, y):
    """Fit y = c + A*exp(-x/lam) without scipy.

    For each lam on a grid the model is linear in (c, A), solved by least
    squares; we keep the lam with the smallest residual. Returns (c, A, lam).
    """
    best = None
    for lam in np.linspace(1.0, 200.0, 20000):
        basis = np.column_stack([np.ones_like(x), np.exp(-x / lam)])
        coef, res, *_ = np.linalg.lstsq(basis, y, rcond=None)
        resid = np.sum((basis @ coef - y) ** 2)
        if best is None or resid < best[0]:
            best = (resid, coef[0], coef[1], lam)
    _, c, A, lam = best
    return c, A, lam

def main():
    sigma_rad = sigma * DEG2RAD

    # equipartition stiffness
    k_kbt   = 1.0 / sigma_rad**2            # k_theta in units of k_B T / rad^2
    k_kcal  = k_kbt * RT_KCAL               # kcal/mol / rad^2
    k_deg   = k_kcal * DEG2RAD**2           # kcal/mol / deg^2

    print(f"Implied physical temperature (T_sim=0.1): {T_KELVIN:6.1f} K")
    print(f"RT at that temperature:                   {RT_KCAL:6.4f} kcal/mol\n")

    hdr = f"{'L(nt)':>6} {'sigma(deg)':>11} {'k [kBT/rad^2]':>15} " \
          f"{'k [kcal/mol/rad^2]':>19} {'k [kcal/mol/deg^2]':>19}"
    print(hdr)
    print("-" * len(hdr))
    for i in range(len(L)):
        print(f"{L[i]:6.0f} {sigma[i]:11.3f} {k_kbt[i]:15.2f} "
              f"{k_kcal[i]:19.3f} {k_deg[i]:19.5f}")

    # relative stiffening across the tested range
    print(f"\nHinge stiffens {k_kbt[-1]/k_kbt[0]:.2f}x going from "
          f"10 -> 40 nt linker (stiffness ~ 1/sigma^2).")

    # exponential decay of sigma(L)
    c, A, lam = exp_decay_fit(L, sigma)
    print(f"\nFit sigma(L) = {c:.3f} + {A:.3f}*exp(-L/{lam:.1f} nt)")
    print(f"  plateau sigma_inf = {c:.3f} deg  (asymptotic angular noise floor)")
    print(f"  decay length      = {lam:.1f} nt  (linker length scale of stiffening)")

    # actuation energy: how much energy to open the scissors by Dtheta
    print("\nActuation energy  E = 1/2 k_theta * Dtheta^2  (kcal/mol):")
    print("  -> the ligand-binding free energy a biosensor must supply to")
    print("     mechanically drive the hinge by Dtheta.")
    print(f"{'L(nt)':>6}" + "".join(f"{d:>10}deg" for d in (10, 20, 30)))
    for i in range(len(L)):
        row = f"{L[i]:6.0f}"
        for dtheta in (10, 20, 30):
            E = 0.5 * k_deg[i] * dtheta**2
            row += f"{E:13.3f}"
        print(row)

    # write a tidy CSV for downstream plotting / the manuscript
    import csv, os
    out = os.path.join(os.path.dirname(__file__), "stiffness_table.csv")
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["L_nt", "sigma_deg", "k_kBT_per_rad2",
                    "k_kcal_per_rad2", "k_kcal_per_deg2"])
        for i in range(len(L)):
            w.writerow([L[i], sigma[i], k_kbt[i], k_kcal[i], k_deg[i]])
    print(f"\nWrote {out}")

    # optional figure if matplotlib is available
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        Lfine = np.linspace(L.min(), L.max(), 200)
        fig, ax = plt.subplots(1, 2, figsize=(10, 4))
        ax[0].plot(L, sigma, "o", color="#2a6f97")
        ax[0].plot(Lfine, c + A*np.exp(-Lfine/lam), "-", color="#2a6f97",
                   label=f"$\\sigma_\\infty$+A$e^{{-L/\\lambda}}$, $\\lambda$={lam:.0f} nt")
        ax[0].set_xlabel("linker length L (nt)")
        ax[0].set_ylabel(r"bond-angle $\sigma$ (deg)")
        ax[0].set_title("Angular fluctuation (thesis)")
        ax[0].legend(fontsize=8)
        ax[1].plot(L, k_deg, "s-", color="#bc4749")
        ax[1].set_xlabel("linker length L (nt)")
        ax[1].set_ylabel(r"hinge stiffness $k_\theta$ (kcal/mol/deg$^2$)")
        ax[1].set_title("Reinterpreted as stiffness")
        for a in ax:
            a.grid(alpha=0.3)
        fig.tight_layout()
        figpath = os.path.join(os.path.dirname(__file__), "stiffness_from_sigma.png")
        fig.savefig(figpath, dpi=150)
        print(f"Wrote {figpath}")
    except ImportError:
        print("(matplotlib not installed -> skipped figure; CSV still written)")

if __name__ == "__main__":
    main()
