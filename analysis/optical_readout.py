#!/usr/bin/env python3
"""
Tier-1 (zero MD): turn the scissor hinge angle into a measurable optical signal.

Chain of models
---------------
    hinge angle  ->  reporter separation  ->  optical observable

1. Geometry.  Two reporters sit a distance rho from the crux, one on each
   beam, on the open-mouth side. With equal arms the law of cosines gives

        r(theta) = 2 * rho * sin(theta / 2)

   (closed scissors theta->0 => r->0; flat-open theta->180 => r->2*rho).
   This is validated against Step A: beam half-length ~50 nm, measured
   tip-tip at 81.1 deg = 65.5 nm, and 2*50*sin(40.55 deg) = 64.7 nm. Match.

2. FRET (short-range ruler, ~1-10 nm).  Forster transfer efficiency

        E(r) = 1 / (1 + (r / R0)^6)

   Useful only when r ~ R0 (5-6 nm), so FRET reporters must be placed NEAR
   the crux (small rho). At the tips (r ~ 65 nm) FRET is identically zero.

3. Plasmon ruler (long-range, tens of nm).  Two Au nanoparticles of diameter
   D; the universal plasmon-ruler equation (Jain et al., Nano Lett. 2007)

        dlambda / lambda0 = A * exp( -(s / D) / tau ),   s = r - D  (edge gap)

   with A ~ 0.18, tau ~ 0.23. Works exactly in the tens-of-nm regime the
   scissor tips live in -- the right transducer for this device.

Biosensor loop
--------------
A binding event that changes the hinge by d(theta) produces a signal change
d(observable); combined with the Tier-0 stiffness (actuation energy per d theta)
this ties ligand-binding free energy to a readable optical shift.
"""

import os, glob, json
import numpy as np

OXDNA_NM = 0.8518

# ---- geometry -------------------------------------------------------------- #
def tip_distance(theta_deg, rho_nm):
    return 2.0 * rho_nm * np.sin(np.radians(theta_deg) / 2.0)

# ---- FRET ------------------------------------------------------------------ #
def fret_efficiency(r_nm, R0_nm=5.4):
    return 1.0 / (1.0 + (r_nm / R0_nm) ** 6)

# ---- plasmon ruler --------------------------------------------------------- #
def plasmon_shift(r_nm, D_nm=30.0, A=0.18, tau=0.23):
    """Fractional resonance shift dlambda/lambda0. NaN where particles overlap
    (r < D, physically impossible)."""
    s = r_nm - D_nm
    out = A * np.exp(-(s / D_nm) / tau)
    out = np.where(r_nm < D_nm, np.nan, out)
    return out

# ---- chiral-plasmon proxy (out-of-plane handedness) ------------------------ #
def chirality_proxy(axis1, axis2, plane_normal):
    """Geometric CD proxy ~ (d1 x d2).n * (d1.d2). Zero for a planar X;
    nonzero only if the two beams acquire an out-of-plane twist. Motivates a
    deliberately non-planar design for circular-dichroism readout."""
    d1 = np.asarray(axis1); d2 = np.asarray(axis2); n = np.asarray(plane_normal)
    return float(np.dot(np.cross(d1, d2), n) * np.dot(d1, d2))

# ---- design helpers -------------------------------------------------------- #
def best_fret_rho(theta_closed, theta_open, R0_nm=5.4, rho_grid=None):
    """Reporter radius (near crux) that maximises FRET contrast |dE| between a
    closed and an open state."""
    if rho_grid is None:
        rho_grid = np.linspace(1.0, 20.0, 400)
    rc = tip_distance(theta_closed, rho_grid)
    ro = tip_distance(theta_open, rho_grid)
    dE = np.abs(fret_efficiency(rc, R0_nm) - fret_efficiency(ro, R0_nm))
    i = int(np.argmax(dE))
    return rho_grid[i], dE[i]


def load_measured():
    """Pull measured (angle, tip distances) from Step A's JSON outputs."""
    here = os.path.dirname(os.path.abspath(__file__))
    pts = []
    for jp in sorted(glob.glob(os.path.join(here, "hinge_*.json"))):
        d = json.load(open(jp))
        pts.append((d["file"], d["hinge_angle_deg"],
                    d.get("nearest_tip_pair_nm", np.nan)))
    return pts


def main():
    here = os.path.dirname(os.path.abspath(__file__))

    RHO_TIP = 50.0          # nm, reporter at the blade tip (beam half-length)
    R0      = 5.4           # nm, Cy3/Cy5-like Forster radius
    D_AUNP  = 30.0          # nm, gold nanoparticle diameter
    TH_CLOSED, TH_OPEN = 40.0, 100.0   # nominal switch end-states

    theta = np.linspace(20, 140, 400)

    # ---- geometry / distances
    r_tip  = tip_distance(theta, RHO_TIP)

    # ---- FRET: optimal near-crux placement
    rho_fret, dE_max = best_fret_rho(TH_CLOSED, TH_OPEN, R0)
    r_fret = tip_distance(theta, rho_fret)
    E_fret = fret_efficiency(r_fret, R0)

    # ---- plasmon at the tips
    dl = plasmon_shift(r_tip, D_AUNP)

    measured = load_measured()

    # -------- printed design summary --------
    print("=== Optical readout design summary ===")
    print(f"Geometry: r(theta) = 2*rho*sin(theta/2)\n")

    print("FRET channel (place dyes NEAR the crux):")
    print(f"  optimal reporter radius rho   = {rho_fret:.1f} nm from pivot")
    print(f"  max FRET contrast |dE| ({TH_CLOSED:.0f}->{TH_OPEN:.0f} deg) = {dE_max:.2f}")
    Ec = fret_efficiency(tip_distance(TH_CLOSED, rho_fret), R0)
    Eo = fret_efficiency(tip_distance(TH_OPEN,   rho_fret), R0)
    print(f"  E(closed {TH_CLOSED:.0f} deg)={Ec:.2f}   E(open {TH_OPEN:.0f} deg)={Eo:.2f}")
    print(f"  (at the tips rho={RHO_TIP:.0f} nm FRET E ~ "
          f"{fret_efficiency(tip_distance(90,RHO_TIP),R0):.1e}  -> dead, as expected)\n")

    print(f"Plasmon-ruler channel (Au D={D_AUNP:.0f} nm at the tips rho={RHO_TIP:.0f} nm):")
    dlc = plasmon_shift(tip_distance(TH_CLOSED, RHO_TIP), D_AUNP)
    dlo = plasmon_shift(tip_distance(TH_OPEN,   RHO_TIP), D_AUNP)
    lam0 = 530.0  # nm, typical AuNP dipole resonance, for an absolute shift
    print(f"  dlambda/lambda0: closed={dlc:.3f}  open={dlo:.3f}  "
          f"contrast={abs(dlc-dlo):.3f}")
    print(f"  ~absolute shift at lambda0={lam0:.0f} nm: "
          f"{abs(dlc-dlo)*lam0:.1f} nm between states\n")

    if measured:
        print("Measured snapshots (from Step A) mapped onto both channels:")
        for name, th, tipnm in measured:
            print(f"  {name:34s} theta={th:5.1f} deg | "
                  f"FRET E(rho={rho_fret:.0f})={fret_efficiency(tip_distance(th,rho_fret),R0):.2f} | "
                  f"plasmon dl/l(tip)={plasmon_shift(tip_distance(th,RHO_TIP),D_AUNP):.3f}")

    # -------- CSV --------
    csv = os.path.join(here, "optical_readout_sweep.csv")
    import csv as _csv
    with open(csv, "w", newline="") as f:
        w = _csv.writer(f)
        w.writerow(["theta_deg", "r_tip_nm", "r_fret_nm",
                    "FRET_E_near_crux", "plasmon_dlambda_over_lambda_tip"])
        for i in range(len(theta)):
            w.writerow([theta[i], r_tip[i], r_fret[i], E_fret[i], dl[i]])
    print(f"\nWrote {csv}")

    # -------- figure --------
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))

        ax[0].plot(theta, r_tip, color="#1d3557", label=f"tips (rho={RHO_TIP:.0f} nm)")
        ax[0].plot(theta, r_fret, color="#e07a5f",
                   label=f"near crux (rho={rho_fret:.0f} nm)")
        ax[0].set_title("Reporter separation r(theta)")
        ax[0].set_xlabel("hinge angle (deg)"); ax[0].set_ylabel("r (nm)")
        ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)

        ax[1].plot(theta, E_fret, color="#e07a5f")
        ax[1].set_title(f"FRET readout (dyes at crux, rho={rho_fret:.0f} nm, R0={R0} nm)")
        ax[1].set_xlabel("hinge angle (deg)"); ax[1].set_ylabel("FRET efficiency E")
        ax[1].set_ylim(-0.02, 1.02); ax[1].grid(alpha=.3)

        ax[2].plot(theta, dl, color="#1d3557")
        ax[2].set_title(f"Plasmon-ruler readout (Au D={D_AUNP:.0f} nm at tips)")
        ax[2].set_xlabel("hinge angle (deg)")
        ax[2].set_ylabel(r"$\Delta\lambda/\lambda_0$")
        ax[2].grid(alpha=.3)

        # overlay measured snapshots
        for name, th, _ in measured:
            for a, y in ((ax[1], fret_efficiency(tip_distance(th, rho_fret), R0)),
                         (ax[2], plasmon_shift(tip_distance(th, RHO_TIP), D_AUNP))):
                a.axvline(th, color="gray", ls=":", lw=.8)
                a.scatter([th], [y], color="k", zorder=5, s=30)
        fig.tight_layout()
        png = os.path.join(here, "optical_readout.png")
        fig.savefig(png, dpi=140)
        print(f"Wrote {png}")
    except ImportError:
        print("(matplotlib missing -> figure skipped)")


if __name__ == "__main__":
    main()
