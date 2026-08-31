#!/usr/bin/env python3
"""
Tier-0/A (zero MD): measure the scissor-hinge geometry of a single DNA-origami
snapshot.

The nanoswitch is a flat "X": two rigid origami beams crossing at a central
crux (the linker). Given only a coordinate snapshot this tool:

  1. loads the point cloud (auto-detects several .xyz-ish formats),
  2. finds the molecular plane by PCA (the structure is ~flat),
  3. separates the two crossing beams with an EM two-line fit
     (global PCA bisects the X, so it cannot give the beam axes directly),
  4. reports the hinge angle, the crux (pivot), each beam's length, the four
     blade tips, and the tip-tip distances -- the quantities a FRET or
     plasmonic reporter reads out (front-end for the Tier-1 optical model).

Optionally, if the file carries base letters (A/T/G/C), it flags the longest
contiguous poly-T run as a candidate engineered linker.

Length unit: oxDNA sim units; 1 unit ~ 0.8518 nm (reported in both).

Usage:
    python hinge_geometry.py [snapshot.xyz] [--json out.json] [--png out.png]
"""

import sys, os, json, argparse
import numpy as np

OXDNA_NM = 0.8518   # nm per oxDNA length unit


# --------------------------------------------------------------------------- #
# Loading
# --------------------------------------------------------------------------- #
def load_xyz(path):
    """Return (coords Nx3 float, bases list[str] or None).

    Handles: standard XYZ with element letters (A/T/G/C or X); MDAnalysis XYZ;
    and LAMMPS-dump-style lines 'id type x y z [ix iy iz ...]'.
    """
    raw = open(path).read().splitlines()
    if not raw:
        raise ValueError("empty file")

    # Standard XYZ: first line is an integer atom count, second is a comment.
    first = raw[0].strip()
    body = raw
    declared = None
    if first.isdigit():
        declared = int(first)
        body = raw[2:2 + declared]

    coords, bases = [], []
    for ln in body:
        p = ln.split()
        if len(p) < 4:
            continue
        # element-letter format: "T x y z"
        if not _is_float(p[0]):
            bases.append(p[0])
            coords.append([float(p[1]), float(p[2]), float(p[3])])
        else:
            # numeric leading columns: could be "type x y z" (4+) or
            # "id type x y z ..." (LAMMPS dump). Detect by column count.
            nums = [float(x) for x in p]
            if len(nums) >= 8:            # id type x y z ix iy iz ...
                bases.append(str(int(nums[1])))
                coords.append(nums[2:5])
            else:                          # type x y z
                bases.append(str(int(nums[0])))
                coords.append(nums[1:4])
    coords = np.asarray(coords, float)
    has_letters = any(b in "ATGC" for b in bases)
    return coords, (bases if has_letters else None)


def _is_float(s):
    try:
        float(s); return True
    except ValueError:
        return False


# --------------------------------------------------------------------------- #
# Geometry
# --------------------------------------------------------------------------- #
def pca_frame(xyz):
    c = xyz - xyz.mean(0)
    w, v = np.linalg.eigh(np.cov(c.T))
    idx = np.argsort(w)[::-1]
    return w[idx], v[:, idx], xyz.mean(0)


def _line_dist(pts, p, d):
    """Perpendicular distance of each point to the line (p, unit d)."""
    q = pts - p
    proj = q @ d
    perp = q - np.outer(proj, d)
    return np.linalg.norm(perp, axis=1)


def em_two_lines(xyz, plane_basis, center, iters=25):
    """Separate two crossing beams. Returns (labels, lines) where each line is
    (point, unit_direction). Initialised from the two dominant in-plane blade
    orientations, then refined by alternating assignment / PCA refit."""
    e1, e2 = plane_basis[:, 0], plane_basis[:, 1]
    d = xyz - center
    u, w = d @ e1, d @ e2
    r = np.hypot(u, w)
    theta = np.degrees(np.arctan2(w, u)) % 180.0     # orientation mod 180

    # distance-weighted orientation histogram -> two peaks = two beam axes
    hist, edges = np.histogram(theta, bins=180, range=(0, 180), weights=r)
    # smooth (circular) so we pick robust maxima
    hist = np.convolve(np.r_[hist[-5:], hist, hist[:5]], np.ones(11)/11, "same")[5:-5]
    order = np.argsort(hist)[::-1]
    peaks = []
    for i in order:
        a = edges[i] + 0.5
        if all(min(abs(a-b), 180-abs(a-b)) > 20 for b in peaks):
            peaks.append(a)
        if len(peaks) == 2:
            break
    dirs = []
    for a in peaks:
        ang = np.radians(a)
        dirs.append(np.cos(ang)*e1 + np.sin(ang)*e2)

    lines = [(center.copy(), dd/np.linalg.norm(dd)) for dd in dirs]
    labels = np.zeros(len(xyz), int)
    for _ in range(iters):
        dists = np.column_stack([_line_dist(xyz, p, dd) for p, dd in lines])
        labels = dists.argmin(1)
        new = []
        for k in (0, 1):
            pts = xyz[labels == k]
            if len(pts) < 3:
                new.append(lines[k]); continue
            cc = pts - pts.mean(0)
            ww, vv = np.linalg.eigh(np.cov(cc.T))
            new.append((pts.mean(0), vv[:, np.argmax(ww)]))
        lines = new
    return labels, lines


def line_intersection(l1, l2):
    """Closest-approach midpoint of two 3D lines = the crux/pivot."""
    p1, d1 = l1; p2, d2 = l2
    r = p1 - p2
    a, b, c = d1 @ d1, d1 @ d2, d2 @ d2
    dd, e = d1 @ r, d2 @ r
    denom = a*c - b*b
    if abs(denom) < 1e-9:
        return (p1 + p2) / 2
    t = (b*e - c*dd) / denom
    s = (a*e - b*dd) / denom
    return ((p1 + t*d1) + (p2 + s*d2)) / 2


def beam_tips(pts, pivot, d):
    """The two extreme points of a beam along its axis, and its length."""
    proj = (pts - pivot) @ d
    lo, hi = pts[proj.argmin()], pts[proj.argmax()]
    return lo, hi, float(proj.max() - proj.min())


def analyze(path):
    xyz, bases = load_xyz(path)
    n = len(xyz)
    evals, basis, centroid = pca_frame(xyz)
    labels, lines = em_two_lines(xyz, basis, centroid)
    pivot = line_intersection(lines[0], lines[1])

    d1, d2 = lines[0][1], lines[1][1]
    cosang = abs(d1 @ d2)
    angle = float(np.degrees(np.arccos(np.clip(cosang, 0, 1))))   # 0..90 acute

    beams = []
    for k in (0, 1):
        pts = xyz[labels == k]
        lo, hi, length = beam_tips(pts, pivot, lines[k][1])
        beams.append(dict(natoms=int(len(pts)),
                          axis=lines[k][1].tolist(),
                          length_units=length, length_nm=length*OXDNA_NM,
                          tip_lo=lo.tolist(), tip_hi=hi.tolist()))

    # tip-tip distances across beams (what a reporter dye/particle pair reads)
    def dist(a, b): return float(np.linalg.norm(np.array(a)-np.array(b)))
    t = [beams[0]["tip_lo"], beams[0]["tip_hi"],
         beams[1]["tip_lo"], beams[1]["tip_hi"]]
    cross = {
        "b0lo-b1lo": dist(t[0], t[2]), "b0lo-b1hi": dist(t[0], t[3]),
        "b0hi-b1lo": dist(t[1], t[2]), "b0hi-b1hi": dist(t[1], t[3]),
    }
    # the two closest cross-pairs are the "open mouth" tip pairs of the scissors
    ordered = sorted(cross.items(), key=lambda kv: kv[1])

    result = dict(
        file=os.path.basename(path), n_atoms=n,
        planarity_pc3_over_pc1=float(evals[2]/evals[0]),
        hinge_angle_deg=angle,
        hinge_opening_deg=180.0 - angle,
        pivot=pivot.tolist(),
        beams=beams,
        tip_tip_distances_units=cross,
        tip_tip_distances_nm={k: v*OXDNA_NM for k, v in cross.items()},
        nearest_tip_pair=ordered[0][0],
        nearest_tip_pair_nm=ordered[0][1]*OXDNA_NM,
    )
    if bases is not None:
        s = "".join(b if b in "ATGC" else "-" for b in bases)
        import re
        runs = list(re.finditer("T+", s))
        if runs:
            m = max(runs, key=lambda mm: mm.end()-mm.start())
            result["longest_polyT_run"] = dict(
                length=m.end()-m.start(), atom_start=m.start(), atom_end=m.end())
    return result, xyz, labels, lines, pivot, basis, centroid


# --------------------------------------------------------------------------- #
# Reporting / figure
# --------------------------------------------------------------------------- #
def print_report(R):
    print(f"\nFile: {R['file']}   atoms: {R['n_atoms']}")
    print(f"Planarity (PC3/PC1 var): {R['planarity_pc3_over_pc1']:.3f}  "
          f"(<0.1 => flat X)")
    print(f"\n>>> HINGE ANGLE: {R['hinge_angle_deg']:.1f} deg  "
          f"(opening supplement {R['hinge_opening_deg']:.1f} deg)")
    print(f"    crux/pivot: [{', '.join(f'{x:.1f}' for x in R['pivot'])}]")
    for i, b in enumerate(R["beams"]):
        print(f"    beam {i}: {b['natoms']} atoms, "
              f"length {b['length_units']:.1f} u = {b['length_nm']:.1f} nm")
    print("\n  Tip-tip distances (nm) -- candidate FRET/plasmon reporter pairs:")
    for k, v in sorted(R["tip_tip_distances_nm"].items(), key=lambda kv: kv[1]):
        tag = "  <- narrowest (open-mouth tips)" if k == R["nearest_tip_pair"] else ""
        print(f"      {k:12s} {v:6.1f} nm{tag}")
    if "longest_polyT_run" in R:
        r = R["longest_polyT_run"]
        print(f"\n  Longest poly-T run: {r['length']} nt "
              f"(atoms {r['atom_start']}-{r['atom_end']}) -- candidate linker")


def make_figure(R, xyz, labels, lines, pivot, basis, centroid, png):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    e1, e2 = basis[:, 0], basis[:, 1]
    P = (xyz - centroid)
    x, y = P @ e1, P @ e2
    piv = (pivot - centroid)
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.scatter(x[labels == 0], y[labels == 0], s=3, c="#2a6f97", label="beam 0")
    ax.scatter(x[labels == 1], y[labels == 1], s=3, c="#bc4749", label="beam 1")
    ax.scatter(piv @ e1, piv @ e2, s=120, c="k", marker="*", zorder=5, label="crux")
    for k, (p, d) in enumerate(lines):
        pp = (p - centroid)
        for sgn in (-1, 1):
            end = pp + sgn * 60 * d
            ax.plot([pp @ e1, end @ e1], [pp @ e2, end @ e2], "k--", lw=1, alpha=.6)
    ax.set_aspect("equal"); ax.legend(loc="upper right", fontsize=9)
    ax.set_xlabel("in-plane PC1"); ax.set_ylabel("in-plane PC2")
    ax.set_title(f"{R['file']}\nhinge angle = {R['hinge_angle_deg']:.1f} deg")
    fig.tight_layout(); fig.savefig(png, dpi=140)
    print(f"\nWrote {png}")


def main():
    ap = argparse.ArgumentParser()
    default = "/Users/josephanthony/Downloads/switch_origami_structure.xyz"
    ap.add_argument("snapshot", nargs="?", default=default)
    ap.add_argument("--json", default=None)
    ap.add_argument("--png", default=None)
    a = ap.parse_args()

    R, xyz, labels, lines, pivot, basis, centroid = analyze(a.snapshot)
    print_report(R)

    here = os.path.dirname(os.path.abspath(__file__))
    stem = os.path.splitext(os.path.basename(a.snapshot))[0]
    jpath = a.json or os.path.join(here, f"hinge_{stem}.json")
    ppath = a.png or os.path.join(here, f"hinge_{stem}.png")
    json.dump(R, open(jpath, "w"), indent=2)
    print(f"\nWrote {jpath}")
    make_figure(R, xyz, labels, lines, pivot, basis, centroid, ppath)


if __name__ == "__main__":
    main()
