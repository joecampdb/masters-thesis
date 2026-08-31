#!/usr/bin/env python3
"""
Assemble an extended edition of the thesis PDF: Part I reproduces the original
M.S. thesis (Campagna & Carnevale, 2025); Part II documents the 2026
zero-additional-simulation extension chain (stiffness, geometry, optical
transduction, binding->signal coupling) with the generated figures.

Output: origami_MNT/extended_thesis.pdf
Greek letters are spelled out (sigma/theta/lambda) so the standard PDF fonts
render cleanly and match the figure labels.
"""

import os
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Image,
                                Table, TableStyle, PageBreak)
from reportlab.lib.utils import ImageReader

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
OUT = os.path.join(REPO, "extended_thesis.pdf")
CONTENT_W = 6.5 * inch

# ------------------------------------------------------------------ styles -- #
ss = getSampleStyleSheet()
H1 = ParagraphStyle("H1", parent=ss["Heading1"], fontSize=14, spaceBefore=14,
                    spaceAfter=6, textColor=colors.HexColor("#1d3557"))
H2 = ParagraphStyle("H2", parent=ss["Heading2"], fontSize=11.5, spaceBefore=10,
                    spaceAfter=4, textColor=colors.HexColor("#2a6f97"))
BODY = ParagraphStyle("BODY", parent=ss["BodyText"], fontSize=9.7, leading=13.5,
                      alignment=TA_JUSTIFY, spaceAfter=6)
CAP = ParagraphStyle("CAP", parent=ss["BodyText"], fontSize=8.3, leading=10.5,
                     textColor=colors.HexColor("#444444"), spaceBefore=3,
                     spaceAfter=10, alignment=TA_CENTER)
TITLE = ParagraphStyle("TITLE", parent=ss["Title"], fontSize=17, leading=21)
SUB = ParagraphStyle("SUB", parent=ss["Normal"], fontSize=10.5,
                     alignment=TA_CENTER, textColor=colors.HexColor("#333333"))
PART = ParagraphStyle("PART", parent=ss["Heading1"], fontSize=15,
                      alignment=TA_CENTER, textColor=colors.HexColor("#bc4749"),
                      spaceBefore=6, spaceAfter=10)
REF = ParagraphStyle("REF", parent=BODY, fontSize=8.4, leading=10.8,
                     alignment=TA_JUSTIFY, spaceAfter=2, leftIndent=14,
                     firstLineIndent=-14)

story = []


def P(t, s=BODY): story.append(Paragraph(t, s))
def gap(h=6): story.append(Spacer(1, h))


def figure(name, caption, max_w=CONTENT_W):
    path = os.path.join(HERE, name)
    iw, ih = ImageReader(path).getSize()
    w = max_w
    h = w * ih / iw
    if h > 4.2 * inch:                      # cap tall/square figures
        h = 4.2 * inch; w = h * iw / ih
    story.append(Image(path, width=w, height=h))
    story.append(Paragraph(caption, CAP))


# =============================================================== title ====== #
P("In silico rational design of DNA origami-based "
  "allosteric metamaterials", TITLE)
gap(4)
P("<b>Extended Edition</b>", SUB)
gap(6)
P("Joseph Campagna &nbsp;&middot;&nbsp; Vincenzo Carnevale", SUB)
gap(2)
P("Original M.S. thesis 2025 &nbsp;&middot;&nbsp; Extended computational "
  "analysis 2026", SUB)
gap(14)

P("<b>Abstract.</b> We characterize the bond-angle distribution across the "
  "linker of DNA-origami scissor nanoswitches with coarse-grained molecular "
  "dynamics (oxDNA2 in LAMMPS), and find the counter-intuitive result that "
  "<i>longer</i> poly-thymine linkers <i>reduce</i> the angular standard "
  "deviation (sigma: 6.24&deg; to 2.27&deg; for 10 to 40 nt), following an "
  "exponential decay. This Extended Edition adds a zero-additional-simulation "
  "analysis chain that (i) reinterprets sigma as an effective hinge stiffness, "
  "(ii) measures the scissor geometry directly from structural snapshots, "
  "(iii) maps the hinge angle onto FRET and plasmonic optical readouts, and "
  "(iv) couples ligand-binding free energy to hinge deflection, optical signal, "
  "and a predicted analyte dose-response with picomolar-scale limits of "
  "detection. Two design principles emerge: the linker's stiffening is a "
  "tunable feature rather than a nuisance, and the mechanical signal-to-noise "
  "of the switch is independent of linker length, so the linker selects "
  "<i>which</i> optical transducer can read the device rather than its "
  "intrinsic detectability.")
gap(6)

# =============================================================== PART I ====== #
story.append(PageBreak())
P("PART I &mdash; Original thesis (2025)", PART)

P("1. Introduction", H1)
P("DNA origami has emerged as a powerful technique in structural DNA "
  "nanotechnology, enabling precisely defined nanoscale structures with "
  "applications in drug delivery, biosensing, and molecular computing. First "
  "introduced by Rothemund in 2006, the method folds a long single-stranded "
  "scaffold into predetermined shapes through hybridization of short staple "
  "strands, yielding devices with nanometer spatial precision that are ideal "
  "candidates for molecular machines and switches.")
P("DNA-origami switches are dynamic nanostructures that undergo conformational "
  "changes in response to molecular triggers. They typically consist of two "
  "rigid origami beams connected by a flexible linker that serves as a hinge. "
  "Despite growing interest, a quantitative understanding of their mechanical "
  "properties remains incomplete; in particular, the relationship between "
  "linker sequence and conformational flexibility has not been systematically "
  "characterized. While persistence length is the traditional descriptor of "
  "DNA rigidity, in complex origami the interplay of sequence-dependent "
  "flexibility, electrostatics, and mechanical strain can deviate from "
  "classical polymer models, motivating a detailed look at bond-angle "
  "distributions that govern the conformational landscape.")
P("We hypothesize that the average standard deviation of the hinge bond-angle "
  "distribution &mdash; a direct measure of conformational flexibility &mdash; "
  "is set by the rigidity of the linker. We test this with oxDNA structural "
  "modeling and coarse-grained molecular-dynamics simulations using the CG-DNA "
  "package in LAMMPS, establishing quantitative relationships between linker "
  "parameters and structural flexibility toward the rational design of "
  "DNA-based nanodevices.")

P("2. Methods", H1)
P("2.1 Structure preparation and modeling.", H2)
P("Initial modeling used oxView; the caDNAno JSON geometry was imported with "
  "the TacoxDNA importer (caDNAno format, hexagonal lattice). The single bond "
  "connecting the two origami beams was nicked (Nick), a thymine linker "
  "generated (Seq), and the beams re-connected (Ligate); the bottom beam was "
  "rotated 90&deg; and ligated to the linker terminus. This produced switches "
  "with 10-, 20-, 30-, and 40-nucleotide poly-T linkers, converted to LAMMPS "
  "data files via the TacoxDNA oxDNA-to-LAMMPS converter.")
P("2.2 Molecular dynamics and analysis.", H2)
P("Simulations used LAMMPS (stable release 29 Aug 2024) with the oxDNA2 force "
  "field, temperature 0.1, timestep 5e-4, for 1e7 timesteps per structure. "
  "Bond angles were measured in VMD (v1.9.4) from three points &mdash; the top "
  "of the linker, the crux, and the bottom beam &mdash; and the mean and "
  "standard deviation of each distribution were computed with custom Python.")

P("3. Results", H1)
P("Contrary to the initial hypothesis that longer linkers would increase "
  "angular variance, an inverse correlation was observed between linker length "
  "and angular deviation. Mean bond angles remained stable, while the "
  "magnitude of angular fluctuation fell sharply with linker length, most "
  "steeply between 10 and 25 nt and plateauing thereafter &mdash; an "
  "exponential-decay trend (Table 1).")

t1 = Table([["Linker length (nt)", "sigma (bond-angle std. dev., deg)"],
            ["10", "6.24"], ["20", "3.74"], ["30", "2.68"], ["40", "2.27"]],
           colWidths=[2.4*inch, 3.0*inch])
t1.setStyle(TableStyle([
    ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#1d3557")),
    ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
    ("FONTSIZE", (0, 0), (-1, -1), 9),
    ("ALIGN", (0, 0), (-1, -1), "CENTER"),
    ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#aaaaaa")),
    ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#eef3f7")]),
]))
story.append(t1)
story.append(Paragraph("Table 1. Linker length vs. sigma, the mean standard "
                       "deviation of the hinge bond-angle distribution.", CAP))

P("4. Discussion", H1)
P("Longer DNA linkers reduce, rather than increase, bond-angle variability "
  "&mdash; challenging the polymer-physics intuition that longer chains grant "
  "more conformational freedom. Beyond a threshold length the linker imparts "
  "structural stability that buffers angular deviation: additional nucleotides "
  "distribute strain and provide a denser network of stacking and "
  "electrostatic interactions, and the entropic landscape stabilizes once a "
  "critical length is reached. This is consistent with regime-specific "
  "behavior of compliant DNA nanostructures &mdash; short linkers as soft "
  "pivots, longer linkers approaching entropic rod-like stability &mdash; and "
  "underscores that persistence length alone cannot capture the localized "
  "angular distributions that dictate hinge mechanics. Limitations include the "
  "coarse-grained abstraction, thymine-only linkers, and static hinge "
  "geometries without applied stimuli.")

P("5. Conclusion", H1)
P("Systematic oxDNA2 simulations establish a clear inverse relationship "
  "between linker length and hinge bond-angle variability (sigma from 6.24&deg; "
  "to 2.27&deg; across 10-40 nt, exponential decay), indicating that longer "
  "linkers provide enhanced conformational stability rather than increased "
  "flexibility &mdash; a design principle for programmable, mechanically "
  "reliable DNA nanodevices.")

# =============================================================== PART II ===== #
story.append(PageBreak())
P("PART II &mdash; Extended analysis (2026)", PART)
P("The following sections extend the study using only the outputs already in "
  "hand (the sigma table and structural snapshots); no additional molecular "
  "dynamics was run. All computations execute on a laptop in seconds. Because "
  "the production trajectory dumps were not retained, geometric quantities "
  "below are measured from individual snapshots and are therefore "
  "instantaneous conformations rather than re-sampled distributions.", BODY)

P("6. From angular variance to hinge stiffness", H1)
P("Near the minimum of its free-energy well the hinge is a harmonic angular "
  "spring, U = 1/2 k_theta (theta - theta0)^2, so equipartition ties the "
  "measured fluctuation directly to a stiffness: k_theta = kB T / sigma^2. "
  "Each sigma in Table 1 therefore becomes an effective torsional stiffness "
  "with no additional simulation. The hinge stiffens 7.5-fold from 10 to 40 nt "
  "(stiffness scales as 1/sigma^2). An exponential fit gives sigma(L) = 1.96 + "
  "10.35 exp(-L/11.3 nt): an ~11 nt stiffening length scale and a ~2&deg; "
  "angular noise floor. The actuation energy 1/2 k_theta (delta theta)^2 sets "
  "how much binding free energy a downstream sensor must supply to drive the "
  "hinge by a given angle (Section 9).")
figure("stiffness_from_sigma.png",
       "Figure 6. Left: the thesis angular fluctuation sigma(L) with the "
       "exponential-decay fit. Right: the same data reinterpreted as hinge "
       "torsional stiffness k_theta, which rises 7.5-fold across 10-40 nt.")

P("7. Direct scissor geometry from single frames", H1)
P("A geometry analyzer loads a structural snapshot, finds the molecular plane "
  "by PCA (the device is nearly flat), and separates the two crossing beams "
  "with an expectation-maximization two-line fit &mdash; necessary because "
  "global PCA bisects the X and cannot yield the beam axes directly. It then "
  "reports the hinge angle, the crux (pivot), each beam's length, and the "
  "tip-tip reporter distances. On the labeled structure the beam split is "
  "essentially exact (8748 vs 8740 nucleotides), the crux sits at the center, "
  "and the hinge opens to 81&deg;; two frames of the open state agree to within "
  "~3&deg;, while a separate closed snapshot reads 32&deg; &mdash; confirming the "
  "tool resolves genuine open/close states. Beams are ~100 nm long and the "
  "open-mouth tip separation is ~65 nm.")
figure("hinge_switch_origami_structure.png",
       "Figure 7. Single-frame geometry of the scissor nanoswitch. The EM "
       "two-line fit cleanly assigns the two beams (blue/red); the crux (star) "
       "is central and the fitted axes give a hinge angle of 81&deg;.",
       max_w=4.4*inch)

P("8. Optical transduction of the hinge angle", H1)
P("With reporters a distance rho from the crux, the reporter separation is "
  "r(theta) = 2 rho sin(theta/2), validated against Section 7 (predicted 64.7 "
  "nm vs measured 65.5 nm at 81&deg;). Two transducers cover complementary "
  "ranges. FRET (efficiency E = 1/(1 + (r/R0)^6), R0 ~ 5.4 nm) is a short-range "
  "ruler: it is essentially dead at the ~65 nm tips but, with dyes placed near "
  "the crux (optimal rho ~ 5.3 nm), swings E from 0.92 (closed) to 0.08 (open) "
  "&mdash; a near-full-scale ratiometric readout. The plasmon ruler "
  "(delta lambda/lambda0 = 0.18 exp(-(r-D)/(0.23 D)) for gold particles of "
  "diameter D) instead operates at the tens-of-nm tips, giving a ~52 nm "
  "spectral shift between states. The device thus supports orthogonal "
  "dual-channel detection: FRET at the crux for single molecules, a plasmon "
  "ruler at the tips for ensembles.")
figure("optical_readout.png",
       "Figure 8. Optical readout of the hinge. Left: reporter separation "
       "r(theta) for tip vs near-crux placement. Center: FRET efficiency vs "
       "hinge angle for near-crux dyes. Right: plasmon-ruler fractional shift "
       "for tip-mounted gold particles. Black points are the measured "
       "snapshots (32/78/81&deg;).")

P("9. Coupling binding energy to signal", H1)
P("Combining the stiffness of Section 6 with the transducers of Section 8 "
  "closes the biosensor loop: a ligand supplying free energy delta G can drive "
  "at most delta theta_max = sqrt(2 delta G / k_theta), producing a computable "
  "signal change. For a typical nanomolar-aptamer budget (delta G ~ 10 "
  "kcal/mol), a soft 10 nt hinge deflects 36&deg; (FRET change 0.42) whereas a "
  "stiff 40 nt hinge deflects only 13&deg; (FRET change 0.24). Crucially, the "
  "thermal noise is sigma = sqrt(kB T / k_theta) &mdash; exactly the thesis "
  "sigma &mdash; so the mechanical signal-to-noise ratio delta theta/sigma = "
  "sqrt(2 delta G / kB T) is <i>independent of linker length</i> (SNR = 5.8 "
  "for every linker at 10 kcal/mol), and reaching a 3-sigma deflection costs a "
  "universal 2.68 kcal/mol (4.5 kB T). The linker therefore sets the absolute "
  "angular scale &mdash; and hence which reporter geometry can transduce the "
  "motion &mdash; not the per-energy detectability.")
figure("binding_to_signal.png",
       "Figure 9. Left: actuation cost 1/2 k_theta (delta theta)^2 per linker "
       "length, with reference ligand affinities (dotted). Right: transduced "
       "FRET change vs deflection; dots mark what a 10 kcal/mol ligand reaches "
       "on each linker.")

P("10. Predicted dose-response and limit of detection", H1)
P("Treating the device as an allosteric sensor with 1:1 (Langmuir) recognition, "
  "the bound fraction f = [A]/([A]+Kd) sets a population-weighted ensemble "
  "signal, S([A]) = S_unbound + (S_bound - S_unbound) f, a sigmoid in log[A] "
  "centred at the affinity Kd. The blank noise follows from the thermal angular "
  "fluctuation propagated through the transducer, sigma_S = |dS/dtheta| sigma, "
  "and the limit of detection is taken at three blank standard deviations "
  "(net signal = 3 sigma_S). For a nanomolar aptamer (Kd = 1 nM), the model "
  "predicts limits of detection of roughly 200-430 pM on the FRET channel and "
  "28-250 pM on the plasmon channel &mdash; i.e. detection well below Kd, the "
  "signature of a well-coupled allosteric switch. The softest hinge gives the "
  "lowest limit (its larger deflection outruns its higher noise), and the "
  "plasmon ruler outperforms FRET at these tip distances (Figure 10).")
figure("dose_response.png",
       "Figure 10. Predicted biosensor performance. Left: sigmoidal dose-"
       "response for a 20 nt linker, FRET and plasmon channels, with Kd (dashed) "
       "and the 3-sigma limit of detection (dotted). Right: predicted limit of "
       "detection versus linker length for both channels at Kd = 1 nM; shorter, "
       "softer hinges reach lower limits.")

P("11. Extended discussion and design principles", H1)
P("<b>Stiffening as a feature.</b> The thesis' counter-intuitive result "
  "becomes a design axis: a short soft hinge throws large angles (ideal for "
  "the tip-mounted plasmon ruler), while a long stiff hinge throws small, "
  "precise angles on a low noise floor (ideal for near-crux FRET, where a few "
  "degrees already crosses R0). One tunes the linker to match the readout.")
P("<b>Universal detectability.</b> Because both the driven deflection and the "
  "thermal noise scale as 1/sqrt(k_theta), stiffness cancels from the "
  "signal-to-noise ratio. Detectability per unit binding energy is set by the "
  "ligand relative to kB T, not by the linker &mdash; a clean, testable "
  "prediction that only emerges once stiffness and transduction are modeled "
  "together.")
P("<b>Limit of detection favours the soft, matched hinge.</b> Because the "
  "detection limit is set by dynamic range relative to blank noise, the softest "
  "linker &mdash; largest deflection per binding event &mdash; reaches the "
  "lowest limit despite its higher angular noise, provided its throw is read by "
  "the transducer whose range it matches (the tip plasmon ruler here). Linker "
  "choice is therefore a joint optimisation of stiffness and reporter geometry, "
  "not of stiffness alone.")
P("<b>Reproducibility and an interactive tool.</b> The four-stage analysis is "
  "released as a small, dependency-light Python chain (stiffness, geometry, "
  "optical readout, binding-to-signal and dose-response), and an interactive "
  "browser calculator exposes the same physics for live what-if design of the "
  "linker, reporter placement, and target affinity.")
P("<b>Limitations.</b> Snapshot geometries are instantaneous, not "
  "distributions (the trajectory dumps were not retained); a circular-dichroism "
  "channel would require an out-of-plane twist the current planar X lacks; and "
  "the analysis inherits oxDNA2's coarse-grained abstractions.")

P("12. Outlook", H1)
P("Two low-to-moderate compute steps follow naturally. First, a metadynamics "
  "reconstruction of the full free-energy profile F(theta) on a truncated "
  "hinge model (linker plus short beam stubs) would capture the global "
  "open/close landscape and its anharmonicity beyond the harmonic stiffness "
  "used here &mdash; and is feasible overnight on a laptop by shrinking the "
  "system 5-10 fold. Second, modeling an aptamer or strand-displacement clamp "
  "bridging the beams would give the two-state allosteric coupling of a "
  "concrete biosensor directly, realizing the functional element envisioned "
  "for the hinge.")

# =============================================================== refs ======= #
story.append(PageBreak())
P("References", H1)
refs = [
 "Seeman, N.C. &amp; Sleiman, H.F. (2018). DNA nanotechnology. Nature Reviews Materials 3, 17068.",
 "Zhang, F. et al. (2014). Structural DNA nanotechnology: state of the art and future perspective. JACS 136, 11198.",
 "Rothemund, P.W.K. (2006). Folding DNA to create nanoscale shapes and patterns. Nature 440, 297.",
 "Marras, A.E. et al. (2015). Programmable motion of DNA origami mechanisms. PNAS 112, 713.",
 "Lauback, S. et al. (2018). Real-time magnetic actuation of DNA nanodevices. Nature Communications 9, 1446.",
 "Ke, Y. et al. (2016). Regulation at a distance using a DNA origami nanoactuator. Nature Communications 7, 10935.",
 "Hudoba, M.W. et al. (2017). Dynamic DNA origami device for measuring compressive depletion forces. ACS Nano 11, 6566.",
 "Gerling, T. et al. (2015). Dynamic DNA devices and assemblies from shape-complementary components. Science 347, 1446.",
 "Zhou, L. et al. (2014). DNA origami compliant nanostructures with tunable mechanical properties. ACS Nano 8, 27.",
 "Schiffels, D. et al. (2013). Nanoscale structure and microscale stiffness of DNA nanotubes. ACS Nano 7, 6700.",
 "Kauert, D.J. et al. (2011). Direct mechanical measurements of 3D DNA origami. Nano Letters 11, 5558.",
 "Shi, Z. et al. (2017). Conformational dynamics of compliant DNA nanostructures from CG-MD. ACS Nano 11, 4617.",
 "Bui, P.T.M. et al. (2019). Thermodynamic and kinetic aspects of DNA nanotechnology. J. Sci. Adv. Mater. Devices 4, 1.",
 "Doye, J.P.K. et al. (2013). Coarse-graining DNA for simulations of DNA nanotechnology. PCCP 15, 20395.",
 "Henrich, O. et al. (2018). Coarse-grained simulation of DNA using LAMMPS. Eur. Phys. J. E 41, 57.",
 "Snodin, B.E.K. et al. (2015). Improved structural properties and salt dependence in a CG DNA model. JCP 142, 234901.",
 "Idili, A., Vallee-Belisle, A. &amp; Ricci, F. (2014). Programmable pH-triggered DNA nanoswitches. JACS 136, 5836.",
 "Douglas, S.M., Bachelet, I. &amp; Church, G.M. (2012). A logic-gated nanorobot for targeted transport. Science 335, 831.",
 "Chatterjee, G. et al. (2017). A spatially localized architecture for fast and modular DNA computing. Nature Nanotechnology 12, 920.",
 "Jain, P.K. et al. (2007). On the universal scaling behavior of the distance decay of plasmon coupling (plasmon ruler equation). Nano Letters 7, 2080. [added, Part II]",
]
for i, r in enumerate(refs, 1):
    P(f"[{i}] {r}", REF)

# ------------------------------------------------------------------ build --- #
def footer(canvas, doc):
    canvas.saveState()
    canvas.setFont("Helvetica", 8)
    canvas.setFillColor(colors.HexColor("#888888"))
    canvas.drawCentredString(letter[0]/2, 0.5*inch,
        f"Campagna & Carnevale - Allosteric metamaterials (Extended Edition) - p.{doc.page}")
    canvas.restoreState()

doc = SimpleDocTemplate(OUT, pagesize=letter,
                        leftMargin=0.9*inch, rightMargin=0.9*inch,
                        topMargin=0.9*inch, bottomMargin=0.8*inch,
                        title="Allosteric metamaterials - Extended Edition",
                        author="Joseph Campagna; Vincenzo Carnevale")
doc.build(story, onFirstPage=footer, onLaterPages=footer)
print("Wrote", OUT)
