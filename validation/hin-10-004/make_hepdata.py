#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11,<3.14"
# dependencies = [
#   "hepdata-lib>=0.20",
#   "hepdata-validator>=0.3",
#   "scipy>=1.11",
#   "numpy>=1.26",
# ]
# ///
"""Build and validate the HEPData submission for the CMS Open Data replication of CMS-HIN-10-004.

Adapted from the hepdata skill template (assets/make_hepdata.py). Reads all
values directly from summary.json and attaches the rendered PDF figures from
plots/. Runs hepdata-validate -d . and prints TABLES=<n>, OUTDIR=<dir>, VALID=1|0.
"""
from __future__ import annotations

import argparse
import json
import math
import numbers
import os
import subprocess
import sys

from hepdata_lib import Submission, Table, Uncertainty, Variable

SQRTS_NN_GEV = 2760.0
REACTION = "PB PB --> JET JET X"
CENTRALITY_ORDER = ["50-100%", "30-50%", "20-30%", "10-20%", "0-10%"]
MERGED_ORDER = ["30-100%", "10-30%", "0-10%"]


def require_files(paths: list[str]) -> None:
    missing = [p for p in paths if not os.path.isfile(p)]
    if missing:
        print("error: missing input(s):\n  " + "\n  ".join(missing), file=sys.stderr)
        sys.exit(2)


def finite(name: str, values):
    """Refuse NaN/inf before hepdata_lib writes `.nan` into YAML."""
    flat = [x for v in values for x in (v if isinstance(v, (tuple, list)) else (v,))]
    bad = [x for x in flat if isinstance(x, numbers.Real) and not math.isfinite(float(x))]
    if bad:
        print(f"error: {name}: {len(bad)} non-finite value(s)", file=sys.stderr)
        sys.exit(2)
    return list(values)


def figure_table(name: str, description: str, location: str, observables: list[str], image: str | None = None, phrases: list[str] | None = None) -> Table:
    t = Table(name)
    t.description = description
    t.location = location
    t.keywords["observables"] = observables
    t.keywords["reactions"] = [REACTION]
    t.keywords["cmenergies"] = [SQRTS_NN_GEV]
    if phrases:
        t.keywords["phrases"] = phrases
    if image:
        t.add_image(os.path.abspath(image))
    return t


def binned_x(name: str, edges: list[float], units: str = "GeV") -> Variable:
    v = Variable(name, is_independent=True, is_binned=True, units=units)
    v.values = finite(name, [(float(lo), float(hi)) for lo, hi in zip(edges[:-1], edges[1:])])
    return v


def category_x(name: str, labels: list[str]) -> Variable:
    v = Variable(name, is_independent=True, is_binned=False)
    v.values = [str(x) for x in labels]
    return v


def y_var(
    name: str,
    values,
    units: str = "",
    qualifiers: dict[str, str] | None = None,
    zero_ok: bool = False,
    digits: int | None = None,
) -> Variable:
    all_int = len(values) > 0 and all(isinstance(x, numbers.Integral) and not isinstance(x, bool) for x in values)
    eff_digits = digits if digits is not None else (9 if all_int else 5)
    v = Variable(
        name,
        is_independent=False,
        is_binned=False,
        units=units,
        zero_uncertainties_warning=not zero_ok,
    )
    v.digits = eff_digits
    v.values = finite(name, [int(x) if isinstance(x, numbers.Integral) and not isinstance(x, bool) else float(x) for x in values])
    for k, q in (qualifiers or {}).items():
        v.add_qualifier(k, q)
    return v


def add_sym(var: Variable, label: str, errors) -> None:
    u = Uncertainty(label, is_symmetric=True)
    u.values = finite(label, [float(e) for e in errors])
    var.add_uncertainty(u)


def add_asym(var: Variable, label: str, down, up) -> None:
    u = Uncertainty(label, is_symmetric=False)
    u.values = finite(label, [(-abs(float(d)), abs(float(p))) for d, p in zip(down, up)])
    var.add_uncertainty(u)


def build(sub: Submission, args) -> None:
    base_dir = os.path.dirname(os.path.abspath(__file__))
    summary_path = os.path.join(base_dir, "summary.json")
    plots_dir = os.path.join(base_dir, "plots")
    for inp in args.input:
        if os.path.isfile(inp) and inp.endswith(".json"):
            summary_path = inp
        elif os.path.isdir(inp):
            plots_dir = inp

    expected_pngs = [
        os.path.join(plots_dir, f"{stem}.png")
        for stem in [
            "fig1a_trigger_turnon",
            "fig3a_hf_centrality",
            "fig5_leading_jet_pt",
            "fig6_dijet_dphi",
            "fig7_rb_dphi_vs_npart",
            "fig8_dijet_asymmetry_aj",
            "fig9_rb_aj_vs_npart",
            "fig10_fractional_imbalance_vs_pt1",
        ]
    ]
    require_files([summary_path] + expected_pngs)
    data = json.loads(open(summary_path, encoding="utf-8").read())

    sub.comment = (
        "CMS Open Data study replicating the 2.76 TeV PbPb calorimeter dijet momentum "
        "asymmetry and jet quenching measurements of CMS-HIN-10-004 (Phys. Rev. C 84, "
        "024906 (2011), arXiv:1102.1957) using 2010 CMS Heavy-Ion Open Data records "
        "14010 (/HICorePhysics/HIRun2010-PromptReco-v3/RECO), 14011 "
        "(/HIAllPhysics/HIRun2010-ZS-v2/RECO), and 14202 (Golden JSON). "
        "Note: This record presents an independent CMS Open Data analysis and is NOT an "
        "official CMS Collaboration result."
    )

    base_qual = {"SQRT(S)/NUCLEON": "2760 GeV", "RE": "PB PB --> JET JET X"}

    # 1. Trigger turn-on (populated bins only so no undefined efficiencies exist)
    trig = data["trigger_turnon"]
    t_edges = trig["pt_edges"]
    populated = [i for i, d in enumerate(trig["n_den"]) if d > 0]
    t1 = figure_table(
        "Trigger efficiency vs leading jet pT",
        "CMS Open Data measurement (not an official CMS result) of the single-jet "
        "HLT_HIJet50U trigger efficiency as a function of corrected leading IC5Calo "
        "jet pT (|eta| < 2.0) in 2.76 TeV PbPb collisions, with 68% Clopper-Pearson "
        "statistical intervals.",
        "Data from Open Data Figure 1a (CMS-HIN-10-004 Fig. 1a equivalent).",
        ["EFF"],
        image=os.path.join(plots_dir, "fig1a_trigger_turnon.png"),
        phrases=["Heavy Ion", "Jet Quenching", "Trigger Efficiency", "CMS Open Data"],
    )
    v_pt = Variable("Leading jet $p_{T,1}$", is_independent=True, is_binned=True, units="GeV")
    v_pt.values = finite("turnon_bins", [(float(t_edges[i]), float(t_edges[i + 1])) for i in populated])
    t1.add_variable(v_pt)
    y_eff = y_var(
        "HLT_HIJet50U efficiency",
        [trig["eff"][i] for i in populated],
        qualifiers={**base_qual, "Trigger": "HLT_HIJet50U"},
    )
    add_asym(y_eff, "stat", [trig["err_low"][i] for i in populated], [trig["err_high"][i] for i in populated])
    t1.add_variable(y_eff)
    t1.add_variable(
        y_var("Denominator events", [int(trig["n_den"][i]) for i in populated], qualifiers=base_qual, zero_ok=True)
    )
    sub.add_table(t1)

    # 2. HF energy distribution
    hf = data["hf_energy_hist"]
    t2 = figure_table(
        "HF total energy distribution",
        "CMS Open Data measurement (not an official CMS result) of the event-normalized "
        "forward hadronic calorimeter total energy sum E_HF in cleaned Minimum-Bias "
        "(HLT_HIMinBiasHfOrBSC) and jet-triggered (HLT_HIJet50U) 2.76 TeV PbPb collisions.",
        "Data from Open Data Figure 3a (CMS-HIN-10-004 Fig. 3a equivalent).",
        ["N"],
        image=os.path.join(plots_dir, "fig3a_hf_centrality.png"),
        phrases=["Heavy Ion", "Centrality", "Forward Calorimeter", "CMS Open Data"],
    )
    t2.add_variable(binned_x("$\\sum E_{\\mathrm{HF}}$", hf["edges_tev"], units="TeV"))
    mb_tot = max(sum(hf["mb_counts"]), 1)
    j50_tot = max(sum(hf["jet50_counts"]), 1)
    y_mb = y_var(
        "Minimum-Bias event fraction",
        hf["mb_fraction"],
        qualifiers={**base_qual, "Trigger": "HLT_HIMinBiasHfOrBSC"},
        zero_ok=True,
    )
    add_sym(y_mb, "stat", [math.sqrt(c) / mb_tot for c in hf["mb_counts"]])
    t2.add_variable(y_mb)
    y_j50 = y_var(
        "HLT_HIJet50U event fraction",
        hf["jet50_fraction"],
        qualifiers={**base_qual, "Trigger": "HLT_HIJet50U"},
        zero_ok=True,
    )
    add_sym(y_j50, "stat", [math.sqrt(c) / j50_tot for c in hf["jet50_counts"]])
    t2.add_variable(y_j50)
    sub.add_table(t2)

    # 3. Leading jet pT spectrum across 5 centrality classes
    pc = data["per_centrality"]
    t3 = figure_table(
        "Leading jet pT spectrum by centrality",
        "CMS Open Data measurement (not an official CMS result) of the dijet-normalized "
        "leading calorimeter jet transverse momentum spectrum (1/N_dijet) dN/dpT1 across "
        "five centrality bins for pT1 > 120 GeV, pT2 > 50 GeV, |eta| < 2.0, dphi12 > 2pi/3.",
        "Data from Open Data Figure 5 (CMS-HIN-10-004 Fig. 5 equivalent).",
        ["1/N DN/DPT"],
        image=os.path.join(plots_dir, "fig5_leading_jet_pt.png"),
        phrases=["Dijet", "Jet Quenching", "Transverse Momentum", "CMS Open Data"],
    )
    t3.add_variable(binned_x("Leading jet $p_{T,1}$", pc["0-10%"]["pt1_hist"]["edges"], units="GeV"))
    for cname in CENTRALITY_ORDER:
        ph = pc[cname]["pt1_hist"]
        yv = y_var(
            f"(1/N_dijet) dN/dpT1 [{cname}]",
            ph["density"],
            units="GeV^-1",
            qualifiers={**base_qual, "Centrality": cname},
            zero_ok=True,
        )
        add_sym(yv, "stat", ph["density_err"])
        t3.add_variable(yv)
    sub.add_table(t3)

    # 4. Dijet azimuthal correlation dphi12 across 5 centrality classes
    t4 = figure_table(
        "Dijet azimuthal correlation dphi12",
        "CMS Open Data measurement (not an official CMS result) of the event-normalized "
        "azimuthal angle correlation (1/N) dN/dDeltaPhi12 between leading (pT1 > 120 GeV) "
        "and subleading (pT2 > 50 GeV) jets in |eta| < 2.0 across five PbPb centrality bins.",
        "Data from Open Data Figure 6 (CMS-HIN-10-004 Fig. 6 equivalent).",
        ["1/N DN/DDPHI"],
        image=os.path.join(plots_dir, "fig6_dijet_dphi.png"),
        phrases=["Dijet", "Azimuthal Correlation", "Jet Quenching", "CMS Open Data"],
    )
    t4.add_variable(binned_x("$\\Delta\\phi_{12}$", pc["0-10%"]["dphi_hist"]["edges"], units="rad"))
    for cname in CENTRALITY_ORDER:
        dh = pc[cname]["dphi_hist"]
        yv = y_var(
            f"(1/N) dN/dDeltaPhi12 [{cname}]",
            dh["density"],
            units="rad^-1",
            qualifiers={**base_qual, "Centrality": cname},
            zero_ok=True,
        )
        add_sym(yv, "stat", dh["density_err"])
        t4.add_variable(yv)
    sub.add_table(t4)

    # 5. Azimuthal back-to-back fraction vs Npart
    t5 = figure_table(
        "Azimuthal back-to-back fraction vs Npart",
        "CMS Open Data measurement (not an official CMS result) of the fraction of "
        "pT1 > 120 GeV, pT2 > 50 GeV dijet candidates satisfying dphi12 > 3.026 rad "
        "as a function of collision centrality and <N_part>, with 68% Clopper-Pearson intervals.",
        "Data from Open Data Figure 7 (CMS-HIN-10-004 Fig. 7 equivalent).",
        ["RATIO"],
        image=os.path.join(plots_dir, "fig7_rb_dphi_vs_npart.png"),
        phrases=["Dijet", "Azimuthal Correlation", "Npart", "CMS Open Data"],
    )
    t5.add_variable(category_x("Centrality", CENTRALITY_ORDER))
    y_np = y_var("$\\langle N_{\\mathrm{part}} \\rangle$", [pc[c]["npart_mean"] for c in CENTRALITY_ORDER], qualifiers=base_qual)
    add_sym(y_np, "syst", [pc[c]["npart_err"] for c in CENTRALITY_ORDER])
    t5.add_variable(y_np)
    y_rb_dp = y_var(
        "$R_{B}(\\Delta\\phi_{12} > 3.026)$",
        [pc[c]["rb_dphi_gt_3026"] for c in CENTRALITY_ORDER],
        qualifiers=base_qual,
    )
    add_asym(
        y_rb_dp,
        "stat",
        [pc[c]["rb_dphi_err_low"] for c in CENTRALITY_ORDER],
        [pc[c]["rb_dphi_err_high"] for c in CENTRALITY_ORDER],
    )
    t5.add_variable(y_rb_dp)
    sub.add_table(t5)

    # 6. Dijet momentum asymmetry AJ across 5 centrality classes
    t6 = figure_table(
        "Dijet asymmetry AJ by centrality",
        "CMS Open Data measurement (not an official CMS result) of the event-normalized "
        "dijet momentum asymmetry distribution AJ = (pT1 - pT2)/(pT1 + pT2) for selected "
        "dijets (pT1 > 120 GeV, pT2 > 50 GeV, |eta| < 2.0, dphi12 > 2pi/3) across five "
        "2.76 TeV PbPb centrality classes. Symmetric Poisson sqrt(N)/N_dijet statistical uncertainties.",
        "Data from Open Data Figure 8 (CMS-HIN-10-004 Fig. 8 equivalent).",
        ["1/N DN/DAJ"],
        image=os.path.join(plots_dir, "fig8_dijet_asymmetry_aj.png"),
        phrases=["Dijet", "Jet Quenching", "Asymmetry", "PbPb", "CMS Open Data"],
    )
    t6.add_variable(binned_x("$A_{J}$", pc["0-10%"]["aj_hist"]["edges"], units=""))
    for cname in CENTRALITY_ORDER:
        ah = pc[cname]["aj_hist"]
        yv = y_var(
            f"Event fraction [{cname}]",
            ah["fraction"],
            qualifiers={**base_qual, "Centrality": cname},
            zero_ok=True,
        )
        add_sym(yv, "stat", ah["fraction_err"])
        t6.add_variable(yv)
    sub.add_table(t6)

    # 7. Balanced dijet fraction R_B(A_J < 0.15) and <A_J> vs Npart
    t7 = figure_table(
        "Balanced dijet fraction and mean AJ vs Npart",
        "CMS Open Data measurement (not an official CMS result) of the balanced dijet "
        "fraction R_B(AJ < 0.15) (normalized to all cleaned HLT_HIJet50U events with "
        "pT1 > 120 GeV, |eta1| < 2.0, including mono-jet events) and mean dijet asymmetry "
        "<AJ> across the five PbPb centrality classes, with statistical and JEC systematic uncertainties.",
        "Data from Open Data Figure 9 (CMS-HIN-10-004 Fig. 9 equivalent).",
        ["RATIO", "MEAN"],
        image=os.path.join(plots_dir, "fig9_rb_aj_vs_npart.png"),
        phrases=["Dijet", "Balanced Fraction", "Jet Quenching", "CMS Open Data"],
    )
    t7.add_variable(category_x("Centrality", CENTRALITY_ORDER))
    y_np2 = y_var("$\\langle N_{\\mathrm{part}} \\rangle$", [pc[c]["npart_mean"] for c in CENTRALITY_ORDER], qualifiers=base_qual)
    add_sym(y_np2, "syst", [pc[c]["npart_err"] for c in CENTRALITY_ORDER])
    t7.add_variable(y_np2)

    y_rb_aj = y_var("$R_{B}(A_{J} < 0.15)$", [pc[c]["rb_aj_lt_015"] for c in CENTRALITY_ORDER], qualifiers=base_qual)
    add_asym(
        y_rb_aj,
        "stat",
        [pc[c]["rb_aj_err_low"] for c in CENTRALITY_ORDER],
        [pc[c]["rb_aj_err_high"] for c in CENTRALITY_ORDER],
    )
    add_sym(y_rb_aj, "syst", [pc[c]["rb_aj_syst"] for c in CENTRALITY_ORDER])
    t7.add_variable(y_rb_aj)

    y_ajm = y_var("$\\langle A_{J} \\rangle$", [pc[c]["aj_mean"] for c in CENTRALITY_ORDER], qualifiers=base_qual)
    add_sym(y_ajm, "stat", [pc[c]["aj_mean_err"] for c in CENTRALITY_ORDER])
    add_sym(y_ajm, "syst", [pc[c]["aj_mean_syst"] for c in CENTRALITY_ORDER])
    t7.add_variable(y_ajm)
    sub.add_table(t7)

    # 8. Mean fractional imbalance <(pT1 - pT2)/pT1> vs pT1
    fi = data["fractional_imbalance_vs_pt1"]
    t8 = figure_table(
        "Mean fractional imbalance vs leading jet pT",
        "CMS Open Data measurement (not an official CMS result) of the profile mean "
        "<(pT1 - pT2)/pT1> as a function of leading jet pT1 in three merged centrality "
        "classes (30-100%, 10-30%, 0-10%) for selected dijets.",
        "Data from Open Data Figure 10 (CMS-HIN-10-004 Fig. 10 equivalent).",
        ["MEAN"],
        image=os.path.join(plots_dir, "fig10_fractional_imbalance_vs_pt1.png"),
        phrases=["Dijet", "Fractional Imbalance", "Jet Quenching", "CMS Open Data"],
    )
    t8.add_variable(binned_x("Leading jet $p_{T,1}$", fi["pt1_edges"], units="GeV"))
    for mname in MERGED_ORDER:
        yv = y_var(
            f"<(pT1 - pT2)/pT1> [{mname}]",
            fi[mname]["mean"],
            qualifiers={**base_qual, "Centrality": mname},
            zero_ok=True,
        )
        add_sym(yv, "stat", fi[mname]["sem"])
        t8.add_variable(yv)
    sub.add_table(t8)

    # 9. Event selection cutflow table
    cf = data["cutflow"]
    cf_order = [
        ("total_events_read", "Total events read in sampled HICorePhysics files"),
        ("golden_json_certified", "Golden JSON certified lumisections"),
        ("hlt_hijet50u", "HLT_HIJet50U single-jet trigger"),
        ("hf_coincidence", "Offline HF coincidence (etHFtowerSumPlus/Minus > 3 GeV)"),
        ("primary_vertex", "Reconstructed heavy-ion primary vertex (|z_vtx| <= 25 cm)"),
        ("noise_and_scraping_cleaned", "HCAL/ECAL noise and beam-scraping rejection"),
        ("two_jets_pt35_eta2", "At least two jets with pT > 35 GeV, |eta| < 2.0"),
        ("leading_jet_pt120", "Leading jet pT1 > 120 GeV, |eta1| < 2.0"),
        ("subleading_jet_pt50", "Subleading jet pT2 > 50 GeV, |eta2| < 2.0"),
        ("dijet_dphi_gt_2pi3", "Azimuthal back-to-back cut dphi12 > 2pi/3"),
    ]
    t9 = figure_table(
        "Event selection cutflow",
        "CMS Open Data event selection cutflow (not an official CMS result) across the "
        "stratified HICorePhysics RECO sample (CMS-HIN-10-004 Table 1 equivalent).",
        "Data from Open Data selection cutflow.",
        ["N"],
        phrases=["Cutflow", "Event Selection", "Heavy Ion", "CMS Open Data"],
    )
    t9.add_variable(category_x("Selection step", [label for _, label in cf_order]))
    t9.add_variable(y_var("Events", [int(cf[k]) for k, _ in cf_order], units="Events", qualifiers=base_qual, zero_ok=True))
    sub.add_table(t9)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--outdir", default="hepdata", help="output directory (replaced)")
    ap.add_argument("--input", action="append", default=[], help="input file or directory for build(), repeatable")
    ap.add_argument("--no-tar", action="store_true", help="delete the submission.tar.gz create_files() writes")
    args = ap.parse_args(argv)
    args.input = [os.path.abspath(x) for x in args.input]

    sub = Submission()
    build(sub, args)
    outdir = os.path.abspath(args.outdir)
    os.makedirs(os.path.dirname(outdir), exist_ok=True)
    cwd = os.getcwd()
    os.chdir(os.path.dirname(outdir))
    try:
        sub.create_files(os.path.basename(outdir), remove_old=True)
    except AssertionError as e:
        print(f"error: hepdata_lib: {e}", file=sys.stderr)
    finally:
        os.chdir(cwd)
    tarball = os.path.join(os.path.dirname(outdir), "submission.tar.gz")
    if args.no_tar and os.path.exists(tarball):
        os.remove(tarball)
    res = subprocess.run(["hepdata-validate", "-d", "."], cwd=outdir, capture_output=True, text=True)
    if res.stdout.strip():
        print(res.stdout.strip(), file=sys.stderr)
    if res.stderr.strip():
        print(res.stderr.strip(), file=sys.stderr)
    print(f"TABLES={len(sub.tables)}")
    print(f"OUTDIR={args.outdir}")
    print(f"VALID={int(res.returncode == 0)}")
    return 0 if res.returncode == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
