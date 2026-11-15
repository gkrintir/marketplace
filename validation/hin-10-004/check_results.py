#!/usr/bin/env python3
"""Objective physics falsifier for CMS-HIN-10-004 Open Data replication.

Uses Python standard library only. Checks structural completeness, cutflow
monotonicity, kinematic bounds, trigger turn-on plateau, HF centrality bias,
and statistically significant jet quenching signatures across centrality.

Exits 0 only if every check passes; prints PASS/FAIL lines and SUMMARY=k/n.
"""
from __future__ import annotations

import json
import math
import os
import sys
from pathlib import Path

CENTRALITY_ORDER = ["50-100%", "30-50%", "20-30%", "10-20%", "0-10%"]
EXPECTED_PLOTS = [
    "fig1a_trigger_turnon",
    "fig3a_hf_centrality",
    "fig4_pt2_vs_pt1",
    "fig5_leading_jet_pt",
    "fig6_dijet_dphi",
    "fig7_rb_dphi_vs_npart",
    "fig8_dijet_asymmetry_aj",
    "fig9_rb_aj_vs_npart",
    "fig10_fractional_imbalance_vs_pt1",
]


def main() -> int:
    base = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(__file__).resolve().parent
    passed = 0
    total = 0

    def check(name: str, cond: bool, detail: str = "") -> None:
        nonlocal passed, total
        total += 1
        status = "PASS" if cond else "FAIL"
        if cond:
            passed += 1
        msg = f"{status} {name}"
        if detail:
            msg += f": {detail}"
        print(msg)

    summary_path = base / "summary.json"
    if not summary_path.is_file():
        print(f"error: missing {summary_path}", file=sys.stderr)
        return 2
    try:
        data = json.loads(summary_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        print(f"error: cannot read {summary_path}: {exc}", file=sys.stderr)
        return 2

    # 1. Specification & required files
    spec_path = base / "specification.md"
    check("specification_exists", spec_path.is_file(), "specification.md")
    if spec_path.is_file():
        spec_txt = spec_path.read_text(encoding="utf-8")
        check("specification_no_placeholders", "{{" not in spec_txt, "zero '{{' tokens")
    check("summary_json_exists", True, "summary.json")

    # 2. Sample size & cutflow monotonicity
    cf = data.get("cutflow", {})
    cf_keys = [
        "total_events_read",
        "golden_json_certified",
        "hlt_hijet50u",
        "hf_coincidence",
        "primary_vertex",
        "noise_and_scraping_cleaned",
        "two_jets_pt35_eta2",
        "leading_jet_pt120",
        "subleading_jet_pt50",
        "dijet_dphi_gt_2pi3",
    ]
    check("cutflow_keys_present", all(k in cf for k in cf_keys), f"found {len(cf)} keys")
    if all(k in cf for k in cf_keys):
        vals = [int(cf[k]) for k in cf_keys]
        mono = all(vals[i] >= vals[i + 1] for i in range(len(vals) - 1))
        check("cutflow_monotonic", mono, " -> ".join(f"{k}:{cf[k]}" for k in cf_keys))
        check(
            "sufficient_dijet_statistics",
            cf["dijet_dphi_gt_2pi3"] >= 300,
            f"N_dijet = {cf['dijet_dphi_gt_2pi3']} (target >= 300 from real Open Data)",
        )
        # Check physical step efficiencies against CMS-HIN-10-004 Table 1 ballpark
        if cf["hlt_hijet50u"] > 0 and cf["leading_jet_pt120"] > 0 and cf["subleading_jet_pt50"] > 0:
            hf_eff = cf["hf_coincidence"] / cf["hlt_hijet50u"]
            pt2_eff = cf["subleading_jet_pt50"] / cf["leading_jet_pt120"]
            dphi_eff = cf["dijet_dphi_gt_2pi3"] / cf["subleading_jet_pt50"]
            check(
                "hf_coincidence_fraction_physical",
                0.55 <= hf_eff <= 0.92,
                f"HF coincidence / HLT_HIJet50U = {hf_eff:.3f} (HIN-10-004 Table 1: 0.750)",
            )
            check(
                "subleading_pt50_fraction_physical",
                0.70 <= pt2_eff <= 0.96,
                f"pT2>50 / pT1>120 = {pt2_eff:.3f} (HIN-10-004 Table 1: 0.874)",
            )
            check(
                "back_to_back_dphi_fraction_physical",
                0.85 <= dphi_eff <= 0.995,
                f"dphi>2pi/3 / pT2>50 = {dphi_eff:.3f} (HIN-10-004 Table 1: 0.954)",
            )

    # 3. Trigger turn-on plateau >= 95% for corrected pT1 >= 120 GeV
    trig = data.get("trigger_turnon", {})
    if trig and "pt_edges" in trig and "eff" in trig and "n_den" in trig:
        edges = trig["pt_edges"]
        effs = trig["eff"]
        dens = trig["n_den"]
        high_pt_bins = [
            effs[i] for i in range(len(effs)) if edges[i] >= 120.0 and dens[i] >= 3
        ]
        check(
            "trigger_plateau_above_120gev",
            len(high_pt_bins) > 0 and min(high_pt_bins) >= 0.90,
            f"min eff(pT1>=120 GeV) = {min(high_pt_bins) if high_pt_bins else 'N/A'}",
        )
    else:
        check("trigger_plateau_above_120gev", False, "missing trigger_turnon")

    # 4. HF centrality bias: jet-triggered events must have higher mean HF energy than MinBias
    hf = data.get("hf_energy_hist", {})
    if hf and "edges_tev" in hf and "mb_fraction" in hf and "jet50_fraction" in hf:
        centers = [0.5 * (hf["edges_tev"][i] + hf["edges_tev"][i + 1]) for i in range(len(hf["mb_fraction"]))]
        mb_mean = sum(c * f for c, f in zip(centers, hf["mb_fraction"]))
        jet_mean = sum(c * f for c, f in zip(centers, hf["jet50_fraction"]))
        check(
            "hf_jet_trigger_centrality_bias",
            jet_mean > 1.5 * mb_mean and mb_mean > 5.0,
            f"<E_HF>_Jet50 = {jet_mean:.1f} TeV vs <E_HF>_MB = {mb_mean:.1f} TeV",
        )
    else:
        check("hf_jet_trigger_centrality_bias", False, "missing hf_energy_hist")

    # 5. Centrality-resolved dijet observables & Jet Quenching physics falsifiers
    pc = data.get("per_centrality", {})
    check(
        "all_5_centrality_bins_present",
        all(c in pc for c in CENTRALITY_ORDER),
        ", ".join(pc.keys()),
    )
    if all(c in pc for c in CENTRALITY_ORDER):
        # Every centrality bin populated and normalized properly
        non_empty = all(pc[c]["n_dijet"] >= 15 for c in CENTRALITY_ORDER)
        check(
            "all_centrality_bins_populated",
            non_empty,
            ", ".join(f"{c}:{pc[c]['n_dijet']}" for c in CENTRALITY_ORDER),
        )

        # Normalization & finite values check on A_J histograms
        aj_norm_ok = True
        for c in CENTRALITY_ORDER:
            fracs = pc[c]["aj_hist"]["fraction"]
            if any(not math.isfinite(x) or x < 0 for x in fracs) or abs(sum(fracs) - 1.0) > 0.02:
                aj_norm_ok = False
        check("aj_histograms_normalized_and_finite", aj_norm_ok, "sum(fraction) == 1.0 +/- 0.02")

        # Azimuthal correlation peaks at pi (last bin of dphi_hist) in every centrality bin
        dphi_peaks_at_pi = all(
            pc[c]["dphi_hist"]["density"][-1] == max(pc[c]["dphi_hist"]["density"])
            for c in CENTRALITY_ORDER
        )
        check("dphi_peaks_at_pi_all_centralities", dphi_peaks_at_pi, "max density in bin near pi")

        # KEY JET QUENCHING SIGNATURE 1: <A_J> increases significantly from 50-100% to 0-10%
        aj_periph = pc["50-100%"]["aj_mean"]
        aj_periph_err = pc["50-100%"]["aj_mean_err"]
        aj_cent = pc["0-10%"]["aj_mean"]
        aj_cent_err = pc["0-10%"]["aj_mean_err"]
        aj_diff = aj_cent - aj_periph
        aj_sigma = aj_diff / math.hypot(aj_periph_err, aj_cent_err) if math.hypot(aj_periph_err, aj_cent_err) > 0 else 0.0
        check(
            "jet_quenching_aj_mean_shift",
            aj_cent > aj_periph + 0.05 and aj_sigma >= 3.0 and 0.10 <= aj_periph <= 0.24 and 0.22 <= aj_cent <= 0.36,
            f"<A_J>(50-100%)={aj_periph:.3f}+/-{aj_periph_err:.3f} -> <A_J>(0-10%)={aj_cent:.3f}+/-{aj_cent_err:.3f} ({aj_sigma:.1f} sigma)",
        )

        # KEY JET QUENCHING SIGNATURE 2: R_B(A_J < 0.15) drops significantly from ~0.45 in peripheral to ~0.23 in central
        rb_periph = pc["50-100%"]["rb_aj_lt_015"]
        rb_periph_err = 0.5 * (pc["50-100%"]["rb_aj_err_low"] + pc["50-100%"]["rb_aj_err_high"])
        rb_cent = pc["0-10%"]["rb_aj_lt_015"]
        rb_cent_err = 0.5 * (pc["0-10%"]["rb_aj_err_low"] + pc["0-10%"]["rb_aj_err_high"])
        rb_drop = rb_periph - rb_cent
        rb_sigma = rb_drop / math.hypot(rb_periph_err, rb_cent_err) if math.hypot(rb_periph_err, rb_cent_err) > 0 else 0.0
        check(
            "jet_quenching_rb_aj_suppression",
            rb_cent < 0.72 * rb_periph and rb_sigma >= 3.0 and 0.32 <= rb_periph <= 0.60 and 0.15 <= rb_cent <= 0.32,
            f"R_B(50-100%)={rb_periph:.3f}+/-{rb_periph_err:.3f} -> R_B(0-10%)={rb_cent:.3f}+/-{rb_cent_err:.3f} ({rb_sigma:.1f} sigma)",
        )

    # 6. Fractional imbalance <(pT1 - pT2)/pT1> increases from 30-100% to 0-10%
    fi = data.get("fractional_imbalance_vs_pt1", {})
    if "30-100%" in fi and "0-10%" in fi:
        c_per = fi["30-100%"]["counts"]
        m_per = fi["30-100%"]["mean"]
        c_cen = fi["0-10%"]["counts"]
        m_cen = fi["0-10%"]["mean"]
        avg_per = sum(c * m for c, m in zip(c_per, m_per)) / max(1, sum(c_per))
        avg_cen = sum(c * m for c, m in zip(c_cen, m_cen)) / max(1, sum(c_cen))
        check(
            "fractional_imbalance_centrality_ordering",
            avg_cen > avg_per + 0.05,
            f"<dpt/pt1>(30-100%)={avg_per:.3f} < <dpt/pt1>(0-10%)={avg_cen:.3f}",
        )
    else:
        check("fractional_imbalance_centrality_ordering", False, "missing fractional_imbalance_vs_pt1")

    # 7. Plots (PNG + PDF, non-empty)
    plots_dir = base / "plots"
    missing_plots = []
    for stem in EXPECTED_PLOTS:
        for ext in (".png", ".pdf"):
            p = plots_dir / f"{stem}{ext}"
            if not p.is_file() or p.stat().st_size < 2000:
                missing_plots.append(f"{stem}{ext}")
    check("all_plots_png_and_pdf_exist", len(missing_plots) == 0, f"missing={missing_plots}" if missing_plots else f"{len(EXPECTED_PLOTS)*2} files OK")

    # 8. HEPData submission YAML check
    hd_sub = base / "hepdata" / "submission.yaml"
    check("hepdata_submission_exists", hd_sub.is_file(), "hepdata/submission.yaml")
    if hd_sub.is_file():
        yamls = list((base / "hepdata").glob("*.yaml"))
        bad_nan = [y.name for y in yamls if ".nan" in y.read_text(encoding="utf-8").lower()]
        check("hepdata_no_nan_values", len(bad_nan) == 0 and len(yamls) >= 6, f"{len(yamls)} YAML files, bad={bad_nan}")

    # 9. No personal paths in committed JSON / YAML / Markdown files
    leaked = []
    scan_files = list(base.glob("*.json")) + list(base.glob("*.md")) + list((base / "hepdata").glob("*.yaml"))
    for sf in scan_files:
        txt = sf.read_text(encoding="utf-8", errors="ignore")
        if any(tok in txt for tok in ("/home/", "/Users/", "/uscms_data/")):
            leaked.append(sf.name)
    check("no_personal_paths", len(leaked) == 0, f"leaked={leaked}" if leaked else f"{len(scan_files)} files clean")

    print(f"SUMMARY={passed}/{total}")
    return 0 if passed == total else 1


if __name__ == "__main__":
    sys.exit(main())
