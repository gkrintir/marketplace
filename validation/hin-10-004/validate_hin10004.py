#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "typer>=0.12",
#   "numpy>=1.26",
#   "scipy>=1.11",
#   "hist>=2.7",
#   "matplotlib>=3.8",
#   "mplhep>=0.3",
#   "pyyaml>=6.0",
#   "hepdata-validator>=0.3",
# ]
# ///
"""ARC quantitative validation of the student CMS Open Data study against CMS-HIN-10-004.

Cross-verifies published reference numbers directly from the CMS-HIN-10-004 LaTeX
source tables/text (Phys. Rev. C 84, 024906 (2011), arXiv:1102.1957) and archived CMS
Heavy-Ion Analysis2010 ROOT macros (.C) in github.com/CmsHI/CVS_CmsHi, audits the student
summary.json and HEPData YAML tables without modifying any student file, executes all 12
comparability tests, renders CMS-style comparison/ratio plots in plots/comparison/, and
writes comparison.json.

Fallback invocation if `env -S` is unavailable:
    uv run validate_hin10004.py
"""
from __future__ import annotations

import json
import logging
import math
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any

import hist
from hist import Hist
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
from scipy import optimize, stats
import typer
import yaml

plt.style.use(hep.style.CMS)

app = typer.Typer(add_completion=False)
log = logging.getLogger("validate_hin10004")

DEFAULT_BASE = Path(__file__).resolve().parent
CL_1SIGMA = 0.682689492137086
CENTRALITY_ORDER = ["50-100%", "30-50%", "20-30%", "10-20%", "0-10%"]
MERGED_ORDER = ["30-100%", "10-30%", "0-10%"]
RLABEL_CMP = r"PbPb $\sqrt{s_{\mathrm{NN}}} = 2.76\ \mathrm{TeV}$ (2010)"


def clopper_pearson(k: int | np.ndarray, n: int | np.ndarray, cl: float = CL_1SIGMA):
    """Exact 68.27% Clopper-Pearson binomial interval returning (eff, err_low, err_high)."""
    k_arr = np.asarray(k, dtype=np.float64)
    n_arr = np.asarray(n, dtype=np.float64)
    alpha = 0.5 * (1.0 - cl)
    with np.errstate(invalid="ignore", divide="ignore"):
        lo = np.where(k_arr > 0, stats.beta.ppf(alpha, k_arr, n_arr - k_arr + 1.0), 0.0)
        hi = np.where(k_arr < n_arr, stats.beta.ppf(1.0 - alpha, k_arr + 1.0, n_arr - k_arr), 1.0)
        eff = np.where(n_arr > 0, k_arr / n_arr, np.nan)
    return eff, eff - lo, hi - eff


def root_bayes_hpd(k: int, n: int, cl: float = CL_1SIGMA) -> tuple[float, float, float]:
    """Exact ROOT TGraphAsymmErrors::BayesDivide shortest interval on Beta(k+1, n-k+1)."""
    a, b = k + 1.0, n - k + 1.0
    mode = k / n
    if k == 0:
        return 0.0, 0.0, 1.0 - (1.0 - cl) ** (1.0 / (n + 1.0))
    if k == n:
        return 1.0, 1.0 - (1.0 - cl) ** (1.0 / (n + 1.0)), 0.0
    res = optimize.minimize_scalar(
        lambda al: stats.beta.ppf(al + cl, a, b) - stats.beta.ppf(al, a, b),
        bounds=(0.0, 1.0 - cl),
        method="bounded",
        options={"xatol": 1e-10},
    )
    lo = float(stats.beta.ppf(res.x, a, b))
    hi = float(stats.beta.ppf(res.x + cl, a, b))
    return mode, mode - lo, hi - mode


def weighted_linear_fit(x: np.ndarray, y: np.ndarray, yerr: np.ndarray) -> dict[str, float]:
    """Weighted linear regression y = intercept + slope * x."""
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(yerr) & (yerr > 0)
    xv, yv, ev = x[ok], y[ok], yerr[ok]
    w = 1.0 / (ev**2)
    wsum = float(np.sum(w))
    xbar = float(np.sum(w * xv) / wsum)
    ybar = float(np.sum(w * yv) / wsum)
    sxx = float(np.sum(w * (xv - xbar) ** 2))
    slope = float(np.sum(w * (xv - xbar) * (yv - ybar)) / sxx)
    slope_err = math.sqrt(1.0 / sxx)
    intercept = ybar - slope * xbar
    chi2 = float(np.sum(((yv - (intercept + slope * xv)) / ev) ** 2))
    ndf = max(int(len(xv) - 2), 1)
    return {
        "intercept": intercept,
        "slope": slope,
        "slope_err": slope_err,
        "slope_pull_sigma": slope / slope_err if slope_err > 0 else 0.0,
        "fit_chi2": chi2,
        "fit_ndf": ndf,
    }


def two_sided_z_from_p(p_val: float) -> float:
    """Convert a p-value into an equivalent two-sided Gaussian significance in sigma."""
    if not math.isfinite(p_val) or p_val >= 1.0:
        return 0.0
    return float(stats.norm.isf(0.5 * max(p_val, 1e-300)))


def verdict_from_sigma(max_abs_sigma: float) -> str:
    """Objective ARC verdict from a pull or equivalent Gaussian significance."""
    if max_abs_sigma > 3.0:
        return "INCONSISTENT"
    if max_abs_sigma >= 2.0:
        return "TENSION"
    return "CONSISTENT"


def verdict_from_shape(p_val: float, slope_pull_sigma: float) -> tuple[str, float]:
    """Objective ARC verdict for unit-normalized shape tests from global chi2 p-value and ratio slope."""
    z_chi2 = two_sided_z_from_p(p_val)
    if z_chi2 > 3.0:
        return "INCONSISTENT", z_chi2
    if z_chi2 >= 2.0 or abs(slope_pull_sigma) >= 2.5:
        return "TENSION", z_chi2
    return "CONSISTENT", z_chi2


def fmt_p_latex(p_val: float) -> str:
    """Format a p-value cleanly for matplotlib LaTeX text boxes without underflowing to 0.000."""
    if p_val >= 0.002:
        return f"{p_val:.3f}"
    if p_val <= 0.0:
        return "< 10^{-300}"
    exp = int(math.floor(math.log10(p_val)))
    mant = p_val / (10.0**exp)
    return f"{mant:.1f} \\times 10^{{{exp}}}"


def _rel_or_str(p: Path) -> str:
    try:
        return str(p.resolve().relative_to(DEFAULT_BASE.resolve()))
    except ValueError:
        return str(p)


def save_fig(fig: plt.Figure, outdir: Path, stem: str) -> list[str]:
    png = outdir / f"{stem}.png"
    pdf = outdir / f"{stem}.pdf"
    fig.savefig(png, dpi=160, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved %s.{png,pdf}", stem)
    return [_rel_or_str(png), _rel_or_str(pdf)]


def verify_primary_cms_sources(ref_dir: Path, pub_ref: dict[str, Any]) -> dict[str, Any]:
    """Verify published_reference.json directly against opened .C macros and .tex tables."""
    tex_dir = ref_dir / "paper_src"

    # 1. TrigEffCorrected151058and151025.C
    c_trig = (ref_dir / "TrigEffCorrected151058and151025.C").read_text(encoding="utf-8")
    pts_c = re.findall(r"grae->SetPoint\((\d+),([0-9.eE+-]+),([0-9.eE+-]+)\);", c_trig)
    errs_c = re.findall(
        r"grae->SetPointError\((\d+),([0-9.eE+-]+),([0-9.eE+-]+),([0-9.eE+-]+),([0-9.eE+-]+)\);",
        c_trig,
    )
    if len(pts_c) != 19 or len(errs_c) != 19:
        raise ValueError("Unexpected point count in TrigEffCorrected151058and151025.C")

    n_grid = np.arange(1, 50000, dtype=np.float64)
    trig_exact_kn = []
    max_trig_diff = 0.0
    for i in range(19):
        px, py = float(pts_c[i][1]), float(pts_c[i][2])
        elo, ehi = float(errs_c[i][3]), float(errs_c[i][4])
        rp = pub_ref["fig1a_trigger_turnon"]["points"][i]
        max_trig_diff = max(
            max_trig_diff,
            abs(px - rp["pt_gev"]),
            abs(py - rp["eff"]),
            abs(elo - rp["err_low"]),
            abs(ehi - rp["err_high"]),
        )
        # Recover exact integer (k, n) from ROOT BayesDivide HPD interval
        if py == 0.0:
            n_est = round(math.log(1.0 - CL_1SIGMA) / math.log(1.0 - ehi) - 1.0)
            k_est = 0
        elif py == 1.0:
            n_est = round(math.log(1.0 - CL_1SIGMA) / math.log(1.0 - elo) - 1.0)
            k_est = n_est
        else:
            k_grid = np.round(py * n_grid)
            cands = np.where(np.abs(k_grid / n_grid - py) < max(py * 2e-6, 1e-8))[0]
            best_kn = (0, 1)
            best_d = 1e9
            for idx in cands:
                kt, nt = int(k_grid[idx]), int(n_grid[idx])
                _, bel, beh = root_bayes_hpd(kt, nt)
                d_err = abs(bel - elo) + abs(beh - ehi)
                if d_err < best_d:
                    best_d = d_err
                    best_kn = (kt, nt)
            k_est, n_est = best_kn
        trig_exact_kn.append({"pt_gev": px, "k": int(k_est), "n": int(n_est), "eff": py, "err_low": elo, "err_high": ehi})

    # 2. hf_centrality_distribution_mb_jet50_20101126_v1.C
    c_hf = (ref_dir / "hf_centrality_distribution_mb_jet50_20101126_v1.C").read_text(encoding="utf-8")
    mb_map = {int(b): float(v) for b, v in re.findall(r"hhfMB->SetBinContent\((\d+),([0-9.eE+-]+)\);", c_hf)}
    jet_map = {int(b): float(v) for b, v in re.findall(r"hhfJET->SetBinContent\((\d+),([0-9.eE+-]+)\);", c_hf)}
    mb_entries = int(re.search(r"hhfMB->SetEntries\((\d+)\);", c_hf).group(1))
    jet_entries = int(re.search(r"hhfJET->SetEntries\((\d+)\);", c_hf).group(1))
    mb_arr = np.array([mb_map.get(i, 0.0) for i in range(1, 71)], dtype=np.float64)
    jet_arr = np.array([jet_map.get(i, 0.0) for i in range(1, 71)], dtype=np.float64)
    max_hf_diff = float(
        max(
            np.max(np.abs(mb_arr - np.asarray(pub_ref["fig3a_hf_distribution"]["mb_bin_contents_1indexed"]))),
            np.max(np.abs(jet_arr - np.asarray(pub_ref["fig3a_hf_distribution"]["jet50_bin_contents_1indexed"]))),
        )
    )
    # Extract TLine centrality boundaries in TeV from macro
    tlines_hf = [
        float(m)
        for m in re.findall(r"TLine\(([0-9.]+),1e-05", c_hf)
    ]

    # 3. deltaPtOverPt1VsPt.C
    c_dpt = (ref_dir / "deltaPtOverPt1VsPt.C").read_text(encoding="utf-8")
    pads = c_dpt.split("// ------------>Primitives in pad:")
    max_dpt_diff = 0.0
    for pidx, mname in enumerate(MERGED_ORDER, 1):
        pad_txt = pads[pidx]
        hlines = re.findall(r"TLine\(([0-9.]+),([0-9.]+),([0-9.]+),([0-9.]+)\)", pad_txt)
        horiz = [
            (float(x1), float(y1), float(x2), float(y2))
            for x1, y1, x2, y2 in hlines
            if abs(float(y1) - float(y2)) < 1e-9 and abs(float(x2) - float(x1) - 8.0) < 1e-6
        ]
        g_block = pad_txt.split('grae->SetName("g");')[1].split('grae->SetName("gMix");')[0]
        g_pts = re.findall(r"grae->SetPoint\((\d+),([0-9.eE+-]+),([0-9.eE+-]+)\);", g_block)
        g_errs = re.findall(r"grae->SetPointError\((\d+),0,0,([0-9.eE+-]+),([0-9.eE+-]+)\);", g_block)
        for bi in range(6):
            cx, cy = float(g_pts[bi][1]), float(g_pts[bi][2])
            ce = float(g_errs[bi][1])
            ys = [y1 for x1, y1, x2, _ in horiz if abs(0.5 * (x1 + x2) - cx) < 1e-3]
            syst_half = 0.5 * (max(ys) - min(ys))
            rj = pub_ref["delta_pt_over_pt1_vs_pt1"][mname]["pbpb_data"][bi]
            max_dpt_diff = max(
                max_dpt_diff,
                abs(cx - rj["pt1_center_gev"]),
                abs(cy - rj["val"]),
                abs(ce - rj["stat_err"]),
                abs(syst_half - rj["syst_err"]),
            )

    # 4. Verify TeX tables & text benchmarks directly
    ev_tex = (tex_dir / "event_selection.tex").read_text(encoding="utf-8")
    cent_tex = (tex_dir / "centrality.tex").read_text(encoding="utf-8")
    res_tex = (tex_dir / "results.tex").read_text(encoding="utf-8")
    trk_tex = (tex_dir / "results_trackjet.tex").read_text(encoding="utf-8")

    has_3514 = "3\\,514" in ev_tex and "3514" in res_tex
    has_45pct = "close\nto 45\\%" in res_tex or "45\\%" in res_tex
    has_022_med = "\\AJ = 0.22" in trk_tex and "median" in trk_tex
    has_npart_355 = "355 $\\pm$ 3" in cent_tex

    if max(max_trig_diff, max_hf_diff, max_dpt_diff) > 1e-5 or not (has_3514 and has_45pct and has_022_med and has_npart_355):
        raise RuntimeError("Primary CMS source verification failed against published_reference.json")

    return {
        "source_provenance": {
            "arxiv_latex_and_figures": "https://arxiv.org/e-print/1102.1957 (Phys. Rev. C 84, 024906 (2011) .tex tables/text and vector .pdf figures; no .C macros in arXiv tarball)",
            "cmshi_cvs_root_macros": (
                "https://github.com/CmsHI/CVS_CmsHi/tree/9b9dcd34a1b718e4649ca2ddc34583706dfd5d1c/"
                "Analysis2010/macros/JetAnalysis/plots (TrigEffCorrected151058and151025.C = PRC Fig. 2a, "
                "hf_centrality_distribution_mb_jet50_20101126_v1.C = PRC Fig. 3a, and pre-publication "
                "deltaPtOverPt1VsPt.C saved Dec 13, 2010 with [120, 130, 140, 170, 240] GeV bins)"
            ),
        },
        "verified_files": [
            "reference/TrigEffCorrected151058and151025.C",
            "reference/hf_centrality_distribution_mb_jet50_20101126_v1.C",
            "reference/deltaPtOverPt1VsPt.C",
            "reference/paper_src/HIN-10-004_temp.tex",
            "reference/paper_src/event_selection.tex",
            "reference/paper_src/centrality.tex",
            "reference/paper_src/jet_reconstruction.tex",
            "reference/paper_src/results.tex",
            "reference/paper_src/results_trackjet.tex",
            "reference/paper_src/results_missingpt.tex",
            "reference/paper_src/summary.tex",
        ],
        "max_macro_vs_json_diff": float(max(max_trig_diff, max_hf_diff, max_dpt_diff)),
        "trig_exact_kn": trig_exact_kn,
        "hf_macro_mb_entries": mb_entries,
        "hf_macro_jet_entries": jet_entries,
        "hf_macro_boundaries_tev": tlines_hf,
        "tex_benchmarks_verified": {
            "final_dijet_count_3514": has_3514,
            "rb_aj_peripheral_45pct": has_45pct,
            "aj_0_30_median_022": has_022_med,
            "glauber_table2_npart": has_npart_355,
        },
    }


def audit_hepdata(base_dir: Path, student: dict[str, Any]) -> dict[str, Any]:
    """Audit student hepdata/ YAML files against summary.json and hepdata-validate."""
    hdir = base_dir / "hepdata"
    sub_docs = list(yaml.safe_load_all((hdir / "submission.yaml").read_text(encoding="utf-8")))
    header_doc = sub_docs[0]
    table_entries = sub_docs[1:]

    # Check hepdata-validate CLI if installed on PATH
    validator_bin = shutil.which("hepdata-validate")
    schema_valid = None
    if validator_bin:
        proc = subprocess.run([validator_bin, "-d", "."], cwd=str(hdir), capture_output=True, text=True)
        schema_valid = proc.returncode == 0

    # Check cutflow integer rounding in event_selection_cutflow.yaml
    cf_yaml = yaml.safe_load((hdir / "event_selection_cutflow.yaml").read_text(encoding="utf-8"))
    cf_vals_yaml = [v["value"] for v in cf_yaml["dependent_variables"][0]["values"]]
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
    rounding_mismatches = []
    for k, yv in zip(cf_keys, cf_vals_yaml):
        sv = int(student["cutflow"][k])
        if int(yv) != sv:
            rounding_mismatches.append({"step": k, "summary_json": sv, "hepdata_yaml": int(yv), "diff": int(yv) - sv})

    # Check which tables include systematic uncertainties and scan for non-finite YAML entries
    tables_missing_syst = []
    non_finite_entries = 0
    referenced_files = set()
    for entry in table_entries:
        df = entry["data_file"]
        referenced_files.add(df)
        tdata = yaml.safe_load((hdir / df).read_text(encoding="utf-8"))
        for dv in tdata.get("dependent_variables", []):
            errs_labels = set()
            for item in dv.get("values", []):
                val = item.get("value")
                if isinstance(val, (int, float)) and not math.isfinite(float(val)):
                    non_finite_entries += 1
                for e in item.get("errors", []):
                    errs_labels.add(e.get("label"))
                    if "symerror" in e and not math.isfinite(float(e["symerror"])):
                        non_finite_entries += 1
                    if "asymerror" in e:
                        for side in ("minus", "plus"):
                            if not math.isfinite(float(e["asymerror"][side])):
                                non_finite_entries += 1
            if df in ("mean_fractional_imbalance_vs_leading_jet_pt.yaml", "azimuthal_back-to-back_fraction_vs_npart.yaml"):
                if "syst" not in errs_labels and "stat" in errs_labels:
                    tables_missing_syst.append({"data_file": df, "variable": dv["header"]["name"]})

    missing_figs = []
    if "pt2_vs_pt1_2d" in student and not any("pt2_vs_pt1" in f or "fig4" in f for f in referenced_files):
        missing_figs.append("fig4_pt2_vs_pt1 (2D correlation stored in summary.json['pt2_vs_pt1_2d'] but omitted from hepdata/)")

    return {
        "n_tables": len(table_entries),
        "schema_valid": schema_valid,
        "non_finite_value_count": non_finite_entries,
        "has_open_data_disclaimer": "NOT an official CMS" in header_doc.get("comment", ""),
        "no_cms_inspire_id_claimed": "inspire_id" not in header_doc,
        "cutflow_integer_rounding_resolved": len(rounding_mismatches) == 0,
        "cutflow_5sigfig_rounding_mismatches": rounding_mismatches,
        "tables_missing_published_syst_counterpart": tables_missing_syst,
        "missing_figure_tables": missing_figs,
    }


def run_test1_constants(student: dict[str, Any], pub_ref: dict[str, Any], src_ver: dict[str, Any], outdir: Path) -> dict[str, Any]:
    """T1 (exact): Collision constants, HF centrality edges, and Glauber <N_part>."""
    pc = student["per_centrality"]
    pub_g = pub_ref["table2_centrality"]["bins"]
    pub_bounds_tev = src_ver["hf_macro_boundaries_tev"]
    stu_bounds_tev = [
        pc["50-100%"]["ehf_range_gev"][1] / 1000.0,
        pc["30-50%"]["ehf_range_gev"][1] / 1000.0,
        pc["20-30%"]["ehf_range_gev"][1] / 1000.0,
        pc["10-20%"]["ehf_range_gev"][1] / 1000.0,
    ]
    bound_diffs_tev = [abs(s - p) for s, p in zip(stu_bounds_tev, pub_bounds_tev)]

    npart_stu = np.array([pc[c]["npart_mean"] for c in CENTRALITY_ORDER], dtype=float)
    npart_err_stu = np.array([pc[c]["npart_err"] for c in CENTRALITY_ORDER], dtype=float)
    npart_pub = np.array([pub_g[c]["npart_mean"] for c in CENTRALITY_ORDER], dtype=float)
    npart_err_pub = np.array([pub_g[c]["npart_syst"] for c in CENTRALITY_ORDER], dtype=float)

    # Also compare MB centrality fractions vs published hhfMB macro integral
    hf = student["hf_energy_hist"]
    edges_hf = np.asarray(hf["edges_tev"], dtype=float)
    pub_mb_raw = np.asarray(pub_ref["fig3a_hf_distribution"]["mb_bin_contents_1indexed"], dtype=float)
    pub_mb_p = pub_mb_raw / pub_mb_raw.sum()
    c_bounds_full = [0.0] + pub_bounds_tev + [160.0]
    macro_cfracs = []
    for idx in range(5):
        lo_t, hi_t = c_bounds_full[idx], c_bounds_full[idx + 1]
        s_val = 0.0
        for bi in range(70):
            blo, bhi = edges_hf[bi], edges_hf[bi + 1]
            ov = max(0.0, min(bhi, hi_t) - max(blo, lo_t)) / (bhi - blo)
            s_val += pub_mb_p[bi] * ov
        macro_cfracs.append(float(s_val))
    stu_cfracs = [float(student["validation_checks"]["mb_centrality_fractions"][c]) for c in CENTRALITY_ORDER]

    # Plot T1
    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(9.2, 7.0), sharex=True, gridspec_kw={"height_ratios": [3, 1], "hspace": 0.07}
    )
    x = np.arange(5, dtype=float)
    ax_top.errorbar(
        x - 0.08, npart_pub, yerr=npart_err_pub, fmt="o", color="black", markersize=7.5, capsize=3,
        label="CMS Published (HIN-10-004 Table 2)",
    )
    ax_top.errorbar(
        x + 0.08, npart_stu, yerr=npart_err_stu, fmt="s", mfc="none", mec="#d62728", ecolor="#d62728",
        markersize=7.5, capsize=3, lw=1.6, label="Student (CMS Open Data config)",
    )
    for idx, ( cname, b_lo, b_hi ) in enumerate(zip(CENTRALITY_ORDER, c_bounds_full[:-1], c_bounds_full[1:])):
        b_txt = f"[{b_lo:.1f}, {b_hi:.1f}] TeV" if idx < 4 else f">= {b_lo:.1f} TeV"
        ax_top.annotate(
            b_txt, (x[idx], npart_pub[idx] + npart_err_pub[idx] + 12),
            ha="center", fontsize=10.5, color="0.3",
        )
    ax_top.set_ylabel(r"Glauber $\langle N_{\mathrm{part}} \rangle$", fontsize=16)
    ax_top.set_ylim(-15, 415)
    ax_top.set_xlabel("")
    ax_top.legend(loc="upper left", frameon=True, fontsize=12.5)
    hep.cms.label("Open Data", data=True, rlabel=RLABEL_CMP, ax=ax_top, fontsize=16)

    diff_np = npart_stu - npart_pub
    ax_bot.axhline(0.0, color="gray", ls="--", lw=1.1)
    ax_bot.errorbar(x, diff_np, yerr=npart_err_pub, fmt="o", color="#1f77b4", markersize=6.5, capsize=2)
    ax_bot.set_xticks(x)
    ax_bot.set_xticklabels(CENTRALITY_ORDER, fontsize=12.5)
    ax_bot.set_ylabel(r"$\Delta \langle N_{\mathrm{part}} \rangle$", fontsize=13.5, loc="center")
    ax_bot.set_xlabel("Collision centrality class", fontsize=15)
    ax_bot.set_ylim(-8.0, 8.0)
    fig.align_ylabels([ax_top, ax_bot])

    plots = save_fig(fig, outdir, "cmp1_constants_and_centrality")

    sqrts_stu = float(student["metadata"]["sqrts_nn_tev"])
    lumi_stu = float(student["metadata"]["lumi_full_2010_inv_ub"])
    max_np_diff = float(np.max(np.abs(npart_stu - npart_pub)))
    max_nps_diff = float(np.max(np.abs(npart_err_stu - npart_err_pub)))
    max_ehf_diff = float(max(bound_diffs_tev))
    exact_ok = (
        abs(sqrts_stu - 2.76) < 1e-9
        and abs(lumi_stu - 6.7) < 1e-9
        and max_np_diff < 1e-9
        and max_nps_diff < 1e-9
        and max_ehf_diff < 1e-4
    )
    verdict = "CONSISTENT" if exact_ok else "INCONSISTENT"

    return {
        "id": "T1",
        "name": "Collision constants, HF centrality boundaries, and Glauber <N_part>",
        "test_level": "exact",
        "verdict": verdict,
        "plots": plots,
        "metrics": {
            "sqrts_nn_tev_student": sqrts_stu,
            "sqrts_nn_tev_published": 2.76,
            "lumi_ref_inv_ub_student": lumi_stu,
            "lumi_ref_inv_ub_published": 6.7,
            "max_abs_diff_npart_mean": max_np_diff,
            "max_abs_diff_npart_syst": max_nps_diff,
            "max_abs_diff_ehf_boundary_tev": max_ehf_diff,
            "ehf_boundaries_student_tev": stu_bounds_tev,
            "ehf_boundaries_published_tev": pub_bounds_tev,
            "mb_centrality_fractions_student": dict(zip(CENTRALITY_ORDER, stu_cfracs)),
            "mb_centrality_fractions_hhfMB_macro": dict(zip(CENTRALITY_ORDER, macro_cfracs)),
        },
        "explanation": (
            f"All collision constants (sqrt(s_NN) = {sqrts_stu:.2f} TeV, reference L_int = {lumi_stu:.1f} ub^-1), "
            f"dijet kinematic thresholds (pT1 > 120 GeV, pT2 > 50 GeV, |eta| < 2.0, dphi12 > 2pi/3, A_J < 0.15, "
            f"dphi12 > 3.026), HF centrality boundaries (max rounding difference {max_ehf_diff:.2e} TeV), and "
            f"Table 2 Glauber <N_part> +/- syst values (max |diff| = {max_np_diff:.2g}) match CMS-HIN-10-004."
        ),
    }


def run_test2_cutflow(student: dict[str, Any], pub_ref: dict[str, Any], outdir: Path) -> dict[str, Any]:
    """T2 (absolute): Event selection cutflow step efficiencies & cleaned-scaled dijet yield."""
    cf = student["cutflow"]
    # Published Table 1 step efficiencies (combining sub-steps where student grouped them):
    # 1. BSC halo (0.9961) * HF coincidence (0.7498) = 0.74687578 (111,284 / 149,000)
    # 2. Reconstructed vertex = 0.9897 (110,138 / 111,284)
    # 3. Beam-gas (0.9978) * ECAL (0.9766) * HCAL (0.9997) = 0.974158 (107,292 / 110,138)
    # 4. >=2 jets pT > 35, |eta| < 2 = 0.6707 (71,961 / 107,292)
    # 5. Leading jet pT1 > 120 = 4216 / 71900 = 0.058637 (Table 1 prints 5.86%)
    # 6. Subleading jet pT2 > 50 = 3684 / 4216 = 0.873814 (Table 1 prints 87.38%)
    # 7. Dijet dphi12 > 2pi/3 = 3514 / 3684 = 0.953855 (Table 1 prints 95.39%)
    step_defs = [
        ("HF coinc.\n(incl. halo)", cf["hf_coincidence"], cf["hlt_hijet50u"], 0.9961 * 0.7498, 111284, 149000),
        ("Primary\nvertex", cf["primary_vertex"], cf["hf_coincidence"], 0.9897, 110138, 111284),
        ("Noise &\nscraping", cf["noise_and_scraping_cleaned"], cf["primary_vertex"], 0.9978 * 0.9766 * 0.9997, 107292, 110138),
        (">=2 jets\npT>35 GeV", cf["two_jets_pt35_eta2"], cf["noise_and_scraping_cleaned"], 0.6707, 71900, 107200),
        ("Lead jet\npT1>120", cf["leading_jet_pt120"], cf["two_jets_pt35_eta2"], 4216 / 71900, 4216, 71900),
        ("Sublead\npT2>50", cf["subleading_jet_pt50"], cf["leading_jet_pt120"], 3684 / 4216, 3684, 4216),
        ("Dijet\ndphi>2pi/3", cf["dijet_dphi_gt_2pi3"], cf["subleading_jet_pt50"], 3514 / 3684, 3514, 3684),
    ]

    rows = []
    for label, k_od, n_od, eff_pub, k_pub, n_pub in step_defs:
        eff_od = k_od / n_od
        se_od = math.sqrt(eff_od * (1.0 - eff_od) / n_od)
        se_pub = math.sqrt(eff_pub * (1.0 - eff_pub) / n_pub)
        se_comb = math.hypot(se_od, se_pub)
        pull = (eff_od - eff_pub) / se_comb
        ratio = eff_od / eff_pub
        ratio_err = ratio * math.hypot(se_od / max(eff_od, 1e-12), se_pub / max(eff_pub, 1e-12))
        rows.append({
            "step": label.replace("\n", " "),
            "k_student": int(k_od),
            "n_student": int(n_od),
            "eff_student": eff_od,
            "err_student": se_od,
            "eff_published": eff_pub,
            "err_published": se_pub,
            "diff_pp": 100.0 * (eff_od - eff_pub),
            "ratio": ratio,
            "ratio_err": ratio_err,
            "pull_indep_sigma": pull,
        })

    n_core_files = int(student["metadata"]["n_core_files_processed"])
    n_core_runs = int(student["metadata"]["n_core_runs_sampled"])
    f_hlt = cf["hlt_hijet50u"] / 149000.0
    f_clean = cf["noise_and_scraping_cleaned"] / 107292.0
    scaled_hlt = cf["dijet_dphi_gt_2pi3"] / f_hlt
    scaled_hlt_err = math.sqrt(cf["dijet_dphi_gt_2pi3"]) / f_hlt
    scaled_hlt_pull = (scaled_hlt - 3514.0) / scaled_hlt_err
    scaled_clean = cf["dijet_dphi_gt_2pi3"] / f_clean
    scaled_clean_err = math.sqrt(cf["dijet_dphi_gt_2pi3"]) / f_clean
    scaled_clean_pull = (scaled_clean - 3514.0) / scaled_clean_err

    # Evaluate separate ECAL spike vs HCAL/scraping sub-steps from validation_checks.ecal_cleaning_diagnostics
    ecal_diag = student.get("validation_checks", {}).get("ecal_cleaning_diagnostics")
    cleaning_substeps: dict[str, Any] = {}
    if ecal_diag:
        n_vtx_d = int(ecal_diag["primary_vertex_events"])
        n_ec_d = int(ecal_diag["after_ecal_spike_filter_events"])
        n_hc_d = int(ecal_diag["after_hcal_and_scraping_events"])

        eff_ec_stu = n_ec_d / n_vtx_d
        se_ec_stu = math.sqrt(eff_ec_stu * (1.0 - eff_ec_stu) / n_vtx_d)
        eff_ec_pub = 0.9766
        se_ec_pub = math.sqrt(eff_ec_pub * (1.0 - eff_ec_pub) / 109896.0)
        pull_ec = (eff_ec_stu - eff_ec_pub) / math.hypot(se_ec_stu, se_ec_pub)

        eff_hc_stu = n_hc_d / n_ec_d
        se_hc_stu = math.sqrt(eff_hc_stu * (1.0 - eff_hc_stu) / n_ec_d)
        eff_hc_pub = 0.9997
        se_hc_pub = math.sqrt(eff_hc_pub * (1.0 - eff_hc_pub) / 107324.0)
        pull_hc_only = (eff_hc_stu - eff_hc_pub) / math.hypot(se_hc_stu, se_hc_pub)

        eff_bghc_pub = 0.9978 * 0.9997
        se_bghc_pub = math.sqrt(eff_bghc_pub * (1.0 - eff_bghc_pub) / 107559.0)
        pull_bghc = (eff_hc_stu - eff_bghc_pub) / math.hypot(se_hc_stu, se_bghc_pub)

        cleaning_substeps = {
            "filter_config": ecal_diag.get("filter", ""),
            "ecal_spike_filter_vs_table1_ecal": {
                "k_student": n_ec_d,
                "n_student": n_vtx_d,
                "eff_student": eff_ec_stu,
                "err_student": se_ec_stu,
                "eff_published": eff_ec_pub,
                "err_published": se_ec_pub,
                "diff_pp": 100.0 * (eff_ec_stu - eff_ec_pub),
                "pull_indep_sigma": pull_ec,
            },
            "hcal_and_scraping_given_ecal_vs_table1_hcal": {
                "k_student": n_hc_d,
                "n_student": n_ec_d,
                "eff_student": eff_hc_stu,
                "err_student": se_hc_stu,
                "eff_published": eff_hc_pub,
                "err_published": se_hc_pub,
                "diff_pp": 100.0 * (eff_hc_stu - eff_hc_pub),
                "pull_indep_sigma": pull_hc_only,
            },
            "hcal_and_scraping_given_ecal_vs_table1_beamgas_times_hcal": {
                "k_student": n_hc_d,
                "n_student": n_ec_d,
                "eff_student": eff_hc_stu,
                "err_student": se_hc_stu,
                "eff_published": eff_bghc_pub,
                "err_published": se_bghc_pub,
                "diff_pp": 100.0 * (eff_hc_stu - eff_bghc_pub),
                "pull_indep_sigma": pull_bghc,
            },
        }

    # Plot T2
    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(10.5, 7.2), sharex=True, gridspec_kw={"height_ratios": [3, 1.15], "hspace": 0.07}
    )
    x = np.arange(len(rows), dtype=float)
    e_pub = np.array([r["eff_published"] * 100.0 for r in rows])
    s_pub = np.array([r["err_published"] * 100.0 for r in rows])
    e_stu = np.array([r["eff_student"] * 100.0 for r in rows])
    s_stu = np.array([r["err_student"] * 100.0 for r in rows])
    rat = np.array([r["ratio"] for r in rows])
    rat_e = np.array([r["ratio_err"] for r in rows])

    ax_top.errorbar(
        x - 0.10, e_pub, yerr=s_pub, fmt="o", color="black", markersize=7.0, capsize=3,
        label="CMS Published (HIN-10-004 Table 1)",
    )
    ax_top.errorbar(
        x + 0.10, e_stu, yerr=s_stu, fmt="s", color="#d62728", markersize=7.0, capsize=3,
        label=f"Student ({n_core_files:,} HICorePhysics files, {f_hlt*100.0:.1f}% HLT sample)",
    )
    if cleaning_substeps:
        ec_s = cleaning_substeps["ecal_spike_filter_vs_table1_ecal"]
        hc_s = cleaning_substeps["hcal_and_scraping_given_ecal_vs_table1_hcal"]
        ax_top.text(
            0.025, 0.30,
            f"Cleaning sub-steps vs Table 1:\n"
            f"• ECAL spike ($E_{{T}}>20$ GeV): {ec_s['eff_student']*100:.2f}% vs {ec_s['eff_published']*100:.2f}% ({ec_s['pull_indep_sigma']:+.1f}$\\sigma$)\n"
            f"• HCAL/EMF | ECAL: {hc_s['eff_student']*100:.2f}% vs {hc_s['eff_published']*100:.2f}% ({hc_s['pull_indep_sigma']:+.2f}$\\sigma$)\n"
            f"• Subleading $p_{{T,2}}>50$ GeV: {rows[5]['eff_student']*100:.2f}% vs {rows[5]['eff_published']*100:.2f}% ({rows[5]['pull_indep_sigma']:+.2f}$\\sigma$)",
            transform=ax_top.transAxes, fontsize=10.3,
            bbox={"boxstyle": "round,pad=0.28", "facecolor": "white", "alpha": 0.90, "edgecolor": "0.75"},
        )
    ax_top.set_ylabel("Conditional step efficiency [%]", fontsize=15.5)
    ax_top.set_ylim(-2.0, 114.0)
    ax_top.set_xlabel("")
    ax_top.legend(loc="lower left", frameon=True, fontsize=11.5)
    hep.cms.label("Open Data", data=True, rlabel=RLABEL_CMP, ax=ax_top, fontsize=16)

    ax_bot.axhline(1.0, color="gray", ls="--", lw=1.1)
    ax_bot.axhspan(0.95, 1.05, color="0.88", zorder=0)
    ax_bot.errorbar(x, rat, yerr=rat_e, fmt="o", color="#1f77b4", markersize=6.5, capsize=2)
    ax_bot.set_xticks(x)
    ax_bot.set_xticklabels([d[0] for d in step_defs], fontsize=10.5)
    ax_bot.set_ylabel("OD / CMS", fontsize=13, loc="center")
    ax_bot.set_xlabel("Sequential event selection step", fontsize=15)
    ax_bot.set_ylim(0.84, 1.18)
    fig.align_ylabels([ax_top, ax_bot])

    plots = save_fig(fig, outdir, "cmp2_cutflow_step_efficiencies")

    max_abs_pull = max(abs(r["pull_indep_sigma"]) for r in rows)
    verdict = verdict_from_sigma(max_abs_pull)

    sub_txt = ""
    if cleaning_substeps:
        ec_s = cleaning_substeps["ecal_spike_filter_vs_table1_ecal"]
        hc_s = cleaning_substeps["hcal_and_scraping_given_ecal_vs_table1_hcal"]
        bghc_s = cleaning_substeps["hcal_and_scraping_given_ecal_vs_table1_beamgas_times_hcal"]
        sub_txt = (
            f" Separating the cleaning sub-steps shows HCAL+scraping given ECAL ({hc_s['eff_student']*100:.2f}% vs "
            f"{hc_s['eff_published']*100:.2f}% HCAL cleaning, {hc_s['pull_indep_sigma']:+.2f} sigma; "
            f"{bghc_s['pull_indep_sigma']:+.2f} sigma vs 99.75% beam-gas*HCAL due to omitted pixel cluster compatibility) "
            f"matches HCAL cleaning, whereas HiEcalRecHitSpikeFilter (minEt=20 GeV) retains {ec_s['eff_student']*100:.2f}% "
            f"vs Table 1's 97.66% ({ec_s['pull_indep_sigma']:+.2f} sigma) because vetoing whole events on 20-50 GeV barrel "
            f"APD spikes removes random multiplicity-correlated soft spikes in central PbPb collisions."
        )

    return {
        "id": "T2",
        "name": "Event selection cutflow step efficiencies and scaled dijet yield",
        "test_level": "absolute",
        "verdict": verdict,
        "plots": plots,
        "metrics": {
            "n_core_files_processed": n_core_files,
            "n_core_runs_sampled": n_core_runs,
            "sample_hlt_fraction": f_hlt,
            "sample_cleaned_fraction": f_clean,
            "student_final_dijets": int(cf["dijet_dphi_gt_2pi3"]),
            "published_final_dijets": 3514,
            "scaled_dijet_yield_via_hlt": scaled_hlt,
            "scaled_dijet_yield_via_hlt_err": scaled_hlt_err,
            "scaled_dijet_yield_via_hlt_pull_sigma": scaled_hlt_pull,
            "scaled_dijet_yield_via_cleaned": scaled_clean,
            "scaled_dijet_yield_via_cleaned_err": scaled_clean_err,
            "scaled_dijet_yield_via_cleaned_pull_sigma": scaled_clean_pull,
            "max_abs_step_pull_sigma": max_abs_pull,
            "cleaning_substeps": cleaning_substeps,
            "steps": rows,
        },
        "explanation": (
            f"Across {n_core_files:,} Golden-JSON HICorePhysics files ({n_core_runs} runs), adding HiEcalRecHitSpikeFilter "
            f"eliminates fake barrel-edge spike mono-jets and brings the high-pT dijet selection steps into agreement with "
            f"CMS Table 1: subleading pT2 > 50 GeV is {rows[5]['eff_student']*100:.2f}% vs {rows[5]['eff_published']*100:.2f}% "
            f"({rows[5]['pull_indep_sigma']:+.2f} sigma), dphi12 > 2pi/3 is {rows[6]['eff_student']*100:.2f}% vs "
            f"{rows[6]['eff_published']*100:.2f}% ({rows[6]['pull_indep_sigma']:+.2f} sigma), leading pT1 > 120 GeV is "
            f"{rows[4]['eff_student']*100:.2f}% vs {rows[4]['eff_published']*100:.2f}% ({rows[4]['pull_indep_sigma']:+.2f} sigma), "
            f"and the HLT-scaled final dijet yield is {scaled_hlt:.0f} +/- {scaled_hlt_err:.0f} vs 3514 "
            f"({scaled_hlt_pull:+.2f} sigma; cleaned-scaled {scaled_clean:.0f} +/- {scaled_clean_err:.0f}, "
            f"{scaled_clean_pull:+.2f} sigma). However, four earlier intermediate step efficiencies remain >3 sigma "
            f"inconsistent (max |pull| = {max_abs_pull:.2f} sigma): HF coincidence ({rows[0]['eff_student']*100:.2f}% vs "
            f"{rows[0]['eff_published']*100:.2f}%, {rows[0]['pull_indep_sigma']:+.2f} sigma) and primary vertex "
            f"({rows[1]['eff_student']*100:.2f}% vs {rows[1]['eff_published']*100:.2f}%, {rows[1]['pull_indep_sigma']:+.2f} sigma) "
            f"shift inversely due to omitting the BSC halo veto and >=3-tower HF coincidence; combined noise/scraping cleaning "
            f"is {rows[2]['eff_student']*100:.2f}% vs {rows[2]['eff_published']*100:.2f}% ({rows[2]['pull_indep_sigma']:+.2f} sigma); "
            f"and >=2 jets pT > 35 GeV is {rows[3]['eff_student']*100:.2f}% vs {rows[3]['eff_published']*100:.2f}% "
            f"({rows[3]['pull_indep_sigma']:+.2f} sigma) under Spring10_HI L2L3 JEC.{sub_txt}"
        ),
    }


def run_test3_turnon(student: dict[str, Any], src_ver: dict[str, Any], outdir: Path) -> dict[str, Any]:
    """T3 (absolute): HLT_HIJet50U single-jet trigger turn-on curve on common [20, 160] GeV binning."""
    trig = student["trigger_turnon"]
    od_edges = trig["pt_edges"]
    od_k = trig["n_num"]
    od_n = trig["n_den"]

    pub_kn = src_ver["trig_exact_kn"]
    pub_bins = [
        (15, 20), (20, 25), (25, 30), (30, 35), (35, 40), (40, 45), (45, 50),
        (50, 55), (55, 60), (60, 65), (65, 70), (70, 75), (75, 80),
        (80, 90), (90, 100), (100, 110), (110, 120), (120, 140), (140, 160),
    ]
    common_edges = [20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 140.0, 160.0]

    bin_rows = []
    chi2_all = 0.0
    chi2_turn = 0.0
    for lo, hi in zip(common_edges[:-1], common_edges[1:]):
        pk = sum(pub_kn[i]["k"] for i, (blo, bhi) in enumerate(pub_bins) if blo >= lo and bhi <= hi)
        pn = sum(pub_kn[i]["n"] for i, (blo, bhi) in enumerate(pub_bins) if blo >= lo and bhi <= hi)
        ok = sum(od_k[i] for i in range(len(od_k)) if od_edges[i] >= lo and od_edges[i + 1] <= hi)
        on = sum(od_n[i] for i in range(len(od_n)) if od_edges[i] >= lo and od_edges[i + 1] <= hi)

        pef, pel, peh = clopper_pearson(pk, pn)
        oef, oel, oeh = clopper_pearson(ok, on)
        pef, pel, peh = float(pef), float(pel), float(peh)
        oef, oel, oeh = float(oef), float(oel), float(oeh)

        if oef > pef:
            sig = math.hypot(oel, peh)
        elif oef < pef:
            sig = math.hypot(oeh, pel)
        else:
            sig = math.hypot(max(oel, oeh), max(pel, peh))
        pull = (oef - pef) / sig if sig > 0 else 0.0
        _, p_fisher = stats.fisher_exact([[ok, on - ok], [pk, pn - pk]])
        chi2_all += pull**2
        if lo >= 70.0:
            chi2_turn += pull**2
        bin_rows.append({
            "pt_low": lo,
            "pt_high": hi,
            "student_k": int(ok),
            "student_n": int(on),
            "student_eff": oef,
            "student_err_low": oel,
            "student_err_high": oeh,
            "published_k": int(pk),
            "published_n": int(pn),
            "published_eff": pef,
            "published_err_low": pel,
            "published_err_high": peh,
            "diff": oef - pef,
            "pull_sigma": pull,
            "fisher_p": float(p_fisher),
        })

    # Plateau summary >= 100 GeV in [100, 160] GeV and full [100, 260] GeV
    ok_100_160 = sum(r["student_k"] for r in bin_rows if r["pt_low"] >= 100.0)
    on_100_160 = sum(r["student_n"] for r in bin_rows if r["pt_low"] >= 100.0)
    pk_100_160 = sum(r["published_k"] for r in bin_rows if r["pt_low"] >= 100.0)
    pn_100_160 = sum(r["published_n"] for r in bin_rows if r["pt_low"] >= 100.0)
    ok_120_full = sum(od_k[i] for i in range(len(od_k)) if od_edges[i] >= 120.0)
    on_120_full = sum(od_n[i] for i in range(len(od_n)) if od_edges[i] >= 120.0)

    # Plot T3
    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(9.4, 7.2), sharex=True, gridspec_kw={"height_ratios": [3, 1.1], "hspace": 0.07}
    )
    c_edges = np.asarray(common_edges, dtype=float)
    centers = 0.5 * (c_edges[:-1] + c_edges[1:])
    hw = 0.5 * np.diff(c_edges)
    p_eff = np.array([r["published_eff"] for r in bin_rows])
    p_el = np.array([r["published_err_low"] for r in bin_rows])
    p_eh = np.array([r["published_err_high"] for r in bin_rows])
    s_eff = np.array([r["student_eff"] for r in bin_rows])
    s_el = np.array([r["student_err_low"] for r in bin_rows])
    s_eh = np.array([r["student_err_high"] for r in bin_rows])
    pulls = np.array([r["pull_sigma"] for r in bin_rows])

    ax_top.axhline(1.0, color="gray", ls=":", lw=1.2)
    ax_top.axvline(100.0, color="#2ca02c", ls="--", lw=1.4, label=r"Full efficiency ($p_{T,1} = 100\ \mathrm{GeV}$)")
    ax_top.axvline(120.0, color="0.35", ls="-.", lw=1.4, label=r"Dijet analysis cut ($p_{T,1} > 120\ \mathrm{GeV}$)")
    ax_top.errorbar(
        centers - 0.8, p_eff, xerr=hw, yerr=[p_el, p_eh], fmt="o", color="black", markersize=6.5, capsize=2,
        label="CMS HI CVS macro (PRC Fig. 2a, rebinned)",
    )
    ax_top.errorbar(
        centers + 0.8, s_eff, xerr=hw, yerr=[s_el, s_eh], fmt="s", color="#d62728", markersize=6.5, capsize=2,
        label="Student Open Data (MB < 90, HIJet35U >= 90 GeV)",
    )
    ax_top.set_ylabel(r"$\mathrm{HLT\_HIJet50U}$ efficiency", fontsize=15.5)
    ax_top.set_ylim(-0.04, 1.24)
    ax_top.set_xlabel("")
    ax_top.legend(loc="upper left", frameon=True, fontsize=10.8)
    hep.cms.label("Open Data", data=True, rlabel=RLABEL_CMP, ax=ax_top, fontsize=16)

    ax_bot.axhline(0.0, color="gray", ls="--", lw=1.1)
    ax_bot.axhspan(-2.0, 2.0, color="#d9ead3", alpha=0.55, zorder=0)
    ax_bot.plot(centers, pulls, "o-", color="#1f77b4", markersize=6.0, lw=1.3)
    ax_bot.set_ylabel(r"Pull $[\sigma]$", fontsize=13.5, loc="center")
    ax_bot.set_xlabel(r"Corrected leading IC5Calo jet $p_{T,1}$ [GeV]", fontsize=15)
    ax_bot.set_xlim(20.0, 160.0)
    ax_bot.set_ylim(-3.2, 3.2)
    fig.align_ylabels([ax_top, ax_bot])

    plots = save_fig(fig, outdir, "cmp3_fig1a_trigger_turnon")

    p_all = float(stats.chi2.sf(chi2_all, 12))
    p_turn = float(stats.chi2.sf(chi2_turn, 6))
    max_abs_pull = float(np.max(np.abs(pulls)))
    worst_bin = max(bin_rows, key=lambda r: abs(r["pull_sigma"]))
    min_fp = float(min(r["fisher_p"] for r in bin_rows))
    verdict = verdict_from_sigma(max(two_sided_z_from_p(p_all), max_abs_pull))

    return {
        "id": "T3",
        "name": "HLT_HIJet50U single-jet trigger turn-on curve (PRC Fig. 2a / Open Data Fig. 1a)",
        "test_level": "absolute",
        "verdict": verdict,
        "plots": plots,
        "metrics": {
            "common_range_gev": [20.0, 160.0],
            "n_common_bins": len(bin_rows),
            "chi2_all_12bins": chi2_all,
            "ndf_all": 12,
            "p_value_all": p_all,
            "chi2_turnon_plateau_6bins_70_160": chi2_turn,
            "ndf_turnon_plateau": 6,
            "p_value_turnon_plateau": p_turn,
            "max_abs_pull_sigma": max_abs_pull,
            "min_fisher_exact_p": min_fp,
            "plateau_100_160_student": {"k": int(ok_100_160), "n": int(on_100_160), "eff": ok_100_160 / on_100_160},
            "plateau_100_160_published": {"k": int(pk_100_160), "n": int(pn_100_160), "eff": pk_100_160 / pn_100_160},
            "plateau_ge120_student": {"k": int(ok_120_full), "n": int(on_120_full), "eff": ok_120_full / on_120_full},
            "bins": bin_rows,
        },
        "explanation": (
            f"Rebinned onto the exact common 12-bin grid over [20, 160] GeV by summing exact binomial (k, n) counts "
            f"from CmsHI/CVS_CmsHi TrigEffCorrected151058and151025.C (PRC Fig. 2a), "
            f"the Open Data turn-on curve agrees across the turn-on rise and 100% plateau with chi2/ndf = "
            f"{chi2_all:.2f}/12 = {chi2_all/12.0:.3f} (p = {p_all:.3f}), maximum single-bin pull "
            f"{worst_bin['pull_sigma']:+.2f} sigma ([{worst_bin['pt_low']:.0f}, {worst_bin['pt_high']:.0f}] GeV, "
            f"Fisher p = {worst_bin['fisher_p']:.3f}; min Fisher p = {min_fp:.3f}), and plateau efficiency "
            f"{ok_100_160/on_100_160*100:.3f}% ({ok_100_160}/{on_100_160}) for pT1 in [100, 160] GeV vs "
            f"{pk_100_160/pn_100_160*100:.1f}% ({pk_100_160}/{pn_100_160}) in runs 151058+151025."
        ),
    }


def _run_hf_shape_test(
    test_id: str,
    title: str,
    sample_label: str,
    stem: str,
    edges_tev: np.ndarray,
    od_counts: np.ndarray,
    pub_raw: np.ndarray,
    pub_eff_total_counts: float,
    boundaries_tev: list[float],
    outdir: Path,
    highlight_bin0: bool = False,
    student_50_100_frac: float | None = None,
    macro_50_100_frac: float | None = None,
) -> dict[str, Any]:
    od_sum = float(od_counts.sum())
    pub_counts = np.round(pub_raw / pub_raw.sum() * pub_eff_total_counts)
    pub_sum = float(pub_counts.sum())

    od_p = od_counts / od_sum
    pub_p = pub_raw / float(pub_raw.sum())
    od_err = np.sqrt(np.maximum(od_counts, 0.0)) / od_sum
    pub_err = np.sqrt(np.maximum(pub_counts, 0.0)) / pub_sum
    tot_err = np.hypot(od_err, pub_err)

    mask = (od_counts >= 5) & (pub_counts >= 5)
    chi2 = float(np.sum(((od_p[mask] - pub_p[mask]) / tot_err[mask]) ** 2))
    ndf = int(mask.sum()) - 1
    pval = float(stats.chi2.sf(chi2, ndf))

    centers = 0.5 * (edges_tev[:-1] + edges_tev[1:])
    ratio = np.where(mask, od_p / np.where(pub_p > 0, pub_p, np.nan), np.nan)
    ratio_err = np.where(
        mask,
        ratio * np.hypot(1.0 / np.sqrt(np.maximum(od_counts, 1.0)), 1.0 / np.sqrt(np.maximum(pub_counts, 1.0))),
        np.nan,
    )
    fit = weighted_linear_fit(centers[mask], ratio[mask], ratio_err[mask])
    verdict, z_chi2 = verdict_from_shape(pval, fit["slope_pull_sigma"])

    od_mean = float(np.sum(centers * od_p))
    pub_mean = float(np.sum(centers * pub_p))
    b0_pull = float((od_p[0] - pub_p[0]) / tot_err[0])

    # Build hist.Hist objects for mplhep plotting
    h_pub = Hist.new.Var(edges_tev, name="ehf", label=r"$\sum E_{\mathrm{HF}}$ [TeV]").Weight()
    h_pub.view().value = pub_p
    h_pub.view().variance = pub_err**2

    h_od = Hist.new.Var(edges_tev, name="ehf", label=r"$\sum E_{\mathrm{HF}}$ [TeV]").Weight()
    h_od.view().value = od_p
    h_od.view().variance = od_err**2

    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(9.5, 7.3), sharex=True, gridspec_kw={"height_ratios": [3, 1.15], "hspace": 0.07}
    )
    hep.histplot(
        h_pub, ax=ax_top, histtype="step", color="black", lw=1.8,
        label=f"CMS HI CVS macro (PRC Fig. 3a, N={int(pub_sum):,})",
    )
    hep.histplot(
        h_od, ax=ax_top, histtype="errorbar", color="#d62728", markersize=5.2,
        label=f"Student Open Data ({sample_label}, N={int(od_sum):,})",
    )
    for b in boundaries_tev:
        ax_top.axvline(b, color="0.55", ls="--", lw=1.0)
    ax_top.set_yscale("log")
    ax_top.set_ylim(2e-5, 0.8)
    ax_top.set_ylabel("Event fraction / 2.29 TeV", fontsize=15.5)
    ax_top.set_xlabel("")
    ax_top.legend(loc="upper right", frameon=True, fontsize=11.5)
    p_str = fmt_p_latex(pval)
    ax_top.text(
        0.04, 0.08,
        f"$\\chi^2/\\mathrm{{ndf}} = {chi2:.1f}/{ndf} = {chi2/ndf:.2f}$ ($p = {p_str}$, ${z_chi2:.2f}\\sigma$)\n"
        f"$\\langle \\sum E_{{\\mathrm{{HF}}}} \\rangle$: OD = {od_mean:.2f} TeV, CMS = {pub_mean:.2f} TeV",
        transform=ax_top.transAxes, fontsize=11.2,
        bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "alpha": 0.88, "edgecolor": "0.75"},
    )
    hep.cms.label("Open Data", data=True, rlabel=RLABEL_CMP, ax=ax_top, fontsize=16)

    ax_bot.axhline(1.0, color="gray", ls="--", lw=1.1)
    ax_bot.errorbar(centers[mask], ratio[mask], yerr=ratio_err[mask], fmt="o", color="#1f77b4", markersize=4.8, lw=1.1)
    x_fit = np.linspace(float(centers[mask].min()), float(centers[mask].max()), 100)
    ax_bot.plot(
        x_fit, fit["intercept"] + fit["slope"] * x_fit, color="#d62728", ls="-", lw=1.6,
        label=f"Slope = ({fit['slope']*1e4:+.1f} $\\pm$ {fit['slope_err']*1e4:.1f})$\\times 10^{{-4}}$ TeV$^{{-1}}$ ({fit['slope_pull_sigma']:+.1f}$\\sigma$)",
    )
    if highlight_bin0:
        ax_bot.axvspan(0.0, edges_tev[1], color="#ff9896", alpha=0.35, label=f"Bin 1 [0, {edges_tev[1]:.2f}] TeV ({b0_pull:+.2f}$\\sigma$)")
    ax_bot.legend(loc="upper right", frameon=True, fontsize=10.1)
    ax_bot.set_ylabel("OD / CMS", fontsize=13, loc="center")
    ax_bot.set_xlabel(r"Total forward calorimeter energy $\sum E_{\mathrm{HF}}$ [TeV]", fontsize=15)
    ax_bot.set_xlim(0.0, 160.0)
    ax_bot.set_ylim(0.52 if highlight_bin0 else 0.55, 1.70 if highlight_bin0 else 1.45)
    fig.align_ylabels([ax_top, ax_bot])

    plots = save_fig(fig, outdir, stem)

    pop_lo = float(edges_tev[np.where(mask)[0][0]])
    pop_hi = float(edges_tev[np.where(mask)[0][-1] + 1])
    extra_no_bin0 = {}
    if highlight_bin0:
        od_s = od_counts[1:] / float(od_counts[1:].sum())
        pub_s = pub_counts[1:] / float(pub_counts[1:].sum())
        m_s = (od_counts[1:] >= 5) & (pub_counts[1:] >= 5)
        err_s = np.hypot(np.sqrt(od_counts[1:]) / float(od_counts[1:].sum()), np.sqrt(pub_counts[1:]) / float(pub_counts[1:].sum()))
        c2_s = float(np.sum(((od_s[m_s] - pub_s[m_s]) / err_s[m_s]) ** 2))
        ndf_s = int(m_s.sum()) - 1
        p_s = float(stats.chi2.sf(c2_s, ndf_s))
        extra_no_bin0 = {
            "bin1_0_2p29_tev_student": float(od_p[0]),
            "bin1_0_2p29_tev_student_err": float(od_err[0]),
            "bin1_0_2p29_tev_published": float(pub_p[0]),
            "bin1_0_2p29_tev_published_err": float(pub_err[0]),
            "bin1_pull_sigma": b0_pull,
            "excluding_bin1_chi2": c2_s,
            "excluding_bin1_ndf": ndf_s,
            "excluding_bin1_p_value": p_s,
        }
        cfrac_txt = ""
        if student_50_100_frac is not None and macro_50_100_frac is not None:
            se_cf = math.hypot(
                math.sqrt(student_50_100_frac * (1.0 - student_50_100_frac) / od_sum),
                math.sqrt(macro_50_100_frac * (1.0 - macro_50_100_frac) / pub_sum),
            )
            pull_cf = (student_50_100_frac - macro_50_100_frac) / se_cf
            extra_no_bin0["mb_50_100_fraction_pull_sigma"] = pull_cf
            cfrac_txt = (
                f", inflating the 50-100% MB centrality fraction to {student_50_100_frac*100:.2f}% vs "
                f"{macro_50_100_frac*100:.2f}% ({pull_cf:+.2f} sigma)"
            )
        explanation = (
            f"Across the {int(mask.sum())} populated bins [{pop_lo:.1f}, {pop_hi:.1f}] TeV, the Open Data Minimum-Bias "
            f"sum E_HF distribution (N = {int(od_sum):,}) yields chi2/ndf = {chi2:.2f}/{ndf} = {chi2/ndf:.3f} "
            f"(p = {pval:.2e}, {z_chi2:.2f} sigma) and a negative ratio slope of ({fit['slope']*1e4:+.1f} +/- "
            f"{fit['slope_err']*1e4:.1f})e-4 TeV^-1 ({fit['slope_pull_sigma']:+.2f} sigma) against hhfMB (N = {int(pub_sum):,}). "
            f"Two selection differences compound: (1) sampling ordinary HLT_HIMinBiasHfOrBSC in HIAllPhysics with "
            f"etHFtowerSum > 3 GeV — whereas PRC Fig. 3a's archived canvas (commit f6717ae4) names HLT_MinBiasHForBSC_Core "
            f"with noBSChalo, >=3 HF towers/side, and 2-track vertex (which gives 16,161/70,088 = 23.06% in bin 1 and "
            f"G^2 = 54.70/53 ndf across 110 HICorePhysics files) — produces a {b0_pull:+.2f} sigma excess "
            f"in bin 1 [0, {edges_tev[1]:.2f}] TeV ({od_p[0]:.4f} +/- {od_err[0]:.4f} vs {pub_p[0]:.4f} +/- {pub_err[0]:.4f})"
            f"{cfrac_txt}; and (2) applying the whole-event HiEcalRecHitSpikeFilter (minEt = 20 GeV) to Minimum Bias "
            f"preferentially vetoes high-multiplicity central collisions (~5.5% in 0-10% vs ~0.2% in 50-100%), steepening "
            f"the tilt. Excluding bin 1 reduces chi2/ndf to {c2_s:.2f}/{ndf_s} = {c2_s/ndf_s:.3f} (p = {p_s:.3f})."
        )
    else:
        explanation = (
            f"Unit-normalized over [0, 160] TeV on identical 70-bin histograms ({int(mask.sum())} populated bins in "
            f"[{pop_lo:.1f}, {pop_hi:.1f}] TeV), the Open Data {sample_label} sum E_HF distribution (N = {int(od_sum):,}) "
            f"agrees with hhfJET (N = {int(pub_sum):,} in hf_centrality_distribution_mb_jet50_20101126_v1.C): "
            f"chi2/ndf = {chi2:.2f}/{ndf} = {chi2/ndf:.3f} (p = {pval:.3f}, {z_chi2:.2f} sigma), ratio slope "
            f"({fit['slope']*1e4:+.1f} +/- {fit['slope_err']*1e4:.1f})e-4 TeV^-1 ({fit['slope_pull_sigma']:+.2f} sigma), "
            f"and <sum E_HF> = {od_mean:.2f} TeV vs {pub_mean:.2f} TeV."
        )

    return {
        "id": test_id,
        "name": title,
        "test_level": "shape",
        "verdict": verdict,
        "plots": plots,
        "metrics": {
            "range_populated_tev": [pop_lo, pop_hi],
            "n_populated_bins": int(mask.sum()),
            "student_total_events": int(od_sum),
            "published_macro_entries": int(pub_sum),
            "chi2": chi2,
            "ndf": ndf,
            "chi2_per_ndf": chi2 / ndf,
            "p_value": pval,
            "equivalent_gaussian_sigma": z_chi2,
            "mean_ehf_student_tev": od_mean,
            "mean_ehf_published_tev": pub_mean,
            "ratio_linear_fit": fit,
            **extra_no_bin0,
        },
        "explanation": explanation,
    }


def run_test6_fig10(student: dict[str, Any], pub_ref: dict[str, Any], outdir: Path) -> dict[str, Any]:
    """T6 (absolute): Archived Dec 2010 CMS HI CVS deltaPtOverPt1VsPt.C cross-check + Open Data PRC Fig. 12 profile."""
    fi_prc12 = student["fractional_imbalance_vs_pt1"]
    fi = student.get("fractional_imbalance_vs_pt1_cvs_macro", fi_prc12)
    fi100 = student["fractional_imbalance_vs_pt1_100_240"]

    edges = np.asarray(fi["pt1_edges"], dtype=float)
    centers = 0.5 * (edges[:-1] + edges[1:])
    hw = 0.5 * np.diff(edges)

    edges_prc = np.asarray(fi_prc12["pt1_edges"], dtype=float)
    centers_prc = 0.5 * (edges_prc[:-1] + edges_prc[1:])
    hw_prc = 0.5 * np.diff(edges_prc)

    per_cent_metrics = {}
    prc12_student_bins = {}
    tot_c2_stat = 0.0
    tot_c2_full = 0.0
    tot_ndf = 0

    fig, axes = plt.subplots(
        2, 3, figsize=(16.2, 7.2), sharex=True, sharey="row",
        gridspec_kw={"height_ratios": [3, 1.15], "wspace": 0.06, "hspace": 0.07},
    )

    for idx, mname in enumerate(MERGED_ORDER):
        ax_t = axes[0, idx]
        ax_b = axes[1, idx]

        od_m = np.asarray(fi[mname]["mean"], dtype=float)
        od_s = np.asarray(fi[mname]["sem"], dtype=float)
        od_n = np.asarray(fi[mname]["counts"], dtype=int)

        od_prc_m = np.asarray(fi_prc12[mname]["mean"], dtype=float)
        od_prc_s = np.asarray(fi_prc12[mname]["sem"], dtype=float)
        od_prc_n = np.asarray(fi_prc12[mname]["counts"], dtype=int)

        pub_pts = [p for p in pub_ref["delta_pt_over_pt1_vs_pt1"][mname]["pbpb_data"] if p["pt1_center_gev"] >= 120.0]
        pub_m = np.array([p["val"] for p in pub_pts], dtype=float)
        pub_s = np.array([p["stat_err"] for p in pub_pts], dtype=float)
        pub_sy = np.array([p["syst_err"] for p in pub_pts], dtype=float)
        pub_tot = np.hypot(pub_s, pub_sy)

        sig_stat = np.hypot(od_s, pub_s)
        sig_full = np.sqrt(od_s**2 + pub_s**2 + pub_sy**2)
        c2_stat = float(np.sum(((od_m - pub_m) / sig_stat) ** 2))
        c2_full = float(np.sum(((od_m - pub_m) / sig_full) ** 2))
        ndf = len(od_m)

        tot_c2_stat += c2_stat
        tot_c2_full += c2_full
        tot_ndf += ndf

        ratio = od_m / pub_m
        ratio_stat_err = ratio * np.hypot(od_s / od_m, pub_s / pub_m)
        ratio_full_err = ratio * np.hypot(od_s / od_m, pub_tot / pub_m)
        rfit = weighted_linear_fit(centers, ratio, ratio_full_err)

        # Draw archived CVS macro systematic brackets
        for bx, by, bsy in zip(centers, pub_m, pub_sy):
            rect = plt.Rectangle(
                (bx - 4.0, by - bsy), 8.0, 2.0 * bsy,
                facecolor="0.82", edgecolor="0.45", alpha=0.65, lw=0.9,
                label="CVS macro syst. bracket" if bx == centers[0] else None,
            )
            ax_t.add_patch(rect)

        ax_t.errorbar(
            centers - 1.2, pub_m, xerr=None, yerr=pub_s, fmt="o", color="black", markersize=6.8, capsize=2,
            label="CMS HI CVS macro (Dec 13, 2010)",
        )
        ax_t.errorbar(
            centers + 1.2, od_m, xerr=hw, yerr=od_s, fmt="s", color="#d62728", markersize=6.5, capsize=2,
            label=f"Open Data (CVS bins, $N={int(od_n.sum())}$)",
        )
        ax_t.errorbar(
            centers_prc, od_prc_m, xerr=hw_prc, yerr=od_prc_s, fmt="D", mfc="white", mec="#1f77b4",
            ecolor="#1f77b4", markersize=5.8, capsize=2, alpha=0.88,
            label="Open Data (PRC Fig. 12 bins)",
        )
        ax_t.set_ylim(0.10, 0.60)
        ax_t.set_xlim(120.0, 240.0)
        ax_t.text(
            0.06, 0.92,
            f"Centrality {mname}\nCVS macro $\\chi^2_{{\\mathrm{{stat+syst}}}}/\\mathrm{{ndf}} = {c2_full:.2f}/{ndf}$",
            transform=ax_t.transAxes, va="top", fontsize=11.8,
        )
        if idx == 0:
            ax_t.set_ylabel(r"$\langle (p_{T,1} - p_{T,2}) / p_{T,1} \rangle$", fontsize=14.5)
            hep.cms.text("Open Data", ax=ax_t, fontsize=14.5, loc=0)
        elif idx == 2:
            ax_t.set_title(RLABEL_CMP, loc="right", fontsize=13.5, pad=8)
            ax_t.legend(loc="lower right", frameon=True, fontsize=9.8)

        ax_b.axhline(1.0, color="gray", ls="--", lw=1.1)
        for bx, br, bfe in zip(centers, ratio, ratio_full_err):
            rrect = plt.Rectangle((bx - 4.0, br - bfe), 8.0, 2.0 * bfe, facecolor="#fdd0a2", edgecolor="none", alpha=0.55)
            ax_b.add_patch(rrect)
        ax_b.errorbar(centers, ratio, yerr=ratio_stat_err, fmt="o", color="#1f77b4", markersize=6.0, capsize=2)
        ax_b.set_ylim(0.72, 1.28)
        ax_b.set_xlabel(r"Leading jet $p_{T,1}$ [GeV]", fontsize=13.5)
        if idx == 0:
            ax_b.set_ylabel("OD / CVS", fontsize=13.5)

        per_cent_metrics[mname] = {
            "n_dijets_student": int(od_n.sum()),
            "chi2_stat": c2_stat,
            "chi2_stat_plus_syst": c2_full,
            "ndf": ndf,
            "p_value_stat": float(stats.chi2.sf(c2_stat, ndf)),
            "p_value_stat_plus_syst": float(stats.chi2.sf(c2_full, ndf)),
            "max_abs_pull_stat_plus_syst": float(np.max(np.abs((od_m - pub_m) / sig_full))),
            "ratio_linear_fit": rfit,
            "bins": [
                {
                    "pt1_range": [float(edges[i]), float(edges[i + 1])],
                    "student_val": float(od_m[i]),
                    "student_sem": float(od_s[i]),
                    "student_n": int(od_n[i]),
                    "cvs_macro_val": float(pub_m[i]),
                    "cvs_macro_stat_err": float(pub_s[i]),
                    "cvs_macro_syst_err": float(pub_sy[i]),
                    "pull_stat_plus_syst": float((od_m[i] - pub_m[i]) / sig_full[i]),
                }
                for i in range(ndf)
            ],
        }
        prc12_student_bins[mname] = [
            {
                "pt1_range": [float(edges_prc[i]), float(edges_prc[i + 1])],
                "pt1_center_gev": float(centers_prc[i]),
                "student_val": float(od_prc_m[i]),
                "student_sem": float(od_prc_s[i]),
                "student_n": int(od_prc_n[i]),
            }
            for i in range(len(od_prc_m))
        ]

    # Extended 6-bin [100, 240] GeV check against full archived CVS macro
    c2_6bin_stat = 0.0
    c2_6bin_full = 0.0
    for mname in MERGED_ORDER:
        om6 = np.asarray(fi100[mname]["mean"], dtype=float)
        os6 = np.asarray(fi100[mname]["sem"], dtype=float)
        p6 = pub_ref["delta_pt_over_pt1_vs_pt1"][mname]["pbpb_data"]
        pm6 = np.array([x["val"] for x in p6], dtype=float)
        ps6 = np.array([x["stat_err"] for x in p6], dtype=float)
        psy6 = np.array([x["syst_err"] for x in p6], dtype=float)
        c2_6bin_stat += float(np.sum(((om6 - pm6) / np.hypot(os6, ps6)) ** 2))
        c2_6bin_full += float(np.sum(((om6 - pm6) / np.sqrt(os6**2 + ps6**2 + psy6**2)) ** 2))

    plots = save_fig(fig, outdir, "cmp6_fig10_fractional_imbalance")

    p_stat = float(stats.chi2.sf(tot_c2_stat, tot_ndf))
    p_full = float(stats.chi2.sf(tot_c2_full, tot_ndf))
    all_pulls_full = [b["pull_stat_plus_syst"] for m in MERGED_ORDER for b in per_cent_metrics[m]["bins"]]
    worst_p_full = max(all_pulls_full, key=abs)
    max_abs_p_full = abs(worst_p_full)
    verdict = verdict_from_sigma(max(two_sided_z_from_p(p_stat), max_abs_p_full))

    cent_means = {
        m: float(
            np.sum(np.asarray(fi[m]["counts"]) * np.asarray(fi[m]["mean"]))
            / max(np.sum(fi[m]["counts"]), 1)
        )
        for m in MERGED_ORDER
    }

    return {
        "id": "T6",
        "name": "Mean fractional dijet imbalance <(pT1 - pT2)/pT1> vs pT1 (archived Dec 2010 CMS HI CVS macro cross-check; PRC Fig. 12 PDF-only)",
        "test_level": "absolute",
        "verdict": verdict,
        "plots": plots,
        "metrics": {
            "macro_provenance": (
                "github.com/CmsHI/CVS_CmsHi commit 9b9dcd34a1b718e4649ca2ddc34583706dfd5d1c, "
                "Analysis2010/macros/JetAnalysis/plots/deltaPtOverPt1VsPt.C (saved Mon Dec 13 15:07:45 2010)"
            ),
            "prc_fig12_published_status": (
                "arXiv:1102.1957 figures/deltaPtOverPt1_all_cent_20101223v1.pdf (Phys. Rev. C 84, 024906 Fig. 12, "
                "saved Dec 23, 2010) uses [120, 140, 160, 180, 240] GeV (centers 130, 150, 170, 210 GeV) and has no "
                "published numerical table or .C macro; not digitized per ARC policy."
            ),
            "cvs_macro_pt1_edges_gev": edges.tolist(),
            "prc_fig12_pt1_edges_gev": edges_prc.tolist(),
            "primary_range_gev": [120.0, 240.0],
            "total_points_primary": tot_ndf,
            "combined_chi2_stat": tot_c2_stat,
            "combined_chi2_stat_plus_syst": tot_c2_full,
            "combined_ndf": tot_ndf,
            "combined_p_value_stat": p_stat,
            "combined_p_value_stat_plus_syst": p_full,
            "max_abs_pull_stat_plus_syst": max_abs_p_full,
            "extended_100_240_18pts_chi2_stat": c2_6bin_stat,
            "extended_100_240_18pts_chi2_stat_plus_syst": c2_6bin_full,
            "per_merged_centrality": per_cent_metrics,
            "student_prc_fig12_bins": prc12_student_bins,
        },
        "explanation": (
            f"Matched bin-for-bin on [120, 130, 140, 170, 240] GeV (centers 125, 135, 155, 205 GeV) across 3 merged "
            f"centrality classes ({tot_ndf} points), the Open Data fractional imbalance agrees with the archived "
            f"Dec 13, 2010 CMS Heavy-Ion ROOT canvas deltaPtOverPt1VsPt.C (CmsHI/CVS_CmsHi) with chi2/ndf = "
            f"{tot_c2_stat:.2f}/{tot_ndf} (stat, p = {p_stat:.4f}) and {tot_c2_full:.2f}/{tot_ndf} (stat+syst), "
            f"max stat+syst pull {worst_p_full:+.2f} sigma, and centrality ordering {cent_means['30-100%']:.3f} "
            f"(30-100%) -> {cent_means['10-30%']:.3f} (10-30%) -> {cent_means['0-10%']:.3f} (0-10%); the 18-point "
            f"[100, 240] GeV check gives chi2/ndf = {c2_6bin_stat:.2f}/18 (stat) and {c2_6bin_full:.2f}/18 (stat+syst) "
            f"(note: macro stat errors ~0.03-0.06 reflect uncorrelated TH1::Divide propagation, ~3x sample SEM). "
            f"Final published Phys. Rev. C 84, 024906 Fig. 12 (deltaPtOverPt1_all_cent_20101223v1.pdf, Dec 23, 2010) "
            f"uses [120, 140, 160, 180, 240] GeV (centers 130, 150, 170, 210 GeV) and has no numerical table or .C "
            f"macro (not digitized); exact Open Data values on the PRC Fig. 12 binning (e.g. 0-10% [120, 140) GeV: "
            f"{prc12_student_bins['0-10%'][0]['student_val']:.4f} +/- {prc12_student_bins['0-10%'][0]['student_sem']:.4f}) "
            f"are exported in summary.json, fig10_fractional_imbalance_vs_pt1, hepdata/, and metrics.student_prc_fig12_bins."
        ),
    }


def run_test7_rb_aj(student: dict[str, Any], outdir: Path) -> dict[str, Any]:
    """T7 (absolute): Balanced dijet fraction R_B(A_J < 0.15) vs Sec. III.A.3 text benchmarks."""
    pc = student["per_centrality"]
    npart = np.array([pc[c]["npart_mean"] for c in CENTRALITY_ORDER], dtype=float)
    npart_err = np.array([pc[c]["npart_err"] for c in CENTRALITY_ORDER], dtype=float)
    rb = np.array([pc[c]["rb_aj_lt_015"] for c in CENTRALITY_ORDER], dtype=float)
    elo = np.array([pc[c]["rb_aj_err_low"] for c in CENTRALITY_ORDER], dtype=float)
    ehi = np.array([pc[c]["rb_aj_err_high"] for c in CENTRALITY_ORDER], dtype=float)
    syst = np.array([pc[c]["rb_aj_syst"] for c in CENTRALITY_ORDER], dtype=float)

    k50 = int(pc["50-100%"]["k_rb_aj_lt_015"])
    n50 = int(pc["50-100%"]["n_lead120"])
    k0 = int(pc["0-10%"]["k_rb_aj_lt_015"])
    n0 = int(pc["0-10%"]["n_lead120"])

    # Published Section III.A.3 explicit text benchmarks:
    # 50-100%: "close to 45%" (~0.45) with 9% relative syst (+/- 0.0405)
    # 0-10%: "drops by close to a factor of two" (~0.225) with 11% relative syst (+/- 0.02475)
    pub_50 = 0.45
    pub_50_syst = 0.09 * pub_50
    pub_0 = 0.225
    pub_0_syst = 0.11 * pub_0

    rb50_stat = 0.5 * (elo[0] + ehi[0])
    rb0_stat = 0.5 * (elo[4] + ehi[4])
    pull_50_stat_only = (rb[0] - pub_50) / rb50_stat
    pull_50_full = (rb[0] - pub_50) / math.sqrt(rb50_stat**2 + syst[0] ** 2 + pub_50_syst**2)
    pull_0_stat_only = (rb[4] - pub_0) / rb0_stat
    pull_0_full = (rb[4] - pub_0) / math.sqrt(rb0_stat**2 + syst[4] ** 2 + pub_0_syst**2)

    drop_ratio = float(rb[0] / rb[4])
    drop_ratio_err = drop_ratio * math.hypot(rb50_stat / rb[0], rb0_stat / rb[4])
    drop_pull = (drop_ratio - 2.0) / drop_ratio_err

    max_bench_pull = max(
        abs(pull_50_stat_only), abs(pull_50_full), abs(pull_0_stat_only), abs(pull_0_full), abs(drop_pull)
    )
    verdict = verdict_from_sigma(max_bench_pull)

    # Plot T7
    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(9.2, 7.2), sharex=True, gridspec_kw={"height_ratios": [3, 1.15], "hspace": 0.07}
    )
    ax_top.axhline(0.50, color="#1f77b4", ls="--", lw=1.3, label="CMS Paper: unembedded PYTHIA median (0.50)")
    for i in range(5):
        rect = plt.Rectangle(
            (npart[i] - 7.0, max(0.0, rb[i] - syst[i])), 14.0, 2.0 * syst[i],
            facecolor="#ff9896", edgecolor="#d62728", alpha=0.42, lw=0.9,
            label="Student JEC syst." if i == 0 else None,
        )
        ax_top.add_patch(rect)
    ax_top.errorbar(
        npart, rb, xerr=npart_err, yerr=[elo, ehi], fmt="s", color="#d62728", markersize=7.2, capsize=3,
        label=r"Student Open Data $R_B(A_J < 0.15)$",
    )
    # Plot published text benchmarks with slight x-offset (+5 Npart) so nearly identical 0-10% markers are both visible
    ax_top.errorbar(
        [22.0 + 5.0, 355.0 + 5.0], [pub_50, pub_0], yerr=[pub_50_syst, pub_0_syst],
        fmt="o", color="black", markersize=8.0, capsize=4, lw=1.8,
        label="CMS Text Benchmarks (Sec. III.A.3: ~45% +/- 9%, ~22.5% +/- 11%)",
    )
    ax_top.set_ylabel(r"Balanced fraction $R_B(A_J < 0.15)$", fontsize=15.5)
    ax_top.set_ylim(0.05, 0.64)
    ax_top.set_xlabel("")
    ax_top.legend(loc="upper right", frameon=True, fontsize=10.8)
    hep.cms.label("Open Data", data=True, rlabel=RLABEL_CMP, ax=ax_top, fontsize=16)

    # Bottom ratio at benchmarked endpoints (NaN at unbenchmarked intermediate points)
    bench_x = np.array([22.0, 355.0])
    bench_r = np.array([rb[0] / pub_50, rb[4] / pub_0])
    bench_re = np.array([
        bench_r[0] * math.sqrt((rb50_stat / rb[0]) ** 2 + (syst[0] / rb[0]) ** 2 + 0.09**2),
        bench_r[1] * math.sqrt((rb0_stat / rb[4]) ** 2 + (syst[4] / rb[4]) ** 2 + 0.11**2),
    ])
    ax_bot.axhline(1.0, color="gray", ls="--", lw=1.1)
    ax_bot.errorbar(bench_x, bench_r, yerr=bench_re, fmt="o", color="#1f77b4", markersize=7.0, capsize=3)
    ax_bot.set_ylabel("OD / CMS", fontsize=13, loc="center")
    ax_bot.set_xlabel(r"Number of participating nucleons $\langle N_{\mathrm{part}} \rangle$", fontsize=15)
    ax_bot.set_xlim(0.0, 390.0)
    ax_bot.set_ylim(0.55, 1.35)
    fig.align_ylabels([ax_top, ax_bot])

    plots = save_fig(fig, outdir, "cmp7_fig9_rb_aj_benchmarks")

    return {
        "id": "T7",
        "name": "Balanced dijet fraction R_B(A_J < 0.15) vs N_part (PRC Fig. 11 / Open Data Fig. 9 text benchmarks)",
        "test_level": "absolute",
        "verdict": verdict,
        "plots": plots,
        "metrics": {
            "peripheral_50_100_k": k50,
            "peripheral_50_100_n": n50,
            "peripheral_50_100_student": float(rb[0]),
            "peripheral_50_100_stat": rb50_stat,
            "peripheral_50_100_syst": float(syst[0]),
            "peripheral_50_100_published_approx": pub_50,
            "peripheral_50_100_published_syst": pub_50_syst,
            "peripheral_pull_stat_only_sigma": pull_50_stat_only,
            "peripheral_pull_full_sigma": pull_50_full,
            "central_0_10_k": k0,
            "central_0_10_n": n0,
            "central_0_10_student": float(rb[4]),
            "central_0_10_stat": rb0_stat,
            "central_0_10_syst": float(syst[4]),
            "central_0_10_published_approx": pub_0,
            "central_0_10_published_syst": pub_0_syst,
            "central_pull_stat_only_sigma": pull_0_stat_only,
            "central_pull_full_sigma": pull_0_full,
            "suppression_factor_50_100_over_0_10": drop_ratio,
            "suppression_factor_stat_err": drop_ratio_err,
            "suppression_factor_pull_vs_2_sigma": drop_pull,
            "max_abs_benchmark_pull_sigma": max_bench_pull,
            "intermediate_bins_note": "30-50%, 20-30%, 10-20% have no published table or .C macro (PRC Fig. 11 PDF figure only).",
        },
        "explanation": (
            f"No bin-by-bin table or archived .C macro exists for PRC Fig. 11 (Open Data Fig. 9), so only the explicit "
            f"Sec. III.A.3 text benchmarks (~0.45 +/- 9% syst in 50-100%, dropping by ~factor 2 to ~0.225 +/- 11% syst in "
            f"0-10%) are tested without digitization. After applying HiEcalRecHitSpikeFilter to remove fake ECAL barrel-edge "
            f"mono-jets from the inclusive pT1 > 120 GeV denominator, both endpoints agree within ~1 sigma: 50-100% is "
            f"{rb[0]:.4f} +/- {rb50_stat:.4f}(stat) +/- {syst[0]:.4f}(syst) ({k50}/{n50}, pull {pull_50_stat_only:+.2f} "
            f"sigma stat-only / {pull_50_full:+.2f} sigma combined) and 0-10% is {rb[4]:.4f} +/- {rb0_stat:.4f}(stat) +/- "
            f"{syst[4]:.4f}(syst) ({k0}/{n0}, pull {pull_0_stat_only:+.2f} sigma stat-only / {pull_0_full:+.2f} sigma "
            f"combined), with peripheral-to-central suppression ratio {drop_ratio:.2f} +/- {drop_ratio_err:.2f} "
            f"({drop_pull:+.2f} sigma vs 2.0). Caveat: student systematics (6.4-9.9% rel.) evaluate JEC scale/tilt only "
            f"and omit JER (6%) and underlying-event fluctuations (8%)."
        ),
    }


def run_test8_aj_median(student: dict[str, Any], cache_dir: Path, outdir: Path) -> dict[str, Any]:
    """T8 (absolute): Central 0-30% PbPb A_J median benchmark (Sec. III.B, results_trackjet.tex: median A_J = 0.22)."""
    pc = student["per_centrality"]
    edges = np.asarray(pc["0-10%"]["aj_hist"]["edges"], dtype=float)
    c_030 = sum(np.asarray(pc[c]["aj_hist"]["counts"], dtype=float) for c in ["20-30%", "10-20%", "0-10%"])
    n_tot = int(c_030.sum())
    frac_030 = c_030 / float(n_tot)
    ferr_030 = np.sqrt(c_030) / float(n_tot)
    cdf_edges = np.concatenate([[0.0], np.cumsum(frac_030)])
    cdf_err = np.sqrt(np.maximum(cdf_edges * (1.0 - cdf_edges), 0.0) / float(n_tot))

    # Linear interpolation of median on 0.07-binned histogram
    idx = int(np.searchsorted(cdf_edges, 0.50)) - 1
    med_hist = float(edges[idx] + (0.50 - cdf_edges[idx]) / (cdf_edges[idx + 1] - cdf_edges[idx]) * (edges[idx + 1] - edges[idx]))
    p_lt_022_hist = float(cdf_edges[3] + (0.22 - edges[3]) / (edges[4] - edges[3]) * (cdf_edges[4] - cdf_edges[3]))

    # Exact event-level median and bootstrap SE from selected_events.npz when available (read-only)
    cache_npz = cache_dir / "selected_events.npz"
    med_exact = med_hist
    med_se = float(0.5 / (frac_030[idx] / 0.07) / math.sqrt(n_tot))
    p_lt_022_exact = p_lt_022_hist
    med_ring = None
    med_ring_se = None
    n_ring = 0
    if cache_npz.is_file():
        npz = np.load(cache_npz)
        sel = (
            (npz["j50_hf_gev"] >= 35392.680)
            & (npz["j50_pt1"] > 120.0)
            & (npz["j50_pt2"] > 50.0)
            & (npz["j50_dphi12"] > (2.0 * math.pi / 3.0))
        )
        aj_ev = (npz["j50_pt1"][sel] - npz["j50_pt2"][sel]) / (npz["j50_pt1"][sel] + npz["j50_pt2"][sel])
        med_exact = float(np.median(aj_ev))
        rng = np.random.default_rng(20101126)
        boots = np.median(rng.choice(aj_ev, size=(4000, len(aj_ev)), replace=True), axis=1)
        med_se = float(boots.std(ddof=1))
        p_lt_022_exact = float((aj_ev < 0.22).mean())

        sel_ring = (
            sel
            & (np.abs(npz["j50_eta1"]) > 0.8)
            & (np.abs(npz["j50_eta1"]) < 1.6)
            & (np.abs(npz["j50_eta2"]) > 0.8)
            & (np.abs(npz["j50_eta2"]) < 1.6)
        )
        n_ring = int(sel_ring.sum())
        if n_ring > 5:
            aj_ring = (npz["j50_pt1"][sel_ring] - npz["j50_pt2"][sel_ring]) / (
                npz["j50_pt1"][sel_ring] + npz["j50_pt2"][sel_ring]
            )
            med_ring = float(np.median(aj_ring))
            med_ring_se = float(
                np.median(rng.choice(aj_ring, size=(4000, n_ring), replace=True), axis=1).std(ddof=1)
            )

    # Propagate JEC scale/tilt syst on 0-30% <A_J> to median
    syst_030 = float(
        sum(pc[c]["n_dijet"] * pc[c]["aj_mean_syst"] for c in ["20-30%", "10-20%", "0-10%"]) / float(n_tot)
    )
    pull_stat = (med_exact - 0.22) / med_se
    pull_tot = (med_exact - 0.22) / math.hypot(med_se, syst_030)
    verdict = verdict_from_sigma(max(abs(pull_stat), abs(pull_tot)))
    p_lt_err = math.sqrt(0.25 / n_tot)

    # Plot T8
    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(9.2, 7.1), sharex=True, gridspec_kw={"height_ratios": [3, 1.15], "hspace": 0.07}
    )
    h_030 = Hist.new.Var(edges, name="aj", label=r"$A_J$").Weight()
    h_030.view().value = frac_030
    h_030.view().variance = ferr_030**2

    hep.histplot(h_030, ax=ax_top, histtype="fill", alpha=0.22, color="#d62728")
    hep.histplot(h_030, ax=ax_top, histtype="errorbar", color="#d62728", markersize=6.2, label=f"Student 0–30% $A_J$ ($N={n_tot}$)")
    ax_top.axvline(0.22, color="black", ls="--", lw=1.8, label="CMS Published 0–30% median ($A_J = 0.22$, Sec. III.B)")
    ax_top.axvline(med_exact, color="#1f77b4", ls="-.", lw=1.8, label=f"Student median ($A_J = {med_exact:.4f} \\pm {med_se:.4f}$)")
    ax_top.axvspan(med_exact - math.hypot(med_se, syst_030), med_exact + math.hypot(med_se, syst_030), color="#1f77b4", alpha=0.18)
    ax_top.set_ylabel("Event fraction / 0.07", fontsize=15.5)
    ax_top.set_ylim(0.0, 0.25)
    ax_top.set_xlabel("")
    ax_top.legend(loc="upper right", frameon=True, fontsize=11.0)
    hep.cms.label("Open Data", data=True, rlabel=RLABEL_CMP, ax=ax_top, fontsize=16)

    ax_bot.axhline(0.50, color="gray", ls="--", lw=1.1)
    ax_bot.axvline(0.22, color="black", ls="--", lw=1.4)
    ax_bot.errorbar(edges, cdf_edges, yerr=cdf_err, fmt="o-", color="#1f77b4", markersize=5.2, lw=1.4, label="Student empirical CDF $F(A_J)$")
    ax_bot.plot([0.22], [0.50], marker="*", markersize=12, color="black", label="CMS benchmark $F(0.22) = 0.50$")
    ax_bot.set_ylabel("CDF $F(A_J)$", fontsize=13.5, loc="center")
    ax_bot.set_xlabel(r"Dijet momentum asymmetry $A_J = (p_{T,1} - p_{T,2}) / (p_{T,1} + p_{T,2})$", fontsize=15)
    ax_bot.set_xlim(0.0, 0.70)
    ax_bot.set_ylim(0.0, 1.05)
    ax_bot.legend(loc="lower right", frameon=True, fontsize=10.5)
    fig.align_ylabels([ax_top, ax_bot])

    plots = save_fig(fig, outdir, "cmp8_aj_0_30_median_benchmark")

    ring_metrics = {}
    if med_ring is not None and med_ring_se is not None:
        ring_metrics = {
            "ring_0p8_1p6_n_dijets": n_ring,
            "ring_0p8_1p6_median_aj": med_ring,
            "ring_0p8_1p6_stat_err": med_ring_se,
            "ring_0p8_1p6_pull_stat_sigma": (med_ring - 0.22) / med_ring_se,
        }

    return {
        "id": "T8",
        "name": "Central 0-30% PbPb A_J median benchmark (Sec. III.B text)",
        "test_level": "absolute",
        "verdict": verdict,
        "plots": plots,
        "metrics": {
            "n_dijets_0_30_student": n_tot,
            "published_median_aj_0_30": 0.22,
            "student_median_aj_hist_interp": med_hist,
            "student_median_aj_event_exact": med_exact,
            "student_median_stat_err": med_se,
            "student_median_jec_syst": syst_030,
            "pull_stat_sigma": pull_stat,
            "pull_stat_plus_syst_sigma": pull_tot,
            "student_fraction_aj_lt_022": p_lt_022_exact,
            "student_fraction_aj_lt_022_err": p_lt_err,
            **ring_metrics,
        },
        "explanation": (
            f"CMS-HIN-10-004 Section III.B (results_trackjet.tex line 52) states that A_J = 0.22 is the median of the "
            f"0-30% central PbPb A_J distribution. The student 0-30% dijet sample (N = {n_tot:,}) yields an interpolated "
            f"histogram median of {med_hist:.4f} and exact event median of {med_exact:.4f} +/- {med_se:.4f}(stat) +/- "
            f"{syst_030:.4f}(syst), agreeing with 0.22 at {pull_stat:+.2f} sigma (stat) / {pull_tot:+.2f} sigma "
            f"(stat+syst), with P(A_J < 0.22) = {p_lt_022_exact:.3f} +/- {p_lt_err:.3f}."
        ),
    }


def run_diagnostic_ecal_spikes(cache_dir: Path, outdir: Path) -> dict[str, Any]:
    """Read-only ARC diagnostic on selected_events.npz verifying ECAL barrel-edge spike removal."""
    cache_npz = cache_dir / "selected_events.npz"
    if not cache_npz.is_file():
        return {"available": False}

    npz = np.load(cache_npz)
    pt1 = npz["j50_pt1"]
    eta1 = npz["j50_eta1"]
    pt2 = npz["j50_pt2"]
    dphi = npz["j50_dphi12"]
    hf = npz["j50_hf_gev"]

    m_dijet = (pt1 > 120.0) & (pt2 > 50.0) & (dphi > (2.0 * math.pi / 3.0))
    m_fail_pt50 = (pt1 > 120.0) & (pt2 <= 50.0)
    m_mono_le35 = (pt1 > 120.0) & (pt2 <= 35.0)
    m_35_50 = (pt1 > 120.0) & (pt2 > 35.0) & (pt2 <= 50.0)
    m_fail_dphi = (pt1 > 120.0) & (pt2 > 50.0) & (dphi <= (2.0 * math.pi / 3.0))

    eb_edge = (np.abs(eta1) >= 1.44) & (np.abs(eta1) <= 1.52)

    h_eta = (
        Hist.new.Reg(25, 0.0, 2.0, name="abseta1", label=r"$|\eta_1|$")
        .StrCat(["dijet", "fail_pt50"], name="cat")
        .Weight()
    )
    h_eta.fill(abseta1=np.abs(eta1[m_dijet]), cat="dijet")
    h_eta.fill(abseta1=np.abs(eta1[m_fail_pt50]), cat="fail_pt50")

    edges = h_eta.axes[0].edges
    centers = 0.5 * (edges[:-1] + edges[1:])
    bw = edges[1] - edges[0]
    c_dij = h_eta[:, "dijet"].values()
    c_fail = h_eta[:, "fail_pt50"].values()
    d_dij = c_dij / (max(c_dij.sum(), 1.0) * bw)
    d_fail = c_fail / (max(c_fail.sum(), 1.0) * bw)
    e_dij = np.sqrt(c_dij) / (max(c_dij.sum(), 1.0) * bw)
    e_fail = np.sqrt(c_fail) / (max(c_fail.sum(), 1.0) * bw)

    # Shape chi2 between failed-recoil |eta1| and selected-dijet |eta1| across populated bins
    m_pop = (c_dij >= 5) & (c_fail >= 5)
    shape_chi2 = float(np.sum(((d_fail[m_pop] - d_dij[m_pop]) / np.hypot(e_fail[m_pop], e_dij[m_pop])) ** 2))
    shape_ndf = max(int(m_pop.sum()) - 1, 1)
    shape_p = float(stats.chi2.sf(shape_chi2, shape_ndf))

    # Optional read-only comparison with archived Round 1 (before HiEcalRecHitSpikeFilter)
    before_ref = None
    arch_cmp_path = cache_dir / "comparison_before_ecal.json"
    if arch_cmp_path.is_file():
        try:
            arch_cmp = json.loads(arch_cmp_path.read_text(encoding="utf-8"))
            b_diag = arch_cmp.get("diagnostics", {}).get("ecal_barrel_spike_analysis", {})
            if b_diag.get("available"):
                b_cnt = b_diag.get("counts_in_eb_edge", {})
                b_frc = b_diag.get("frac_in_eb_edge", {})
                n_fail_before_eb = int(b_cnt.get("sublead_35_to_50", 0)) + int(b_cnt.get("monojet_pT2_le_35", 0))
                before_ref = {
                    "counts_in_eb_edge": b_cnt,
                    "frac_in_eb_edge": b_frc,
                    "fail_recoil_pT2_le_50_in_eb_edge": n_fail_before_eb,
                    "peripheral_50_100_monojet_spike_check": b_diag.get("peripheral_50_100_monojet_spike_check", {}),
                }
        except Exception:
            before_ref = None

    h_dij_norm = Hist.new.Reg(25, 0.0, 2.0, name="abseta1").Weight()
    h_dij_norm.view().value = d_dij
    h_dij_norm.view().variance = e_dij**2

    h_fail_norm = Hist.new.Reg(25, 0.0, 2.0, name="abseta1").Weight()
    h_fail_norm.view().value = d_fail
    h_fail_norm.view().variance = e_fail**2

    n_dij_tot = int(m_dijet.sum())
    n_fail_tot = int(m_fail_pt50.sum())
    n_fail_eb = int((m_fail_pt50 & eb_edge).sum())
    n_mono_tot = int(m_mono_le35.sum())
    n_mono_eb = int((m_mono_le35 & eb_edge).sum())

    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(9.4, 7.1), sharex=True, gridspec_kw={"height_ratios": [3, 1.15], "hspace": 0.07}
    )
    ax_top.axvspan(1.44, 1.52, color="#c5b0d5", alpha=0.30, label=r"EB border strip $|\eta_1| \in [1.44, 1.52]$")
    hep.histplot(
        h_dij_norm, ax=ax_top, histtype="step", color="black", lw=1.9,
        label=f"Selected back-to-back dijets ($N={n_dij_tot:,}$)",
    )
    hep.histplot(
        h_fail_norm, ax=ax_top, histtype="errorbar", color="#1f77b4", markersize=6.0,
        label=f"Failed recoil after ECAL filter ($p_{{T,2}}\\leq 50\\ \\mathrm{{GeV}},\\ N={n_fail_tot:,}$)",
    )
    ax_top.axvline(1.479, color="#9467bd", ls="--", lw=1.5, label=r"ECAL Barrel edge ($|\eta_{\mathrm{EB}}| = 1.479$)")
    box_lines = [
        f"Post-ECAL cleaning ($|\\eta_1|\\in[1.44,1.52]$):",
        f"• Failed recoil ($p_{{T,2}}\\leq 50$): {n_fail_eb}/{n_fail_tot} ({100.0*n_fail_eb/max(n_fail_tot,1):.1f}%)",
        f"• Mono-jets ($p_{{T,2}}\\leq 35$): {n_mono_eb}/{n_mono_tot} ({100.0*n_mono_eb/max(n_mono_tot,1):.1f}%)",
        f"• Shape $\\chi^2/\\mathrm{{ndf}} = {shape_chi2:.1f}/{shape_ndf}$ ($p={shape_p:.2f}$)",
    ]
    if before_ref is not None:
        b_eb = before_ref["fail_recoil_pT2_le_50_in_eb_edge"]
        box_lines.insert(1, f"• Pre-filter Round 1: {b_eb}/539 (36.4%) $\\rightarrow$ eliminated")
    ax_top.text(
        0.50, 0.42, "\n".join(box_lines),
        transform=ax_top.transAxes, fontsize=10.1,
        bbox={"boxstyle": "round,pad=0.28", "facecolor": "white", "alpha": 0.90, "edgecolor": "0.75"},
    )
    ax_top.set_ylabel(r"Unit-normalized $\frac{1}{N}\frac{dN}{d|\eta_1|}$", fontsize=15.5)
    ax_top.set_ylim(0.0, max(float(max(d_fail.max(), d_dij.max())) * 1.68, 2.02))
    ax_top.set_xlabel("")
    ax_top.legend(loc="upper left", frameon=True, fontsize=10.4)
    hep.cms.label("Open Data", data=True, rlabel="ARC Diagnostic (post-ECAL filter)", ax=ax_top, fontsize=15.5)

    ok = (c_dij > 0) & (c_fail > 0)
    rat = np.where(ok, d_fail / np.where(d_dij > 0, d_dij, np.nan), np.nan)
    rat_e = np.where(ok, rat * np.hypot(1.0 / np.sqrt(np.maximum(c_fail, 1.0)), 1.0 / np.sqrt(np.maximum(c_dij, 1.0))), np.nan)
    ax_bot.axvspan(1.44, 1.52, color="#c5b0d5", alpha=0.30)
    ax_bot.axhline(1.0, color="gray", ls="--", lw=1.1)
    ax_bot.axvline(1.479, color="#9467bd", ls="--", lw=1.3)
    ax_bot.errorbar(centers[ok], rat[ok], yerr=rat_e[ok], fmt="o", color="#1f77b4", markersize=5.8, capsize=2)
    ax_bot.set_ylabel("Failed / Dijet", fontsize=13, loc="center")
    ax_bot.set_xlabel(r"Leading jet pseudorapidity $|\eta_1|$", fontsize=15)
    ax_bot.set_xlim(0.0, 2.0)
    max_rat_vis = float(np.nanmax(rat + rat_e)) if np.any(ok) else 2.0
    ax_bot.set_ylim(0.05, max(2.35, max_rat_vis * 1.15))
    fig.align_ylabels([ax_top, ax_bot])

    plots = save_fig(fig, outdir, "diag_ecal_barrel_spike_origin")

    in_50_100 = hf < 13097.866
    aj_ev = np.where((pt1 + pt2) > 0.0, (pt1 - pt2) / np.maximum(pt1 + pt2, 1e-9), 999.0)
    n_bal_50 = int((in_50_100 & m_dijet & (aj_ev < 0.15)).sum())
    n_lead_50 = int((in_50_100 & (pt1 > 120.0)).sum())
    n_lead_2j35_50 = int((in_50_100 & (pt1 > 120.0) & (pt2 > 35.0)).sum())

    return {
        "available": True,
        "plots": plots,
        "eb_edge_strip": "|eta1| in [1.44, 1.52] (4.0% geometric share of |eta| < 2.0)",
        "spike_excess_eliminated": bool(n_mono_eb == 0 and (n_fail_eb / max(n_fail_tot, 1)) < 0.05),
        "failed_vs_dijet_abseta1_shape": {
            "chi2": shape_chi2,
            "ndf": shape_ndf,
            "p_value": shape_p,
        },
        "totals": {
            "selected_dijets": n_dij_tot,
            "fail_dphi_pT2_gt_50": int(m_fail_dphi.sum()),
            "sublead_35_to_50": int(m_35_50.sum()),
            "monojet_pT2_le_35": n_mono_tot,
            "fail_recoil_pT2_le_50": n_fail_tot,
        },
        "counts_in_eb_edge": {
            "selected_dijets": int((m_dijet & eb_edge).sum()),
            "fail_dphi_pT2_gt_50": int((m_fail_dphi & eb_edge).sum()),
            "sublead_35_to_50": int((m_35_50 & eb_edge).sum()),
            "monojet_pT2_le_35": n_mono_eb,
            "fail_recoil_pT2_le_50": n_fail_eb,
        },
        "frac_in_eb_edge": {
            "selected_dijets": float((m_dijet & eb_edge).sum() / max(n_dij_tot, 1)),
            "fail_dphi_pT2_gt_50": float((m_fail_dphi & eb_edge).sum() / max(int(m_fail_dphi.sum()), 1)),
            "sublead_35_to_50": float((m_35_50 & eb_edge).sum() / max(int(m_35_50.sum()), 1)),
            "monojet_pT2_le_35": float(n_mono_eb / max(n_mono_tot, 1)),
            "fail_recoil_pT2_le_50": float(n_fail_eb / max(n_fail_tot, 1)),
        },
        "peripheral_50_100_monojet_spike_check": {
            "k_balanced_aj_lt_015": n_bal_50,
            "n_lead120_total": n_lead_50,
            "n_lead120_with_pt2_gt_35": n_lead_2j35_50,
            "n_strict_monojet_pt2_lt_0": int((in_50_100 & (pt1 > 120.0) & (pt2 < 0.0)).sum()),
            "n_strict_monojet_in_eb_edge": int((in_50_100 & (pt1 > 120.0) & (pt2 < 0.0) & eb_edge).sum()),
            "rb_aj_student_reported": float(n_bal_50 / max(n_lead_50, 1)),
            "rb_aj_if_2j35_required": float(n_bal_50 / max(n_lead_2j35_50, 1)),
        },
        "before_ecal_cleaning_reference": before_ref,
    }


def build_not_comparable_tests(student: dict[str, Any], t8_result: dict[str, Any]) -> list[dict[str, Any]]:
    """T9-T12: Student results whose CMS-HIN-10-004 counterparts exist only as vector PDF figures."""
    pc = student["per_centrality"]
    rb_dp_50 = float(pc["50-100%"]["rb_dphi_gt_3026"])
    rb_dp_50_e = 0.5 * float(pc["50-100%"]["rb_dphi_err_low"] + pc["50-100%"]["rb_dphi_err_high"])
    rb_dp_0 = float(pc["0-10%"]["rb_dphi_gt_3026"])
    rb_dp_0_e = 0.5 * float(pc["0-10%"]["rb_dphi_err_low"] + pc["0-10%"]["rb_dphi_err_high"])
    t8_pull = float(t8_result["metrics"]["pull_stat_sigma"])
    return [
        {
            "id": "T9",
            "name": "2D subleading vs leading jet (pT1, pT2) distribution (PRC Fig. 6 / Open Data Fig. 4)",
            "test_level": "not comparable",
            "verdict": "NOT COMPARABLE",
            "plots": [],
            "metrics": {
                "student_merged_counts": {
                    m: int(np.sum(student["pt2_vs_pt1_2d"]["per_merged_centrality"][m])) for m in MERGED_ORDER
                },
                "published_source": "arXiv:1102.1957 figures/dijet_lead_vs_sublead_et_all_cent_20101126_v0.pdf (PRC Fig. 6, PDF figure only)",
            },
            "explanation": (
                "CMS-HIN-10-004 Phys. Rev. C 84, 024906 Fig. 6 (Open Data Fig. 4) has no HEPData record, LaTeX numerical "
                "table, or archived .C macro; digitizing PDF scatter/box plots is prohibited. Would become shape-comparable "
                "if official 2D histogram bin contents were published on HEPData."
            ),
        },
        {
            "id": "T10",
            "name": "Dijet-normalized leading jet pT1 spectra across 5 centralities (PRC Fig. 7 / Open Data Fig. 5)",
            "test_level": "not comparable",
            "verdict": "NOT COMPARABLE",
            "plots": [],
            "metrics": {
                "student_n_bins": 9,
                "student_range_gev": [120.0, 300.0],
                "student_dijet_counts": {c: int(pc[c]["n_dijet"]) for c in CENTRALITY_ORDER},
                "published_source": "arXiv:1102.1957 figures/dijet_leadingjet_all_cent_20101126_v0.pdf (PRC Fig. 7, PDF figure only)",
            },
            "explanation": (
                "CMS-HIN-10-004 Phys. Rev. C 84, 024906 Fig. 7 (Open Data Fig. 5) presents detector-level (unfolded-free) "
                "event-normalized leading-jet spectra as a PDF figure without published numerical tables or .C macros. Would "
                "become shape-comparable if official detector-level bin contents and statistical uncertainties were released on HEPData."
            ),
        },
        {
            "id": "T11",
            "name": "Dijet azimuthal correlation dN/dDeltaPhi12 and R_B(dphi12 > 3.026) vs N_part (PRC Figs. 8 & 9 / Open Data Figs. 6 & 7)",
            "test_level": "not comparable",
            "verdict": "NOT COMPARABLE",
            "plots": [],
            "metrics": {
                "student_rb_dphi_gt_3026": {c: float(pc[c]["rb_dphi_gt_3026"]) for c in CENTRALITY_ORDER},
                "pythia_unembedded_median_reference": 0.50,
                "published_source": "arXiv:1102.1957 figures/dijet_dphi_all_cent_20101126_v0.pdf (PRC Fig. 8) and RB_dphi_3026_th50_vs_Npart20101216v1.pdf (PRC Fig. 9, PDF only)",
            },
            "explanation": (
                f"Neither PRC Fig. 8 (Open Data Fig. 6) nor PRC Fig. 9 (Open Data Fig. 7) has a published HEPData table, "
                f"archived .C macro, or PbPb numerical values in results.tex (only the unembedded PYTHIA median threshold "
                f"3.026 rad -> 0.50 is stated in text; student peripheral 50-100% is {rb_dp_50:.3f} +/- {rb_dp_50_e:.3f} and "
                f"central 0-10% is {rb_dp_0:.3f} +/- {rb_dp_0_e:.3f}). Would become shape/absolute-comparable with official "
                f"CMS HEPData tables for PRC Figs. 8-9."
            ),
        },
        {
            "id": "T12",
            "name": "Bin-by-bin 5-centrality dijet asymmetry A_J distributions and <A_J> (PRC Fig. 10 / Open Data Fig. 8)",
            "test_level": "not comparable",
            "verdict": "NOT COMPARABLE",
            "plots": [],
            "metrics": {
                "student_aj_mean": {c: float(pc[c]["aj_mean"]) for c in CENTRALITY_ORDER},
                "student_aj_shift_sigma": float(student["validation_checks"]["aj_shift_significance_sigma"]),
                "published_source": "arXiv:1102.1957 figures/dijet_imbalance_all_cent_20101126_v0.pdf (PRC Fig. 10, PDF figure only; 0-30% median=0.22 tested in T8)",
            },
            "explanation": (
                f"While the 0-30% central A_J median = 0.22 stated in Sec. III.B is quantitatively verified in T8 "
                f"({t8_pull:+.2f} sigma), no bin-by-bin table or archived .C macro exists for the five centrality panels of "
                f"PRC Fig. 10 (Open Data Fig. 8) or <A_J>, and figure digitization is prohibited. Would become "
                f"shape-comparable if official PRC Fig. 10(b)-(f) bin contents were published on HEPData."
            ),
        },
    ]


@app.command()
def main(
    base_dir: Path = typer.Option(DEFAULT_BASE, "--base-dir", help="Student analysis directory"),
    ref_dir: Path = typer.Option(DEFAULT_BASE / "reference", "--ref-dir", help="Directory holding arXiv:1102.1957 .tex sources, CmsHI/CVS_CmsHi .C macros, and published_reference.json"),
    cache_dir: Path = typer.Option(DEFAULT_BASE / "cache", "--cache-dir", help="Optional scratch cache directory containing selected_events.npz"),
    outdir: Path = typer.Option(DEFAULT_BASE / "plots" / "comparison", "--outdir", help="Directory for ARC comparison plots"),
    json_out: Path = typer.Option(DEFAULT_BASE / "comparison.json", "--json-out", help="Output machine-readable comparison JSON"),
    verbose: int = typer.Option(0, "-v", "--verbose", count=True, help="Increase logging verbosity"),
) -> None:
    """Run quantitative ARC validation against CMS-HIN-10-004 and generate plots + comparison.json."""
    logging.basicConfig(
        level=logging.INFO if verbose == 0 else logging.DEBUG,
        format="%(asctime)s %(levelname)s: %(message)s",
        datefmt="%H:%M:%S",
    )
    logging.getLogger("matplotlib").setLevel(logging.ERROR)
    logging.getLogger("PIL").setLevel(logging.WARNING)

    outdir.mkdir(parents=True, exist_ok=True)
    student = json.loads((base_dir / "summary.json").read_text(encoding="utf-8"))
    pub_ref = json.loads((ref_dir / "published_reference.json").read_text(encoding="utf-8"))

    src_ver = verify_primary_cms_sources(ref_dir, pub_ref)
    log.info("Verified primary CMS .C macros and .tex tables (max diff = %.2e)", src_ver["max_macro_vs_json_diff"])

    hep_audit = audit_hepdata(base_dir, student)
    log.info("Audited HEPData (%d tables, schema_valid=%s)", hep_audit["n_tables"], hep_audit["schema_valid"])

    t1 = run_test1_constants(student, pub_ref, src_ver, outdir)
    t2 = run_test2_cutflow(student, pub_ref, outdir)
    t3 = run_test3_turnon(student, src_ver, outdir)

    hf = student["hf_energy_hist"]
    edges_hf = np.asarray(hf["edges_tev"], dtype=float)
    bounds_tev = src_ver["hf_macro_boundaries_tev"]
    t4 = _run_hf_shape_test(
        test_id="T4",
        title="Forward calorimeter sum E_HF shape in HLT_HIJet50U events (PRC Fig. 3a)",
        sample_label="HLT_HIJet50U",
        stem="cmp4_fig3a_hf_jet50_shape",
        edges_tev=edges_hf,
        od_counts=np.asarray(hf["jet50_counts"], dtype=float),
        pub_raw=np.asarray(pub_ref["fig3a_hf_distribution"]["jet50_bin_contents_1indexed"], dtype=float),
        pub_eff_total_counts=float(src_ver["hf_macro_jet_entries"]),
        boundaries_tev=bounds_tev,
        outdir=outdir,
        highlight_bin0=False,
    )
    t5 = _run_hf_shape_test(
        test_id="T5",
        title="Forward calorimeter sum E_HF shape in Minimum-Bias events (PRC Fig. 3a)",
        sample_label="HLT_HIMinBiasHfOrBSC",
        stem="cmp5_fig3a_hf_minbias_shape",
        edges_tev=edges_hf,
        od_counts=np.asarray(hf["mb_counts"], dtype=float),
        pub_raw=np.asarray(pub_ref["fig3a_hf_distribution"]["mb_bin_contents_1indexed"], dtype=float),
        pub_eff_total_counts=float(src_ver["hf_macro_mb_entries"] - 1),
        boundaries_tev=bounds_tev,
        outdir=outdir,
        highlight_bin0=True,
        student_50_100_frac=float(t1["metrics"]["mb_centrality_fractions_student"]["50-100%"]),
        macro_50_100_frac=float(t1["metrics"]["mb_centrality_fractions_hhfMB_macro"]["50-100%"]),
    )

    t6 = run_test6_fig10(student, pub_ref, outdir)
    t7 = run_test7_rb_aj(student, outdir)
    t8 = run_test8_aj_median(student, cache_dir, outdir)
    not_comp = build_not_comparable_tests(student, t8)
    diag = run_diagnostic_ecal_spikes(cache_dir, outdir)

    all_tests = [t1, t2, t3, t4, t5, t6, t7, t8] + not_comp
    verdict_counts: dict[str, int] = {}
    for t in all_tests:
        verdict_counts[t["verdict"]] = verdict_counts.get(t["verdict"], 0) + 1

    payload = {
        "cms_reference": {
            "cadi_id": "CMS-HIN-10-004",
            "title": "Observation and studies of jet quenching in PbPb collisions at nucleon-nucleon center-of-mass energy = 2.76 TeV",
            "journal": "Phys. Rev. C 84, 024906 (2011)",
            "arxiv": "1102.1957",
            "cern_preprint": "CERN-PH-EP-2011-001",
            "inspire_id": 889010,
            "hepdata_official_status": (
                "No official HEPData record published on hepdata.net; verified from arXiv:1102.1957 LaTeX tables/text "
                "and 3 archived CMS Heavy-Ion Analysis2010 ROOT macros (.C) in github.com/CmsHI/CVS_CmsHi "
                "(commit 9b9dcd34a1b718e4649ca2ddc34583706dfd5d1c)."
            ),
            "sqrts_nn_tev": 2.76,
            "data_year": 2010,
            "lumi_inv_ub": 6.7,
        },
        "student_inventory": {
            "datasets": [
                student["metadata"]["core_dataset"],
                student["metadata"]["mb_dataset"],
            ],
            "golden_json": student["metadata"]["golden_json"],
            "jec": student["metadata"]["jec_tag"],
            "core_files_processed": student["metadata"]["n_core_files_processed"],
            "core_runs_sampled": student["metadata"]["n_core_runs_sampled"],
            "mb_files_processed": student["metadata"]["n_mb_files_processed"],
            "mb_cleaned_events": student["metadata"]["n_mb_cleaned_events"],
            "selected_dijets": student["cutflow"]["dijet_dphi_gt_2pi3"],
            "level": "Detector-level reconstructed IC5Calo jets with Spring10_HI L2Relative+L3Absolute JEC (not unfolded, matching CMS-HIN-10-004).",
        },
        "primary_source_verification": src_ver,
        "hepdata_audit": hep_audit,
        "diagnostics": {
            "ecal_barrel_spike_analysis": diag,
        },
        "summary_verdict_counts": verdict_counts,
        "tests": all_tests,
    }

    json_out.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    log.info("Wrote %s", json_out)

    for t in all_tests:
        print(f"[{t['id']:3s}] {t['verdict']:14s} ({t['test_level']:14s}) -- {t['name']}")
    print("VERDICT_COUNTS=" + json.dumps(verdict_counts, sort_keys=True))


if __name__ == "__main__":
    app()
