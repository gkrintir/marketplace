#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11,<3.14"
# dependencies = [
#   "typer>=0.12",
#   "matplotlib>=3.8",
#   "mplhep>=0.3",
#   "hist>=2.7",
#   "numpy>=1.26",
#   "scipy>=1.11",
# ]
# ///
"""Publication-quality CMS Open Data figures for CMS-HIN-10-004 (2.76 TeV PbPb).

Reads summary.json (and optional cache/selected_events.npz) and renders all 9
figures in PNG and PDF formats under plots/, strictly labeled "CMS Open Data"
with sqrt(s_NN) = 2.76 TeV, PbPb 2010.
"""
from __future__ import annotations

import json
import logging
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import typer

plt.style.use(hep.style.CMS)

app = typer.Typer(add_completion=False)
log = logging.getLogger("plot_hin10004")

BASE_DIR = Path(__file__).resolve().parent
CENTRALITY_ORDER = ["50-100%", "30-50%", "20-30%", "10-20%", "0-10%"]
MERGED_ORDER = ["30-100%", "10-30%", "0-10%"]
RLABEL = r"PbPb 2010 ($\sqrt{s_{\mathrm{NN}}} = 2.76\ \mathrm{TeV}$)"

CENT_COLORS = {
    "50-100%": "#1f77b4",
    "30-50%": "#2ca02c",
    "20-30%": "#ff7f0e",
    "10-20%": "#9467bd",
    "0-10%": "#d62728",
}
CENT_MARKERS = {
    "50-100%": "o",
    "30-50%": "s",
    "20-30%": "^",
    "10-20%": "D",
    "0-10%": "p",
}


def save_both(fig: plt.Figure, plots_dir: Path, stem: str) -> None:
    png_path = plots_dir / f"{stem}.png"
    pdf_path = plots_dir / f"{stem}.pdf"
    fig.savefig(png_path, dpi=160, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved %s.{png,pdf}", stem)


def add_open_data_label(ax: plt.Axes, fontsize: int | None = None, loc: int = 0) -> None:
    kw = {"data": True, "rlabel": RLABEL, "ax": ax, "loc": loc}
    if fontsize is not None:
        kw["fontsize"] = fontsize
    hep.cms.label("Open Data", **kw)


def plot_fig1a_trigger_turnon(data: dict, plots_dir: Path) -> None:
    trig = data["trigger_turnon"]
    edges = np.asarray(trig["pt_edges"], dtype=float)
    centers = 0.5 * (edges[:-1] + edges[1:])
    half_w = 0.5 * np.diff(edges)

    fig, ax = plt.subplots(figsize=(8.5, 6.5))

    # Plot MinBias direct turn-on
    mb = trig.get("mb_only", trig)
    mb_den = np.asarray(mb["n_den"], dtype=int)
    mb_eff = np.asarray(mb["eff"], dtype=float)
    mb_elo = np.asarray(mb["err_low"], dtype=float)
    mb_ehi = np.asarray(mb["err_high"], dtype=float)
    ok_mb = mb_den > 0
    if ok_mb.any():
        ax.errorbar(
            centers[ok_mb],
            mb_eff[ok_mb],
            xerr=half_w[ok_mb],
            yerr=[mb_elo[ok_mb], mb_ehi[ok_mb]],
            fmt="o",
            color="#1f77b4",
            markersize=6.5,
            capsize=2,
            label=r"Ref: $\mathrm{HLT\_HIMinBiasHfOrBSC}$",
        )

    # Plot lower-threshold jet reference HLT_HIJet35U on its plateau (pT1 >= 80 GeV)
    j35 = trig.get("jet35_ref")
    if j35:
        j35_den = np.asarray(j35["n_den"], dtype=int)
        j35_eff = np.asarray(j35["eff"], dtype=float)
        j35_elo = np.asarray(j35["err_low"], dtype=float)
        j35_ehi = np.asarray(j35["err_high"], dtype=float)
        ok_j35 = (j35_den > 0) & (edges[:-1] >= 80.0)
        if ok_j35.any():
            ax.errorbar(
                centers[ok_j35] + 1.2,
                j35_eff[ok_j35],
                xerr=None,
                yerr=[j35_elo[ok_j35], j35_ehi[ok_j35]],
                fmt="s",
                color="black",
                markersize=6.0,
                capsize=2,
                label=r"Ref: $\mathrm{HLT\_HIJet35U}$ ($p_{T,1} \geq 80\ \mathrm{GeV}$)",
            )

    ax.axhline(1.0, ls=":", color="gray", lw=1.2)
    ax.axvline(100.0, ls="--", color="#2ca02c", lw=1.5, label=r"Full efficiency ($p_{T,1} > 100\ \mathrm{GeV}$)")
    ax.axvline(120.0, ls="-.", color="#d62728", lw=1.6, label=r"Analysis cut ($p_{T,1} > 120\ \mathrm{GeV}$)")

    ax.set_xlim(20.0, 260.0)
    ax.set_ylim(-0.03, 1.24)
    ax.set_xlabel(r"Corrected leading jet $p_{T,1}$ [GeV]")
    ax.set_ylabel(r"$\mathrm{HLT\_HIJet50U}$ trigger efficiency")
    ax.legend(loc="lower right", frameon=True, fontsize=13)
    ax.text(
        28.0,
        1.10,
        "IC5Calo jets (Spring10_HI L2L3), " + r"$|\eta_1| < 2.0$",
        fontsize=13,
    )
    add_open_data_label(ax, fontsize=18)
    save_both(fig, plots_dir, "fig1a_trigger_turnon")


def plot_fig3a_hf_centrality(data: dict, plots_dir: Path) -> None:
    hf = data["hf_energy_hist"]
    edges = np.asarray(hf["edges_tev"], dtype=float)
    centers = 0.5 * (edges[:-1] + edges[1:])
    mb_f = np.asarray(hf["mb_fraction"], dtype=float)
    j50_f = np.asarray(hf["jet50_fraction"], dtype=float)
    j50_c = np.asarray(hf.get("jet50_counts", np.ones_like(j50_f)), dtype=float)
    j50_tot = max(float(j50_c.sum()), 1.0)
    j50_err = np.sqrt(j50_c) / j50_tot

    fig, ax = plt.subplots(figsize=(9.0, 6.8))
    hep.histplot(
        mb_f,
        bins=edges,
        ax=ax,
        histtype="fill",
        alpha=0.28,
        color="#1f77b4",
        label=r"Minimum Bias ($\mathrm{HLT\_HIMinBiasHfOrBSC}$)",
    )
    hep.histplot(mb_f, bins=edges, ax=ax, histtype="step", lw=1.8, color="#1f77b4")

    ok = j50_f > 0
    ax.errorbar(
        centers[ok],
        j50_f[ok],
        yerr=j50_err[ok],
        fmt="o",
        color="#d62728",
        markersize=5.5,
        label=r"Jet-triggered ($\mathrm{HLT\_HIJet50U}$)",
    )

    # Centrality boundary vertical lines in TeV
    boundaries = [13.097866, 35.392680, 53.479020, 79.370000]
    for b in boundaries:
        ax.axvline(b, ls="--", color="0.45", lw=1.1)

    ax.set_yscale("log")
    ax.set_xlim(0.0, 160.0)
    ax.set_ylim(3e-5, 1.2)
    ax.set_xlabel(r"Total forward calorimeter energy $\sum E_{\mathrm{HF}}$ [TeV]")
    ax.set_ylabel("Event fraction / 2.29 TeV")

    label_pos = [
        (6.5, "50–100%"),
        (24.0, "30–50%"),
        (44.4, "20–30%"),
        (66.4, "10–20%"),
        (105.0, "0–10%"),
    ]
    for x_p, txt in label_pos:
        ax.text(x_p, 2.2e-4, txt, ha="center", va="bottom", fontsize=11, rotation=90, color="0.25")

    ax.legend(loc="upper right", frameon=True, fontsize=13)
    add_open_data_label(ax, fontsize=18)
    save_both(fig, plots_dir, "fig3a_hf_centrality")


def plot_fig4_pt2_vs_pt1(data: dict, plots_dir: Path) -> None:
    d2 = data["pt2_vs_pt1_2d"]
    xedges = np.asarray(d2["pt1_edges"], dtype=float)
    yedges = np.asarray(d2["pt2_edges"], dtype=float)

    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.6), sharey=True)
    fig.subplots_adjust(wspace=0.08)

    vmax = max(
        max(1.0, float(np.max(d2["per_merged_centrality"][m]))) for m in MERGED_ORDER
    )
    mesh = None
    for idx, mname in enumerate(MERGED_ORDER):
        ax = axes[idx]
        mat = np.asarray(d2["per_merged_centrality"][mname], dtype=float)
        mat_masked = np.ma.masked_where(mat <= 0, mat)
        mesh = ax.pcolormesh(xedges, yedges, mat_masked.T, cmap="viridis", vmin=1.0, vmax=vmax)
        ax.plot([120.0, 250.0], [120.0, 250.0], ls="--", color="#d62728", lw=1.4, label=r"$p_{T,2} = p_{T,1}$")
        ax.set_xlim(120.0, 250.0)
        ax.set_ylim(50.0, 250.0)
        ax.set_xlabel(r"Leading jet $p_{T,1}$ [GeV]", fontsize=15)
        if idx == 0:
            ax.set_ylabel(r"Subleading jet $p_{T,2}$ [GeV]", fontsize=15)
        ax.text(
            0.06,
            0.90,
            f"Centrality {mname}\n$N_{{\\mathrm{{dijet}}}} = {int(mat.sum())}$",
            transform=ax.transAxes,
            va="top",
            fontsize=13,
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "alpha": 0.85, "edgecolor": "0.7"},
        )
        if idx == 0:
            hep.cms.text(" Open Data", ax=ax, fontsize=14, loc=0)
        elif idx == 2:
            ax.set_title(RLABEL, loc="right", fontsize=14, pad=8)

    cbar = fig.colorbar(mesh, ax=axes.ravel().tolist(), pad=0.02, aspect=25)
    cbar.set_label("Selected dijet events / bin", fontsize=14)
    save_both(fig, plots_dir, "fig4_pt2_vs_pt1")


def plot_fig5_leading_jet_pt(data: dict, plots_dir: Path) -> None:
    pc = data["per_centrality"]
    fig, ax = plt.subplots(figsize=(8.8, 7.0))

    scale_factors = [1.0, 2.0, 4.0, 8.0, 16.0]
    for idx, cname in enumerate(CENTRALITY_ORDER):
        h = pc[cname]["pt1_hist"]
        edges = np.asarray(h["edges"], dtype=float)
        centers = 0.5 * (edges[:-1] + edges[1:])
        half_w = 0.5 * np.diff(edges)
        dens = np.asarray(h["density"], dtype=float) * scale_factors[idx]
        derr = np.asarray(h["density_err"], dtype=float) * scale_factors[idx]
        ok = dens > 0
        sf_txt = f" ($\\times {int(scale_factors[idx])}$)" if scale_factors[idx] > 1.0 else ""
        if ok.any():
            ax.errorbar(
                centers[ok],
                dens[ok],
                xerr=half_w[ok],
                yerr=derr[ok],
                fmt=CENT_MARKERS[cname],
                color=CENT_COLORS[cname],
                markersize=6.5,
                label=f"{cname}{sf_txt}",
            )

    ax.set_yscale("log")
    ax.set_xlim(120.0, 300.0)
    ax.set_ylim(2e-5, 2.5)
    ax.set_xlabel(r"Leading jet $p_{T,1}$ [GeV]")
    ax.set_ylabel(r"$\frac{1}{N_{\mathrm{dijet}}} \frac{dN}{dp_{T,1}}$ [GeV$^{-1}$]")
    ax.text(
        125.0,
        8e-5,
        r"$p_{T,1} > 120\ \mathrm{GeV},\ p_{T,2} > 50\ \mathrm{GeV}$" + "\n" + r"$|\eta_{1,2}| < 2.0,\ \Delta\phi_{12} > 2\pi/3$",
        fontsize=12.5,
    )
    ax.legend(loc="upper right", frameon=True, fontsize=12.5, title="Centrality")
    add_open_data_label(ax, fontsize=18)
    save_both(fig, plots_dir, "fig5_leading_jet_pt")


def plot_fig6_dijet_dphi(data: dict, plots_dir: Path) -> None:
    pc = data["per_centrality"]
    fig, axes = plt.subplots(2, 3, figsize=(15.5, 9.5), sharex=True, sharey=True)
    fig.subplots_adjust(wspace=0.08, hspace=0.10)

    # Panel (0, 0): info & legend panel matching CMS-HIN-10-004 multi-panel convention
    ax_info = axes[0, 0]
    ax_info.axis("off")
    ax_info.text(
        0.05,
        0.88,
        "CMS Open Data\n"
        + r"PbPb $\sqrt{s_{\mathrm{NN}}} = 2.76\ \mathrm{TeV}$ (2010)"
        + "\n\nIC5Calo Jets (L2L3 JEC)\n"
        + r"$p_{T,1} > 120\ \mathrm{GeV}$"
        + "\n"
        + r"$p_{T,2} > 50\ \mathrm{GeV}$"
        + "\n"
        + r"$|\eta_{1,2}| < 2.0$"
        + "\n\nDashed line: $\\Delta\\phi_{12} = 2\\pi/3$",
        transform=ax_info.transAxes,
        va="top",
        fontsize=13.5,
    )

    periph_dens = np.asarray(pc["50-100%"]["dphi_hist"]["density"], dtype=float)
    panel_slots = [(0, 1), (0, 2), (1, 0), (1, 1), (1, 2)]
    for idx, cname in enumerate(CENTRALITY_ORDER):
        r, c = panel_slots[idx]
        ax = axes[r, c]
        h = pc[cname]["dphi_hist"]
        edges = np.asarray(h["edges"], dtype=float)
        centers = 0.5 * (edges[:-1] + edges[1:])
        dens = np.asarray(h["density"], dtype=float)
        derr = np.asarray(h["density_err"], dtype=float)

        if cname != "50-100%":
            hep.histplot(
                np.maximum(periph_dens, 1e-4),
                bins=edges,
                ax=ax,
                histtype="step",
                color="#1f77b4",
                ls="--",
                lw=1.4,
                label="50–100% ref.",
            )

        ok = dens > 0
        if ok.any():
            ax.errorbar(
                centers[ok],
                dens[ok],
                yerr=derr[ok],
                fmt="o",
                color=CENT_COLORS[cname],
                markersize=6.0,
                label=f"PbPb {cname}",
            )
        ax.axvline(2.0 * math.pi / 3.0, ls=":", color="0.35", lw=1.4)
        ax.set_yscale("log")
        ax.set_xlim(0.0, math.pi)
        ax.set_ylim(2e-3, 18.0)
        ax.text(
            0.07,
            0.88,
            f"Centrality {cname}\n$N = {pc[cname]['n_dphi_sample']}$",
            transform=ax.transAxes,
            va="top",
            fontsize=13,
        )
        if ax.get_legend_handles_labels()[0]:
            ax.legend(loc="lower left", frameon=False, fontsize=11.5)
        if r == 1:
            ax.set_xlabel(r"$\Delta\phi_{12}$ [rad]", fontsize=15)
        if c == 0 or (r == 0 and c == 1):
            ax.set_ylabel(r"$\frac{1}{N} \frac{dN}{d\Delta\phi_{12}}$ [rad$^{-1}$]", fontsize=14)
        if r == 0 and c == 2:
            ax.set_title(RLABEL, loc="right", fontsize=14, pad=8)
        if r == 0 and c == 1:
            ax.tick_params(labelleft=True)
            hep.cms.text("Open Data", ax=ax, fontsize=15, loc=0)

    save_both(fig, plots_dir, "fig6_dijet_dphi")


def plot_fig7_rb_dphi_vs_npart(data: dict, plots_dir: Path) -> None:
    pc = data["per_centrality"]
    npart = np.array([pc[c]["npart_mean"] for c in CENTRALITY_ORDER], dtype=float)
    npart_err = np.array([pc[c]["npart_err"] for c in CENTRALITY_ORDER], dtype=float)
    rb = np.array([pc[c]["rb_dphi_gt_3026"] for c in CENTRALITY_ORDER], dtype=float)
    elo = np.array([pc[c]["rb_dphi_err_low"] for c in CENTRALITY_ORDER], dtype=float)
    ehi = np.array([pc[c]["rb_dphi_err_high"] for c in CENTRALITY_ORDER], dtype=float)

    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    ax.errorbar(
        npart,
        rb,
        xerr=npart_err,
        yerr=[elo, ehi],
        fmt="o",
        color="#1f77b4",
        markersize=8.0,
        capsize=3,
        lw=1.8,
        label=r"PbPb Open Data ($\Delta\phi_{12} > 3.026\ \mathrm{rad}$)",
    )
    for x_p, y_p, dy_hi, cname in zip(npart, rb, ehi, CENTRALITY_ORDER):
        ax.annotate(cname, (x_p, y_p + dy_hi), textcoords="offset points", xytext=(0, 6), ha="center", fontsize=11, color="0.25")

    ax.set_xlim(0.0, 400.0)
    ax.set_ylim(0.0, 1.05)
    ax.set_xlabel(r"Number of participating nucleons $N_{\mathrm{part}}$")
    ax.set_ylabel(r"Azimuthal back-to-back fraction $R_{B}(\Delta\phi_{12} > 3.026)$")
    ax.text(
        20.0,
        0.12,
        r"$p_{T,1} > 120\ \mathrm{GeV},\ p_{T,2} > 50\ \mathrm{GeV},\ |\eta_{1,2}| < 2.0$"
        + "\n68% Clopper–Pearson statistical intervals",
        fontsize=12.5,
    )
    ax.legend(loc="upper right", frameon=True, fontsize=13)
    add_open_data_label(ax, fontsize=18)
    save_both(fig, plots_dir, "fig7_rb_dphi_vs_npart")


def plot_fig8_dijet_asymmetry_aj(data: dict, plots_dir: Path) -> None:
    pc = data["per_centrality"]
    fig, axes = plt.subplots(2, 3, figsize=(15.5, 9.5), sharex=True, sharey=True)
    fig.subplots_adjust(wspace=0.08, hspace=0.10)

    ax_info = axes[0, 0]
    ax_info.axis("off")
    ax_info.text(
        0.05,
        0.90,
        "CMS Open Data\n"
        + r"PbPb $\sqrt{s_{\mathrm{NN}}} = 2.76\ \mathrm{TeV}$ (2010)"
        + "\n\nIC5Calo Jets (L2L3 JEC)\n"
        + r"$p_{T,1} > 120\ \mathrm{GeV}$"
        + "\n"
        + r"$p_{T,2} > 50\ \mathrm{GeV}$"
        + "\n"
        + r"$|\eta_{1,2}| < 2.0$"
        + "\n"
        + r"$\Delta\phi_{12} > 2\pi/3$"
        + "\n\nBlue shaded: 50–100% ref.",
        transform=ax_info.transAxes,
        va="top",
        fontsize=13.5,
    )

    periph_h = pc["50-100%"]["aj_hist"]
    periph_edges = np.asarray(periph_h["edges"], dtype=float)
    periph_frac = np.asarray(periph_h["fraction"], dtype=float)

    ymax = max(
        0.42,
        1.32 * max(max(pc[c]["aj_hist"]["fraction"]) for c in CENTRALITY_ORDER),
    )

    panel_slots = [(0, 1), (0, 2), (1, 0), (1, 1), (1, 2)]
    for idx, cname in enumerate(CENTRALITY_ORDER):
        r, c = panel_slots[idx]
        ax = axes[r, c]
        h = pc[cname]["aj_hist"]
        edges = np.asarray(h["edges"], dtype=float)
        centers = 0.5 * (edges[:-1] + edges[1:])
        half_w = 0.5 * np.diff(edges)
        frac = np.asarray(h["fraction"], dtype=float)
        ferr = np.asarray(h["fraction_err"], dtype=float)

        # Reference 50-100% peripheral PbPb distribution
        hep.histplot(
            periph_frac,
            bins=periph_edges,
            ax=ax,
            histtype="fill",
            alpha=0.22,
            color="#1f77b4",
            label="50–100% peripheral ref." if cname != "50-100%" else None,
        )
        hep.histplot(periph_frac, bins=periph_edges, ax=ax, histtype="step", color="#1f77b4", lw=1.5)

        ok = (frac > 0) | (centers < 0.70)
        ax.errorbar(
            centers[ok],
            frac[ok],
            xerr=half_w[ok],
            yerr=ferr[ok],
            fmt="o",
            color="black" if cname != "50-100%" else "#1f77b4",
            markersize=6.2,
            capsize=2,
            label=f"PbPb {cname}",
        )

        ax.set_xlim(0.0, 0.70)
        ax.set_ylim(0.0, ymax)
        aj_m = pc[cname]["aj_mean"]
        aj_e = pc[cname]["aj_mean_err"]
        ax.text(
            0.95,
            0.90,
            f"Centrality {cname}\n$\\langle A_J \\rangle = {aj_m:.3f} \\pm {aj_e:.3f}$\n$N_{{\\mathrm{{dijet}}}} = {pc[cname]['n_dijet']}$",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=12.5,
        )
        ax.legend(loc="upper left", frameon=False, fontsize=11)
        if r == 1:
            ax.set_xlabel(r"Dijet asymmetry $A_{J} = (p_{T,1} - p_{T,2}) / (p_{T,1} + p_{T,2})$", fontsize=13.5)
        if c == 0 or (r == 0 and c == 1):
            ax.set_ylabel("Event fraction", fontsize=14.5)
        if r == 0 and c == 2:
            ax.set_title(RLABEL, loc="right", fontsize=14, pad=8)
        if r == 0 and c == 1:
            ax.tick_params(labelleft=True)
            hep.cms.text("Open Data", ax=ax, fontsize=15, loc=0)

    save_both(fig, plots_dir, "fig8_dijet_asymmetry_aj")


def plot_fig9_rb_aj_vs_npart(data: dict, plots_dir: Path) -> None:
    pc = data["per_centrality"]
    npart = np.array([pc[c]["npart_mean"] for c in CENTRALITY_ORDER], dtype=float)
    npart_err = np.array([pc[c]["npart_err"] for c in CENTRALITY_ORDER], dtype=float)
    rb = np.array([pc[c]["rb_aj_lt_015"] for c in CENTRALITY_ORDER], dtype=float)
    elo = np.array([pc[c]["rb_aj_err_low"] for c in CENTRALITY_ORDER], dtype=float)
    ehi = np.array([pc[c]["rb_aj_err_high"] for c in CENTRALITY_ORDER], dtype=float)
    syst = np.array([pc[c].get("rb_aj_syst", 0.0) for c in CENTRALITY_ORDER], dtype=float)

    fig, ax = plt.subplots(figsize=(8.6, 6.6))
    # Draw JEC systematic boxes
    for i in range(len(npart)):
        w_box = 12.0
        rect = plt.Rectangle(
            (npart[i] - w_box / 2.0, max(0.0, rb[i] - syst[i])),
            w_box,
            2.0 * syst[i],
            facecolor="#ff9896",
            edgecolor="#d62728",
            alpha=0.45,
            lw=1.0,
            label="JEC systematic uncertainty" if i == 0 else None,
        )
        ax.add_patch(rect)

    ax.errorbar(
        npart,
        rb,
        xerr=npart_err,
        yerr=[elo, ehi],
        fmt="s",
        color="#d62728",
        markersize=8.0,
        capsize=3,
        lw=1.8,
        label=r"PbPb Open Data ($A_{J} < 0.15$)",
    )
    for x_p, y_p, dy_hi, sy, cname in zip(npart, rb, ehi, syst, CENTRALITY_ORDER):
        ax.annotate(cname, (x_p, y_p + max(dy_hi, sy)), textcoords="offset points", xytext=(0, 6), ha="center", fontsize=11, color="0.25")

    ax.set_xlim(0.0, 400.0)
    ax.set_ylim(0.0, 0.72)
    ax.set_xlabel(r"Number of participating nucleons $N_{\mathrm{part}}$")
    ax.set_ylabel(r"Balanced dijet fraction $R_{B}(A_{J} < 0.15)$")
    ax.text(
        18.0,
        0.05,
        r"Denominator: all events with $p_{T,1} > 120\ \mathrm{GeV},\ |\eta_1| < 2.0$"
        + "\n"
        + r"Numerator: $p_{T,2} > 50\ \mathrm{GeV},\ |\eta_2| < 2.0,\ \Delta\phi_{12} > 2\pi/3,\ A_J < 0.15$",
        fontsize=12,
    )
    ax.legend(loc="upper right", frameon=True, fontsize=13)
    add_open_data_label(ax, fontsize=18)
    save_both(fig, plots_dir, "fig9_rb_aj_vs_npart")


def plot_fig10_fractional_imbalance_vs_pt1(data: dict, plots_dir: Path) -> None:
    fi = data["fractional_imbalance_vs_pt1"]
    edges = np.asarray(fi["pt1_edges"], dtype=float)
    centers = 0.5 * (edges[:-1] + edges[1:])
    half_w = 0.5 * np.diff(edges)

    fig, axes = plt.subplots(1, 3, figsize=(16.0, 5.6), sharey=True)
    fig.subplots_adjust(wspace=0.06)

    periph_cnt = np.asarray(fi["30-100%"]["counts"], dtype=int)
    periph_mean = np.asarray(fi["30-100%"]["mean"], dtype=float)
    periph_sem = np.asarray(fi["30-100%"]["sem"], dtype=float)
    ok_per = periph_cnt > 0

    m_colors = {"30-100%": "#1f77b4", "10-30%": "#ff7f0e", "0-10%": "#d62728"}
    for idx, mname in enumerate(MERGED_ORDER):
        ax = axes[idx]
        sub = fi[mname]
        cnt = np.asarray(sub["counts"], dtype=int)
        mean = np.asarray(sub["mean"], dtype=float)
        sem = np.asarray(sub["sem"], dtype=float)
        ok = cnt > 0

        if mname != "30-100%" and ok_per.any():
            ax.errorbar(
                centers[ok_per],
                periph_mean[ok_per],
                xerr=half_w[ok_per],
                yerr=periph_sem[ok_per],
                fmt="o",
                mfc="white",
                mec="#1f77b4",
                ecolor="#1f77b4",
                markersize=6.5,
                alpha=0.75,
                label="PbPb 30–100% ref.",
            )

        if ok.any():
            ax.errorbar(
                centers[ok],
                mean[ok],
                xerr=half_w[ok],
                yerr=sem[ok],
                fmt="s",
                color=m_colors[mname],
                markersize=7.0,
                capsize=2,
                label=f"PbPb {mname}",
            )

        ax.set_xlim(120.0, 240.0)
        ax.set_ylim(0.0, 0.68)
        ax.set_xlabel(r"Leading jet $p_{T,1}$ [GeV]", fontsize=15)
        if idx == 0:
            ax.set_ylabel(r"$\langle (p_{T,1} - p_{T,2}) / p_{T,1} \rangle$", fontsize=15)
        ax.text(
            0.07,
            0.90,
            f"Centrality {mname}\n$N_{{\\mathrm{{dijet}}}} = {int(cnt.sum())}$",
            transform=ax.transAxes,
            va="top",
            fontsize=13,
        )
        if ax.get_legend_handles_labels()[0]:
            ax.legend(loc="lower right", frameon=True, fontsize=12)
        if idx == 0:
            hep.cms.text("Open Data", ax=ax, fontsize=14.5, loc=0)
        elif idx == 2:
            ax.set_title(RLABEL, loc="right", fontsize=14, pad=8)

    save_both(fig, plots_dir, "fig10_fractional_imbalance_vs_pt1")


@app.command()
def main(
    summary_path: Path = typer.Option(BASE_DIR / "summary.json", "--summary", help="Input summary.json path"),
    plots_dir: Path = typer.Option(BASE_DIR / "plots", "--outdir", help="Output plots directory"),
    verbose: int = typer.Option(0, "-v", "--verbose", count=True),
) -> None:
    """Render all 9 CMS Open Data figures (PNG + PDF) for CMS-HIN-10-004."""
    logging.basicConfig(
        level=logging.INFO if verbose == 0 else logging.DEBUG,
        format="%(asctime)s %(levelname)s: %(message)s",
        datefmt="%H:%M:%S",
    )
    logging.getLogger("matplotlib").setLevel(logging.ERROR)
    logging.getLogger("PIL").setLevel(logging.WARNING)
    if not summary_path.is_file():
        raise typer.BadParameter(f"Missing {summary_path}; run analyze_hin10004.py first.")
    plots_dir.mkdir(parents=True, exist_ok=True)
    data = json.loads(summary_path.read_text(encoding="utf-8"))

    plot_fig1a_trigger_turnon(data, plots_dir)
    plot_fig3a_hf_centrality(data, plots_dir)
    plot_fig4_pt2_vs_pt1(data, plots_dir)
    plot_fig5_leading_jet_pt(data, plots_dir)
    plot_fig6_dijet_dphi(data, plots_dir)
    plot_fig7_rb_dphi_vs_npart(data, plots_dir)
    plot_fig8_dijet_asymmetry_aj(data, plots_dir)
    plot_fig9_rb_aj_vs_npart(data, plots_dir)
    plot_fig10_fractional_imbalance_vs_pt1(data, plots_dir)
    log.info("Successfully generated all 9 figures (18 files) in %s", plots_dir)


if __name__ == "__main__":
    app()
