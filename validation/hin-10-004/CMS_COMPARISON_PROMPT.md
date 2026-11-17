# Task for the cms-results-validator agent: compare with published CMS-HIN-10-004

Run headless from the repository root with Claude Code in a fresh session
(clean context, separate from the producer session):

```bash
tmux new -s hin10004-arc
{ sed -n '/^## Task/,$p' validation/hin-10-004/CMS_COMPARISON_PROMPT.md; \
  echo; echo "REF_DIR=<scratch>/hin10004_ref"; \
  echo "CACHE_DIR=<scratch>/hin10004_cache"; } | \
  claude -p --agent cms:cms-results-validator \
    --plugin-dir plugins/cms --plugin-dir plugins/hep-python-tools \
    --permission-mode acceptEdits --allowedTools "Bash Read Write Edit Glob Grep Skill TodoWrite" \
    > <scratch>/hin10004_validator_report.md 2>&1
```

Or invoke `cms-results-validator` as an independent clean-context subagent in
Google Antigravity after `./scripts/install_antigravity.sh`.

## Task

Act as an independent CMS Analysis Review Committee (ARC) validator using
`cms-results-validator`. Audit the Open Data replication in
`validation/hin-10-004/` against the published **CMS-HIN-10-004** result
(_"Observation and studies of jet quenching in PbPb collisions at
nucleon-nucleon center-of-mass energy = 2.76 TeV"_, Phys. Rev. C 84 (2011)
024906, arXiv:1102.1957, INSPIRE `889010`).

Do **not** modify any producer output (`summary.json`, `specification.md`,
`hepdata/`, `analyze_hin10004.py`, `plot_hin10004.py`, `make_hepdata.py`).
Keep downloaded reference tarballs and macros in `$REF_DIR` outside the
repository.

### Primary CMS reference sources to verify directly

`CMS-HIN-10-004` predates mandatory HEPData submission and has no official
record on `hepdata.net`. Verify that fact via the INSPIRE and HEPData APIs,
then extract all quantitative CMS reference numbers directly from two primary
sources (without any figure digitization):

1. Published `arXiv:1102.1957` LaTeX source archive (`https://arxiv.org/e-print/1102.1957`,
   `Phys. Rev. C 84, 024906 (2011)` — contains `.tex` files and vector `.pdf`
   figures `PRC Figs. 1–15`, with **no** `.C` macros in the arXiv tarball):
   - `event_selection.tex` (Table I): `HLT_HIJet50U` event selection cutflow
     (`149k` triggered -> `99.61%` BSC halo -> `74.98%` HF coincidence ->
     `98.97%` vertex -> `99.78%` pixel compatibility -> `97.66%` ECAL cleaning ->
     `99.97%` HCAL cleaning -> `71.9k` `>=2` jets `p_T > 35 GeV` -> `4,216`
     `p_{T,1} > 120 GeV` -> `3,684` `p_{T,2} > 50 GeV` -> `3,514` `Δφ_{12} > 2π/3`)
   - `centrality.tex` (Table II): Glauber `<N_part> ± syst` across the 5 centrality classes
   - `results.tex` (Section III.A.3) & `results_trackjet.tex` (Section III.B):
     balanced dijet fraction `R_B(A_J < 0.15)` benchmarks (`~0.45 ± 9%` in
     `50-100%`, `~0.225 ± 11%` in `0-10%`; `PRC Fig. 11`) and central `0-30%`
     `A_J` median (`0.22`; `PRC Fig. 10`)
2. Archived CMS Heavy-Ion `Analysis2010` ROOT `.C` figure macros in
   `github.com/CmsHI/CVS_CmsHi`
   (`https://github.com/CmsHI/CVS_CmsHi/tree/9b9dcd34a1b718e4649ca2ddc34583706dfd5d1c/Analysis2010/macros/JetAnalysis/plots`):
   - `TrigEffCorrected151058and151025.C` (`PRC Fig. 2a` / Open Data Fig. 1a
     `HLT_HIJet50U` turn-on; invert ROOT 5 `TGraphAsymmErrors::BayesDivide`
     shortest Bayesian HPD intervals to recover exact binomial `(k, n)` counts)
   - `hf_centrality_distribution_mb_jet50_20101126_v1.C` (`PRC Fig. 3a` / Open
     Data Fig. 3a `E_HF` histograms `hhfMB` and `hhfJET` plus exact centrality
     boundary `TLine` coordinates)
   - `deltaPtOverPt1VsPt.C` (saved Dec 13, 2010: pre-publication `<(p_{T,1} - p_{T,2})/p_{T,1}>`
     vs `p_{T,1}` canvas with `[120, 130, 140, 170, 240] GeV` bins at centers
     `125, 135, 155, 205 GeV`, `TH1::Divide` statistical errors, and `8.2%`
     systematic boxes across `30-100%`, `10-30%`, and `0-10%`; note that final
     published `PRC Fig. 12`, `figures/deltaPtOverPt1_all_cent_20101223v1.pdf`
     saved Dec 23, 2010, uses `[120, 140, 160, 180, 240] GeV` at centers
     `130, 150, 170, 210 GeV` and has no `.C` macro or numerical table)

Never digitize PDF figures (`PRC Figs. 6–12` have no published numerical tables;
mark bin-by-bin comparisons without `.C` macros `NOT COMPARABLE` and state what
official HEPData table would make them comparable).

### Tests to execute

State the comparability level first (`exact`, `shape`, `ratio`, `absolute`, or
`not comparable`), compute the test statistic (`χ²/ndf`, `p`-value, ratio, or
pull in `σ`), and assign a verdict (`CONSISTENT`, `TENSION`, `INCONSISTENT`, or
`NOT COMPARABLE`):

- **T1 (exact):** Collision constants, HF centrality boundaries, and Table II
  Glauber `<N_part>` values.
- **T2 (absolute):** Cutflow step efficiencies and luminosity-scaled dijet
  yield vs Table I (separating high-`p_T` dijet steps from intermediate event
  cleaning sub-steps).
- **T3 (absolute):** `HLT_HIJet50U` turn-on curve (`PRC Fig. 2a`) rebinned onto
  the common `[20, 160] GeV` grid vs `TrigEffCorrected151058and151025.C`.
- **T4 (shape):** Unit-normalized `E_HF` distribution in `HLT_HIJet50U` events
  (`PRC Fig. 3a`) vs `hhfJET` (`χ²/ndf` and linear ratio slope vs `E_HF`).
- **T5 (shape):** Unit-normalized `E_HF` distribution in Minimum-Bias events
  (`PRC Fig. 3a`) vs `hhfMB` (`χ²/ndf`, ratio slope, and ultra-peripheral bin 1
  sensitivity).
- **T6 (absolute):** `<(p_{T,1} - p_{T,2})/p_{T,1}>` vs `p_{T,1}` in `30-100%`,
  `10-30%`, and `0-10%` matched bin-for-bin on `[120, 130, 140, 170, 240] GeV`
  against the archived Dec 13, 2010 `CmsHI/CVS_CmsHi` `deltaPtOverPt1VsPt.C`
  canvas, alongside recording and overlaying the Open Data `[120, 140, 160, 180, 240] GeV`
  profile corresponding to published `PRC Fig. 12` (PDF-only, not digitized).
- **T7 (absolute):** `R_B(A_J < 0.15)` (`PRC Fig. 11`) in `50-100%` and `0-10%`
  vs Section III.A.3 text benchmarks (`~0.45 ± 9%` and `~0.225 ± 11%`).
- **T8 (absolute):** Unbinned `0-30%` `A_J` median (bootstrap SE from
  `$CACHE_DIR/selected_events.npz` if available) vs `0.22` (`PRC Fig. 10` text).
- **T9–T12 (not comparable):** Document `PRC Figs. 6, 7, 8–9, and 10` (Open Data
  `Figs. 4, 5, 6–7, and 8`, `NOT COMPARABLE` without digitization), audit
  `hepdata/` with `hepdata-validate`, and run barrel ECAL spike diagnostics
  (`EMF_1` and `(η_1, φ_1)` distributions).

### Deliverables in `validation/hin-10-004/`

- `validate_hin10004.py`: rerunnable PEP 723 script accepting `--base-dir`,
  `--ref-dir`, `--cache-dir`, `--outdir`, and `--json-out`
- `comparison.json`: machine-readable record of all 12 tests (using relative
  paths only; zero personal paths or usernames)
- `COMPARISON.md`: full ARC report with comparability matrix, per-test plots
  and statistics, HEPData audit, and Round-1 vs Round-2 ECAL spike diagnosis
- `plots/comparison/*.{png,pdf}`: CMS-style comparison + ratio/pull panels
