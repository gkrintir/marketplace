# Validation: Observation of jet quenching in 2010 PbPb collisions (CMS-HIN-10-004)

This validation accompanies the `heavy-ion-reco` skill PR and exercises both
`cms-analysis-coder` and `cms-results-validator` in a two-agent
producer/clean-context-reviewer loop on Run-1 CMS Heavy-Ion EDM `RECO` Open
Data. Together with `heavy-ion-reco`, it exercises `analysis-spec-builder`,
`awkward-array`, `vector-awkward`, `trigger-efficiency`, `hist`, `hepdata`,
`standalone-script`, and `cli-creator`.

## Reference analysis

- **CMS Collaboration**, "Observation and studies of jet quenching in PbPb
  collisions at nucleon-nucleon center-of-mass energy = 2.76 TeV",
  `CMS-HIN-10-004`, Phys. Rev. C 84 (2011) 024906, arXiv:1102.1957
  (CERN-PH-EP-2011-001; INSPIRE 889010; 6.7 μb⁻¹ of 2010 PbPb data at
  √s_NN = 2.76 TeV).
- Primary reference artifacts: `CMS-HIN-10-004` has no official `hepdata.net`
  entry, so quantitative reference numbers come directly from the published
  LaTeX tables in `arXiv:1102.1957` (Tables I–II and Section III benchmarks)
  and the three official CMS ROOT `.C` figure macros in the `arXiv:1102.1957`
  source archive (`TrigEffCorrected151058and151025.C`,
  `hf_centrality_distribution_mb_jet50_20101126_v1.C`, and
  `deltaPtOverPt1VsPt.C`).

This is an **Open Data reproduction for validating tools**. It is not a CMS
result and makes no physics claim. Plots are labeled "CMS Open Data".

## Differences from the reference

| Aspect                      | Reference (`CMS-HIN-10-004`)                                                                      | This validation (2010 CMS Open Data)                                                                                                                                                          |
| --------------------------- | ------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Primary jet dataset         | Full 2010 `HICorePhysics` (`PromptReco`, 6.7 μb⁻¹, ~149k `HLT_HIJet50U`, 3,514 final dijets)      | Public CERN Open Data Record `14010` (`/HICorePhysics/HIRun2010-PromptReco-v3/RECO`, all 2,154 certified files across 66 runs; 77,987 `HLT_HIJet50U` ≈ 3.51 μb⁻¹ [52.3%], 1,879 final dijets) |
| Minimum-bias dataset        | `HLT_HIMinBiasHfOrBSC` (`PromptReco`, ≥3 HF towers/side + BSC halo + pixel compatibility)         | Public CERN Open Data Record `14011` (`/HIAllPhysics/HIRun2010-ZS-v2/RECO`, 11 certified files across 10 runs, 14,255 cleaned MB events with `etHFtowerSum > 3 GeV`)                          |
| Offline ECAL spike cleaning | `HiEcalRecHitSpikeFilter` (Table I ECAL step passes 97.66%, corresponding to jet-scale threshold) | `HiEcalRecHitSpikeFilter` on `EcalRecHitsEB` at the `CMSSW_3_9_2_patch5` default `minEt = 20 GeV` (`94.31%` pass rate; documented in `heavy-ion-reco` vs `minEt = 50 GeV`)                    |
| Jet energy corrections      | `Spring10_HI_IC5Calo` L2Relative + L3Absolute                                                     | `Spring10_HI_IC5Calo` L2Relative + L3Absolute on UE-subtracted `iterativeConePu5CaloJets`                                                                                                     |
| Simulation baselines        | Unquenched `PYTHIA Z2` and `PYTHIA + HYDJET` embedded MC                                          | Data-only Open Data validation (published text/macro benchmarks cited where tabulated)                                                                                                        |

## Results (2,154 `HICorePhysics` files + 11 `HIAllPhysics` files)

| Quantity                                                 | Round 1 (before `HiEcalRecHitSpikeFilter`)          | Round 2 (with `heavy-ion-reco` skill)                                  | Published `CMS-HIN-10-004`         |
| -------------------------------------------------------- | --------------------------------------------------- | ---------------------------------------------------------------------- | ---------------------------------- |
| Certified `HLT_HIJet50U` events                          | 77,987                                              | 77,987                                                                 | ~149,000                           |
| Offline noise + scraping survival                        | 99.96% (missing EB Swiss-cross filter)              | 94.28% (`HiEcalRecHitSpikeFilter` removes 3,324 spike events)          | 97.42% (ECAL × HCAL × scraping)    |
| Subleading-jet step (`pT,2 > 50` / `pT,1 > 120 GeV`)     | 82.22% (−5.62σ deficit from fake EB-edge mono-jets) | **87.16% ± 0.70%** (1,982 / 2,274; −0.26σ)                             | **87.38% ± 0.51%** (3,684 / 4,216) |
| Azimuthal dijet step (`Δφ_12 > 2π/3` / `pT,2 > 50 GeV`)  | 93.81% (−2.49σ)                                     | **94.80% ± 0.50%** (1,879 / 1,982; −0.96σ)                             | **95.39% ± 0.35%** (3,514 / 3,684) |
| HLT-scaled final dijet yield (`6.7 μb⁻¹` equiv.)         | 3,789                                               | **3,590 ± 83** (+0.92σ)                                                | **3,514**                          |
| `HLT_HIJet50U` plateau (`pT,1 in [100, 160] GeV`)        | 99.95%                                              | **99.98%** (5,829 / 5,830)                                             | **100.0%** (54 / 54)               |
| Peripheral `50-100%` → central `0-10%` `<A_J>`           | 0.181 → 0.240                                       | **0.184 ± 0.013 → 0.237 ± 0.005** (+3.90σ shift)                       | —                                  |
| Central `0-30%` `A_J` median                             | 0.226                                               | **0.223 ± 0.005(stat) ± 0.006(syst)** (+0.39σ)                         | **0.22**                           |
| Peripheral `50-100%` → central `0-10%` `R_B(A_J < 0.15)` | 0.361 → 0.206 (`TENSION` at −1.88σ)                 | **0.398 ± 0.052 → 0.225 ± 0.015** (3.23σ drop; −0.70σ / +0.01σ vs pub) | **~0.45 → ~0.225**                 |

`check_results.py` runs 22 objective stdlib-only physics and sanitization checks
on `summary.json`, `plots/`, `hepdata/`, and Markdown/JSON outputs. CI runs it.

## How it was made (two-agent producer + clean-context ARC review)

1. **Round 1 (before `heavy-ion-reco`):** `cms-analysis-coder` executed
   [PROMPT.md](PROMPT.md) over the 2,154 `HICorePhysics` EDM `RECO` files on
   EOS Open Data. Because `HcalNoiseSummary.nproblemRBXs_` uses pp occupancy
   thresholds and flags >99.8% of central PbPb events, the producer disabled
   `nproblemRBXs_` without adding `HiEcalRecHitSpikeFilter` on `EcalRecHitsEB`.
2. **Independent ARC diagnosis:** A clean-context `cms-results-validator` agent
   ran [CMS_COMPARISON_PROMPT.md](CMS_COMPARISON_PROMPT.md) against the primary
   `CMS-HIN-10-004` `.tex` tables and `.C` macros. It caught a −5.62σ Table I
   subleading-jet step deficit (`82.22%` vs `87.38%`) and `TENSION` in
   peripheral `R_B(A_J < 0.15)` (`0.361` vs `0.45`), proving from `(η_1, φ_1)`
   and `EMF_1` that `36.4%` (`196 / 539`) of failed-recoil events sat in the
   ECAL Barrel border strip `|η_1| in [1.44, 1.52]` with `EMF_1 ~ 1.0`
   (uncleaned barrel APD spikes), and also caught `hepdata_lib` 5-sig-fig
   rounding on 7-digit cutflow integers (`2,452,052 -> 2,452,100`).
3. **Round 2 & skill updates:** Implementing `HiEcalRecHitSpikeFilter` on
   `EcalRecHitsEB` reduced barrel-edge failed-recoil spikes from `36.36%` to
   `1.19%` (`0` mono-jets with `pT,2 <= 35 GeV`), upgraded `T7` (`R_B(A_J < 0.15)`)
   from `TENSION` to `CONSISTENT`, brought the high-`p_T` dijet cutflow steps to
   `−0.26σ` and `−0.96σ`, and motivated `plugins/cms/skills/heavy-ion-reco/` plus
   fixes in `hepdata`, `awkward-array`, `hist`, and `trigger-efficiency`.

## Comparison with published CMS results

`validate_hin10004.py` runs 12 quantitative tests against `CMS-HIN-10-004`
(`6 CONSISTENT`, `0 TENSION`, `2 INCONSISTENT`, `4 NOT COMPARABLE`), documented
in [COMPARISON.md](COMPARISON.md) and `comparison.json`:

| ID     | Test                                                     | CMS counterpart                                                | Result                                                                                                                                                               | Verdict        |
| ------ | -------------------------------------------------------- | -------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------- |
| T1     | Collision constants, HF boundaries & Glauber `<N_part>`  | Table II & `hf_centrality_distribution_mb_jet50_20101126_v1.C` | ΔN_part = 0; max HF edge diff 4.7e−7 TeV                                                                                                                             | CONSISTENT     |
| T2     | Event selection cutflow step efficiencies & scaled yield | Table I (`event_selection.tex`)                                | High-`pT` dijet steps −0.26σ, −0.96σ, scaled yield +0.92σ; intermediate cleaning steps differ by 20–39σ (`minEt = 20 GeV` whole-event veto & `etHFtowerSum > 3 GeV`) | INCONSISTENT   |
| T3     | `HLT_HIJet50U` trigger turn-on vs `pT,1` (Fig. 1a)       | `TrigEffCorrected151058and151025.C`                            | Common 12-bin χ²/ndf = 6.34/12 (p = 0.90); plateau 99.98% vs 100%                                                                                                    | CONSISTENT     |
| T4     | Forward `E_HF` shape in `HLT_HIJet50U` (Fig. 3a)         | `hf_centrality_distribution_mb_jet50_20101126_v1.C` (`hhfJET`) | χ²/ndf = 57.93/54 (p = 0.33, 0.97σ)                                                                                                                                  | CONSISTENT     |
| T5     | Forward `E_HF` shape in Minimum-Bias (Fig. 3a)           | `hf_centrality_distribution_mb_jet50_20101126_v1.C` (`hhfMB`)  | Full: χ²/ndf = 102.84/54 (3.98σ, ultra-peripheral bin 1 + `minEt = 20 GeV` central MB veto); bins ≥2: χ²/ndf = 64.20/53 (p = 0.14)                                   | INCONSISTENT   |
| T6     | `<(pT,1 - pT,2)/pT,1>` vs `pT,1` (Fig. 10)               | `deltaPtOverPt1VsPt.C`                                         | 12-point χ²/ndf = 1.80/12 (stat), 1.13/12 (stat+syst)                                                                                                                | CONSISTENT     |
| T7     | Balanced dijet fraction `R_B(A_J < 0.15)` (Fig. 9)       | Section III.A.3 text benchmarks (~0.45, ~0.225)                | `50-100%`: 0.398 (−0.70σ); `0-10%`: 0.225 (+0.01σ); ratio 1.77 ± 0.26                                                                                                | CONSISTENT     |
| T8     | Central `0-30%` `A_J` median                             | Section III.B text benchmark (0.22)                            | Event median 0.2230 ± 0.0051(stat) ± 0.0058(syst) (+0.39σ)                                                                                                           | CONSISTENT     |
| T9–T12 | 2D `(pT,1, pT,2)`, `pT,1`, `Δφ_12`, 5-bin `A_J`          | Figs. 4–8 (PDF figures only; no `.C` / HEPData)                | Digitization prohibited; physics trends verified in `check_results.py`                                                                                               | NOT COMPARABLE |

## Files

| File                                                 | What it is                                                                                                                     |
| ---------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------ |
| [PROMPT.md](PROMPT.md)                               | Task given to the `cms-analysis-coder` producer agent                                                                          |
| [CMS_COMPARISON_PROMPT.md](CMS_COMPARISON_PROMPT.md) | Task given to the clean-context `cms-results-validator` reviewer agent                                                         |
| [specification.md](specification.md)                 | Complete analysis specification (`analysis-spec-builder` + `heavy-ion-reco`)                                                   |
| `analyze_hin10004.py`                                | PEP 723 streaming EDM `RECO` analysis (`HiEcalRecHitSpikeFilter`, IC5Calo JEC, histograms, `summary.json`)                     |
| `plot_hin10004.py`                                   | Renders the 9 CMS-style `CMS Open Data` figures (`plots/*.{png,pdf}`)                                                          |
| `make_hepdata.py`, `hepdata/`                        | Validated 9-table HEPData record (`hepdata-validate` passes)                                                                   |
| `validate_hin10004.py`                               | Quantitative ARC comparison script producing `comparison.json` and `plots/comparison/`                                         |
| `data/`                                              | Certified 2010 HI Golden JSON (`Cert_150436-152957_...`, Record `14202`) and `Spring10_HI_IC5Calo` L2/L3 JEC tables            |
| `summary.json`, `comparison.json`                    | Machine-readable analysis outputs and 12-test ARC comparison record                                                            |
| [COMPARISON.md](COMPARISON.md)                       | Detailed ARC comparability matrix, diagnostic findings, and per-test verdicts                                                  |
| [SKILLS_USED.md](SKILLS_USED.md)                     | Skills exercised and bugs/gotchas fed back into `heavy-ion-reco`, `hepdata`, `awkward-array`, `hist`, and `trigger-efficiency` |
| `check_results.py`                                   | 22 stdlib-only physics, completeness, and path-sanitization checks (exit 0 = all pass)                                         |

## Reproduce

Per-file streaming caches (`<scratch>/hin10004_cache`, ~150 MB) and downloaded
reference sources (`<scratch>/hin10004_ref`, `arXiv:1102.1957` tarball) are kept
in a scratch directory outside the repository:

```bash
cd validation/hin-10-004
uv run analyze_hin10004.py --cache-dir <scratch>/hin10004_cache --workers 8 --summary summary.json
uv run plot_hin10004.py --summary summary.json --cache <scratch>/hin10004_cache/selected_events.npz --outdir plots
uv run make_hepdata.py --summary summary.json --plots-dir plots --outdir hepdata --no-tar
uv run validate_hin10004.py --base-dir . --ref-dir <scratch>/hin10004_ref --cache-dir <scratch>/hin10004_cache --outdir plots/comparison --json-out comparison.json
python3 check_results.py .   # SUMMARY=22/22
```
