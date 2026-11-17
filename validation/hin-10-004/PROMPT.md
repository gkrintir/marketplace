# Task for the cms-analysis-coder agent

Run headless from the repository root with Claude Code:

```bash
tmux new -s hin10004
sed -n '/^## Task/,$p' validation/hin-10-004/PROMPT.md | \
  claude -p --agent cms:cms-analysis-coder \
    --plugin-dir plugins/cms --plugin-dir plugins/hep-python-tools \
    --permission-mode acceptEdits --allowedTools "Bash Read Write Edit Glob Grep Skill TodoWrite" \
    > <scratch>/hin10004_agent_report.md 2>&1
```

Or run in Google Antigravity after `./scripts/install_antigravity.sh` by
invoking the `cms-analysis-coder` subagent with the `## Task` section below.

## Task

Replicate, using 2010 CMS Heavy-Ion Open Data (`√s_NN = 2.76 TeV` PbPb), the
detector-level dijet asymmetry and jet-quenching analysis of **CMS-HIN-10-004**
(_"Observation and studies of jet quenching in PbPb collisions at
nucleon-nucleon center-of-mass energy = 2.76 TeV"_, Phys. Rev. C 84 (2011)
024906, arXiv:1102.1957). Read `validation/hin-10-004/README.md` first for the
Open Data dataset scope and differences from the published reference.

Work in `validation/hin-10-004/`, keeping per-file streaming caches in a
scratch directory outside the repository (`--cache-dir <scratch>/hin10004_cache`).
First write `specification.md` using `analysis-spec-builder` and
`heavy-ion-reco`, then implement the end-to-end streaming pipeline.

**Inputs (CERN Open Data, `https://opendata.cern.ch/eos/opendata/...`, no proxy required):**

- Primary jet-triggered PbPb data: CERN Open Data Record `14010`
  (`/HICorePhysics/HIRun2010-PromptReco-v3/RECO`, 2,154 certified EDM `RECO` files).
- Minimum-bias PbPb reference data: CERN Open Data Record `14011`
  (`/HIAllPhysics/HIRun2010-ZS-v2/RECO`, sample >=10 certified files across the
  run range for minimum-bias `E_HF` centrality calibration and low-`p_T`
  trigger turn-on).
- Luminosity certification: 2010 HI Golden JSON
  `data/Cert_150436-152957_HI7TeV_StreamExpress_Collisions10_JSON_v2.txt`
  (CERN Open Data Record `14202`).
- Jet energy corrections: `data/Spring10_HI_L2Relative_IC5Calo.txt` and
  `data/Spring10_HI_L3Absolute_IC5Calo.txt` applied to raw UE-subtracted
  `recoCaloJets_iterativeConePu5CaloJets__RECO`.

**Event, centrality, and dijet selection (follow `heavy-ion-reco`):**

1. Certified `(run, lumisection)` in Golden JSON; decode exact per-run `HLT`
   indices from `edmTriggerResults_TriggerResults__HLT` via the in-file
   `MetaData` `ParameterSet` registry (`HLT_HIJet50U` signal, `HLT_HIJet35U`
   and `HLT_HIMinBiasHfOrBSC` references).
2. Collision selection: HF forward coincidence (`etHFtowerSumPlus >= 3 GeV` and
   `etHFtowerSumMinus >= 3 GeV`) and reconstructed `hiSelectedVertex`
   (`|z_vtx| <= 25 cm`, `validity`, `ndof > 0`).
3. Heavy-ion background & spike rejection: apply `HiEcalRecHitSpikeFilter` on
   `EcalRecHitsEB` (Swiss-cross `1 - E4/E1 > 0.95`, `|time| <= 4 ns`,
   `avoidIeta85 = False`), `HcalNoiseSummary.maxzeros_ < 100` (do **not** cut
   on pp-tuned `nproblemRBXs_ == 0`), pixel-vs-HF monster compatibility, and
   leading-jet `emEnergyFraction >= 0.01`.
4. Centrality: classify `E_HF = etHFhitSumPlus + etHFhitSumMinus` into
   `50-100%`, `30-50%`, `20-30%`, `10-20%`, and `0-10%` using the official
   2010 boundaries (`[0, 13.0979, 35.3927, 53.4790, 79.3700, 160.0]` TeV) and
   Table II `<N_part>` values (`22, 108, 187, 261, 355`).
5. Dijet selection: sort JEC-corrected `iterativeConePu5CaloJets` with
   `|η| < 2.0` by `p_T`; require `p_{T,1} > 120 GeV`, `p_{T,2} > 50 GeV`, and
   `Δφ_{12} > 2π/3`, and compute `A_J = (p_{T,1} - p_{T,2}) / (p_{T,1} + p_{T,2})`
   plus inclusive `R_B(A_J < 0.15)` normalized to all `p_{T,1} > 120 GeV` events.

**Plots (`CMS Open Data` style, PNG + PDF, in `validation/hin-10-004/plots/`):**

1. `fig1a_trigger_turnon` (`PRC Fig. 2a`): `HLT_HIJet50U` turn-on efficiency vs corrected `p_{T,1}`
2. `fig3a_hf_centrality` (`PRC Fig. 3a`): total `E_HF` distribution in MinBias vs `HLT_HIJet50U` with centrality boundaries
3. `fig4_pt2_vs_pt1` (`PRC Fig. 6`): 2D `(p_{T,1}, p_{T,2})` distribution in peripheral (`30-100%`) vs central (`0-30%`)
4. `fig5_leading_jet_pt` (`PRC Fig. 7`): leading-jet `p_{T,1}` spectrum across centrality
5. `fig6_dijet_dphi` (`PRC Fig. 8`): normalized `(1/N) dN/dΔφ_{12}` in all 5 centrality bins
6. `fig7_rb_dphi_vs_npart` (`PRC Fig. 9`): azimuthal back-to-back fraction `R_B(Δφ_{12} > 3.026)` vs `<N_part>`
7. `fig8_dijet_asymmetry_aj` (`PRC Fig. 10`): normalized `A_J` distributions across all 5 centrality bins
8. `fig9_rb_aj_vs_npart` (`PRC Fig. 11`): balanced dijet fraction `R_B(A_J < 0.15)` and `<A_J>` vs `<N_part>`
9. `fig10_fractional_imbalance_vs_pt1` (`PRC Fig. 12`): `<(p_{T,1} - p_{T,2}) / p_{T,1}>` vs `p_{T,1}` (`30-100%`, `10-30%`, `0-10%`) in the published `PRC Fig. 12` binning `[120, 140, 160, 180, 240] GeV` (also saving `[120, 130, 140, 170, 240] GeV` in `summary.json` for cross-checking the archived Dec 2010 `CmsHI/CVS_CmsHi` `deltaPtOverPt1VsPt.C` macro)

**Deliverables in `validation/hin-10-004/`:**

- `specification.md`
- `analyze_hin10004.py`, `plot_hin10004.py`, `make_hepdata.py`, `check_results.py`
- `plots/*.{png,pdf}`
- `summary.json` (no personal paths)
- `hepdata/` validated with `hepdata-validate`
- `SKILLS_USED.md`
