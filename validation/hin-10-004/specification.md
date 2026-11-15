# Replication of CMS-HIN-10-004: Dijet Quenching in PbPb Collisions at 2.76 TeV

Data-driven CMS Open Data replication of the calorimeter dijet momentum imbalance, azimuthal correlation, and balanced-jet fraction measurements in PbPb collisions at `$\sqrt{s_{\mathrm{NN}}} = 2.76\text{ TeV}$` (CMS-HIN-10-004, Phys. Rev. C 84, 024906 (2011), arXiv:1102.1957).

## Data Samples

| Dataset (DAS name or Open Data path)          | Tier   | Era / campaign                                          | Role(s) in analysis                                                                                  |
| --------------------------------------------- | ------ | ------------------------------------------------------- | ---------------------------------------------------------------------------------------------------- |
| `/HICorePhysics/HIRun2010-PromptReco-v3/RECO` | `RECO` | `HIRun2010-PromptReco-v3` (CERN Open Data Record 14010) | Primary jet-triggered PbPb collision data (5,436 files total; 4,640 files in validated runs)         |
| `/HIAllPhysics/HIRun2010-ZS-v2/RECO`          | `RECO` | `HIRun2010-ZS-v2` (CERN Open Data Record 14011)         | Minimum-bias PbPb collision reference for HF centrality distribution and trigger turn-on cross-check |

- Luminosity / certification: `Cert_150436-152957_HI7TeV_StreamExpress_Collisions10_JSON_v2.txt` (CERN Open Data Record 14202), corresponding to `$L_{\mathrm{int}} = 6.7\ \mu\mathrm{b}^{-1}$` across the full 2010 PbPb dataset (runs `150436` to `152957`, global tag `GR_R_39X_V6B` in `CMSSW_3_9_2_patch5`). Because streaming all 4,640 `RECO` files over HTTPS is I/O-bound, the workflow processes a stratified multi-run subset of **2,154 Golden-JSON validated `HICorePhysics` files across 66 early, mid, and late 2010 runs** (`52.34%` of the full `HLT_HIJet50U` dataset, `$\approx 3.51\ \mu\mathrm{b}^{-1}$`, scaling up from an initial 2-file yardstick check) alongside `11` Golden-JSON `HIAllPhysics` Minimum-Bias files, and reports both raw yields and per-event normalized observables.
- Triggers:
  - Signal dijet analysis (`/HICorePhysics/HIRun2010-PromptReco-v3/RECO`): inclusive single-jet trigger `HLT_HIJet50U` (requiring L1 uncorrected jet `$p_{T} > 30\text{ GeV}$` and HLT uncorrected jet `$p_{T} > 50\text{ GeV}$`, identified per event via active filter tag `hltHI1jet50U::HLT` with non-empty `filterKeys_` in `triggerTriggerEvent_hltTriggerSummaryAOD__HLT.obj`). Fully efficient for corrected leading jet `$p_{T,1} > 100\text{ GeV}$`.
  - Minimum-bias centrality reference (`/HIAllPhysics/HIRun2010-ZS-v2/RECO` or minimum-bias triggers in `HICorePhysics`): `HLT_HIMinBiasHfOrBSC` (active L1/HLT filter tag `hltL1sHIMinBiasHfOrBSC::HLT` in `triggerTriggerEvent_hltTriggerSummaryAOD__HLT.obj`).
- MC normalisation: Data-only Open Data measurement (shape-normalized per event or per dijet pair, matching CMS-HIN-10-004 Section III.A). The private 2010 `PYTHIA+DATA` data-mixed embedded simulation used in CMS-HIN-10-004 is not publicly released on CERN Open Data; comparisons are performed directly against the published CMS-HIN-10-004 PbPb data figures, cutflow efficiencies, and official `Analysis2010` reference macros.

## Histograms

All centrality-resolved histograms use the 5 CMS-HIN-10-004 collision centrality categories: `50-100%` (`$\langle N_{\mathrm{part}} \rangle = 22 \pm 2$`), `30-50%` (`$108 \pm 5$`), `20-30%` (`$187 \pm 5$`), `10-20%` (`$261 \pm 4$`), and `0-10%` (`$355 \pm 3$`), or the 3 merged categories `30-100%`, `10-30%`, and `0-10%` where specified.

1. Total HF Energy Centrality Distribution (Figure 3a equivalent)

   - x Axis
     - Total forward hadronic calorimeter hit energy sum `$E_{\mathrm{HF,tot}} = (\text{etHFhitSumPlus} + \text{etHFhitSumMinus}) / 1000$` in cleaned minimum-bias events vs `HLT_HIJet50U` events
     - `$\sum E_{\mathrm{HF}}\text{ [TeV]}$`
     - 70 uniform bins from `0` to `160` TeV
     - Entries per event: one per cleaned event
   - Other axes: trigger category (`MinBias`, `HLT_HIJet50U`)
   - Weights: unweighted (normalized to unit integral per trigger category)

2. Single-Jet Trigger Turn-On Curve (Figure 1a equivalent)

   - x Axis
     - Corrected leading calorimeter jet transverse momentum `$p_{T,1}$` (`$|\eta_{1}| < 2.0$`) in minimum-bias cleaned events (denominator) and events additionally passing `HLT_HIJet50U` (numerator)
     - `$\text{Corrected leading jet } p_{T,1}\text{ [GeV]}$`
     - 24 uniform bins from `20` to `260` GeV
     - Entries per event: one per cleaned event with at least one jet in `$|\eta| < 2.0$`
   - Weights: unweighted (efficiency with Clopper-Pearson 68% confidence intervals)

3. Subleading vs Leading Jet Transverse Momentum 2D Correlation (Figure 4 equivalent)

   - x Axis: Leading jet `$p_{T,1}\text{ [GeV]}$`, 26 uniform bins from `120` to `250` GeV
   - y Axis: Subleading jet `$p_{T,2}\text{ [GeV]}$`, 40 uniform bins from `50` to `250` GeV
   - Entries per event: one per selected dijet event (`$p_{T,1} > 120\text{ GeV}$`, `$p_{T,2} > 50\text{ GeV}$`, `$|\eta_{1,2}| < 2.0$`, `$\Delta\phi_{12} > 2\pi/3$`)
   - Other axes: 3 centrality categories (`30-100%`, `10-30%`, `0-10%`)
   - Weights: unweighted

4. Leading Jet Transverse Momentum Spectrum (Figure 5 equivalent)

   - x Axis
     - Corrected leading jet transverse momentum `$p_{T,1}$` in selected dijet events (`$p_{T,1} > 120\text{ GeV}$`, `$p_{T,2} > 50\text{ GeV}$`, `$\Delta\phi_{12} > 2\pi/3$`)
     - `$\text{Leading jet } p_{T,1}\text{ [GeV]}$`
     - 9 uniform bins from `120` to `300` GeV (`20` GeV bin width)
     - Entries per event: one per selected dijet event
   - Other axes: 5 centrality categories (`50-100%`, `30-50%`, `20-30%`, `10-20%`, `0-10%`)
   - Weights: unweighted (normalized to event fraction `$(1/N_{\mathrm{dijet}})\, \Delta N / \Delta p_{T,1}$` or unit sum per centrality bin)

5. Dijet Azimuthal Correlation `$\Delta\phi_{12}$` (Figure 6 equivalent)

   - x Axis
     - Azimuthal angle difference `$\Delta\phi_{12} = |\phi_{1} - \phi_{2}|$` wrapped into `$[0, \pi]$` for events with `$p_{T,1} > 120\text{ GeV}$` and `$p_{T,2} > 50\text{ GeV}$` (`$|\eta_{1,2}| < 2.0$`), **without** applying the `$\Delta\phi_{12} > 2\pi/3$` cut
     - `$\Delta\phi_{12}\text{ [rad]}$`
     - 20 uniform bins from `0` to `$\pi$`
     - Entries per event: one per event passing `$p_{T,1} > 120\text{ GeV}$` and `$p_{T,2} > 50\text{ GeV}$`
   - Other axes: 5 centrality categories (`50-100%`, `30-50%`, `20-30%`, `10-20%`, `0-10%`)
   - Weights: unweighted (normalized to `$(1/N)\, dN/d\Delta\phi_{12}$`)

6. Azimuthal Back-to-Back Fraction `$R_{B}(\Delta\phi_{12} > 3.026)$` vs `$N_{\mathrm{part}}$` (Figure 7 equivalent)

   - x Axis
     - Mean number of participating nucleons `$\langle N_{\mathrm{part}} \rangle$` for each of the 5 centrality classes (`22`, `108`, `187`, `261`, `355`)
     - `$N_{\mathrm{part}}$`
     - 5 centrality categories
     - Entries per event: fraction of events with `$p_{T,1} > 120\text{ GeV}$` and `$p_{T,2} > 50\text{ GeV}$` (`$|\eta_{1,2}| < 2.0$`, no `$\Delta\phi_{12} > 2\pi/3$` pre-cut) that satisfy `$\Delta\phi_{12} > 3.026\text{ rad}$`
   - Weights: unweighted (binomial / Clopper-Pearson statistical uncertainties)

7. Dijet Momentum Asymmetry Ratio `$A_{J}$` (Figure 8 equivalent, primary jet-quenching observable)

   - x Axis
     - Dijet transverse momentum asymmetry `$A_{J} = (p_{T,1} - p_{T,2}) / (p_{T,1} + p_{T,2})$` for selected dijet events (`$p_{T,1} > 120\text{ GeV}$`, `$p_{T,2} > 50\text{ GeV}$`, `$|\eta_{1,2}| < 2.0$`, `$\Delta\phi_{12} > 2\pi/3$`)
     - `$A_{J} = (p_{T,1} - p_{T,2}) / (p_{T,1} + p_{T,2})$`
     - 10 uniform bins from `0.0` to `0.7` (bin width `0.07`, matching the kinematic range of CMS-HIN-10-004 Figure 8) plus fine/full `[0.0, 0.9]` summary check (`13` bins of width `~0.069` or `10` bins `[0, 0.7]`; note CMS-HIN-10-004 Figure 8 uses 13 uniform bins across `[0, 0.91]` with bin width `0.07`, where `$A_{J} < 0.7$` is populated)
     - Entries per event: one per selected dijet event
   - Other axes: 5 centrality categories (`50-100%`, `30-50%`, `20-30%`, `10-20%`, `0-10%`)
   - Weights: unweighted (normalized to event fraction `$\Delta N / N_{\mathrm{dijet}}$` summing to unity in each centrality bin)

8. Balanced Dijet Fraction `$R_{B}(A_{J} < 0.15)$` vs `$N_{\mathrm{part}}$` (Figure 9 equivalent)

   - x Axis
     - Mean number of participating nucleons `$\langle N_{\mathrm{part}} \rangle$` (`22`, `108`, `187`, `261`, `355`)
     - `$N_{\mathrm{part}}$`
     - 5 centrality categories
     - Entries per event: denominator is **all** cleaned `HLT_HIJet50U` events having a leading jet `$p_{T,1} > 120\text{ GeV}$` within `$|\eta_{1}| < 2.0$` (explicitly including mono-jet events with no partner above `50` GeV or with `$\Delta\phi_{12} \le 2\pi/3$`); numerator requires a subleading jet (`$|\eta_{2}| < 2.0$`) with `$p_{T,2} > 50\text{ GeV}$`, `$\Delta\phi_{12} > 2\pi/3$`, and `$A_{J} < 0.15$`
   - Weights: unweighted (Clopper-Pearson 68% statistical intervals)

9. Mean Fractional Imbalance `$\langle (p_{T,1} - p_{T,2}) / p_{T,1} \rangle$` vs `$p_{T,1}$` (Figure 10 equivalent)
   - x Axis
     - Leading jet `$p_{T,1}\text{ [GeV]}$` in variable bins `[120, 130, 140, 170, 240]` GeV (bin centers `125, 135, 155, 205` GeV, matching `deltaPtOverPt1VsPt.C` for `$p_{T,1} > 120\text{ GeV}$`)
     - `$\text{Leading jet } p_{T,1}\text{ [GeV]}$`
     - Profile mean and standard error of the mean of `$(p_{T,1} - p_{T,2}) / p_{T,1}$` per `$p_{T,1}$` bin for selected dijet events (`$p_{T,1} > 120\text{ GeV}$`, `$p_{T,2} > 50\text{ GeV}$`, `$|\eta_{1,2}| < 2.0$`, `$\Delta\phi_{12} > 2\pi/3$`)
   - Other axes: 3 centrality categories (`30-100%`, `10-30%`, `0-10%`)
   - Weights: unweighted

## Analysis Steps

1. Luminosity Certification and Trigger Selection

   - Decode `(run, lumi, event)` per entry from `EventAuxiliary` (big-endian `uint32` triplet at byte offset 43 in `CMSSW_3_9_2_patch5` EDM `RECO` records) and retain only luminosity sections listed in `Cert_150436-152957_HI7TeV_StreamExpress_Collisions10_JSON_v2.txt`.
   - Require `HLT_HIJet50U` for the jet sample (`hltHI1jet50U::HLT` in `triggerTriggerEvent_hltTriggerSummaryAOD__HLT.obj.triggerFilters_.filterTag_` with `len(filterKeys_) > 0`).

2. Offline Heavy-Ion Collision Event Cleaning (CMS-HIN-10-004 Table 1)

   - Offline HF coincidence (UPC and beam-gas rejection): require both forward calorimeters to record above-threshold tower energy, `etHFtowerSumPlus > 3.0` GeV and `etHFtowerSumMinus > 3.0` GeV in `recoCentrality_hiCentrality__RECO.obj` (verified event-by-event on `HLT_HIJet50U` data to match the exact `CaloTowers` requirement of at least 3 towers with `$E > 3\text{ GeV}$` on each side).
   - Primary vertex requirement: require at least one valid reconstructed heavy-ion primary vertex (`recoVertexs_hiSelectedVertex__RECO.obj`) with `validity == True`, `ndof > 0` (`$\ge 2$` tracks), and `$|z_{\mathrm{vtx}}| \le 25\text{ cm}$`.
   - ECAL barrel APD spike rejection (`HiEcalRecHitSpikeFilter` from 2010 Heavy-Ion `CmsHi/PhotonAnalysis` & `EventSelection_cff.py`): inspect `EcalRecHitsSorted_ecalRecHit_EcalRecHitsEB_RECO.obj` (`EcalRecHitsEB`) and reject any event containing at least one barrel crystal with `$E_{T,\mathrm{rh}} = E / \cosh\eta > 20.0\text{ GeV}$` that either has out-of-time arrival `$|t| > 4.0\text{ ns}$` (`timeThreshold = 4.0`) or fails the topological Swiss-cross criterion `$1 - E_{4}/E_{1} > 0.95$` (`swissThreshold = 0.95`, summing raw 4-neighbor `ieta`/`iphi` crystal energies across `$z = 0$` and wrapping `$\mathrm{iphi}\in[1,360]$` with **`avoidIeta85 = False`** so barrel-endcap border `$|\mathrm{ieta}| = 85$` (`$|\eta| \approx 1.475$`) spikes unflagged in `CMSSW_3_9_2_patch5` `PromptReco-v3` are rejected).
   - Beam-scraping / pixel-HF compatibility and HCAL noise rejection: require `HcalNoiseSummary_hcalnoise__RECO.obj.maxzeros_ < 100` (matching 2010 Heavy-Ion `EventSelection_cff.py`; note `nproblemRBXs_` in `CMSSW_3_9_2_patch5` uses `pp` occupancy thresholds that flag all 36 HB RBXs in central PbPb, so only `0 <= nproblemRBXs_ <= 72` sanity bounds apply), reject extreme beam-scraping outliers where `pixelMultiplicity < 50` while `etHFhitSum > 5000` GeV, and reject HCAL spike/noise artifacts where a barrel/endcap leading jet (`$|\eta| < 2.0$`, `$p_{T} > 50\text{ GeV}$`) has `mEnergyFractionEm < 0.01`.

3. Centrality Classification

   - Compute total HF hit sum `$E_{\mathrm{HF,hit}} = \text{etHFhitSumPlus} + \text{etHFhitSumMinus}$` (in GeV) from `recoCentrality_hiCentrality__RECO.obj`.
   - Assign each cleaned collision event to one of the 5 official CMS 2010 (`GR_R_39X_V6B` / `hf_centrality_distribution_mb_jet50_20101126_v1.C`) centrality classes:
     - `50-100%`: `$E_{\mathrm{HF,hit}} < 13097.866\text{ GeV}$` (`$< 13.098\text{ TeV}$`, `$\langle N_{\mathrm{part}} \rangle = 22$`)
     - `30-50%`: `$13097.866 \le E_{\mathrm{HF,hit}} < 35392.680\text{ GeV}$` (`$13.098\text{--}35.393\text{ TeV}$`, `$\langle N_{\mathrm{part}} \rangle = 108$`)
     - `20-30%`: `$35392.680 \le E_{\mathrm{HF,hit}} < 53479.020\text{ GeV}$` (`$35.393\text{--}53.479\text{ TeV}$`, `$\langle N_{\mathrm{part}} \rangle = 187$`)
     - `10-20%`: `$53479.020 \le E_{\mathrm{HF,hit}} < 79370.000\text{ GeV}$` (`$53.479\text{--}79.370\text{ TeV}$`, `$\langle N_{\mathrm{part}} \rangle = 261$`)
     - `0-10%`: `$E_{\mathrm{HF,hit}} \ge 79370.000\text{ GeV}$` (`$\ge 79.370\text{ TeV}$`, `$\langle N_{\mathrm{part}} \rangle = 355$`)

4. Jet Energy Calibration and Object Selection

   - Read underlying-event-subtracted iterative cone `$R = 0.5$` calorimeter jets from `recoCaloJets_iterativeConePu5CaloJets__RECO.obj`.
   - Apply the official CMS `Spring10_HI` `IC5Calo` factorized jet energy corrections (`data/Spring10_HI_L2Relative_IC5Calo.txt` as a function of `$\eta$` and raw UE-subtracted `$p_{T}$`, followed by `data/Spring10_HI_L3Absolute_IC5Calo.txt` as a function of `$p_{T,\mathrm{L2}}$`) to obtain `pt_corr` and `mass_corr`.
   - Select calorimeter jets within `$|\eta| < 2.0$` and sort them in descending order of corrected `$p_{T}$` per event.

5. Dijet Combinatorics and Observable Calculation
   - Define `jet1` as the highest-`$p_{T}$` (`pt_corr`) jet within `$|\eta| < 2.0$` and `jet2` as the second-highest-`$p_{T}$` jet within `$|\eta| < 2.0$` (explicitly making no veto or requirement on a third jet, matching CMS-HIN-10-004 Section II.B).
   - Leading-jet inclusive sample (denominator of `$R_{B}(A_{J} < 0.15)$`): events with `$p_{T,1} > 120\text{ GeV}$`.
   - Azimuthal correlation sample (Figures 6 and 7): events with `$p_{T,1} > 120\text{ GeV}$` and `$p_{T,2} > 50\text{ GeV}$` (no `$\Delta\phi_{12}$` cut), computing `$\Delta\phi_{12} = |\text{jet1.deltaphi(jet2)}| \in [0, \pi]$` via `vector` `Momentum4D`.
   - Standard dijet sample (Figures 4, 5, 8, and 10): events with `$p_{T,1} > 120\text{ GeV}$`, `$p_{T,2} > 50\text{ GeV}$`, and `$\Delta\phi_{12} > 2\pi/3$`. Compute:
     - `$A_{J} = (p_{T,1} - p_{T,2}) / (p_{T,1} + p_{T,2})$`
     - `$\Delta p_{T,\mathrm{rel}} = (p_{T,1} - p_{T,2}) / p_{T,1}$`

## Corrections and Systematics

- Jet Energy Corrections (JEC): Factorized `L2Relative` (`$\eta$` uniformity) and `L3Absolute` (`$p_{T}$` response) corrections from `Spring10_HI_L2Relative_IC5Calo.txt` and `Spring10_HI_L3Absolute_IC5Calo.txt` (corresponding to `ic5CaloL2L3` in `CMSSW_3_9_2_patch5` / `GR_R_39X_V6B`, without data/MC residual L2L3 corrections, matching `PatAna_cff.py` and CMS-HIN-10-004 Section II.D).
- Systematic cross-check variations recorded in `summary.json`:
  - JEC scale shift `$\pm 5\%$` (verifying exact cancellation in `$A_{J}$` for constant shifts apart from the `$p_{T,1} > 120\text{ GeV}$` and `$p_{T,2} > 50\text{ GeV}$` threshold migration).
  - Centrality-dependent residual JEC tilt (`$-5\%$` in `50-100%` to `$+5\%$` at low `$p_{T}$` in `0-10%`, motivated by CMS-HIN-10-004 Section II.D.3 Figure 2) to bracket systematic impact on `$R_{B}(A_{J} < 0.15)$` and `$\langle \Delta p_{T,\mathrm{rel}} \rangle$`.

## Workflow

1. Data Extraction and Skimming
   - Stream validated 2010 `RECO` files over HTTPS Range requests (`https://opendata.cern.ch//eos/opendata/cms/hidata/HIRun2010/...`, rewritten from `root://eospublic.cern.ch/`).
   - Read only lightweight primitive leaf branches (avoiding the 23 MB/file `CaloTowers` branch):
     - `EventAuxiliary` (raw basket bytes for `(run, lumi, event)`)
     - `triggerTriggerEvent_hltTriggerSummaryAOD__HLT./triggerTriggerEvent_hltTriggerSummaryAOD__HLT.obj`
     - `recoCaloJets_iterativeConePu5CaloJets__RECO./recoCaloJets_iterativeConePu5CaloJets__RECO.obj/recoCaloJets_iterativeConePu5CaloJets__RECO.obj.{pt_,eta_,phi_,mass_,m_specific.mEnergyFractionEm}`
     - `recoCentrality_hiCentrality__RECO./recoCentrality_hiCentrality__RECO.obj/recoCentrality_hiCentrality__RECO.obj.{etHFhitSumPlus_,etHFhitSumMinus_,etHFtowerSumPlus_,etHFtowerSumMinus_,pixelMultiplicity_,trackMultiplicity_}`
     - `recoVertexs_hiSelectedVertex__RECO./recoVertexs_hiSelectedVertex__RECO.obj/recoVertexs_hiSelectedVertex__RECO.obj.{position_.fCoordinates.fZ,ndof_,validity_}`
     - `HcalNoiseSummary_hcalnoise__RECO./HcalNoiseSummary_hcalnoise__RECO.obj/HcalNoiseSummary_hcalnoise__RECO.obj.{maxzeros_,nproblemRBXs_}`
     - `recoCaloClusters_uncleanedHybridSuperClusters_hybridBarrelBasicClusters_RECO` (`energy_` and `seedId_.id_`) as an exact algebraic lossless pre-filter (`$E_{\mathrm{BC}} / \cosh\eta_{\mathrm{seed}} > 19.5\text{ GeV}$`) for candidate events, followed by targeted `fBasketSeek`/`fBasketBytes` HTTP Range reads and binary decoding of `EcalRecHitsSorted_ecalRecHit_EcalRecHitsEB_RECO.obj` to evaluate `HiEcalRecHitSpikeFilter` bit-for-bit without downloading unneeded 850 KB/event `EcalRecHitsEB` baskets.
   - Cache compact per-event selected skim arrays to disk (`cache/selected_events.parquet` or `.npz`) so plotting and HEPData packaging run instantaneously and reproducibly without re-streaming `RECO` files.
2. Analysis and Histogramming
   - Construct `Momentum4D` jet collections with `vector.register_awkward()`, apply `Spring10_HI` `IC5Calo` L2L3 JEC, sort by `pt`, apply event and dijet cuts, fill `hist.Hist` objects with `.Weight()` storage, render CMS Open Data figures (`mplhep.style.CMS`, `hep.cms.label("Open Data", ...)`), write `summary.json`, and build validated HEPData YAML tables in `hepdata/`.

## Statistical Analysis

None (model-independent differential cross-section ratios, event-normalized kinematic distributions, and Clopper-Pearson binomial fractions with statistical and JEC-variation uncertainties; compared quantitatively against CMS-HIN-10-004 via `check_results.py` and `cms-results-validator`).

## Tools

- `uv` + PEP 723 inline script metadata (`standalone-script` skill): reproducible isolated execution without manual virtualenvs.
- `uproot` + `aiohttp` / `requests`: HTTP byte-range streaming of CMS Run-1 EDM `RECO` files from `https://opendata.cern.ch/`.
- `awkward` (`awkward-array` skill): jagged calorimeter jet masking, `pt` sorting (`ak.argsort`), leading/subleading extraction (`ak.pad_none`), and per-event cutflow tracking.
- `vector` (`vector-awkward` skill): `Momentum4D` four-vector construction and wrapped azimuthal separation `deltaphi`.
- `hist` + `mplhep` + `matplotlib` (`hist` skill): N-dimensional categorical + regular/variable histograms and CMS-style publication figures (`PNG` and `PDF`, labeled `"CMS Open Data"`).
- `hepdata_lib` + `hepdata-validator` (`hepdata` skill): packaging Open Data distributions and summary fractions into validated HEPData YAML tables (`VALID=1`).
