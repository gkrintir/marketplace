# 2010 CMS Heavy-Ion PbPb RECO physics reference

Read this reference when analysing 2010 CMS $\sqrt{s_{\mathrm{NN}}} = 2.76\ \mathrm{TeV}$
PbPb Open Data (`HIRun2010`, `CMSSW_3_9_2_patch5` EDM `RECO`): looking up Open
Data record IDs, HF centrality boundaries and Glauber $\langle N_{\mathrm{part}} \rangle$,
offline Heavy-Ion collision and noise filters, ECAL Barrel APD spike thresholds,
or `Spring10_HI` `IC5Calo` jet energy corrections. Tested on CERN Open Data
records 14010 and 14011 (Python 3.11, uproot 5.7.7, awkward 2.15.0), 2026-11-14.

## Contents

1. [2010 PbPb Open Data records](#2010-pbpb-open-data-records)
2. [Centrality and Glauber N_part](#centrality-and-glauber-n_part)
3. [Offline Heavy-Ion event cleaning](#offline-heavy-ion-event-cleaning)
4. [ECAL Barrel APD spike cleaning: photon vs jet minEt](#ecal-barrel-apd-spike-cleaning-photon-vs-jet-minet)
5. [Iterative-cone UE-subtracted jets and Spring10_HI JEC](#iterative-cone-ue-subtracted-jets-and-spring10_hi-jec)
6. [Inclusive vs dijet-conditional quench observables](#inclusive-vs-dijet-conditional-quench-observables)

## 2010 PbPb Open Data records

| Dataset / artefact                                                 | CERN Open Data record                                                                     | Split level & size                                                         | Purpose                                                                                                                   |
| ------------------------------------------------------------------ | ----------------------------------------------------------------------------------------- | -------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------- |
| `/HICorePhysics/HIRun2010-PromptReco-v3/RECO`                      | [14010](https://opendata.cern.ch/record/14010)                                            | `splitLevel=99` (2,276 files total; 2,154 certified across 66 runs)        | Unprescaled high-$p_T$ jet (`HLT_HIJet50U`), photon, dimuon, and core Minimum-Bias (`HLT_HIMinBiasHfOrBSC_Core`) triggers |
| `/HIAllPhysics/HIRun2010-ZS-v2/RECO`                               | [14011](https://opendata.cern.ch/record/14011)                                            | `splitLevel=0` (12,243 files)                                              | Bulk Minimum-Bias (`HLT_HIMinBiasHfOrBSC`) and lower-threshold jet (`HLT_HIJet35U`) triggers                              |
| `Cert_150436-152957_HI7TeV_StreamExpress_Collisions10_JSON_v2.txt` | [14202](https://opendata.cern.ch/record/14202)                                            | Golden JSON (runs `150436`–`152957`)                                       | Official validated PbPb run/lumisection mask ($\mathcal{L}_{\mathrm{int}} \approx 6.7\ \mu\mathrm{b}^{-1}$)               |
| `HeavyIonDijetAsymmetry2010` / `GR_R_39X_V6B` `IC5Calo` JEC        | [470](https://opendata.cern.ch/record/470) / [1811](https://opendata.cern.ch/record/1811) | Factorized ASCII JEC tables (`Spring10_HI_IC5Calo`, `Spring10_V8_IC5Calo`) | `L2Relative` ($\eta$) and `L3Absolute` ($p_T$) corrections for `iterativeConePu5CaloJets`                                 |

Convert `root://eospublic.cern.ch//eos/opendata/...` file-index URLs to
`https://opendata.cern.ch/eos/opendata/...` when streaming from environments
without an XRootD client.

## Centrality and Glauber N_part

Event centrality in 2010 PbPb `RECO` is defined from the total electromagnetic
plus hadronic energy deposited in both Forward Hadron calorimeters
($2.9 < |\eta| < 5.2$), stored in `recoCentrality_hiCentrality__RECO`:

```python
hf_sum_gev = etHFhitSumPlus_ + etHFhitSumMinus_   # total HF hit energy in GeV (NOT transverse ET!)
```

Despite the C++ member name `etHFhitSum*`, `CentralityProducer` in
`CMSSW_3_9_2_patch5` stores **total HF hit energy $\sum E_{\mathrm{HF}}$ in
GeV** (peaking near $110\text{--}125\ \mathrm{TeV}$ in $0\text{--}10\%$ PbPb),
whereas `etHFtowerSumPlus_ + etHFtowerSumMinus_` stores transverse tower
$\sum E_{T,\mathrm{HF}}$ ($\le 4\ \mathrm{TeV}$) used for the 3 GeV tower
coincidence check.

Five-bin centrality boundaries (`hf_centrality_distribution_mb_jet50` in
CMS-HIN-10-004, `Phys. Rev. C 84 (2011) 024906`) and Monte Carlo Glauber
$\langle N_{\mathrm{part}} \rangle$:

| Centrality | $\sum E_{\mathrm{HF}}$ low [GeV] | $\sum E_{\mathrm{HF}}$ high [GeV] | $\langle N_{\mathrm{part}} \rangle$ |
| ---------- | -------------------------------- | --------------------------------- | ----------------------------------- |
| `50-100%`  | `0.0`                            | `13097.866`                       | $22 \pm 2$                          |
| `30-50%`   | `13097.866`                      | `35392.680`                       | $108 \pm 5$                         |
| `20-30%`   | `35392.680`                      | `53479.020`                       | $187 \pm 5$                         |
| `10-20%`   | `53479.020`                      | `79370.000`                       | $261 \pm 4$                         |
| `0-10%`    | `79370.000`                      | $+\infty$                         | $355 \pm 3$                         |

On cleaned `HLT_HIMinBiasHfOrBSC` events from `HIAllPhysics`, these fixed
boundaries give `0.506 / 0.197 / 0.097 / 0.098 / 0.101` (matching the target
`0.50 / 0.20 / 0.10 / 0.10 / 0.10` fractions within $0.6\%$). Note that
`CMS-HIN-10-004` Fig. 3a (`hhfMB` in `hf_centrality_distribution_mb_jet50_20101126_v1.C`,
commit `f6717ae4`) names `HLT_MinBiasHForBSC_Core` (`HLT_HIMinBiasHfOrBSC_Core`
in `HICorePhysics`): with `noBSChalo`, $\ge 3$ HF `CaloTowers`/side, and the
2-track vertex across 110 `HICorePhysics` files, `_Core` gives `23.06%` in the
lowest `[0, 2.286] TeV` bin vs `23.23%` in `hhfMB` ($G^2 = 54.70 / 53\ \mathrm{ndf}$),
whereas ordinary `HLT_HIMinBiasHfOrBSC` in `HIAllPhysics` differs in low-HF
trigger sampling and primary-vertex track attachment (`5` of `591` shared `_Core`
events have `2–11` attached vertex tracks in `HIAllPhysics` vs `0` in
`HICorePhysics`). Real Minimum-Bias data is $\sim 97\pm 3\%$ efficient in
low-$\sum E_{\mathrm{HF}}$ hadronic collisions ($\sum E_{\mathrm{HF}} < 2.75\ \mathrm{TeV}$).

## Offline Heavy-Ion event cleaning

Apply the 2010 Heavy-Ion collision and monster-event filters
(`HeavyIonsAnalysis/Configuration/python/EventSelection_cff.py`) in order:

1. **Golden JSON lumisection mask** on `(run, lumi)` from `EventAuxiliary`,
   followed by the **BSC beam-halo veto (`noBSChalo`)** rejecting L1 technical
   trigger bits `36, 37, 38, 39` (`99.60%` pass rate on `HLT_HIJet50U` vs
   `99.61%` in `CMS-HIN-10-004` Table I).
2. **HF 3-tower coincidence (`hfCoincFilter3`: `hfPosFilter3` + `hfNegFilter3`)**:
   at least **3 HF `CaloTowers` with $E \ge 3.0\ \mathrm{GeV}$ in
   $3.0 \le |\eta| \le 6.0$ on each side** (`75.19%` pass rate after
   `noBSChalo`, vs `74.98%` in Table I). Using the summary proxy
   `(etHFtowerSumPlus_ >= 3.0) & (etHFtowerSumMinus_ >= 3.0)` is faster when
   skipping `CaloTowersSorted`, and admits ~5% extra low-multiplicity jet
   triggers that subsequently fail the 2-track primary vertex cut (`0 / 148`
   pass `tracksSize >= 2`).
3. **Reconstructed Heavy-Ion primary vertex (`primaryVertexFilter` on
   `hiSelectedVertex`)**: at least one vertex with `!isFake` (`validity_ == True`),
   `tracksSize >= 2` (`ndof_ > 0`), `|fZ| <= 25.0` cm, and transverse radius
   `Rho <= 2.0` cm (`98.78%` pass rate after `hfCoincFilter3`, vs `98.97%` in
   Table I).
4. **Beam-halo / scraping pixel-HF compatibility**: reject events where pixel
   hit multiplicity exceeds the HF-scaled envelope, e.g.
   `pixelMultiplicity_ <= 150.0 + 7.5 * hf_sum_gev` (or cluster-shape vertex
   compatibility `hltPixelClusterShapeFilter` when reading `siPixelRecHits`).
5. **HCAL RBX noise (`HcalNoiseSummary_hcalnoise__RECO`)**: require
   `maxzeros_ < 100`. **Never** cut on `nproblemRBXs_ == 0` in Run-1 PbPb:
   `HcalNoiseSummary` in `CMSSW_3_9_2_patch5` evaluates proton-proton RBX
   multiplicity thresholds and flags all 36 HB readout boxes
   (`nproblemRBXs_ = 36..38`) in $>99.8\%$ of central PbPb events.
6. **ECAL Barrel APD spike rejection (`HiEcalRecHitSpikeFilter`)** plus leading
   jet `mEnergyFractionEm >= 0.01` (rejects residual pure-HCAL noise).

## ECAL Barrel APD spike cleaning: photon vs jet minEt

Single-crystal silicon avalanche photodiode (APD) direct-ionization spikes in
`EcalRecHitsEB` reconstruct as narrow jets with $\mathrm{EMF} \approx 1.0$, so
`EMF >= 0.01` does **not** remove them. Because barrel-edge crystals sit at
$|\eta_{\mathrm{EB}}| = 1.479$ (`|ieta| = 85`, $\cosh\eta \approx 2.31$), APD
spikes accumulate at $|\eta_1| \in [1.44, 1.52]$ and inflate mono-jet
denominators ($p_{T,1} > 120\ \mathrm{GeV}$ with no recoil jet).

`RecoHI/HiEgammaAlgos/HiEcalRecHitSpikeFilter` in `CMSSW_3_9_2_patch5` flags any
`EcalRecHitsEB` crystal with $E_T > E_{T,\min}$ (`avoidIeta85 = False`) that has:

- Swiss-cross $1 - E_4 / E_1 > 0.95$, where $E_4$ sums raw crystal energies
  (including negative electronics noise, wrapping in `iphi` $1 \leftrightarrow 360$,
  with 3 neighbours at `|ieta| = 85` and `|ieta| = 1`), **or**
- $|t| > 4.0\ \mathrm{ns}$.

**Whole-event `minEt` threshold sensitivity in central PbPb:**

| `minEt` setting                                 | `HLT_HIJet50U` pass rate                  | `0-10%` MinBias veto rate | `50-100%` MinBias veto rate | Behaviour & historical usage                                                                                                                                                         |
| ----------------------------------------------- | ----------------------------------------- | ------------------------- | --------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `20.0 GeV` (`hiEcalSpikeFilter` module default) | `94.25%–94.31%`                           | `5.49%`                   | `0.18%`                     | Isolated photon analyses, or conservative whole-event cleaning (appears in `CmsHI/CVS_CmsHi` `37ffd5e4`; also vetoes `20-50 GeV` isolated EB crystal deposits outside jet cones)     |
| `50.0 GeV` (jet-scale whole-event veto)         | `97.46%` (`97.66%` in HIN-10-004 Table I) | `< 0.5%`                  | `< 0.05%`                   | High-$p_T$ jet analyses ($p_{T,1} > 120\ \mathrm{GeV}$, $p_{T,2} > 50\ \mathrm{GeV}$): restores out-of-cone dijets (`4` of `6` restored in 110 files); exact Table I job unconfirmed |

Use `scripts/hi_ecal_spike_filter.py` to evaluate either threshold on real EDM
`RECO` files.

## Iterative-cone UE-subtracted jets and Spring10_HI JEC

`recoCaloJets_iterativeConePu5CaloJets__RECO` stores $R = 0.5$ iterative-cone
calorimeter jets **after** variant-B $\eta$-ring underlying-event
pedestal subtraction, at raw (uncorrected) calorimeter scale. Apply factorized
`Spring10_HI` `L2Relative` $\times$ `L3Absolute` corrections before sorting in
$p_T$:

```python
# L2Relative formula in Spring10_HI_L2Relative_IC5Calo.txt:
# max(0.0001, p0 + x*(p1 + x*(p2 + x*(p3 + x*(p4 + x*p5))))), x = clip(pt_raw, pt_min, pt_max)
# L3Absolute formula in Spring10_HI_L3Absolute_IC5Calo.txt:
# p0 + p1 / ((log10(x))**p2 + p3), x = clip(pt_l2, pt_min, pt_max)
scale = l2_factor(eta_raw, pt_raw) * l3_factor(pt_raw * l2_factor(eta_raw, pt_raw))
pt_corr = pt_raw * scale
mass_corr = mass_raw * scale
order = ak.argsort(pt_corr, axis=1, ascending=False)   # ALWAYS re-sort after JEC!
```

Typical JEC systematic uncertainty on `IC5Calo` in 2010 PbPb is $5\%$ absolute
energy scale plus a $2\%$ central-to-peripheral underlying-event residual scale
shift. Note that two public factorized `IC5Calo` payloads exist on CERN Open
Data: `Spring10_HI_IC5Calo` (Record [470](https://opendata.cern.ch/record/470))
and the released global-tag payload `GR_R_39X_V6B` (`Spring10_V8_IC5Calo`,
Record [1811](https://opendata.cern.ch/record/1811)), which reproduces jet
$p_T$ slightly more closely in archived 2010 CMS dijet ntuples ($p_{T,2} > 50\ \mathrm{GeV}$)
while selecting nearly identical final back-to-back dijet samples (`96` vs `98`
in a 110-file check, `93` in common).

## Inclusive vs dijet-conditional quench observables

Do not confuse the denominators of CMS-HIN-10-004's balanced fractions:

- **Azimuthal correlation & asymmetry shapes (`Figs. 6, 8`)**: normalized per
  selected back-to-back dijet ($p_{T,1} > 120\ \mathrm{GeV}$,
  $p_{T,2} > 50\ \mathrm{GeV}$, $|\eta_{1,2}| < 2$, $\Delta\phi_{12} > 2\pi/3$).
- **Balanced fractions $R_B(\Delta\phi_{12} > 3.026)$ (`Fig. 7`) and
  $R_B(A_J < 0.15)$ (`Fig. 9`)**: denominator is **all events with a leading jet
  $p_{T,1} > 120\ \mathrm{GeV}$, $|\eta_1| < 2$** (including events where the
  recoil jet was quenched below $50\ \mathrm{GeV}$ and mono-jets!), NOT the
  selected dijet count. Dividing by $N_{\mathrm{dijet}}$ artificially inflates
  $R_B$ by $\sim 1 / 0.84 \approx 1.19\times$.
- **Mean fractional imbalance $\langle (p_{T,1} - p_{T,2}) / p_{T,1} \rangle$
  (`Fig. 10`)**: evaluated over leading-jet events ($p_{T,1} > 120\ \mathrm{GeV}$,
  $|\eta_1| < 2$) with the highest-$p_T$ away-side ($\Delta\phi_{12} > 2\pi/3$,
  $|\eta_2| < 2$) jet above $p_{T,2} > 25\ \mathrm{GeV}$, setting $p_{T,2} = 0$
  (ratio $= 1$) when no away-side jet exceeds $25\ \mathrm{GeV}$ (note also
  whether MC jet-finding efficiency reweighting is applied; unweighted
  detector-level Open Data sits $\sim 0.04\text{--}0.08$ higher due to
  $25\text{--}50\ \mathrm{GeV}$ calo-jet inefficiency).
