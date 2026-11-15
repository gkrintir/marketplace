---
name: heavy-ion-reco
description: >-
  Use when analysing CMS Heavy-Ion (PbPb, XeXe, pPb) collisions or Run-1 CMS
  EDM RECO / AOD Open Data in Python with uproot: computing HF forward
  calorimeter centrality (hiCentrality, sum E_HF, Glauber N_part), applying
  Heavy-Ion offline collision and noise filters (HF tower coincidence,
  hiSelectedVertex, HcalNoiseSummary, ECAL Barrel APD spike cleaning with
  HiEcalRecHitSpikeFilter), applying UE-subtracted calorimeter jet corrections
  (iterativeConePu5CaloJets, akPu*CaloJets), decoding EventAuxiliary or
  edm::TriggerResults from EDM RECO, or debugging errors such as
  nproblemRBXs_ rejecting central PbPb events, mono-jet spikes at |eta| ~ 1.48,
  "SyntaxError: invalid syntax" on "./" EDM branch names, "cannot produce
  Awkward Arrays for interpretation AsObjects", or "ClientResponseError object
  is not subscriptable" when streaming EOS HTTP files.
metadata:
  summary: Heavy-Ion PbPb centrality, APD/HCAL noise cleaning, UE-subtracted jets, and Run-1 EDM RECO streaming with uproot
---

# Heavy-Ion RECO Analysis

## Overview

CMS Heavy-Ion (`PbPb`, `XeXe`, `pPb`) analyses differ from proton-proton
NanoAOD workflows in three critical ways:

1. **Event centrality & underlying event:** observables are binned in
   collision centrality from Forward Hadron (`HF`, $2.9 < |\eta| < 5.2$) total
   energy $\sum E_{\mathrm{HF}}$ mapped to Glauber $\langle N_{\mathrm{part}} \rangle$,
   and jets require $\eta$-ring underlying-event subtraction before factorized `L2Relative`
   $\times$ `L3Absolute` JEC.
2. **High-multiplicity detector cleaning:** `pp` noise flags such as
   `HcalNoiseSummary.nproblemRBXs_ == 0` fail in $>99\%$ of central PbPb
   collisions, while single-crystal ECAL Barrel APD spikes (`EMF ~ 1.0`) pile up
   at $|\eta_{\mathrm{EB}}| = 1.479$ unless `HiEcalRecHitSpikeFilter` is applied.
3. **Run-1 EDM `RECO` format:** 2010–2011 Heavy-Ion Open Data is distributed as
   EDM `RECO` (`CMSSW_3_9_2_patch5`) rather than NanoAOD, requiring `uproot`
   patterns for `./` branch names, `EventAuxiliary`, `edm::TriggerResults`, and
   batched HTTP range reads.

[scripts/hi_ecal_spike_filter.py](scripts/hi_ecal_spike_filter.py) runs
crystal-level `HiEcalRecHitSpikeFilter` checks (`--self-test` or on a local /
remote EDM `RECO` file) without CMSSW. Offline kinematics and $\sum E_{\mathrm{HF}}$
in `CMSSW_3_9_2_patch5` are in **GeV**.

## When to Use

- Analysing CMS Heavy-Ion Open Data (`HICorePhysics`, `HIAllPhysics`, PbPb/pPb)
- Classifying centrality from `hiCentrality` and mapping to Glauber $\langle N_{\mathrm{part}} \rangle$
- Cleaning Heavy-Ion collisions (HF coincidence, vertex, HCAL RBX noise, ECAL APD spikes)
- Applying UE-subtracted `IC5Calo` / `akPu*` JEC and building dijet quenching observables ($A_J$, $R_B$)
- Reading Run-1 EDM `RECO` / `AOD` trees in Python with `uproot` over XRootD or EOS HTTPS

## Key Concepts

| Concept             | Branch / object in 2010 PbPb `RECO`                                        | Key rule                                                                                                      |
| ------------------- | -------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------- |
| Centrality variable | `recoCentrality_hiCentrality__RECO` (`etHFhitSumPlus_ + etHFhitSumMinus_`) | Stores **total HF energy $\sum E_{\mathrm{HF}}$ in GeV** (up to ~150 TeV), despite the `et*` prefix           |
| HF coincidence      | `etHFtowerSumPlus_ >= 3.0` and `etHFtowerSumMinus_ >= 3.0`                 | Suppresses ultra-peripheral EM dissociation and beam-gas (~75% of `HLT_HIJet50U` pass)                        |
| Heavy-Ion vertex    | `recoVertexs_hiSelectedVertex__RECO`                                       | Require `validity_`, `ndof_ > 0`, and `\|fZ\| <= 25` cm                                                       |
| HCAL RBX filter     | `HcalNoiseSummary_hcalnoise__RECO` (`maxzeros_ < 100`)                     | **Never** require `nproblemRBXs_ == 0` in PbPb (flags all 36 HB RBXs in central collisions)                   |
| ECAL APD spikes     | `EcalRecHitsEB` (`HiEcalRecHitSpikeFilter`)                                | Reject crystal $E_T > E_{T,\min}$ with $1 - E_4/E_1 > 0.95$ or $\|t\| > 4\ \mathrm{ns}$ (`avoidIeta85=False`) |
| UE-subtracted jets  | `recoCaloJets_iterativeConePu5CaloJets__RECO`                              | Raw UE-subtracted scale; apply `L2Relative * L3Absolute` JEC, **then** sort by `pt_corr`                      |

## Canonical Patterns

### 1. Read EDM `RECO` branches safely and compute HF centrality

```python
import awkward as ak, numpy as np, uproot

# Always pass filter_name= and how=dict on EDM "./" branch names; use HTTPSource over HTTP
f = uproot.open(url, handler=uproot.source.http.HTTPSource)
ev = f["Events"]
LEAF_C = "recoCentrality_hiCentrality__RECO.obj."
PFX_C = f"recoCentrality_hiCentrality__RECO./recoCentrality_hiCentrality__RECO.obj/{LEAF_C}"
arr = ev.arrays(
    filter_name=[PFX_C + s for s in ("etHFhitSumPlus_", "etHFhitSumMinus_", "etHFtowerSumPlus_", "etHFtowerSumMinus_", "pixelMultiplicity_")],
    how=dict,  # keys in arr are leaf branch names: LEAF_C + s
)
hf_sum_gev = ak.to_numpy(arr[LEAF_C + "etHFhitSumPlus_"] + arr[LEAF_C + "etHFhitSumMinus_"])
hf_coinc = ak.to_numpy((arr[LEAF_C + "etHFtowerSumPlus_"] >= 3.0) & (arr[LEAF_C + "etHFtowerSumMinus_"] >= 3.0))

# 2010 2.76 TeV PbPb 5-bin sum(E_HF) [GeV] edges: 50-100%, 30-50%, 20-30%, 10-20%, 0-10%
HF_EDGES_GEV = np.array([0.0, 13097.866, 35392.680, 53479.020, 79370.000, np.inf])
cent_idx = np.digitize(hf_sum_gev, HF_EDGES_GEV) - 1
```

### 2. Heavy-Ion event cleaning and ECAL Barrel APD spike rejection

Run the bundled script to verify or measure `HiEcalRecHitSpikeFilter` rejection
on an EDM `RECO` file:

```bash
python3 <skill>/scripts/hi_ecal_spike_filter.py --self-test
# SELF_TEST=PASS  CHECKS=6/6

python3 <skill>/scripts/hi_ecal_spike_filter.py https://opendata.cern.ch/eos/opendata/cms/hidata/HIRun2010/HICorePhysics/RECO/PromptReco-v3/000/151/238/92004464-88F1-DF11-A6EF-001D09F23A20.root --max-events 200 --min-et 20.0
# EVENTS=200  PREFILTER_CANDIDATES=40  PASSED=162  REJECTED=38  PASS_FRACTION=0.810000  MIN_ET_GEV=20
```

Choose `--min-et` to match the physics objects:

- **Isolated photon / $Z \to ee$ analyses (`min_et = 20.0` GeV):** matches the
  default `hiEcalSpikeFilter` module in `CMSSW_3_9_2_patch5`. Note that as a
  whole-event veto, `20 GeV` rejects ~5.5% of `0-10%` central Minimum-Bias
  collisions due to multiplicity-scaling slow-neutron APD hits of `20-50 GeV`.
- **High-$p_T$ jet analyses (`min_et = 50.0` GeV, plus leading-jet
  `mEnergyFractionEm >= 0.01`):** removes every jet-scale APD spike (including
  at the `|ieta| = 85` barrel boundary) without biasing central PbPb event counts.

### 3. Apply factorized JEC and compute inclusive balanced fractions ($R_B$)

```python
# Flatten jagged raw jets, evaluate L2Relative(eta_raw, pt_raw) * L3Absolute(pt_l2), unflatten, and re-sort
counts = ak.num(pt_raw, axis=1)
scale = ak.unflatten(jec_l2l3(ak.to_numpy(ak.flatten(eta_raw)), ak.to_numpy(ak.flatten(pt_raw))), counts)
jets = ak.zip({"pt": pt_raw * scale, "eta": eta_raw, "phi": phi_raw, "mass": mass_raw * scale, "emf": emf_raw})
jets = jets[ak.argsort(jets.pt, axis=1, ascending=False)]

in_acc = jets[abs(jets.eta) < 2.0]
padded = ak.pad_none(in_acc, 2, axis=1)
j1, j2 = padded[:, 0], padded[:, 1]

# Inclusive leading-jet denominator (includes quenched recoil pT2 <= 50 GeV and mono-jets!)
lead_sel = (ak.fill_none(j1.pt, -1.0) > 120.0) & (ak.fill_none(j1.emf, 0.0) >= 0.01)
dphi12 = np.abs((ak.fill_none(j1.phi, 0.0) - ak.fill_none(j2.phi, 0.0) + np.pi) % (2 * np.pi) - np.pi)
dijet_sel = lead_sel & (ak.fill_none(j2.pt, -1.0) > 50.0) & (dphi12 > 2.0 * np.pi / 3.0)
aj = (j1.pt - j2.pt) / (j1.pt + j2.pt)
# R_B(A_J < 0.15) = count(dijet_sel & (aj < 0.15)) / count(lead_sel)   <-- NOT divided by count(dijet_sel)!
```

## Gotchas

- **`HcalNoiseSummary.nproblemRBXs_ == 0` wipes out central PbPb events.** In
  `CMSSW_3_9_2_patch5`, `HcalNoiseSummary` evaluates `pp` occupancy thresholds
  and sets `nproblemRBXs_ = 36..38` in $>99.8\%$ of central PbPb events
  (collapsing $\langle \sum E_{\mathrm{HF}} \rangle_{\mathrm{MB}}$ from
  `27.5 TeV` to `0.96 TeV`). Use `maxzeros_ < 100` instead, matching
  `HeavyIonsAnalysis/Configuration/python/EventSelection_cff.py`.
- **Jet `EMF >= 0.01` does NOT remove ECAL Barrel APD spikes.** APD direct
  ionisation hits have `EMF ~ 1.0` and pile up at $|\eta_1| \in [1.44, 1.52]$
  (`|ieta| = 85`, where `avoidIeta85` is `False` in 2010 `RECO`), inflating
  mono-jet counts by $>10\times$. Apply `HiEcalRecHitSpikeFilter` on
  `EcalRecHitsEB`.
- **Dividing $R_B(A_J < 0.15)$ or $R_B(\Delta\phi_{12} > 3.026)$ by
  $N_{\mathrm{dijet}}$ gives values ~19% too high.** Published CMS Heavy-Ion
  balanced fractions (`CMS-HIN-10-004` Figs. 7, 9) normalize to **all** events
  with $p_{T,1} > 120\ \mathrm{GeV}$ ($|\eta_1| < 2$), including events whose
  recoil jet was quenched below $50\ \mathrm{GeV}$.
- **`SyntaxError: invalid syntax` in `ev.arrays([...])` on EDM branches.** When
  any branch in a positional list is absent, `uproot` parses `./` in EDM branch
  names via `ast.parse` and crashes. Always call
  `ev.arrays(filter_name=[...], how=dict)`, and probe branch existence with
  `try: ev[name] except KeyError:` (never `"a/b" in ev`).
- **`ValueError: cannot produce Awkward Arrays for interpretation AsObjects`**
  on `HIAllPhysics/HIRun2010-ZS-v2/RECO`: `ZS-v2` was written with
  `splitLevel=0`. Unpack `TBasket` records directly or use `library="np"` on
  small slices; see [references/edm-reco-streaming.md](references/edm-reco-streaming.md).
- **`HLT_HIMinBiasHfOrBSC` is missing from `hltTriggerSummaryAOD`.** L1-seeded
  Minimum-Bias paths produce no HLT filter objects (`filterKeys_` has length 0).
  Decode `edmTriggerResults_TriggerResults__HLT.obj` matched against
  `@trigger_paths` in `ParameterSets`.
- **`TypeError: 'ClientResponseError' object is not subscriptable`** or zlib
  decompression errors over EOS HTTPS: `uproot`'s default `FSSpecSource`
  (`aiohttp`) leaks exception objects on `503` retries, and EOS
  `multipart/byteranges` uses `\n\n` headers that an unbounded `\r\n\r\n` scan
  confuses with bytes inside compressed `TBasket` payloads. Bound header scans
  to `[after, after + 512]`; see [references/edm-reco-streaming.md](references/edm-reco-streaming.md).

## Reference Files

- [scripts/hi_ecal_spike_filter.py](scripts/hi_ecal_spike_filter.py): standalone
  CLI (`--self-test`, `--max-events`, `--min-et`, `--swiss-threshold`,
  `--time-threshold`) implementing lossless `seedId` pre-filtering and
  `HiEcalRecHitSpikeFilter` (`KEY=VALUE` stdout, exit 0/1/2).
- [references/tested-2010-pbpb.md](references/tested-2010-pbpb.md): read for
  2010 PbPb Open Data records (14010, 14011, 14202), $\sum E_{\mathrm{HF}}$
  centrality boundaries and $\langle N_{\mathrm{part}} \rangle$, offline
  cleaning cuts, photon vs jet `minEt` spike behaviour, and `Spring10_HI` JEC.
- [references/edm-reco-streaming.md](references/edm-reco-streaming.md): read
  when streaming Run-1 EDM `RECO` files with `uproot` (`EventAuxiliary`,
  `edm::TriggerResults` + `ParameterSets`, `splitLevel=0` binary unpacking, and
  EOS HTTP `multipart/byteranges` coalescing).

## Docs

- CMS Collaboration, "Observation and studies of jet quenching in PbPb
  collisions at nucleon-nucleon center-of-mass energy = 2.76 TeV",
  CMS-HIN-10-004, Phys. Rev. C 84 (2011) 024906, <https://arxiv.org/abs/1102.1957>
- CERN Open Data 2010 PbPb datasets: `HICorePhysics`
  <https://opendata.cern.ch/record/14010>, `HIAllPhysics`
  <https://opendata.cern.ch/record/14011>, validated runs JSON
  <https://opendata.cern.ch/record/14202>
- CMS Open Data 2010 Heavy-Ion dijet example (`Spring10_HI` IC5Calo JEC):
  <https://github.com/cms-opendata-analyses/HeavyIonDijetAsymmetry2010>
- CMSSW Heavy-Ion reconstruction and event selection (`CMSSW_3_9_X`):
  <https://github.com/cms-sw/cmssw> (`HeavyIonsAnalysis/Configuration`,
  `RecoHI/HiEgammaAlgos`, `RecoHI/HiJetAlgos`)
