# Streaming Run-1 CMS EDM RECO files in Python with uproot

Read this reference when reading Run-1 (2010–2012) CMS EDM `RECO` or `AOD`
files in Python with `uproot` (without CMSSW): avoiding `./` branch-name
`SyntaxError` crashes, decoding `EventAuxiliary` `(run, lumi, event)` and
`edm::TriggerResults` + `ParameterSets`, reading unsplit (`splitLevel=0`)
collections, and streaming thousands of multi-GB EOS HTTP files without `503`
throttling or `multipart/byteranges` decompression corruption. Tested with
Python 3.11.2, uproot 5.7.7, awkward 2.15.0 and requests 2.32 on 2,154
`HICorePhysics` (`splitLevel=99`) and 120 `HIAllPhysics` (`splitLevel=0`) files,
2026-11-14.

## Contents

1. [EDM branch lookup and `./` AST fallback gotchas](#edm-branch-lookup-and--ast-fallback-gotchas)
2. [Decoding `EventAuxiliary` and `edm::TriggerResults`](#decoding-eventauxiliary-and-edmtriggerresults)
3. [Split (`PromptReco`) vs unsplit (`ZS-v2`) collections](#split-promptreco-vs-unsplit-zs-v2-collections)
4. [High-throughput EOS HTTP `multipart/byteranges` streaming](#high-throughput-eos-http-multipartbyteranges-streaming)
5. [Lossless `seedId` pre-filtering and surgical `TBasket` reads](#lossless-seedid-pre-filtering-and-surgical-tbasket-reads)

## EDM branch lookup and `./` AST fallback gotchas

EDM `RECO` branch names contain dots and slashes, such as
`recoCaloJets_iterativeConePu5CaloJets__RECO./recoCaloJets_iterativeConePu5CaloJets__RECO.obj/recoCaloJets_iterativeConePu5CaloJets__RECO.obj.pt_`.
Two `uproot` traps trigger on these names:

1. **`SyntaxError: invalid syntax` in `ev.arrays([...])`**: When a positional
   list is passed to `ev.arrays([br1, br2])` and any branch is missing (for
   example in an unsplit `ZS-v2` file), `uproot` falls back to parsing the
   string as a Python expression with `ast.parse`, crashing on `./`. Always pass
   `filter_name=` and `how=dict`:

   ```python
   arr = ev.arrays(filter_name=[br1, br2], how=dict)   # never triggers ast.parse
   # Note: dict keys are leaf branch names (e.g. "recoCaloJets_...RECO.obj.pt_"), without parent "/" prefixes
   ```

2. **`"Parent/SubBranch" in ev` returns `False`**: Membership checks (`in`) on
   `uproot` `TTree` objects do not reliably resolve `/` subbranch paths across
   all split levels. Probe subbranches with `try: br = ev[name] except KeyError:`.

## Decoding `EventAuxiliary` and `edm::TriggerResults`

In Run-1 EDM `RECO`, `(run, lumi, event)` live in the unsplit `EventAuxiliary`
branch, and L1-seeded Minimum-Bias paths such as `HLT_HIMinBiasHfOrBSC` leave
`filterKeys_` empty in `triggerTriggerEvent_hltTriggerSummaryAOD` on every
event. Read `EventAuxiliary` and `edm::TriggerResults` directly from coalesced
`TBasket` bytes:

```python
import queue, re, struct, numpy as np, uproot

def read_all_baskets(branch):
    nb = int(np.sum(np.asarray(branch.member("fBasketSeek"), dtype=np.int64) > 0))
    seeks = branch.member("fBasketSeek")[:nb]
    nbytes = branch.member("fBasketBytes")[:nb]
    ranges = [(int(s), int(s + b)) for s, b in zip(seeks, nbytes)]
    chunks = branch.file.source.chunks(ranges, notifications=queue.Queue())
    return [
        uproot.models.TBasket.Model_TBasket.read(
            c, uproot.source.cursor.Cursor(c.start), {"basket_num": i}, branch.file, branch.file, branch
        )
        for i, c in enumerate(chunks)
    ]

def decode_event_aux(ev):
    runs, lumis, events = [], [], []
    for bsk in read_all_baskets(ev["EventAuxiliary"]):
        raw, off = bytes(bsk.data), bsk.byte_offsets
        for i in range(len(off) - 1):
            r, l, e = struct.unpack(">III", raw[off[i] + 43 : off[i] + 55])
            runs.append(r); lumis.append(l); events.append(e)
    return np.asarray(runs), np.asarray(lumis), np.asarray(events)
```

For `edmTriggerResults_TriggerResults__HLT.obj`:

- `npaths = struct.unpack(">I", rec[26:30])[0]` gives the number of HLT paths.
- Each path's 2-bit state (`1 == Pass`) is `rec[30 + 12 * p + 11] & 0x3`.
- The 16-byte `ParameterSetID` hex hash follows `40000013000810` in the record
  tail (`rec[30 + 12 * npaths :].hex()`), matching `blob[17:33].hex()` in
  `f["ParameterSets"]["IdToParameterSetsBlobs"]`, where `@trigger_paths=+s({...})`
  stores hex-encoded ASCII path names. Cache the `pset_hash -> [path_names]`
  dictionary to JSON: only a handful of HLT menus exist across an entire run era,
  eliminating a full `ParameterSets` read on every subsequent file.
- **Avoid `branch.entries_to_ranges_or_baskets(0, branch.num_entries)` over HTTP:**
  on `ParameterSets/IdToParameterSetsBlobs` (and other branches with unpopulated
  embedded basket-entry tables), `entries_to_ranges_or_baskets` falls back to
  scanning 80-byte `TBasket` headers one HTTP round-trip at a time (~31 s/file),
  whereas indexing `fBasketSeek[:nb]` and `fBasketBytes[:nb]` directly coalesces
  the read into ~0.7 s.

## Split (`PromptReco`) vs unsplit (`ZS-v2`) collections

- **`PromptReco-v3` (`splitLevel=99`)**: primitive subleaves (`pt_`, `eta_`,
  `fCoordinates.fZ`, `mEnergyFractionEm`) vectorize cleanly via `AwkwardForth`
  in `ev.arrays(filter_name=..., how=dict)`.
- **`ZS-v2` (`splitLevel=0`)**: `reco::Vertex` and `reco::CaloJet` are stored as
  opaque C++ objects (`AsObjects`) containing `ROOT::Math::PositionVector3D<...>`
  and fail in `ev.arrays()` with
  `ValueError: cannot produce Awkward Arrays for interpretation AsObjects`.
  Either read with `library="np"` on small samples or unpack the fixed big-endian
  `TBasket` byte offsets directly (verified bit-exact with `0` mismatches against
  `library="np"`; see `validation/hin-10-004/analyze_hin10004.py`).

## High-throughput EOS HTTP `multipart/byteranges` streaming

Three transport issues dominate when streaming multi-GB EDM `RECO` files over
`https://opendata.cern.ch`:

1. **`FSSpecSource` (`aiohttp`) chunk-exception bug**: When `aiohttp` is
   installed, `uproot.open("https://...")` defaults to `FSSpecSource`. On
   transient range timeouts or `503` responses, `fsspec`'s `_cat_ranges` returns
   `ClientResponseError` objects inside the result list instead of raising them,
   crashing `uproot/source/coalesce.py` with
   `TypeError: 'ClientResponseError' object is not subscriptable`. Pass
   `handler=uproot.source.http.HTTPSource` (for a few files) or a custom
   `requests.Session` `Source` subclass with `Retry(status_forcelist=[429, 500, 502, 503, 504])`.
2. **OpenShift sticky-session cookies & range coalescing**: Reject cookies on
   `requests.Session` (`DefaultCookiePolicy(allowed_domains=[])`) so concurrent
   workers load-balance across backend EOS HTTP pods instead of pinning to one
   pod. Batch up to 64 byte ranges per GET (`Range: bytes=s1-e1,s2-e2,...`):
   this reduces `PromptReco-v3` reads from ~215 HTTP GETs/file (~36 s) down to
   ~8 GETs/file (~5.5 s).
3. **EOS `multipart/byteranges` `\n\n` header delimiter trap**: CERN EOS HTTP
   (`EOSMULTIPARTBOUNDARY`) separates multipart headers from binary body bytes
   using `\n\n` (`LF LF`), **not** `\r\n\r\n` (`CRLF CRLF`). An unbounded
   `raw.find(b"\r\n\r\n", after)` scans past the header and matches accidental
   `0d 0a 0d 0a` 4-byte sequences inside compressed zlib `TBasket` streams,
   truncating payloads and failing decompression. Always bound the header search
   to `[after, after + 512]` and pick whichever of `\r\n\r\n` or `\n\n` appears
   earliest:

   ```python
   win = raw[after : after + 512]
   r_rn, r_nn = win.find(b"\r\n\r\n"), win.find(b"\n\n")
   if r_rn != -1 and (r_nn == -1 or r_rn <= r_nn):
       hdr_end, body_start = after + r_rn, after + r_rn + 4
   elif r_nn != -1:
       hdr_end, body_start = after + r_nn, after + r_nn + 2
   ```

4. **Multi-process sharding for EDM `TTree` metadata**: Deserializing ~4,000 EDM
   `TBranch` headers in `f["Events"]` takes ~1.0 s of pure-Python CPU time per
   file under the GIL, capping a single Python process with threads at
   ~0.65 files/s. Shard files across `os.cpu_count()` worker processes (each
   running 6–12 HTTP threads) and write atomic `.tmp -> .pkl` per-file caches so
   interrupted runs resume in milliseconds.

## Lossless `seedId` pre-filtering and surgical `TBasket` reads

Reading all 1,500–4,000 `EcalRecHitsEB` crystals in every PbPb event downloads
tens of megabytes per file. Because `hybridSuperClusters` in `CMSSW_3_9_2_patch5`
seeds a `reco::BasicCluster` around every barrel crystal with
$E_T > 1.0\ \mathrm{GeV}$, any event containing a crystal with
$E_T > 20\ \mathrm{GeV}$ is guaranteed to have a `hybridBarrelBasicClusters`
entry whose `energy_ / cosh(min_eta_from_seedId)` exceeds $20\ \mathrm{GeV}$
(verified lossless: `0 / 5,000` false negatives).

Read the tiny split `hybridBarrelBasicClusters` (`energy_`, `seedId_.id_`) leaves
first, identify the ~2–5% of candidate event indices, map them via
`np.searchsorted(fBasketEntry, entry_indices, side="right") - 1`, and fetch only
those exact `fBasketSeek` byte ranges (passing `notifications=queue.Queue()` to
`source.chunks(...)`; `notifications=None` raises
`AttributeError: 'NoneType' object has no attribute 'put'` in `HTTPSource`).
See `scripts/hi_ecal_spike_filter.py` for the complete implementation.
