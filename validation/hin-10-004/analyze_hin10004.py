#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11,<3.14"
# dependencies = [
#   "typer>=0.12",
#   "uproot>=5.3",
#   "awkward>=2.6",
#   "vector>=1.4",
#   "hist>=2.7",
#   "scipy>=1.11",
#   "numpy>=1.26",
#   "requests>=2.31",
#   "aiohttp>=3.9",
# ]
# ///
"""End-to-end CMS Open Data analysis for CMS-HIN-10-004 (2.76 TeV PbPb dijet quenching).

Streams validated 2010 CMS HICorePhysics and HIAllPhysics RECO files over HTTPS,
decodes EventAuxiliary and edm::TriggerResults, applies Heavy-Ion offline event
cleaning, HF centrality classification, and Spring10_HI IC5Calo L2L3 JEC, builds
hist.Hist objects and Clopper-Pearson intervals, caches compact event arrays in
cache/, and writes summary.json.
"""
from __future__ import annotations

import hashlib
import json
import logging
import math
import os
import pickle
import queue
import re
import struct
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Optional

import awkward as ak
import hist
from hist import Hist
import numpy as np
import requests
from requests.adapters import HTTPAdapter
import scipy
from scipy import stats
import typer
import uproot
from urllib3.util.retry import Retry
import vector

vector.register_awkward()

app = typer.Typer(add_completion=False)
log = logging.getLogger("hin10004")

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"

CENTRALITY_BINS = [
    ("50-100%", 0.0, 13097.866, 22.0, 2.0),
    ("30-50%", 13097.866, 35392.680, 108.0, 5.0),
    ("20-30%", 35392.680, 53479.020, 187.0, 5.0),
    ("10-20%", 53479.020, 79370.000, 261.0, 4.0),
    ("0-10%", 79370.000, float("inf"), 355.0, 3.0),
]
CENTRALITY_NAMES = [c[0] for c in CENTRALITY_BINS]
MERGED_CENTRALITY_MAP = {
    "50-100%": "30-100%",
    "30-50%": "30-100%",
    "20-30%": "10-30%",
    "10-20%": "10-30%",
    "0-10%": "0-10%",
}
MERGED_CENTRALITY_NAMES = ["30-100%", "10-30%", "0-10%"]

CORE_SPLIT_BRANCHES = [
    "recoCaloJets_iterativeConePu5CaloJets__RECO./recoCaloJets_iterativeConePu5CaloJets__RECO.obj/recoCaloJets_iterativeConePu5CaloJets__RECO.obj.pt_",
    "recoCaloJets_iterativeConePu5CaloJets__RECO./recoCaloJets_iterativeConePu5CaloJets__RECO.obj/recoCaloJets_iterativeConePu5CaloJets__RECO.obj.eta_",
    "recoCaloJets_iterativeConePu5CaloJets__RECO./recoCaloJets_iterativeConePu5CaloJets__RECO.obj/recoCaloJets_iterativeConePu5CaloJets__RECO.obj.phi_",
    "recoCaloJets_iterativeConePu5CaloJets__RECO./recoCaloJets_iterativeConePu5CaloJets__RECO.obj/recoCaloJets_iterativeConePu5CaloJets__RECO.obj.mass_",
    "recoCaloJets_iterativeConePu5CaloJets__RECO./recoCaloJets_iterativeConePu5CaloJets__RECO.obj/recoCaloJets_iterativeConePu5CaloJets__RECO.obj.m_specific.mEnergyFractionEm",
    "recoCentrality_hiCentrality__RECO./recoCentrality_hiCentrality__RECO.obj/recoCentrality_hiCentrality__RECO.obj.etHFhitSumPlus_",
    "recoCentrality_hiCentrality__RECO./recoCentrality_hiCentrality__RECO.obj/recoCentrality_hiCentrality__RECO.obj.etHFhitSumMinus_",
    "recoCentrality_hiCentrality__RECO./recoCentrality_hiCentrality__RECO.obj/recoCentrality_hiCentrality__RECO.obj.etHFtowerSumPlus_",
    "recoCentrality_hiCentrality__RECO./recoCentrality_hiCentrality__RECO.obj/recoCentrality_hiCentrality__RECO.obj.etHFtowerSumMinus_",
    "recoCentrality_hiCentrality__RECO./recoCentrality_hiCentrality__RECO.obj/recoCentrality_hiCentrality__RECO.obj.pixelMultiplicity_",
    "recoVertexs_hiSelectedVertex__RECO./recoVertexs_hiSelectedVertex__RECO.obj/recoVertexs_hiSelectedVertex__RECO.obj.position_.fCoordinates.fZ",
    "recoVertexs_hiSelectedVertex__RECO./recoVertexs_hiSelectedVertex__RECO.obj/recoVertexs_hiSelectedVertex__RECO.obj.ndof_",
    "recoVertexs_hiSelectedVertex__RECO./recoVertexs_hiSelectedVertex__RECO.obj/recoVertexs_hiSelectedVertex__RECO.obj.validity_",
    "HcalNoiseSummary_hcalnoise__RECO./HcalNoiseSummary_hcalnoise__RECO.obj/HcalNoiseSummary_hcalnoise__RECO.obj.maxzeros_",
    "HcalNoiseSummary_hcalnoise__RECO./HcalNoiseSummary_hcalnoise__RECO.obj/HcalNoiseSummary_hcalnoise__RECO.obj.nproblemRBXs_",
    "recoCaloClusters_uncleanedHybridSuperClusters_hybridBarrelBasicClusters_RECO./recoCaloClusters_uncleanedHybridSuperClusters_hybridBarrelBasicClusters_RECO.obj/recoCaloClusters_uncleanedHybridSuperClusters_hybridBarrelBasicClusters_RECO.obj.energy_",
    "recoCaloClusters_uncleanedHybridSuperClusters_hybridBarrelBasicClusters_RECO./recoCaloClusters_uncleanedHybridSuperClusters_hybridBarrelBasicClusters_RECO.obj/recoCaloClusters_uncleanedHybridSuperClusters_hybridBarrelBasicClusters_RECO.obj.seedId_.id_",
]

ECAL_EB_BRANCH = "EcalRecHitsSorted_ecalRecHit_EcalRecHitsEB_RECO./EcalRecHitsSorted_ecalRecHit_EcalRecHitsEB_RECO.obj"
CACHE_CORE_PREFIX = "core_v2"
CACHE_MB_PREFIX = "mb_v2"

import http.cookiejar as cookiejar

_PSET_LOCK = threading.Lock()
_PSET_TRIGGER_PATHS: dict[str, list[str]] = {}
PSET_CACHE_JSON = BASE_DIR / "cache" / "pset_trigger_paths.json"
_TLS_SESSION = threading.local()


def _get_http_session() -> requests.Session:
    sess = getattr(_TLS_SESSION, "session", None)
    if sess is None:
        sess = requests.Session()
        # Reject OpenShift sticky-session cookies so each worker load-balances across backend pods
        sess.cookies.set_policy(cookiejar.DefaultCookiePolicy(allowed_domains=[]))
        adapter = HTTPAdapter(
            pool_connections=6,
            pool_maxsize=6,
            max_retries=Retry(total=5, status_forcelist=[429, 500, 502, 503, 504], backoff_factor=0.5),
        )
        sess.mount("https://", adapter)
        sess.mount("http://", adapter)
        _TLS_SESSION.session = sess
    return sess


class PooledHTTPSource(uproot.source.chunk.Source):
    """Keep-alive HTTPS range reader with multipart/byteranges batching and 503 retry."""

    MAX_RANGES_PER_REQ = 64

    def __init__(self, file_path: str, **options):
        self._file_path = file_path
        self._num_bytes: Optional[int] = None
        self._timeout = options.get("timeout", 45)
        self._num_requests = 0
        self._num_requested_chunks = 0
        self._num_requested_bytes = 0
        self._closed = False

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self._closed = True

    @property
    def closed(self) -> bool:
        return self._closed

    @property
    def num_bytes(self) -> int:
        if self._num_bytes is None:
            r = _get_http_session().get(self._file_path, headers={"Range": "bytes=0-0"}, timeout=self._timeout)
            cr = r.headers.get("Content-Range", "")
            self._num_bytes = int(cr.split("/")[-1]) if "/" in cr else int(r.headers["Content-Length"])
        return self._num_bytes

    def chunk(self, start: int, stop: int) -> uproot.source.chunk.Chunk:
        self._num_requests += 1
        self._num_requested_chunks += 1
        self._num_requested_bytes += stop - start
        r = _get_http_session().get(self._file_path, headers={"Range": f"bytes={start}-{stop - 1}"}, timeout=self._timeout)
        if r.status_code not in (200, 206):
            raise OSError(f"HTTP {r.status_code} for {self._file_path}")
        data = r.content if r.status_code == 206 else r.content[start:stop]
        fut = uproot.source.futures.TrivialFuture(np.frombuffer(data, dtype=np.uint8))
        return uproot.source.chunk.Chunk(self, start, stop, fut)

    @staticmethod
    def _parse_multipart_ranges(raw: bytes, ctype: str) -> list[tuple[int, int, memoryview]]:
        m_b = re.search(r'boundary="?([^"\s;]+)"?', ctype, re.I)
        if not m_b:
            raise OSError(f"Missing multipart boundary in Content-Type: {ctype}")
        b_delim = b"--" + m_b.group(1).encode("ascii")
        mv = memoryview(raw)
        parts: list[tuple[int, int, memoryview]] = []
        pos = 0
        while True:
            idx = raw.find(b_delim, pos)
            if idx < 0:
                break
            after = idx + len(b_delim)
            if raw[after : after + 2] == b"--":
                break
            hdr_rn = raw.find(b"\r\n\r\n", after, after + 512)
            hdr_nn = raw.find(b"\n\n", after, after + 512)
            if hdr_rn >= 0 and (hdr_nn < 0 or hdr_rn < hdr_nn):
                hdr_end, sep_len = hdr_rn, 4
            elif hdr_nn >= 0:
                hdr_end, sep_len = hdr_nn, 2
            else:
                break
            m_cr = re.search(rb"Content-Range:\s*bytes\s+(\d+)-(\d+)/", raw[after:hdr_end], re.I)
            if not m_cr:
                raise OSError("Missing Content-Range in multipart/byteranges part")
            p_start = int(m_cr.group(1))
            p_stop = int(m_cr.group(2)) + 1
            data_start = hdr_end + sep_len
            p_len = p_stop - p_start
            parts.append((p_start, p_stop, mv[data_start : data_start + p_len]))
            pos = data_start + p_len
        return parts

    def chunks(self, ranges: list[tuple[int, int]], notifications: queue.Queue) -> list[uproot.source.chunk.Chunk]:
        if not ranges:
            return []
        indexed = sorted(enumerate(ranges), key=lambda x: x[1][0])
        groups = []
        cur_start, cur_stop = indexed[0][1]
        cur_items = [indexed[0]]
        for item in indexed[1:]:
            _, (s, e) = item
            if s - cur_stop <= 4096 and (max(cur_stop, e) - cur_start) <= 512 * 1024:
                cur_stop = max(cur_stop, e)
                cur_items.append(item)
            else:
                groups.append((cur_start, cur_stop, cur_items))
                cur_start, cur_stop = s, e
                cur_items = [item]
        groups.append((cur_start, cur_stop, cur_items))

        out_chunks: list[Any] = [None] * len(ranges)
        for b_idx in range(0, len(groups), self.MAX_RANGES_PER_REQ):
            batch = groups[b_idx : b_idx + self.MAX_RANGES_PER_REQ]
            self._num_requests += 1
            self._num_requested_chunks += sum(len(it) for _, _, it in batch)
            self._num_requested_bytes += sum(ge - gs for gs, ge, _ in batch)
            range_hdr = "bytes=" + ",".join(f"{gs}-{ge - 1}" for gs, ge, _ in batch)
            r = _get_http_session().get(self._file_path, headers={"Range": range_hdr}, timeout=self._timeout)
            if r.status_code not in (200, 206):
                raise OSError(f"HTTP {r.status_code} for {self._file_path}")
            ctype = r.headers.get("Content-Type", "")
            raw = r.content
            mv = memoryview(raw)

            if r.status_code == 200:
                for _, _, items in batch:
                    for orig_idx, (s, e) in items:
                        arr = np.frombuffer(bytes(mv[s:e]), dtype=np.uint8)
                        chk = uproot.source.chunk.Chunk(self, s, e, uproot.source.futures.TrivialFuture(arr))
                        out_chunks[orig_idx] = chk
                        notifications.put(chk)
            elif "multipart/byteranges" in ctype.lower():
                parts = self._parse_multipart_ranges(raw, ctype)
                exact_map = {(ps, pe): pmv for ps, pe, pmv in parts}
                for g_start, g_stop, items in batch:
                    buf = exact_map.get((g_start, g_stop))
                    base_off = g_start
                    if buf is None:
                        for ps, pe, pmv in parts:
                            if ps <= g_start and g_stop <= pe:
                                buf = pmv
                                base_off = ps
                                break
                    if buf is None:
                        single_chk = self.chunk(g_start, g_stop)
                        buf = memoryview(single_chk.raw_data)
                        base_off = g_start
                    for orig_idx, (s, e) in items:
                        arr = np.frombuffer(bytes(buf[s - base_off : e - base_off]), dtype=np.uint8)
                        chk = uproot.source.chunk.Chunk(self, s, e, uproot.source.futures.TrivialFuture(arr))
                        out_chunks[orig_idx] = chk
                        notifications.put(chk)
            else:
                cr = r.headers.get("Content-Range", "")
                m_cr = re.search(r"bytes\s+(\d+)-(\d+)/", cr, re.I)
                base_off = int(m_cr.group(1)) if m_cr else batch[0][0]
                for _, _, items in batch:
                    for orig_idx, (s, e) in items:
                        arr = np.frombuffer(bytes(mv[s - base_off : e - base_off]), dtype=np.uint8)
                        chk = uproot.source.chunk.Chunk(self, s, e, uproot.source.futures.TrivialFuture(arr))
                        out_chunks[orig_idx] = chk
                        notifications.put(chk)
        return out_chunks



def clopper_pearson(k: np.ndarray | float, n: np.ndarray | float, cl: float = 0.6827):
    """Exact Clopper-Pearson binomial confidence interval (68.27% default)."""
    k_arr = np.asarray(k, dtype=np.float64)
    n_arr = np.asarray(n, dtype=np.float64)
    alpha = (1.0 - cl) / 2.0
    with np.errstate(invalid="ignore", divide="ignore"):
        lo = np.where(k_arr > 0, stats.beta.ppf(alpha, k_arr, n_arr - k_arr + 1.0), 0.0)
        hi = np.where(k_arr < n_arr, stats.beta.ppf(1.0 - alpha, k_arr + 1.0, n_arr - k_arr), 1.0)
        eff = np.where(n_arr > 0, k_arr / n_arr, np.nan)
    return eff, lo, hi


def load_l2_table(path: Path) -> np.ndarray:
    lines = path.read_text(encoding="utf-8").splitlines()[1:]
    rows = [[float(x) for x in ln.split()] for ln in lines if ln.strip()]
    return np.asarray(rows, dtype=np.float64)


def compute_jec_factor(pt_raw: np.ndarray, eta: np.ndarray, l2_arr: np.ndarray) -> np.ndarray:
    """Evaluate Spring10_HI IC5Calo L2Relative * L3Absolute JEC factor."""
    idx = np.clip(np.searchsorted(l2_arr[:, 1], eta, side="right"), 0, len(l2_arr) - 1)
    x_l2 = np.clip(pt_raw, l2_arr[idx, 3], l2_arr[idx, 4])
    lx = np.log10(np.maximum(x_l2, 1e-6))
    c_l2 = (
        l2_arr[idx, 5]
        + l2_arr[idx, 6] * lx
        + l2_arr[idx, 7] * lx**2
        + l2_arr[idx, 8] * lx**3
        + l2_arr[idx, 9] * lx**4
    )
    x_l3 = np.clip(pt_raw * c_l2, 4.0, 5000.0)
    c_l3 = 1.00298 + 5.95184 / (np.log10(x_l3) ** 3.42297 + 4.15457)
    return c_l2 * c_l3


def apply_jec_jagged(pt_raw_jag: ak.Array, eta_jag: ak.Array, mass_raw_jag: ak.Array, l2_arr: np.ndarray):
    counts = ak.num(pt_raw_jag, axis=1)
    pt_flat = ak.to_numpy(ak.flatten(pt_raw_jag, axis=1)).astype(np.float64)
    eta_flat = ak.to_numpy(ak.flatten(eta_jag, axis=1)).astype(np.float64)
    mass_flat = ak.to_numpy(ak.flatten(mass_raw_jag, axis=1)).astype(np.float64)
    if len(pt_flat) == 0:
        return pt_raw_jag, mass_raw_jag
    scale = compute_jec_factor(pt_flat, eta_flat, l2_arr)
    return ak.unflatten(pt_flat * scale, counts), ak.unflatten(mass_flat * scale, counts)


def classify_centrality_idx(hf_hit_sum_gev: np.ndarray) -> np.ndarray:
    """Return integer centrality index in 0..4 matching CENTRALITY_BINS."""
    edges = np.array([13097.866, 35392.680, 53479.020, 79370.000], dtype=np.float64)
    return np.searchsorted(edges, hf_hit_sum_gev, side="right").astype(np.int8)


def read_all_baskets_coalesced(branch) -> list:
    """Fetch all TBaskets of a branch in one coalesced HTTP chunk request."""
    if branch.num_entries == 0 or branch.num_baskets == 0:
        return []
    rob = branch.entries_to_ranges_or_baskets(0, branch.num_entries)
    ranges = [r for _, r in rob if isinstance(r, tuple)]
    notif = queue.Queue()
    chunks = branch.file.source.chunks(ranges, notifications=notif) if ranges else []
    chunk_iter = iter(chunks)
    baskets = []
    for bi, r in rob:
        if isinstance(r, tuple):
            chk = next(chunk_iter)
            cur = uproot.source.cursor.Cursor(chk.start)
            bsk = uproot.models.TBasket.Model_TBasket.read(
                chk, cur, {"basket_num": bi}, branch.file, branch.file, branch
            )
            baskets.append(bsk)
        else:
            baskets.append(r)
    return baskets


def read_selected_baskets(branch, entry_indices: np.ndarray) -> dict[int, bytes]:
    """Fetch only the TBaskets needed for `entry_indices` via in-memory TBranch tables (zero header scans)."""
    if len(entry_indices) == 0:
        return {}
    nb = int(branch.num_baskets)
    offsets = np.asarray(branch.member("fBasketEntry")[: nb + 1], dtype=np.int64).copy()
    if offsets[-1] == 0 or offsets[-1] < branch.num_entries:
        offsets[-1] = int(branch.num_entries)
    seeks = np.asarray(branch.member("fBasketSeek")[:nb], dtype=np.int64)
    nbytes = np.asarray(branch.member("fBasketBytes")[:nb], dtype=np.int64)

    bsk_ids = np.searchsorted(offsets[1:], entry_indices, side="right")
    uniq_bsks = np.unique(bsk_ids)

    ranges = [(int(seeks[bi]), int(seeks[bi] + nbytes[bi])) for bi in uniq_bsks]
    notif = queue.Queue()
    chunks = []
    # Batch at most 8 ~850 KB EcalRecHitsEB baskets per multipart HTTP GET
    for s in range(0, len(ranges), 8):
        chunks.extend(branch.file.source.chunks(ranges[s : s + 8], notifications=notif))

    bsk_map = {}
    for bi, chk in zip(uniq_bsks, chunks):
        bi_int = int(bi)
        cur = uproot.source.cursor.Cursor(chk.start)
        bsk_map[bi_int] = uproot.models.TBasket.Model_TBasket.read(
            chk, cur, {"basket_num": bi_int}, branch.file, branch.file, branch
        )

    out: dict[int, bytes] = {}
    for eidx, bi in zip(entry_indices, bsk_ids):
        bi_int = int(bi)
        bsk = bsk_map[bi_int]
        local_i = int(eidx - offsets[bi_int])
        off = bsk.byte_offsets
        if off is None:
            out[int(eidx)] = bytes(bsk.data)
        else:
            out[int(eidx)] = bytes(bsk.data[off[local_i] : off[local_i + 1]])
    return out


def parse_eb_rechits_binary(raw_rec: bytes) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Decode EcalRecHitsSorted_ecalRecHit_EcalRecHitsEB_RECO.obj record into (ids, energies, times)."""
    header_len = 20
    n = struct.unpack(">I", raw_rec[header_len - 4 : header_len])[0]
    if n == 0:
        return np.empty(0, dtype=np.uint32), np.empty(0, dtype=np.float32), np.empty(0, dtype=np.float32)
    rem = len(raw_rec) - header_len
    stride = 54 if rem == n * 54 else (rem // n)
    dt = np.dtype(
        [
            ("pad", f"V{stride - 24}"),
            ("id", ">u4"),
            ("energy", ">f4"),
            ("time", ">f4"),
            ("flags", ">u4"),
            ("aux", ">u4"),
            ("flagBits", ">u4"),
        ]
    )
    arr = np.frombuffer(raw_rec, dtype=dt, count=n, offset=header_len)
    return arr["id"], arr["energy"], arr["time"]


def hi_ecal_spike_filter_exact(
    ids: np.ndarray,
    energies: np.ndarray,
    times: np.ndarray,
    min_et: float = 20.0,
    swiss_thresh: float = 0.95,
    time_thresh: float = 4.0,
) -> bool:
    """Exact 2010 CMS Heavy-Ion HiEcalRecHitSpikeFilter (EcalTools::swissCross, avoidIeta85=False)."""
    if len(ids) == 0:
        return True
    cand_e = np.flatnonzero(energies > min_et)
    if len(cand_e) == 0:
        return True

    cand_ids = ids[cand_e]
    zside_c = np.where((cand_ids & 0x10000) != 0, 1, -1)
    ieta_abs_c = ((cand_ids >> 9) & 0x7F).astype(np.int32)
    eta_c = zside_c * (ieta_abs_c - 0.5) * 0.017453292519943295
    et_c = energies[cand_e] / np.cosh(eta_c)

    hi_mask = et_c > min_et
    if not np.any(hi_mask):
        return True

    hi_indices = cand_e[hi_mask]
    hi_times = times[hi_indices]
    if np.any(np.abs(hi_times) > time_thresh):
        return False

    zside = np.where((ids & 0x10000) != 0, 1, -1)
    ieta_abs = ((ids >> 9) & 0x7F).astype(np.int32)
    row = np.where(zside < 0, -ieta_abs + 85, ieta_abs + 84)  # -85->0 .. -1->84, +1->85 .. +85->169
    col = ((ids & 0x1FF) - 1).astype(np.int32)  # 0..359

    grid = np.zeros((170, 360), dtype=np.float32)
    grid[row, col] = energies  # raw energy including negative noise per EcalTools::recHitE

    for idx in hi_indices:
        r = int(row[idx])
        c = int(col[idx])
        e1 = float(energies[idx])
        s4_sum = float(grid[r, (c - 1) % 360]) + float(grid[r, (c + 1) % 360])
        if r > 0:
            s4_sum += float(grid[r - 1, c])
        if r < 169:
            s4_sum += float(grid[r + 1, c])
        swiss = 1.0 - (s4_sum / e1) if e1 != 0.0 else 0.0
        if swiss > swiss_thresh:
            return False
    return True


def decode_unsplit_max_sc_et(ev) -> np.ndarray:
    """Fast binary pre-filter on ZS-v2 uncleanedHybridSuperClusters returning upper-bound max barrel ET per event."""
    br_name = "recoSuperClusters_uncleanedHybridSuperClusters__RECO./recoSuperClusters_uncleanedHybridSuperClusters__RECO.obj"
    try:
        sc_br = ev[br_name]
    except KeyError:
        return np.full(ev.num_entries, 1e6, dtype=np.float64)
    bsks = read_all_baskets_coalesced(sc_br)
    max_ets: list[float] = []
    for bsk in bsks:
        raw, off = bytes(bsk.data), bsk.byte_offsets
        for i in range(len(off) - 1):
            rec = raw[off[i] : off[i + 1]]
            nsc = struct.unpack(">I", rec[6:10])[0]
            pos = 10
            m_et = 0.0
            for _ in range(nsc):
                bc = struct.unpack(">I", rec[pos : pos + 4])[0] & 0x3FFFFFFF
                if bc >= 68:
                    ene = struct.unpack(">d", rec[pos + 20 : pos + 28])[0]
                    fx, fy, fz = struct.unpack(">ddd", rec[pos + 48 : pos + 72])
                    r_perp = math.hypot(fx, fy)
                    r_3d = math.sqrt(fx * fx + fy * fy + fz * fz)
                    et = (ene * r_perp / r_3d) if (r_perp > 50.0 and r_3d > 0.0) else ene
                    if et > m_et:
                        m_et = et
                pos += 4 + bc
            max_ets.append(m_et)
    return np.asarray(max_ets, dtype=np.float64)


def decode_event_aux(ev) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    bsks = read_all_baskets_coalesced(ev["EventAuxiliary"])
    runs, lumis, events = [], [], []
    for bsk in bsks:
        raw, off = bsk.data, bsk.byte_offsets
        for i in range(len(off) - 1):
            r, l, e = struct.unpack(">III", bytes(raw[off[i] + 43 : off[i] + 55]))
            runs.append(r)
            lumis.append(l)
            events.append(e)
    return (
        np.asarray(runs, dtype=np.int64),
        np.asarray(lumis, dtype=np.int64),
        np.asarray(events, dtype=np.int64),
    )


def ensure_pset_paths(f, pset_hash: str) -> Optional[list[str]]:
    with _PSET_LOCK:
        if not _PSET_TRIGGER_PATHS and PSET_CACHE_JSON.is_file():
            try:
                _PSET_TRIGGER_PATHS.update(json.loads(PSET_CACHE_JSON.read_text(encoding="utf-8")))
            except Exception:
                pass
        if pset_hash in _PSET_TRIGGER_PATHS:
            return _PSET_TRIGGER_PATHS[pset_hash]
    ps_bsks = read_all_baskets_coalesced(f["ParameterSets"]["IdToParameterSetsBlobs"])
    found_local: dict[str, list[str]] = {}
    for bsk in ps_bsks:
        raw, off = bytes(bsk.data), bsk.byte_offsets
        for i in range(len(off) - 1):
            blob = raw[off[i] : off[i + 1]]
            m = re.search(rb"@trigger_paths=\+s\(\{([0-9A-Fa-f,]+)\}\)", blob)
            if m:
                pid = blob[17:33].hex()
                paths = [bytes.fromhex(h.decode()).decode() for h in m.group(1).split(b",") if h]
                found_local[pid] = paths
    with _PSET_LOCK:
        _PSET_TRIGGER_PATHS.update(found_local)
        try:
            PSET_CACHE_JSON.parent.mkdir(parents=True, exist_ok=True)
            tmp_ps = PSET_CACHE_JSON.with_suffix(f".{os.getpid()}.{threading.get_ident()}.tmp")
            tmp_ps.write_text(json.dumps(_PSET_TRIGGER_PATHS), encoding="utf-8")
            tmp_ps.replace(PSET_CACHE_JSON)
        except Exception:
            pass
        return _PSET_TRIGGER_PATHS.get(pset_hash)


def decode_hlt_trigger_results(f, ev) -> dict[str, np.ndarray]:
    """Decode edm::TriggerResults HLT pass masks by path name."""
    tr_bsks = read_all_baskets_coalesced(
        ev["edmTriggerResults_TriggerResults__HLT./edmTriggerResults_TriggerResults__HLT.obj"]
    )
    states_list = []
    tail_sample = ""
    for bsk in tr_bsks:
        raw, off = bytes(bsk.data), bsk.byte_offsets
        for i in range(len(off) - 1):
            rec = raw[off[i] : off[i + 1]]
            npaths = struct.unpack(">I", rec[26:30])[0]
            if not tail_sample:
                tail_sample = rec[30 + 12 * npaths :].hex()
            states_list.append([rec[30 + 12 * p + 11] & 0x3 for p in range(npaths)])
    if not states_list:
        return {}
    states = np.asarray(states_list, dtype=np.uint8)
    # Extract 16-byte ParameterSetID hash from tail_sample (follows 40000013000810)
    m_hash = re.search(r"40000013000810([0-9a-f]{32})", tail_sample)
    pset_hash = m_hash.group(1) if m_hash else ""
    paths = ensure_pset_paths(f, pset_hash) if pset_hash else None
    if not paths or len(paths) != states.shape[1]:
        return {}
    return {pname: (states[:, idx] == 1) for idx, pname in enumerate(paths)}


def decode_unsplit_vertices_and_jets(ev):
    """Fast binary parser for unsplit ZS-v2 reco::Vertex and reco::CaloJet collections."""
    v_bsks = read_all_baskets_coalesced(
        ev["recoVertexs_hiSelectedVertex__RECO./recoVertexs_hiSelectedVertex__RECO.obj"]
    )
    vtx_ok = []
    for bsk in v_bsks:
        raw, off = bytes(bsk.data), bsk.byte_offsets
        for i in range(len(off) - 1):
            rec = raw[off[i] : off[i + 1]]
            nv = struct.unpack(">I", rec[6:10])[0]
            pos = 10
            ok = False
            for _ in range(nv):
                bc = struct.unpack(">I", rec[pos : pos + 4])[0] & 0x3FFFFFFF
                v_rec = rec[pos : pos + 4 + bc]
                ndof = struct.unpack(">d", v_rec[18:26])[0]
                vz = struct.unpack(">d", v_rec[62:70])[0]
                valid = bool(v_rec[-1])
                if valid and ndof > 0.0 and abs(vz) <= 25.0:
                    ok = True
                    break
                pos += 4 + bc
            vtx_ok.append(ok)

    j_bsks = read_all_baskets_coalesced(
        ev["recoCaloJets_iterativeConePu5CaloJets__RECO./recoCaloJets_iterativeConePu5CaloJets__RECO.obj"]
    )
    pts, etas, phis, masses, emfs = [], [], [], [], []
    for bsk in j_bsks:
        raw, off = bytes(bsk.data), bsk.byte_offsets
        for i in range(len(off) - 1):
            rec = raw[off[i] : off[i + 1]]
            nj = struct.unpack(">I", rec[6:10])[0]
            pos = 10
            ev_pt, ev_eta, ev_phi, ev_m, ev_emf = [], [], [], [], []
            for _ in range(nj):
                bc = struct.unpack(">I", rec[pos : pos + 4])[0] & 0x3FFFFFFF
                jrec = rec[pos : pos + 4 + bc]
                pt, eta, phi, mass = struct.unpack(">ffff", jrec[54:70])
                emf = struct.unpack(">f", jrec[-8:-4])[0]
                ev_pt.append(pt)
                ev_eta.append(eta)
                ev_phi.append(phi)
                ev_m.append(mass)
                ev_emf.append(emf)
                pos += 4 + bc
            pts.append(ev_pt)
            etas.append(ev_eta)
            phis.append(ev_phi)
            masses.append(ev_m)
            emfs.append(ev_emf)

    return (
        np.asarray(vtx_ok, dtype=bool),
        ak.Array(pts),
        ak.Array(etas),
        ak.Array(phis),
        ak.Array(masses),
        ak.Array(emfs),
    )


def in_golden_mask(runs: np.ndarray, lumis: np.ndarray, golden: dict[str, list[list[int]]]) -> np.ndarray:
    out = np.zeros(len(runs), dtype=bool)
    for r_val in np.unique(runs):
        intervals = golden.get(str(int(r_val)))
        if not intervals:
            continue
        r_mask = runs == r_val
        l_sub = lumis[r_mask]
        ok = np.zeros(len(l_sub), dtype=bool)
        for lo, hi in intervals:
            ok |= (l_sub >= lo) & (l_sub <= hi)
        out[r_mask] = ok
    return out


def process_core_file(raw_url: str, golden: dict, l2_arr: np.ndarray, retries: int = 4) -> Optional[dict[str, Any]]:
    """Process one HICorePhysics PromptReco-v3 split-branch RECO file."""
    url = raw_url.replace("root://eospublic.cern.ch/", "https://opendata.cern.ch/")
    for attempt in range(retries):
        try:
            with uproot.open(url, timeout=45, handler=PooledHTTPSource) as f:
                ev = f["Events"]
                runs, lumis, _ = decode_event_aux(ev)
                n_total = len(runs)
                if n_total == 0:
                    return None
                gold = in_golden_mask(runs, lumis, golden)
                hlt_paths = decode_hlt_trigger_results(f, ev)
                pass_j50 = hlt_paths.get("HLT_HIJet50U", np.zeros(n_total, dtype=bool))
                pass_j35 = hlt_paths.get("HLT_HIJet35U", np.zeros(n_total, dtype=bool))
                if not np.any(gold & (pass_j50 | pass_j35)):
                    empty_f = np.array([], dtype=np.float64)
                    return {
                        "runs": np.unique(runs[gold]).tolist(),
                        "cutflow": {
                            "total_events_read": int(n_total),
                            "golden_json_certified": int(gold.sum()),
                            "hlt_hijet50u": 0,
                            "hf_coincidence": 0,
                            "primary_vertex": 0,
                            "ecal_spike_cleaned": 0,
                            "noise_and_scraping_cleaned": 0,
                            "two_jets_pt35_eta2": 0,
                            "leading_jet_pt120": 0,
                            "subleading_jet_pt50": 0,
                            "dijet_dphi_gt_2pi3": 0,
                            "nprob_zero_among_clean": 0,
                        },
                        "j50_hf_gev": empty_f,
                        "j50_pt1": empty_f,
                        "j50_eta1": empty_f,
                        "j50_phi1": empty_f,
                        "j50_pt1_raw": empty_f,
                        "j50_pt2": empty_f,
                        "j50_eta2": empty_f,
                        "j50_phi2": empty_f,
                        "j50_pt2_raw": empty_f,
                        "j50_dphi12": empty_f,
                        "j35_pt1": empty_f,
                        "j35_pass_j50": np.array([], dtype=bool),
                    }
                arrs = ev.arrays(filter_name=CORE_SPLIT_BRANCHES, how=dict)

                def get_leaf(suffix: str):
                    for k, v in arrs.items():
                        if k.endswith(suffix):
                            return v
                    raise KeyError(suffix)

                hf_plus = ak.to_numpy(get_leaf("etHFhitSumPlus_"))
                hf_minus = ak.to_numpy(get_leaf("etHFhitSumMinus_"))
                tw_plus = ak.to_numpy(get_leaf("etHFtowerSumPlus_"))
                tw_minus = ak.to_numpy(get_leaf("etHFtowerSumMinus_"))
                pix = ak.to_numpy(get_leaf("pixelMultiplicity_"))
                maxzeros = ak.to_numpy(get_leaf("maxzeros_"))
                nprob = ak.to_numpy(get_leaf("nproblemRBXs_"))

                vz = get_leaf("position_.fCoordinates.fZ")
                ndof = get_leaf("ndof_")
                valid = get_leaf("validity_")
                vtx_ok = ak.to_numpy(ak.any(valid & (ndof > 0.0) & (np.abs(vz) <= 25.0), axis=1))

                j_pt_raw = get_leaf("pt_")
                j_eta = get_leaf("eta_")
                j_phi = get_leaf("phi_")
                j_mass_raw = get_leaf("mass_")
                j_emf = get_leaf("m_specific.mEnergyFractionEm")

                j_pt_corr, j_mass_corr = apply_jec_jagged(j_pt_raw, j_eta, j_mass_raw, l2_arr)

                # Filter jets to |eta| < 2.0 and sort descending by corrected pT
                in_eta = np.abs(j_eta) < 2.0
                jets_all = ak.zip(
                    {
                        "pt": j_pt_corr[in_eta],
                        "eta": j_eta[in_eta],
                        "phi": j_phi[in_eta],
                        "mass": j_mass_corr[in_eta],
                        "pt_raw": j_pt_raw[in_eta],
                        "emf": j_emf[in_eta],
                    },
                    with_name="Momentum4D",
                )
                ord_idx = ak.argsort(jets_all.pt, axis=1, ascending=False)
                jets = jets_all[ord_idx]

                padded = ak.pad_none(jets, 2, axis=1)
                j1 = padded[:, 0]
                j2 = padded[:, 1]

                pt1 = ak.to_numpy(ak.fill_none(j1.pt, -1.0)).astype(np.float64)
                eta1 = ak.to_numpy(ak.fill_none(j1.eta, 0.0)).astype(np.float64)
                phi1 = ak.to_numpy(ak.fill_none(j1.phi, 0.0)).astype(np.float64)
                pt1_raw = ak.to_numpy(ak.fill_none(j1.pt_raw, -1.0)).astype(np.float64)
                emf1 = ak.to_numpy(ak.fill_none(j1.emf, 1.0)).astype(np.float64)

                pt2 = ak.to_numpy(ak.fill_none(j2.pt, -1.0)).astype(np.float64)
                eta2 = ak.to_numpy(ak.fill_none(j2.eta, 0.0)).astype(np.float64)
                phi2 = ak.to_numpy(ak.fill_none(j2.phi, 0.0)).astype(np.float64)
                pt2_raw = ak.to_numpy(ak.fill_none(j2.pt_raw, -1.0)).astype(np.float64)

                # Compute wrapped |deltaphi| via vector Momentum4D
                dphi12 = ak.to_numpy(ak.fill_none(np.abs(j1.deltaphi(j2)), -1.0)).astype(np.float64)
                n_pt35 = ak.to_numpy(ak.num(jets[jets.pt > 35.0], axis=1))

                hf_tot = hf_plus + hf_minus
                hf_coinc_ok = (tw_plus > 3.0) & (tw_minus > 3.0)
                spike_bad = (pt1 > 50.0) & (emf1 < 0.01)
                scraping_bad = (pix < 50.0) & (hf_tot > 5000.0)

                # Exact algebraic lossless BasicCluster pre-filter for HiEcalRecHitSpikeFilter (minEt = 20 GeV):
                # HybridClusterAlgo seeds on all unexcluded EB crystals (E_T > 1 GeV) in descending E_T order and
                # sums non-negative dominoes, so the highest-E_T EB crystal s* always forms a cluster with
                # E_BC / cosh(eta(seedId)) >= rhEt(s*).
                bc_e = get_leaf("hybridBarrelBasicClusters_RECO.obj.energy_")
                bc_id = get_leaf("hybridBarrelBasicClusters_RECO.obj.seedId_.id_")
                bc_ieta_abs = (bc_id >> 9) & 0x7F
                bc_et = bc_e / np.cosh((bc_ieta_abs - 0.5) * 0.017453292519943295)
                max_bc_et = ak.to_numpy(ak.fill_none(ak.max(bc_et, axis=1), 0.0))

                s_gold = gold
                s_hlt = s_gold & pass_j50
                s_hf = s_hlt & hf_coinc_ok
                s_vtx = s_hf & vtx_ok
                j35_hi_pre = s_gold & pass_j35 & hf_coinc_ok & vtx_ok & (pt1 >= 80.0)
                cand_idx = np.flatnonzero((s_vtx | j35_hi_pre) & (max_bc_et > 19.5))

                ecal_ok = np.ones(n_total, dtype=bool)
                if len(cand_idx) > 0:
                    eb_br = ev[ECAL_EB_BRANCH]
                    raw_map = read_selected_baskets(eb_br, cand_idx)
                    for eidx in cand_idx:
                        ids, ene, tim = parse_eb_rechits_binary(raw_map[int(eidx)])
                        ecal_ok[eidx] = hi_ecal_spike_filter_exact(ids, ene, tim)

                # Note: in CMSSW_3_9_2_patch5, HcalNoiseSummary nproblemRBXs_ counts HBHE RBXs above pp
                # occupancy (36-38 out of 72 RBXs in every central PbPb collision); HI hcalNoiseFilter cuts on
                # maxzeros_ < 100 and sanity-checks 0 <= nproblemRBXs_ <= 72.
                hcal_scraping_ok = (maxzeros < 100) & (nprob >= 0) & (nprob <= 72) & (~scraping_bad) & (~spike_bad)
                noise_ok = ecal_ok & hcal_scraping_ok

                # Sequential cutflow for HLT_HIJet50U
                s_ecal = s_vtx & ecal_ok
                s_clean = s_vtx & noise_ok
                s_2j35 = s_clean & (n_pt35 >= 2)
                s_pt120 = s_2j35 & (pt1 > 120.0)
                s_pt50 = s_pt120 & (pt2 > 50.0)
                s_dphi = s_pt50 & (dphi12 > (2.0 * math.pi / 3.0))

                cf = {
                    "total_events_read": int(n_total),
                    "golden_json_certified": int(s_gold.sum()),
                    "hlt_hijet50u": int(s_hlt.sum()),
                    "hf_coincidence": int(s_hf.sum()),
                    "primary_vertex": int(s_vtx.sum()),
                    "ecal_spike_cleaned": int(s_ecal.sum()),
                    "noise_and_scraping_cleaned": int(s_clean.sum()),
                    "two_jets_pt35_eta2": int(s_2j35.sum()),
                    "leading_jet_pt120": int(s_pt120.sum()),
                    "subleading_jet_pt50": int(s_pt50.sum()),
                    "dijet_dphi_gt_2pi3": int(s_dphi.sum()),
                    "nprob_zero_among_clean": int((s_clean & (nprob == 0)).sum()),
                }

                # Also extract cleaned HLT_HIJet35U events for high-pT turn-on reference
                clean_any = s_gold & hf_coinc_ok & vtx_ok & noise_ok
                j35_mask = clean_any & pass_j35 & (pt1 > 0.0)

                # Keep all cleaned HLT_HIJet50U events for HF centrality + R_B(A_J<0.15) + dijet spectra
                return {
                    "runs": np.unique(runs[s_gold]).tolist(),
                    "cutflow": cf,
                    "j50_hf_gev": hf_tot[s_clean],
                    "j50_pt1": pt1[s_clean],
                    "j50_eta1": eta1[s_clean],
                    "j50_phi1": phi1[s_clean],
                    "j50_pt1_raw": pt1_raw[s_clean],
                    "j50_pt2": pt2[s_clean],
                    "j50_eta2": eta2[s_clean],
                    "j50_phi2": phi2[s_clean],
                    "j50_pt2_raw": pt2_raw[s_clean],
                    "j50_dphi12": dphi12[s_clean],
                    "j35_pt1": pt1[j35_mask],
                    "j35_pass_j50": pass_j50[j35_mask],
                }
        except Exception as exc:
            log.warning("Retry %d/%d on Core file %s: %s", attempt + 1, retries, url.split("/")[-1], exc)
            time.sleep(1.5 * (attempt + 1))
    log.error("Failed Core file after %d retries: %s", retries, url)
    return None


def _stream_core_shard(args: tuple) -> int:
    """Worker process executing multi-threaded HTTPS range reads on its own GIL."""
    global _TLS_SESSION, _PSET_LOCK
    _TLS_SESSION = threading.local()
    _PSET_LOCK = threading.Lock()
    shard_id, urls, per_file_dir_str, golden, l2_arr, threads, deadline = args
    pdir = Path(per_file_dir_str)
    saved = 0
    t0 = time.perf_counter()

    def run_one(raw_url: str) -> int:
        key = hashlib.sha1(raw_url.encode("utf-8")).hexdigest()[:20]
        dst = pdir / f"{CACHE_CORE_PREFIX}_{key}.pkl"
        if dst.is_file():
            return 0
        if deadline > 0.0 and time.perf_counter() > deadline:
            return 0
        res = process_core_file(raw_url, golden, l2_arr)
        if res is not None:
            tmp = pdir / f".{CACHE_CORE_PREFIX}_{key}.{os.getpid()}.{threading.get_ident()}.tmp"
            tmp.write_bytes(pickle.dumps(res, protocol=pickle.HIGHEST_PROTOCOL))
            tmp.replace(dst)
            return 1
        return 0

    with ThreadPoolExecutor(max_workers=max(1, threads)) as pool:
        futs = [pool.submit(run_one, u) for u in urls]
        for idx, fut in enumerate(as_completed(futs), 1):
            saved += fut.result()
            if saved > 0 and idx % 20 == 0:
                dt = time.perf_counter() - t0
                log.info("[Shard %d] %d/%d checked, new_saved=%d (%.1fs, %.2f new/s)", shard_id, idx, len(urls), saved, dt, saved / max(dt, 1e-3))
    return saved


def process_mb_file(raw_url: str, golden: dict, l2_arr: np.ndarray, retries: int = 4) -> Optional[dict[str, Any]]:
    """Process one HIAllPhysics ZS-v2 file; returns gold_events=0 if outside Golden JSON."""
    url = raw_url.replace("root://eospublic.cern.ch/", "https://opendata.cern.ch/")
    for attempt in range(retries):
        try:
            with uproot.open(url, timeout=60, handler=PooledHTTPSource) as f:
                ev = f["Events"]
                runs, lumis, _ = decode_event_aux(ev)
                if len(runs) == 0:
                    return {"gold_events": 0}
                gold = in_golden_mask(runs, lumis, golden)
                if not gold.any():
                    return {"gold_events": 0}

                hlt_paths = decode_hlt_trigger_results(f, ev)
                pass_mb = hlt_paths.get("HLT_HIMinBiasHfOrBSC", np.zeros(len(runs), dtype=bool))
                if not (gold & pass_mb).any():
                    return {"gold_events": 0}
                pass_j50 = hlt_paths.get("HLT_HIJet50U", np.zeros(len(runs), dtype=bool))
                pass_j35 = hlt_paths.get("HLT_HIJet35U", np.zeros(len(runs), dtype=bool))

                arrs = ev.arrays(
                    filter_name=[
                        "recoCentrality_hiCentrality__RECO./recoCentrality_hiCentrality__RECO.obj",
                        "HcalNoiseSummary_hcalnoise__RECO./HcalNoiseSummary_hcalnoise__RECO.obj",
                    ],
                    how=dict,
                )
                cobj = next(v for k, v in arrs.items() if "hiCentrality" in k)
                nobj = next(v for k, v in arrs.items() if "hcalnoise" in k)

                hf_tot = ak.to_numpy(cobj["etHFhitSumPlus_"] + cobj["etHFhitSumMinus_"])
                tw_plus = ak.to_numpy(cobj["etHFtowerSumPlus_"])
                tw_minus = ak.to_numpy(cobj["etHFtowerSumMinus_"])
                pix = ak.to_numpy(cobj["pixelMultiplicity_"])
                maxzeros = ak.to_numpy(nobj["maxzeros_"])
                nprob = ak.to_numpy(nobj["nproblemRBXs_"])

                vtx_ok, j_pt_raw, j_eta, j_phi, j_mass_raw, j_emf = decode_unsplit_vertices_and_jets(ev)
                j_pt_corr, _ = apply_jec_jagged(j_pt_raw, j_eta, j_mass_raw, l2_arr)

                in_eta = np.abs(j_eta) < 2.0
                j_pt_in = j_pt_corr[in_eta]
                j_emf_in = j_emf[in_eta]
                ord_idx = ak.argsort(j_pt_in, axis=1, ascending=False)
                pt1 = ak.to_numpy(ak.fill_none(ak.firsts(j_pt_in[ord_idx]), -1.0)).astype(np.float64)
                emf1 = ak.to_numpy(ak.fill_none(ak.firsts(j_emf_in[ord_idx]), 1.0)).astype(np.float64)

                spike_bad = (pt1 > 50.0) & (emf1 < 0.01)
                scraping_bad = (pix < 50.0) & (hf_tot > 5000.0)
                pre_mb = (
                    gold
                    & pass_mb
                    & (tw_plus > 3.0)
                    & (tw_minus > 3.0)
                    & vtx_ok
                    & (maxzeros < 100)
                    & (nprob >= 0)
                    & (nprob <= 72)
                    & (~scraping_bad)
                    & (~spike_bad)
                )
                max_sc_et = decode_unsplit_max_sc_et(ev)
                cand_mb_idx = np.flatnonzero(pre_mb & (max_sc_et > 18.5))
                ecal_mb_ok = np.ones(len(runs), dtype=bool)
                if len(cand_mb_idx) > 0:
                    eb_br = ev[ECAL_EB_BRANCH]
                    raw_map = read_selected_baskets(eb_br, cand_mb_idx)
                    for eidx in cand_mb_idx:
                        ids, ene, tim = parse_eb_rechits_binary(raw_map[int(eidx)])
                        ecal_mb_ok[eidx] = hi_ecal_spike_filter_exact(ids, ene, tim)
                clean_mb = pre_mb & ecal_mb_ok
                return {
                    "gold_events": int(gold.sum()),
                    "runs": np.unique(runs[gold]).tolist(),
                    "mb_hf_gev": hf_tot[clean_mb],
                    "mb_pt1": pt1[clean_mb],
                    "mb_pass_j50": pass_j50[clean_mb],
                    "mb_pass_j35": pass_j35[clean_mb],
                }
        except Exception as exc:
            log.warning("Retry %d/%d on MB file %s: %s", attempt + 1, retries, url.split("/")[-1], exc)
            time.sleep(1.5 * (attempt + 1))
    log.error("Failed MB file after %d retries: %s", retries, url)
    return None


def build_summary_and_histograms(
    cache_arrays: dict[str, Any],
    cutflow_totals: dict[str, int],
    meta: dict[str, Any],
) -> dict[str, Any]:
    """Compute all histograms, per-centrality observables, and systematics using hist.Hist."""
    j50_hf_gev = cache_arrays["j50_hf_gev"]
    j50_pt1 = cache_arrays["j50_pt1"]
    j50_eta1 = cache_arrays["j50_eta1"]
    j50_pt2 = cache_arrays["j50_pt2"]
    j50_eta2 = cache_arrays["j50_eta2"]
    j50_dphi12 = cache_arrays["j50_dphi12"]

    mb_hf_gev = cache_arrays["mb_hf_gev"]
    mb_pt1 = cache_arrays["mb_pt1"]
    mb_pass_j50 = cache_arrays["mb_pass_j50"]
    mb_pass_j35 = cache_arrays["mb_pass_j35"]

    j35_pt1 = cache_arrays["j35_pt1"]
    j35_pass_j50 = cache_arrays["j35_pass_j50"]

    # 1. HF Energy Centrality Distribution (Figure 3a): 70 uniform bins [0, 160] TeV
    h_hf = (
        Hist.new.Reg(70, 0.0, 160.0, name="ehf", label=r"$\sum E_{\mathrm{HF}}$ [TeV]")
        .StrCat(["MinBias", "HLT_HIJet50U"], name="trig")
        .Weight()
    )
    h_hf.fill(ehf=mb_hf_gev / 1000.0, trig="MinBias")
    h_hf.fill(ehf=j50_hf_gev / 1000.0, trig="HLT_HIJet50U")

    mb_hf_vals = h_hf[:, "MinBias"].values()
    j50_hf_vals = h_hf[:, "HLT_HIJet50U"].values()
    mb_hf_sum = max(float(mb_hf_vals.sum()), 1.0)
    j50_hf_sum = max(float(j50_hf_vals.sum()), 1.0)
    hf_edges_tev = h_hf.axes[0].edges.tolist()

    # Check MB 5-centrality fractions
    mb_cidx = classify_centrality_idx(mb_hf_gev)
    mb_cent_fracs = {
        CENTRALITY_NAMES[i]: float((mb_cidx == i).mean()) if len(mb_cidx) > 0 else 0.0 for i in range(5)
    }

    # 2. Single-Jet Trigger Turn-On Curve (Figure 1a): 24 uniform bins [20, 260] GeV
    turnon_edges = np.linspace(20.0, 260.0, 25)
    mb_has_jet = mb_pt1 > 0.0
    den_mb, _ = np.histogram(mb_pt1[mb_has_jet], bins=turnon_edges)
    num_mb, _ = np.histogram(mb_pt1[mb_has_jet & mb_pass_j50], bins=turnon_edges)

    # Check HIJet35U efficiency in MB events to verify plateau at pT1 >= 85-90 GeV
    num_mb_j35, _ = np.histogram(mb_pt1[mb_has_jet & mb_pass_j35], bins=turnon_edges)

    # Orthogonal lower-threshold jet reference HLT_HIJet35U on HICorePhysics above 90 GeV
    den_j35, _ = np.histogram(j35_pt1, bins=turnon_edges)
    num_j35, _ = np.histogram(j35_pt1[j35_pass_j50], bins=turnon_edges)

    # Combined unbiased turn-on: pure MinBias below 90 GeV, HLT_HIJet35U (100% plateau) + MinBias >= 90 GeV
    # Note: to avoid double counting any identical event, HIAllPhysics and HICorePhysics are distinct datasets,
    # and we use MinBias for bins < 90 GeV and HLT_HIJet35U (from HICorePhysics) for bins >= 90 GeV!
    den_comb = np.where(turnon_edges[:-1] < 90.0, den_mb, den_j35 + den_mb)
    num_comb = np.where(turnon_edges[:-1] < 90.0, num_mb, num_j35 + num_mb)
    # Wait: could an event exist in BOTH HIAllPhysics MB and HICorePhysics HIJet35U if same run/lumi?
    # Using pure MinBias for pT1 < 90 GeV and pure HLT_HIJet35U for pT1 >= 90 GeV guarantees strictly disjoint
    # independent binomial trials in every single bin! Let's check:
    den_stitch = np.where(turnon_edges[:-1] < 90.0, den_mb, np.where(den_j35 > 0, den_j35, den_mb))
    num_stitch = np.where(turnon_edges[:-1] < 90.0, num_mb, np.where(den_j35 > 0, num_j35, num_mb))

    def pack_turnon(num_arr: np.ndarray, den_arr: np.ndarray) -> dict[str, Any]:
        eff, lo, hi = clopper_pearson(num_arr, den_arr)
        return {
            "n_num": [int(x) for x in num_arr],
            "n_den": [int(x) for x in den_arr],
            "eff": [float(eff[i]) if den_arr[i] > 0 else 0.0 for i in range(len(den_arr))],
            "lo": [float(lo[i]) if den_arr[i] > 0 else 0.0 for i in range(len(den_arr))],
            "hi": [float(hi[i]) if den_arr[i] > 0 else 1.0 for i in range(len(den_arr))],
            "err_low": [float(eff[i] - lo[i]) if den_arr[i] > 0 else 0.0 for i in range(len(den_arr))],
            "err_high": [float(hi[i] - eff[i]) if den_arr[i] > 0 else 0.0 for i in range(len(den_arr))],
        }

    # Check overlap in [80, 120] GeV between MB and HLT_HIJet35U
    ov_mask_bins = (turnon_edges[:-1] >= 80.0) & (turnon_edges[1:] <= 120.0)
    mb_ov_den = int(den_mb[ov_mask_bins].sum())
    mb_ov_num = int(num_mb[ov_mask_bins].sum())
    j35_ov_den = int(den_j35[ov_mask_bins].sum())
    j35_ov_num = int(num_j35[ov_mask_bins].sum())

    trig_summary = {
        "pt_edges": turnon_edges.tolist(),
        **pack_turnon(num_stitch, den_stitch),
        "stitch_gev": 90.0,
        "mb_only": pack_turnon(num_mb, den_mb),
        "jet35_ref": pack_turnon(num_j35, den_j35),
        "mb_jet35_turnon": pack_turnon(num_mb_j35, den_mb),
        "overlap_80_120_gev": {
            "mb_num": mb_ov_num,
            "mb_den": mb_ov_den,
            "mb_eff": float(mb_ov_num / mb_ov_den) if mb_ov_den > 0 else None,
            "jet35_num": j35_ov_num,
            "jet35_den": j35_ov_den,
            "jet35_eff": float(j35_ov_num / j35_ov_den) if j35_ov_den > 0 else None,
        },
    }

    # 3-9. Centrality-resolved dijet observables with hist.Hist
    cidx = classify_centrality_idx(j50_hf_gev)
    cent_labels = np.array(CENTRALITY_NAMES, dtype=object)[cidx]
    merged_labels = np.array([MERGED_CENTRALITY_MAP[c] for c in cent_labels], dtype=object)

    sel_lead120 = j50_pt1 > 120.0
    sel_dphi_sample = sel_lead120 & (j50_pt2 > 50.0)
    sel_dijet = sel_dphi_sample & (j50_dphi12 > (2.0 * math.pi / 3.0))

    aj_all = np.where(sel_dijet, (j50_pt1 - j50_pt2) / np.maximum(j50_pt1 + j50_pt2, 1e-9), -1.0)
    dpt_rel_all = np.where(sel_dijet, (j50_pt1 - j50_pt2) / np.maximum(j50_pt1, 1e-9), -1.0)

    # Book hist.Hist objects
    h_pt2_vs_pt1 = (
        Hist.new.Reg(26, 120.0, 250.0, name="pt1", label=r"Leading jet $p_{T,1}$ [GeV]")
        .Reg(40, 50.0, 250.0, name="pt2", label=r"Subleading jet $p_{T,2}$ [GeV]")
        .StrCat(MERGED_CENTRALITY_NAMES, name="cent")
        .Weight()
    )
    h_pt2_vs_pt1.fill(
        pt1=j50_pt1[sel_dijet],
        pt2=j50_pt2[sel_dijet],
        cent=merged_labels[sel_dijet],
    )

    h_pt1 = (
        Hist.new.Reg(9, 120.0, 300.0, name="pt1", label=r"Leading jet $p_{T,1}$ [GeV]")
        .StrCat(CENTRALITY_NAMES, name="cent")
        .Weight()
    )
    h_pt1.fill(pt1=j50_pt1[sel_dijet], cent=cent_labels[sel_dijet])

    h_dphi = (
        Hist.new.Reg(20, 0.0, math.pi, name="dphi", label=r"$\Delta\phi_{12}$ [rad]")
        .StrCat(CENTRALITY_NAMES, name="cent")
        .Weight()
    )
    h_dphi.fill(dphi=j50_dphi12[sel_dphi_sample], cent=cent_labels[sel_dphi_sample])

    # 13 uniform bins [0.0, 0.91] (bin width 0.07, exact CMS-HIN-10-004 Figure 8 binning; zero overflow)
    h_aj = (
        Hist.new.Reg(13, 0.0, 0.91, name="aj", label=r"$A_{J} = (p_{T,1} - p_{T,2}) / (p_{T,1} + p_{T,2})$")
        .StrCat(CENTRALITY_NAMES, name="cent")
        .Weight()
    )
    h_aj.fill(aj=aj_all[sel_dijet], cent=cent_labels[sel_dijet])

    per_centrality: dict[str, Any] = {}
    dphi_bw = math.pi / 20.0
    pt1_bw = (300.0 - 120.0) / 9.0

    for i, (cname, ehf_lo, ehf_hi, npart, npart_err) in enumerate(CENTRALITY_BINS):
        in_c = cidx == i
        m_lead = in_c & sel_lead120
        m_dphi = in_c & sel_dphi_sample
        m_dij = in_c & sel_dijet

        n_lead = int(m_lead.sum())
        n_dphi = int(m_dphi.sum())
        n_dij = int(m_dij.sum())

        # A_J mean & SEM
        aj_c = aj_all[m_dij]
        aj_mean = float(aj_c.mean()) if n_dij > 0 else 0.0
        aj_sem = float(aj_c.std(ddof=1) / math.sqrt(n_dij)) if n_dij > 1 else 0.0

        # R_B(A_J < 0.15): denominator = all cleaned HLT_HIJet50U with pT1 > 120 GeV (including mono-jets!)
        k_rb_aj = int((m_dij & (aj_all < 0.15)).sum())
        rb_aj, rb_aj_lo, rb_aj_hi = clopper_pearson(k_rb_aj, n_lead)
        rb_aj_val = float(rb_aj) if n_lead > 0 else 0.0
        rb_aj_elo = float(rb_aj - rb_aj_lo) if n_lead > 0 else 0.0
        rb_aj_ehi = float(rb_aj_hi - rb_aj) if n_lead > 0 else 0.0

        # R_B(dphi12 > 3.026): denominator = pT1 > 120, pT2 > 50 (no dphi pre-cut)
        k_rb_dp = int((m_dphi & (j50_dphi12 > 3.026)).sum())
        rb_dp, rb_dp_lo, rb_dp_hi = clopper_pearson(k_rb_dp, n_dphi)
        rb_dp_val = float(rb_dp) if n_dphi > 0 else 0.0
        rb_dp_elo = float(rb_dp - rb_dp_lo) if n_dphi > 0 else 0.0
        rb_dp_ehi = float(rb_dp_hi - rb_dp) if n_dphi > 0 else 0.0

        aj_counts = h_aj[:, cname].values()
        aj_tot = max(float(aj_counts.sum()), 1.0)
        aj_frac = (aj_counts / aj_tot).tolist()
        aj_ferr = (np.sqrt(aj_counts) / aj_tot).tolist()

        dp_counts = h_dphi[:, cname].values()
        dp_tot = max(float(dp_counts.sum()), 1.0)
        dp_dens = (dp_counts / (dp_tot * dphi_bw)).tolist()
        dp_derr = (np.sqrt(dp_counts) / (dp_tot * dphi_bw)).tolist()

        p1_counts = h_pt1[:, cname].values()
        p1_tot = max(float(p1_counts.sum()), 1.0)
        p1_dens = (p1_counts / (p1_tot * pt1_bw)).tolist()
        p1_derr = (np.sqrt(p1_counts) / (p1_tot * pt1_bw)).tolist()

        per_centrality[cname] = {
            "ehf_range_gev": [ehf_lo, None if math.isinf(ehf_hi) else ehf_hi],
            "npart_mean": npart,
            "npart_err": npart_err,
            "n_cleaned_jet50": int(in_c.sum()),
            "n_lead120": n_lead,
            "n_dphi_sample": n_dphi,
            "n_dijet": n_dij,
            "aj_mean": aj_mean,
            "aj_mean_err": aj_sem,
            "k_rb_aj_lt_015": k_rb_aj,
            "rb_aj_lt_015": rb_aj_val,
            "rb_aj_err_low": rb_aj_elo,
            "rb_aj_err_high": rb_aj_ehi,
            "k_rb_dphi_gt_3026": k_rb_dp,
            "rb_dphi_gt_3026": rb_dp_val,
            "rb_dphi_err_low": rb_dp_elo,
            "rb_dphi_err_high": rb_dp_ehi,
            "aj_hist": {
                "edges": [round(float(x), 6) for x in h_aj.axes[0].edges],
                "counts": [int(x) for x in aj_counts],
                "fraction": [float(x) for x in aj_frac],
                "fraction_err": [float(x) for x in aj_ferr],
            },
            "dphi_hist": {
                "edges": [float(x) for x in h_dphi.axes[0].edges],
                "counts": [int(x) for x in dp_counts],
                "density": [float(x) for x in dp_dens],
                "density_err": [float(x) for x in dp_derr],
            },
            "pt1_hist": {
                "edges": [float(x) for x in h_pt1.axes[0].edges],
                "counts": [int(x) for x in p1_counts],
                "density": [float(x) for x in p1_dens],
                "density_err": [float(x) for x in p1_derr],
                "fraction": [float(x / p1_tot) for x in p1_counts],
                "fraction_err": [float(math.sqrt(x) / p1_tot) for x in p1_counts],
            },
        }

    # Figure 10: Mean fractional imbalance <(pT1 - pT2)/pT1> vs pT1 in 3 merged centrality bins
    fi_edges = np.array([120.0, 130.0, 140.0, 170.0, 240.0], dtype=np.float64)
    fi_summary: dict[str, Any] = {"pt1_edges": fi_edges.tolist()}
    for mname in MERGED_CENTRALITY_NAMES:
        in_m = (merged_labels == mname) & sel_dijet
        p1_m = j50_pt1[in_m]
        dpt_m = dpt_rel_all[in_m]
        counts_m, means_m, sems_m = [], [], []
        for bi in range(len(fi_edges) - 1):
            in_b = (p1_m >= fi_edges[bi]) & (p1_m < fi_edges[bi + 1])
            vals = dpt_m[in_b]
            cnt = int(len(vals))
            counts_m.append(cnt)
            means_m.append(float(vals.mean()) if cnt > 0 else 0.0)
            sems_m.append(float(vals.std(ddof=1) / math.sqrt(cnt)) if cnt > 1 else 0.0)
        fi_summary[mname] = {
            "pt1_edges": fi_edges.tolist(),
            "counts": counts_m,
            "mean": means_m,
            "sem": sems_m,
        }

    # Full 6-bin [100, 240] GeV profile matching deltaPtOverPt1VsPt.C down to pT1 = 100 GeV
    fi100_edges = np.array([100.0, 110.0, 120.0, 130.0, 140.0, 170.0, 240.0], dtype=np.float64)
    sel_dijet_100 = (j50_pt1 >= 100.0) & (j50_pt2 > 50.0) & (j50_dphi12 > (2.0 * math.pi / 3.0))
    dpt_rel_100 = np.where(sel_dijet_100, (j50_pt1 - j50_pt2) / np.maximum(j50_pt1, 1e-9), -1.0)
    fi_summary_100: dict[str, Any] = {"pt1_edges": fi100_edges.tolist()}
    for mname in MERGED_CENTRALITY_NAMES:
        in_m = (merged_labels == mname) & sel_dijet_100
        p1_m = j50_pt1[in_m]
        dpt_m = dpt_rel_100[in_m]
        counts_m, means_m, sems_m = [], [], []
        for bi in range(len(fi100_edges) - 1):
            in_b = (p1_m >= fi100_edges[bi]) & (p1_m < fi100_edges[bi + 1])
            vals = dpt_m[in_b]
            cnt = int(len(vals))
            counts_m.append(cnt)
            means_m.append(float(vals.mean()) if cnt > 0 else 0.0)
            sems_m.append(float(vals.std(ddof=1) / math.sqrt(cnt)) if cnt > 1 else 0.0)
        fi_summary_100[mname] = {
            "pt1_edges": fi100_edges.tolist(),
            "counts": counts_m,
            "mean": means_m,
            "sem": sems_m,
        }

    # Figure 4: 2D subleading vs leading jet pT correlation in 3 merged centrality bins
    pt2_vs_pt1_summary: dict[str, Any] = {
        "pt1_edges": h_pt2_vs_pt1.axes[0].edges.tolist(),
        "pt2_edges": h_pt2_vs_pt1.axes[1].edges.tolist(),
        "per_merged_centrality": {
            mname: h_pt2_vs_pt1[:, :, mname].values().astype(int).tolist() for mname in MERGED_CENTRALITY_NAMES
        },
    }

    # Systematics: JEC +/- 5% scale shift & centrality-dependent residual JEC tilt
    # Tilt model (CMS-HIN-10-004 Sec. II.D.3 Fig. 2): -5% in 50-100% up to +5%*(100/pT) at low pT in 0-10%
    syst_summary: dict[str, Any] = {}
    cent_scale_offset = np.array([-0.05, -0.025, 0.0, +0.025, +0.05], dtype=np.float64)[cidx]
    for var_name, s1, s2 in [
        ("jec_up_5pct", np.full_like(j50_pt1, 1.05), np.full_like(j50_pt2, 1.05)),
        ("jec_down_5pct", np.full_like(j50_pt1, 0.95), np.full_like(j50_pt2, 0.95)),
        (
            "jec_centrality_pt_tilt",
            1.0 + cent_scale_offset * np.clip(100.0 / np.maximum(j50_pt1, 30.0), 0.3, 1.5),
            1.0 + cent_scale_offset * np.clip(100.0 / np.maximum(j50_pt2, 30.0), 0.3, 1.5),
        ),
    ]:
        p1_v = np.where(j50_pt1 > 0, j50_pt1 * s1, -1.0)
        p2_v = np.where(j50_pt2 > 0, j50_pt2 * s2, -1.0)
        # Re-order if p2_v > p1_v under non-uniform tilt
        swap = (p2_v > p1_v) & (p2_v > 0)
        p1_s = np.where(swap, p2_v, p1_v)
        p2_s = np.where(swap, p1_v, p2_v)
        m_l120 = p1_s > 120.0
        m_dij = m_l120 & (p2_s > 50.0) & (j50_dphi12 > (2.0 * math.pi / 3.0))
        aj_v = np.where(m_dij, (p1_s - p2_s) / np.maximum(p1_s + p2_s, 1e-9), -1.0)
        var_pc = {}
        for i, cname in enumerate(CENTRALITY_NAMES):
            ic = cidx == i
            nl = int((ic & m_l120).sum())
            nd = int((ic & m_dij).sum())
            aj_sub = aj_v[ic & m_dij]
            k_rb = int((ic & m_dij & (aj_v < 0.15)).sum())
            var_pc[cname] = {
                "n_dijet": nd,
                "aj_mean": float(aj_sub.mean()) if nd > 0 else 0.0,
                "rb_aj_lt_015": float(k_rb / nl) if nl > 0 else 0.0,
            }
        syst_summary[var_name] = var_pc

    # Attach total systematic uncertainty envelope per centrality onto per_centrality
    for cname in CENTRALITY_NAMES:
        nom_aj = per_centrality[cname]["aj_mean"]
        nom_rb = per_centrality[cname]["rb_aj_lt_015"]
        aj_shifts = [abs(syst_summary[vn][cname]["aj_mean"] - nom_aj) for vn in syst_summary]
        rb_shifts = [abs(syst_summary[vn][cname]["rb_aj_lt_015"] - nom_rb) for vn in syst_summary]
        per_centrality[cname]["aj_mean_syst"] = float(max(aj_shifts)) if aj_shifts else 0.0
        per_centrality[cname]["rb_aj_syst"] = float(max(rb_shifts)) if rb_shifts else 0.0

    # Key physics headline checks
    aj_per = per_centrality["50-100%"]["aj_mean"]
    aj_per_e = per_centrality["50-100%"]["aj_mean_err"]
    aj_cen = per_centrality["0-10%"]["aj_mean"]
    aj_cen_e = per_centrality["0-10%"]["aj_mean_err"]
    aj_sig = (
        (aj_cen - aj_per) / math.hypot(aj_per_e, aj_cen_e)
        if math.hypot(aj_per_e, aj_cen_e) > 0
        else 0.0
    )

    rb_per = per_centrality["50-100%"]["rb_aj_lt_015"]
    rb_per_e = 0.5 * (per_centrality["50-100%"]["rb_aj_err_low"] + per_centrality["50-100%"]["rb_aj_err_high"])
    rb_cen = per_centrality["0-10%"]["rb_aj_lt_015"]
    rb_cen_e = 0.5 * (per_centrality["0-10%"]["rb_aj_err_low"] + per_centrality["0-10%"]["rb_aj_err_high"])
    rb_sig = (
        (rb_per - rb_cen) / math.hypot(rb_per_e, rb_cen_e)
        if math.hypot(rb_per_e, rb_cen_e) > 0
        else 0.0
    )

    centers_tev = 0.5 * (np.array(hf_edges_tev[:-1]) + np.array(hf_edges_tev[1:]))
    mb_frac = (mb_hf_vals / mb_hf_sum).tolist()
    j50_frac = (j50_hf_vals / j50_hf_sum).tolist()
    ehf_mb_mean = float(np.sum(centers_tev * np.array(mb_frac)))
    ehf_j50_mean = float(np.sum(centers_tev * np.array(j50_frac)))

    n_vtx_tot = int(cutflow_totals.get("primary_vertex", 0))
    n_ecal_tot = int(cutflow_totals.get("ecal_spike_cleaned", cutflow_totals.get("noise_and_scraping_cleaned", 0)))
    n_clean_tot = int(cutflow_totals.get("noise_and_scraping_cleaned", 0))

    return {
        "metadata": meta,
        "cutflow": {
            k: v for k, v in cutflow_totals.items() if k not in ("nprob_zero_among_clean", "ecal_spike_cleaned")
        },
        "validation_checks": {
            "mb_ehf_mean_tev": ehf_mb_mean,
            "jet50_ehf_mean_tev": ehf_j50_mean,
            "mb_centrality_fractions": mb_cent_fracs,
            "aj_shift_50_100_to_0_10": float(aj_cen - aj_per),
            "aj_shift_significance_sigma": float(aj_sig),
            "rb_aj_drop_50_100_to_0_10": float(rb_per - rb_cen),
            "rb_aj_drop_significance_sigma": float(rb_sig),
            "ecal_cleaning_diagnostics": {
                "filter": "HiEcalRecHitSpikeFilter (EcalRecHitsEB, minEt=20.0 GeV, swissThreshold=0.95, timeThreshold=4.0 ns, avoidIeta85=False)",
                "primary_vertex_events": n_vtx_tot,
                "after_ecal_spike_filter_events": n_ecal_tot,
                "ecal_spike_filter_pass_fraction": float(n_ecal_tot / n_vtx_tot) if n_vtx_tot > 0 else 1.0,
                "after_hcal_and_scraping_events": n_clean_tot,
                "hcal_and_scraping_pass_fraction_given_ecal": float(n_clean_tot / n_ecal_tot) if n_ecal_tot > 0 else 1.0,
                "combined_cleaning_pass_fraction": float(n_clean_tot / n_vtx_tot) if n_vtx_tot > 0 else 1.0,
            },
            "nprob_rbx_note": (
                f"In CMSSW_3_9_2_patch5 RECO, HcalNoiseSummary.nproblemRBXs_ uses pp occupancy thresholds "
                f"and flags 36-38 HB RBXs in central PbPb collisions ({cutflow_totals['nprob_zero_among_clean']} / "
                f"{cutflow_totals['noise_and_scraping_cleaned']} cleaned HLT_HIJet50U events have nproblemRBXs_ == 0, "
                f"all ultra-peripheral). Following 2010 HI EventSelection_cff.py, offline cleaning combines "
                f"HiEcalRecHitSpikeFilter on EcalRecHitsEB (minEt=20 GeV, swissThreshold=0.95, avoidIeta85=False, "
                f"|time|<=4 ns) with HCAL maxzeros_ < 100, pixel-HF scraping rejection, and jet EMF >= 0.01."
            ),
        },
        "trigger_turnon": trig_summary,
        "hf_energy_hist": {
            "edges_tev": hf_edges_tev,
            "mb_counts": [int(x) for x in mb_hf_vals],
            "jet50_counts": [int(x) for x in j50_hf_vals],
            "mb_fraction": [float(x) for x in mb_frac],
            "jet50_fraction": [float(x) for x in j50_frac],
        },
        "per_centrality": per_centrality,
        "fractional_imbalance_vs_pt1": fi_summary,
        "fractional_imbalance_vs_pt1_100_240": fi_summary_100,
        "pt2_vs_pt1_2d": pt2_vs_pt1_summary,
        "systematics": syst_summary,
    }


@app.command()
def main(
    max_files: int = typer.Option(560, "--max-files", help="Stratified Golden-JSON HICorePhysics files to stream"),
    max_mb_files: Optional[int] = typer.Option(
        None, "--max-mb-files", help="Validated HIAllPhysics Minimum-Bias files (default: min(max_files, 25))"
    ),
    workers: int = typer.Option(10, "--workers", help="Concurrent HTTPS worker threads"),
    max_runtime_sec: float = typer.Option(0.0, "--max-runtime-sec", help="Optional wall-clock budget per wave in seconds (0=unlimited)"),
    cache_dir: Path = typer.Option(BASE_DIR / "cache", "--cache-dir", help="Directory for compact event cache"),
    from_cache: bool = typer.Option(False, "--from-cache", help="Recompute summary.json from existing cache"),
    out_summary: Path = typer.Option(BASE_DIR / "summary.json", "--out-summary", help="Output summary JSON path"),
    verbose: int = typer.Option(0, "-v", "--verbose", count=True, help="Increase logging verbosity"),
) -> None:
    """Run the CMS-HIN-10-004 2.76 TeV PbPb dijet asymmetry Open Data analysis."""
    lvl = logging.INFO if verbose == 0 else logging.DEBUG
    logging.basicConfig(level=lvl, format="%(asctime)s %(levelname)s: %(message)s", datefmt="%H:%M:%S")

    t_start = time.perf_counter()
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_npz = cache_dir / "selected_events.npz"
    cache_meta_json = cache_dir / "cache_meta.json"

    if from_cache and cache_npz.is_file() and cache_meta_json.is_file():
        log.info("Loading cached skim from %s", cache_npz)
        loaded = np.load(cache_npz)
        cache_arrays = {k: loaded[k] for k in loaded.files}
        cmeta = json.loads(cache_meta_json.read_text(encoding="utf-8"))
        summary = build_summary_and_histograms(cache_arrays, cmeta["cutflow_totals"], cmeta["metadata"])
        out_summary.write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        log.info("Wrote %s from cache in %.2fs", out_summary, time.perf_counter() - t_start)
        return

    golden = json.loads(
        (DATA_DIR / "Cert_150436-152957_HI7TeV_StreamExpress_Collisions10_JSON_v2.txt").read_text(encoding="utf-8")
    )
    l2_arr = load_l2_table(DATA_DIR / "Spring10_HI_L2Relative_IC5Calo.txt")

    # Sanity-check L2L3 JEC at reference points
    jec_test = compute_jec_factor(np.array([50.0, 100.0, 150.0]), np.array([0.05, 0.05, 1.2]), l2_arr)
    log.info(
        "JEC check (Spring10_HI IC5Calo L2L3): scale(50 GeV, eta=0.05)=%.4f, scale(100 GeV, eta=0.05)=%.4f, scale(150 GeV, eta=1.2)=%.4f",
        jec_test[0],
        jec_test[1],
        jec_test[2],
    )

    core_all = [
        ln.strip()
        for ln in (DATA_DIR / "CMS_HIRun2010_HICorePhysics_PromptReco-v3_RECO_file_index.txt")
        .read_text(encoding="utf-8")
        .splitlines()
        if ln.strip()
    ]
    valid_core_pairs = []
    for ln in core_all:
        parts = ln.split("/")
        run_str = str(int(parts[-3] + parts[-2]))
        if run_str in golden:
            valid_core_pairs.append((int(run_str), ln))
    valid_core_pairs.sort()

    per_file_dir = cache_dir / "per_file"
    per_file_dir.mkdir(parents=True, exist_ok=True)

    def url_key(raw_url: str) -> str:
        return hashlib.sha1(raw_url.encode("utf-8")).hexdigest()[:20]

    n_core = min(max_files, len(valid_core_pairs))
    cached_core_pairs = [p for p in valid_core_pairs if (per_file_dir / f"{CACHE_CORE_PREFIX}_{url_key(p[1])}.pkl").is_file()]
    v1_core_pairs = [p for p in valid_core_pairs if (per_file_dir / f"core_v1_{url_key(p[1])}.pkl").is_file()]
    if len(v1_core_pairs) >= n_core:
        stride_idx = np.linspace(0, len(v1_core_pairs) - 1, n_core, dtype=int)
        chosen_pairs = [v1_core_pairs[i] for i in stride_idx]
    elif len(cached_core_pairs) >= n_core:
        stride_idx = np.linspace(0, len(cached_core_pairs) - 1, n_core, dtype=int)
        chosen_pairs = [cached_core_pairs[i] for i in stride_idx]
    else:
        cached_urls_set = {p[1] for p in cached_core_pairs}
        rem_pairs = [p for p in valid_core_pairs if p[1] not in cached_urls_set]
        n_rem = n_core - len(cached_core_pairs)
        rem_idx = np.linspace(0, len(rem_pairs) - 1, n_rem, dtype=int)
        chosen_pairs = sorted(cached_core_pairs + [rem_pairs[i] for i in rem_idx])

    selected_core_urls = [p[1] for p in chosen_pairs]
    selected_core_runs = sorted({p[0] for p in chosen_pairs})
    n_selected_cached = sum(1 for u in selected_core_urls if (per_file_dir / f"{CACHE_CORE_PREFIX}_{url_key(u)}.pkl").is_file())
    log.info(
        "Selected %d / %d Golden-JSON HICorePhysics files across %d runs (%d already cached, workers=%d)",
        len(selected_core_urls),
        len(valid_core_pairs),
        len(selected_core_runs),
        n_selected_cached,
        workers,
    )

    n_mb_target = min(max_files, 25) if max_mb_files is None else max_mb_files
    mb_all = [
        ln.strip()
        for ln in (DATA_DIR / "CMS_HIRun2010_HIAllPhysics_ZS-v2_RECO_file_index.txt")
        .read_text(encoding="utf-8")
        .splitlines()
        if ln.strip()
    ]
    cached_mb_urls = [
        u for u in mb_all if (per_file_dir / f"{CACHE_MB_PREFIX}_{url_key(u)}.pkl").is_file()
    ]
    v1_good_mb_urls: list[str] = []
    seen_v2_mb = set(cached_mb_urls)
    for u in mb_all:
        if u in seen_v2_mb:
            continue
        p_v1 = per_file_dir / f"mb_v1_{url_key(u)}.pkl"
        if p_v1.is_file():
            try:
                d_v1 = pickle.loads(p_v1.read_bytes())
                if d_v1.get("gold_events", 0) > 0 and len(d_v1.get("mb_hf_gev", [])) > 0:
                    v1_good_mb_urls.append(u)
            except Exception:
                pass
    n_mb_cand = min(max(n_mb_target * 3, n_mb_target + 10), len(mb_all))
    mb_cand_idx = np.linspace(0, len(mb_all) - 1, n_mb_cand, dtype=int)
    seen_mb = set(cached_mb_urls) | set(v1_good_mb_urls)
    mb_cand_urls = cached_mb_urls + v1_good_mb_urls + [mb_all[i] for i in mb_cand_idx if mb_all[i] not in seen_mb]

    def cached_file_exec(prefix: str, raw_url: str, fn):
        key = url_key(raw_url)
        dst = per_file_dir / f"{prefix}_{key}.pkl"
        if dst.is_file():
            try:
                return pickle.loads(dst.read_bytes())
            except Exception:
                pass
        if max_runtime_sec > 0.0 and (time.perf_counter() - t_start) > max_runtime_sec:
            return None
        res = fn(raw_url)
        if res is not None:
            tmp = per_file_dir / f".{prefix}_{key}.{os.getpid()}.{threading.get_ident()}.tmp"
            tmp.write_bytes(pickle.dumps(res, protocol=pickle.HIGHEST_PROTOCOL))
            tmp.replace(dst)
        return res

    # Load cached URLs first, then shuffle uncached URLs across runs
    cached_url_lookup = {p[1] for p in cached_core_pairs}
    core_cached_part = [u for u in selected_core_urls if u in cached_url_lookup]
    core_uncached_part = [u for u in selected_core_urls if u not in cached_url_lookup]
    perm = np.random.default_rng(42).permutation(len(core_uncached_part))
    shuffled_core_urls = core_cached_part + [core_uncached_part[i] for i in perm]

    # Pre-warm ParameterSets HLT path cache on 1 file only if uncached files will be streamed and disk cache is empty
    if PSET_CACHE_JSON.is_file() and not _PSET_TRIGGER_PATHS:
        try:
            _PSET_TRIGGER_PATHS.update(json.loads(PSET_CACHE_JSON.read_text(encoding="utf-8")))
        except Exception:
            pass
    has_uncached_core = len(core_uncached_part) > 0
    if core_uncached_part and has_uncached_core and not _PSET_TRIGGER_PATHS:
        try:
            u_warm = core_uncached_part[0].replace("root://eospublic.cern.ch/", "https://opendata.cern.ch/")
            with uproot.open(u_warm, timeout=45, handler=PooledHTTPSource) as fw:
                decode_hlt_trigger_results(fw, fw["Events"])
        except Exception:
            pass

    # Stream HICorePhysics in parallel
    cf_totals = {
        "total_events_read": 0,
        "golden_json_certified": 0,
        "hlt_hijet50u": 0,
        "hf_coincidence": 0,
        "primary_vertex": 0,
        "ecal_spike_cleaned": 0,
        "noise_and_scraping_cleaned": 0,
        "two_jets_pt35_eta2": 0,
        "leading_jet_pt120": 0,
        "subleading_jet_pt50": 0,
        "dijet_dphi_gt_2pi3": 0,
        "nprob_zero_among_clean": 0,
    }
    core_results = []
    runs_observed: set[int] = set()

    t_core0 = time.perf_counter()
    shuffled_uncached = [core_uncached_part[i] for i in perm]
    n_cpus = os.cpu_count() or 1
    if shuffled_uncached and n_cpus > 1 and workers >= 4:
        import multiprocessing as mp

        n_procs = min(n_cpus, max(1, workers // 4), len(shuffled_uncached))
        threads_per_proc = max(2, (workers + n_procs - 1) // n_procs)
        deadline = (t_start + max_runtime_sec) if max_runtime_sec > 0.0 else 0.0
        shard_args = [
            (
                p,
                shuffled_uncached[p::n_procs],
                str(per_file_dir),
                golden,
                l2_arr,
                threads_per_proc,
                deadline,
            )
            for p in range(n_procs)
        ]
        log.info(
            "Streaming %d uncached HICorePhysics files across %d processes x %d threads",
            len(shuffled_uncached),
            n_procs,
            threads_per_proc,
        )
        ctx = mp.get_context("fork")
        with ctx.Pool(processes=n_procs) as mpool:
            new_counts = mpool.map(_stream_core_shard, shard_args)
        dt_stream = time.perf_counter() - t_core0
        log.info("Multi-process wave saved %d new Core caches in %.1fs (%.2f new files/s)", sum(new_counts), dt_stream, sum(new_counts) / max(dt_stream, 1e-3))
        # Aggregate all cached files in selected_core_urls (plus any newly cached files if wave was bounded)
        all_now_cached = [p[1] for p in valid_core_pairs if (per_file_dir / f"{CACHE_CORE_PREFIX}_{url_key(p[1])}.pkl").is_file()]
        agg_urls = [u for u in selected_core_urls if (per_file_dir / f"{CACHE_CORE_PREFIX}_{url_key(u)}.pkl").is_file()]
        if len(agg_urls) < min(n_core, len(all_now_cached)):
            agg_urls = all_now_cached[:n_core]
        for u in agg_urls:
            dst = per_file_dir / f"{CACHE_CORE_PREFIX}_{url_key(u)}.pkl"
            try:
                res = pickle.loads(dst.read_bytes())
            except Exception:
                continue
            core_results.append(res)
            runs_observed.update(res["runs"])
            for k in cf_totals:
                cf_totals[k] += res["cutflow"].get(k, 0)
        log.info(
            "[Core aggregated] ok=%d runs=%d N_jet50=%d N_dijet=%d (%.1fs total)",
            len(core_results),
            len(runs_observed),
            cf_totals["noise_and_scraping_cleaned"],
            cf_totals["dijet_dphi_gt_2pi3"],
            time.perf_counter() - t_core0,
        )
    else:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futs = {
                pool.submit(cached_file_exec, CACHE_CORE_PREFIX, u, lambda url: process_core_file(url, golden, l2_arr)): idx
                for idx, u in enumerate(shuffled_core_urls)
            }
            done_cnt = 0
            for fut in as_completed(futs):
                done_cnt += 1
                res = fut.result()
                if res is not None:
                    core_results.append(res)
                    runs_observed.update(res["runs"])
                    for k in cf_totals:
                        cf_totals[k] += res["cutflow"].get(k, 0)
                if done_cnt % 25 == 0 or done_cnt == len(shuffled_core_urls):
                    elapsed = time.perf_counter() - t_core0
                    log.info(
                        "[Core %d/%d] ok=%d runs=%d N_jet50=%d N_dijet=%d (%.1fs, %.2f files/s)",
                        done_cnt,
                        len(shuffled_core_urls),
                        len(core_results),
                        len(runs_observed),
                        cf_totals["noise_and_scraping_cleaned"],
                        cf_totals["dijet_dphi_gt_2pi3"],
                        elapsed,
                        len(core_results) / max(elapsed, 1e-3),
                    )

    # Stream HIAllPhysics until n_mb_target Golden-JSON files are collected (drain cached files first!)
    mb_results = []
    mb_runs_observed: set[int] = set()
    t_mb0 = time.perf_counter()
    cursor = 0
    while len(mb_results) < n_mb_target and cursor < len(mb_cand_urls):
        need = n_mb_target - len(mb_results)
        if cursor < len(cached_mb_urls):
            batch = mb_cand_urls[cursor : len(cached_mb_urls)]
        elif cursor < len(cached_mb_urls) + len(v1_good_mb_urls):
            batch = mb_cand_urls[cursor : min(len(cached_mb_urls) + len(v1_good_mb_urls), cursor + need)]
        else:
            batch = mb_cand_urls[cursor : cursor + max(need + 4, need)]
        cursor += len(batch)
        before_n = len(mb_results)
        with ThreadPoolExecutor(max_workers=min(workers, max(1, len(batch)))) as pool:
            for res in pool.map(
                lambda u: cached_file_exec(CACHE_MB_PREFIX, u, lambda url: process_mb_file(url, golden, l2_arr)),
                batch,
            ):
                if res is not None and res.get("gold_events", 0) > 0 and len(res["mb_hf_gev"]) > 0:
                    mb_results.append(res)
                    mb_runs_observed.update(res["runs"])
                    if len(mb_results) >= n_mb_target:
                        break
        log.info(
            "[MB %d/%d] cursor=%d/%d runs=%d N_mb_clean=%d (%.1fs)",
            len(mb_results),
            n_mb_target,
            cursor,
            len(mb_cand_urls),
            len(mb_runs_observed),
            sum(len(r["mb_hf_gev"]) for r in mb_results),
            time.perf_counter() - t_mb0,
        )
        if max_runtime_sec > 0.0 and (time.perf_counter() - t_start) > max_runtime_sec and len(mb_results) == before_n:
            break

    def cat(res_list: list[dict], key: str, dtype=np.float64) -> np.ndarray:
        non_empty = [r[key] for r in res_list if len(r[key]) > 0]
        return np.concatenate(non_empty).astype(dtype) if non_empty else np.array([], dtype=dtype)

    # Combine HIAllPhysics HLT_HIMinBiasHfOrBSC with HICorePhysics prescaled HLT_HIMinBiasHfOrBSC_Core if MB files == 0,
    # otherwise use pure HIAllPhysics HLT_HIMinBiasHfOrBSC for HF centrality & combine disjoint MB for turn-on!
    mb_hf_all = cat(mb_results, "mb_hf_gev", np.float64)
    mb_pt1_all = cat(mb_results, "mb_pt1", np.float64)
    mb_j50_all = cat(mb_results, "mb_pass_j50", bool)
    mb_j35_all = cat(mb_results, "mb_pass_j35", bool)

    cache_arrays = {
        "j50_hf_gev": cat(core_results, "j50_hf_gev", np.float64),
        "j50_pt1": cat(core_results, "j50_pt1", np.float64),
        "j50_eta1": cat(core_results, "j50_eta1", np.float64),
        "j50_phi1": cat(core_results, "j50_phi1", np.float64),
        "j50_pt1_raw": cat(core_results, "j50_pt1_raw", np.float64),
        "j50_pt2": cat(core_results, "j50_pt2", np.float64),
        "j50_eta2": cat(core_results, "j50_eta2", np.float64),
        "j50_phi2": cat(core_results, "j50_phi2", np.float64),
        "j50_pt2_raw": cat(core_results, "j50_pt2_raw", np.float64),
        "j50_dphi12": cat(core_results, "j50_dphi12", np.float64),
        "j35_pt1": cat(core_results, "j35_pt1", np.float64),
        "j35_pass_j50": cat(core_results, "j35_pass_j50", bool),
        "mb_hf_gev": mb_hf_all,
        "mb_pt1": mb_pt1_all,
        "mb_pass_j50": mb_j50_all,
        "mb_pass_j35": mb_j35_all,
    }

    total_runtime_s = time.perf_counter() - t_start
    meta = {
        "analysis": "CMS-HIN-10-004 Open Data Replication",
        "collision_system": "PbPb",
        "sqrts_nn_tev": 2.76,
        "lumi_full_2010_inv_ub": 6.7,
        "core_dataset": "/HICorePhysics/HIRun2010-PromptReco-v3/RECO (CERN Open Data Record 14010)",
        "mb_dataset": "/HIAllPhysics/HIRun2010-ZS-v2/RECO (CERN Open Data Record 14011)",
        "golden_json": "Cert_150436-152957_HI7TeV_StreamExpress_Collisions10_JSON_v2.txt",
        "jec_tag": "Spring10_HI_IC5Calo (L2Relative + L3Absolute)",
        "n_core_files_requested": n_core,
        "n_core_files_processed": len(core_results),
        "n_core_runs_sampled": len(runs_observed),
        "core_runs_sampled": sorted(runs_observed),
        "n_mb_files_processed": len(mb_results),
        "n_mb_runs_sampled": len(mb_runs_observed),
        "n_mb_cleaned_events": int(len(mb_hf_all)),
        "runtime_seconds": round(total_runtime_s, 2),
        "software_versions": {
            "python": sys.version.split()[0],
            "uproot": uproot.__version__,
            "awkward": ak.__version__,
            "vector": vector.__version__,
            "hist": hist.__version__,
            "scipy": scipy.__version__,
            "numpy": np.__version__,
        },
    }

    np.savez_compressed(cache_npz, **cache_arrays)
    cache_meta_json.write_text(
        json.dumps({"cutflow_totals": cf_totals, "metadata": meta}, indent=2) + "\n",
        encoding="utf-8",
    )
    log.info("Saved compact event cache to %s (%.1f KB)", cache_npz, cache_npz.stat().st_size / 1024.0)

    summary = build_summary_and_histograms(cache_arrays, cf_totals, meta)
    out_summary.write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    log.info(
        "Done in %.1fs! N_dijet=%d | <A_J>(50-100%%)=%.3f -> <A_J>(0-10%%)=%.3f (%.1f sigma) | R_B(50-100%%)=%.3f -> R_B(0-10%%)=%.3f (%.1f sigma)",
        total_runtime_s,
        cf_totals["dijet_dphi_gt_2pi3"],
        summary["per_centrality"]["50-100%"]["aj_mean"],
        summary["per_centrality"]["0-10%"]["aj_mean"],
        summary["validation_checks"]["aj_shift_significance_sigma"],
        summary["per_centrality"]["50-100%"]["rb_aj_lt_015"],
        summary["per_centrality"]["0-10%"]["rb_aj_lt_015"],
        summary["validation_checks"]["rb_aj_drop_significance_sigma"],
    )


if __name__ == "__main__":
    app()
