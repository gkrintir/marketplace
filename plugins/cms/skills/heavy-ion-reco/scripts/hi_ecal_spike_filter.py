#!/usr/bin/env python3
"""Evaluate the Run-1 CMS Heavy-Ion ECAL Barrel APD spike filter (HiEcalRecHitSpikeFilter).

Implements the exact crystal-level Swiss-cross (1 - E_4/E_1) and timing (|t|)
rejection of RecoHI/HiEgammaAlgos/HiEcalRecHitSpikeFilter (CMSSW_3_9_2_patch5,
avoidIeta85=False) on Run-1 CMS EDM RECO files without requiring CMSSW, using a
lossless BasicCluster/SuperCluster seed-crystal E_T pre-filter and surgical
uproot TBasket reads.

Usage:
    hi_ecal_spike_filter.py --self-test
    hi_ecal_spike_filter.py <file_or_url.root> [--max-events 500] [--min-et 20.0]

Stdout: KEY=VALUE lines (SELF_TEST=PASS or EVENTS, PREFILTER_CANDIDATES, PASSED,
        REJECTED, PASS_FRACTION, MIN_ET_GEV). Diagnostics on stderr.
Exit codes: 0 success, 1 self-test failed, 2 invalid arguments or unreadable input.
Tested with Python 3.11, numpy 2.4, uproot 5.7 on 2010 PbPb HICorePhysics RECO
(CERN Open Data record 14010).
"""
from __future__ import annotations

import argparse
import math
import struct
import sys
from typing import Iterable

import numpy as np

ECAL_EB_BRANCH = "EcalRecHitsSorted_ecalRecHit_EcalRecHitsEB_RECO./EcalRecHitsSorted_ecalRecHit_EcalRecHitsEB_RECO.obj"
BC_ENERGY_BRANCH = (
    "recoCaloClusters_uncleanedHybridSuperClusters_hybridBarrelBasicClusters_RECO."
    "/recoCaloClusters_uncleanedHybridSuperClusters_hybridBarrelBasicClusters_RECO.obj"
    "/recoCaloClusters_uncleanedHybridSuperClusters_hybridBarrelBasicClusters_RECO.obj.energy_"
)
BC_SEED_BRANCH = (
    "recoCaloClusters_uncleanedHybridSuperClusters_hybridBarrelBasicClusters_RECO."
    "/recoCaloClusters_uncleanedHybridSuperClusters_hybridBarrelBasicClusters_RECO.obj"
    "/recoCaloClusters_uncleanedHybridSuperClusters_hybridBarrelBasicClusters_RECO.obj.seedId_.id_"
)
BC_ENERGY_LEAF = "recoCaloClusters_uncleanedHybridSuperClusters_hybridBarrelBasicClusters_RECO.obj.energy_"
BC_SEED_LEAF = "recoCaloClusters_uncleanedHybridSuperClusters_hybridBarrelBasicClusters_RECO.obj.seedId_.id_"
SC_UNSPLIT_BRANCH = (
    "recoSuperClusters_uncleanedHybridSuperClusters__RECO."
    "/recoSuperClusters_uncleanedHybridSuperClusters__RECO.obj"
)

DEG_TO_RAD = 0.017453292519943295  # pi / 180 (crystal pitch = 1 deg in eta/phi)


def pack_eb_detid(ieta: int, iphi: int) -> int:
    """Pack signed crystal indices ieta in [-85..-1, +1..+85], iphi in [1..360] into an EBDetId uint32."""
    if ieta == 0 or abs(ieta) > 85 or not (1 <= iphi <= 360):
        raise ValueError(f"invalid EBDetId indices: ieta={ieta}, iphi={iphi}")
    det_ecal_eb = (3 << 28) | (1 << 25)
    zside_bit = 0x10000 if ieta > 0 else 0
    return det_ecal_eb | zside_bit | ((abs(ieta) & 0x7F) << 9) | (iphi & 0x1FF)


def unpack_eb_detid(raw_ids: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (zside, ieta_abs, iphi) arrays from packed uint32 EBDetId values."""
    ids = np.asarray(raw_ids, dtype=np.uint32)
    zside = np.where((ids & 0x10000) != 0, 1, -1)
    ieta_abs = ((ids >> 9) & 0x7F).astype(np.int32)
    iphi = (ids & 0x1FF).astype(np.int32)
    return zside, ieta_abs, iphi


def eb_crystal_eta(raw_ids: np.ndarray) -> np.ndarray:
    """Approximate crystal pseudorapidity eta = zside * (|ieta| - 0.5) * (pi / 180)."""
    zside, ieta_abs, _ = unpack_eb_detid(raw_ids)
    return zside * (ieta_abs - 0.5) * DEG_TO_RAD


def eb_seed_et_upper_bound(cluster_energies: np.ndarray, seed_ids: np.ndarray) -> np.ndarray:
    """Lossless upper bound on seed-crystal E_T from BasicCluster energy and seed EBDetId.

    Because E_cluster >= E_seed and |eta_crystal| >= max(0, (|ieta| - 2) * deg),
    E_cluster / cosh(min_eta) strictly bounds any constituent crystal E_T from above.
    """
    if len(cluster_energies) == 0:
        return np.empty(0, dtype=np.float64)
    _, ieta_abs, _ = unpack_eb_detid(seed_ids)
    min_abs_eta = np.maximum(0.0, (ieta_abs - 2.0) * DEG_TO_RAD)
    return np.asarray(cluster_energies, dtype=np.float64) / np.cosh(min_abs_eta)


def hi_ecal_spike_filter(
    ids: np.ndarray,
    energies: np.ndarray,
    times: np.ndarray,
    min_et: float = 20.0,
    swiss_threshold: float = 0.95,
    time_threshold: float = 4.0,
    avoid_ieta85: bool = False,
) -> bool:
    """Return True if event passes HiEcalRecHitSpikeFilter, False if an APD spike is found."""
    if len(ids) == 0:
        return True
    cand = np.flatnonzero(energies > min_et)
    if len(cand) == 0:
        return True

    cand_ids = ids[cand]
    _, ieta_abs_c, _ = unpack_eb_detid(cand_ids)
    if avoid_ieta85:
        keep = ieta_abs_c < 85
        cand = cand[keep]
        cand_ids = cand_ids[keep]
        if len(cand) == 0:
            return True

    et_c = energies[cand] / np.cosh(eb_crystal_eta(cand_ids))
    hi_mask = et_c > min_et
    if not np.any(hi_mask):
        return True

    hi_indices = cand[hi_mask]
    if np.any(np.abs(times[hi_indices]) > time_threshold):
        return False

    zside, ieta_abs, iphi = unpack_eb_detid(ids)
    row = np.where(zside < 0, -ieta_abs + 85, ieta_abs + 84)  # -85->0 .. -1->84, +1->85 .. +85->169
    col = iphi - 1  # 0..359

    grid = np.zeros((170, 360), dtype=np.float32)
    grid[row, col] = energies  # EcalTools::recHitE sums raw crystal energies (including negative noise)

    for idx in hi_indices:
        r = int(row[idx])
        c = int(col[idx])
        e1 = float(energies[idx])
        s4 = float(grid[r, (c - 1) % 360]) + float(grid[r, (c + 1) % 360])
        if r > 0:
            s4 += float(grid[r - 1, c])
        if r < 169:
            s4 += float(grid[r + 1, c])
        swiss = 1.0 - (s4 / e1) if e1 != 0.0 else 0.0
        if swiss > swiss_threshold:
            return False
    return True


def parse_eb_rechits_record(raw_rec: bytes) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Unpack a big-endian ROOT EcalRecHitsSorted record into (ids, energies, times)."""
    header_len = 20
    if len(raw_rec) < header_len:
        return np.empty(0, dtype=np.uint32), np.empty(0, dtype=np.float32), np.empty(0, dtype=np.float32)
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


def _read_selected_basket_records(branch, entry_indices: Iterable[int]) -> dict[int, bytes]:
    """Read only the TBaskets covering entry_indices via fBasketSeek/fBasketBytes."""
    import queue
    import uproot

    indices = sorted({int(i) for i in entry_indices})
    if not indices:
        return {}
    offsets = np.asarray(branch.member("fBasketEntry"), dtype=np.int64)
    seeks = np.asarray(branch.member("fBasketSeek"), dtype=np.int64)
    nbytes = np.asarray(branch.member("fBasketBytes"), dtype=np.int64)
    valid_nb = int(np.sum(seeks > 0))
    offsets = offsets[: valid_nb + 1]

    bsk_ids = np.searchsorted(offsets, indices, side="right") - 1
    uniq_bsks = sorted({int(b) for b in bsk_ids if 0 <= int(b) < valid_nb})
    ranges = [(int(seeks[b]), int(seeks[b] + nbytes[b])) for b in uniq_bsks]
    chunks = branch.file.source.chunks(ranges, notifications=queue.Queue())

    bsk_map = {}
    for b, chk in zip(uniq_bsks, chunks):
        cur = uproot.source.cursor.Cursor(chk.start)
        bsk_map[b] = uproot.models.TBasket.Model_TBasket.read(
            chk, cur, {"basket_num": b}, branch.file, branch.file, branch
        )

    out: dict[int, bytes] = {}
    for idx, b in zip(indices, bsk_ids):
        b_int = int(b)
        bsk = bsk_map[b_int]
        local_i = int(idx - offsets[b_int])
        off = bsk.byte_offsets
        out[idx] = bytes(bsk.data) if off is None else bytes(bsk.data[off[local_i] : off[local_i + 1]])
    return out


def run_self_test() -> int:
    """Deterministic unit test of EBDetId geometry, Swiss-cross, timing, and ieta=85 boundary."""
    # 1. Multi-crystal EM shower at ieta=+40, iphi=100: E1=45 GeV, neighbors=18 GeV -> swiss=0.60 -> PASS
    ids_em = np.array(
        [pack_eb_detid(40, 100), pack_eb_detid(40, 101), pack_eb_detid(41, 100)],
        dtype=np.uint32,
    )
    ene_em = np.array([45.0, 10.0, 8.0], dtype=np.float32)
    tim_em = np.array([0.3, 0.2, 0.4], dtype=np.float32)
    ok_em = hi_ecal_spike_filter(ids_em, ene_em, tim_em, min_et=20.0)

    # 2. Barrel-edge |ieta|=85 (|eta|=1.479) single-crystal APD spike: E=285 GeV (ET=123.7 GeV), E4=0.1 GeV -> REJECT
    ids_spk = np.array([pack_eb_detid(-85, 360), pack_eb_detid(-85, 1)], dtype=np.uint32)
    ene_spk = np.array([285.0, 0.1], dtype=np.float32)
    tim_spk = np.array([-0.5, 0.0], dtype=np.float32)
    rej_ieta85 = not hi_ecal_spike_filter(ids_spk, ene_spk, tim_spk, min_et=20.0, avoid_ieta85=False)
    pass_if_avoid85 = hi_ecal_spike_filter(ids_spk, ene_spk, tim_spk, min_et=20.0, avoid_ieta85=True)

    # 3. Out-of-time crystal hit (|t| = 6.5 ns > 4.0 ns) even with broad Swiss-cross -> REJECT
    tim_oot = np.array([6.5, 0.2, 0.4], dtype=np.float32)
    rej_oot = not hi_ecal_spike_filter(ids_em, ene_em, tim_oot, min_et=20.0)

    # 4. Moderate 32 GeV slow-neutron APD hit: rejected at photon minEt=20, accepted at jet minEt=50
    ids_mod = np.array([pack_eb_detid(5, 180)], dtype=np.uint32)
    ene_mod = np.array([32.0], dtype=np.float32)
    tim_mod = np.array([0.0], dtype=np.float32)
    rej_at_20 = not hi_ecal_spike_filter(ids_mod, ene_mod, tim_mod, min_et=20.0)
    pass_at_50 = hi_ecal_spike_filter(ids_mod, ene_mod, tim_mod, min_et=50.0)

    # 5. Seed upper-bound E_T strictly >= true crystal E_T for all |ieta| in 1..85
    all_ids = np.array([pack_eb_detid(ie, 1) for ie in list(range(-85, 0)) + list(range(1, 86))], dtype=np.uint32)
    enes = np.full(len(all_ids), 100.0, dtype=np.float64)
    true_et = enes / np.cosh(eb_crystal_eta(all_ids))
    bound_et = eb_seed_et_upper_bound(enes, all_ids)
    bound_ok = bool(np.all(bound_et >= true_et - 1e-12))

    checks = {
        "em_shower_passes": ok_em,
        "ieta85_spike_rejected": rej_ieta85,
        "avoid_ieta85_flag_respected": pass_if_avoid85,
        "out_of_time_rejected": rej_oot,
        "min_et_threshold_scaling": rej_at_20 and pass_at_50,
        "seed_upper_bound_lossless": bound_ok,
    }
    failed = [k for k, v in checks.items() if not v]
    if failed:
        print(f"error: failed self-test checks: {failed}", file=sys.stderr)
        print("SELF_TEST=FAIL")
        return 1
    print("SELF_TEST=PASS")
    print(f"CHECKS={len(checks)}/{len(checks)}")
    return 0


def inspect_reco_file(path: str, max_events: int, min_et: float, swiss_thr: float, time_thr: float) -> int:
    try:
        import awkward as ak
        import uproot
    except ImportError as exc:
        print(f"error: inspecting RECO files requires uproot and awkward ({exc})", file=sys.stderr)
        return 2

    try:
        open_kw = (
            {"handler": uproot.source.http.HTTPSource}
            if path.startswith(("http://", "https://"))
            else {}
        )
        f = uproot.open(path, **open_kw)
        ev = f["Events"]
    except Exception as exc:
        print(f"error: cannot open Events tree in {path}: {exc}", file=sys.stderr)
        return 2

    n_ev = min(int(ev.num_entries), int(max_events)) if max_events > 0 else int(ev.num_entries)
    if n_ev == 0:
        print("EVENTS=0\nPREFILTER_CANDIDATES=0\nPASSED=0\nREJECTED=0\nPASS_FRACTION=1.000000")
        return 0

    try:
        eb_branch = ev[ECAL_EB_BRANCH]
    except KeyError:
        print(f"error: branch {ECAL_EB_BRANCH} not found in {path}", file=sys.stderr)
        return 2

    # Lossless pre-filter via split BasicClusters when present; otherwise inspect all n_ev entries
    try:
        arr = ev.arrays(filter_name=[BC_ENERGY_BRANCH, BC_SEED_BRANCH], entry_stop=n_ev, how=dict)
        if BC_ENERGY_LEAF in arr and BC_SEED_LEAF in arr:
            bc_e = arr[BC_ENERGY_LEAF]
            bc_id = arr[BC_SEED_LEAF]
            counts = ak.to_numpy(ak.num(bc_e, axis=1))
            flat_ub = eb_seed_et_upper_bound(
                ak.to_numpy(ak.flatten(bc_e)),
                ak.to_numpy(ak.flatten(bc_id)),
            )
            ub_jagged = ak.unflatten(flat_ub, counts)
            max_ub = ak.to_numpy(ak.fill_none(ak.max(ub_jagged, axis=1), 0.0))
            candidates = np.flatnonzero(max_ub > (min_et - 0.5)).tolist()
        else:
            candidates = list(range(n_ev))
    except Exception:
        candidates = list(range(n_ev))

    recs = _read_selected_basket_records(eb_branch, candidates)
    rejected = 0
    for idx in candidates:
        ids, enes, times = parse_eb_rechits_record(recs[idx])
        if not hi_ecal_spike_filter(
            ids,
            enes,
            times,
            min_et=min_et,
            swiss_threshold=swiss_thr,
            time_threshold=time_thr,
            avoid_ieta85=False,
        ):
            rejected += 1

    passed = n_ev - rejected
    frac = passed / n_ev
    print(f"EVENTS={n_ev}")
    print(f"PREFILTER_CANDIDATES={len(candidates)}")
    print(f"PASSED={passed}")
    print(f"REJECTED={rejected}")
    print(f"PASS_FRACTION={frac:.6f}")
    print(f"MIN_ET_GEV={min_et:g}")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", nargs="?", help="local path or root:///https:// URL of a Run-1 EDM RECO file")
    ap.add_argument("--self-test", action="store_true", help="run deterministic crystal/geometry unit tests and exit")
    ap.add_argument("--max-events", type=int, default=500, help="max events to inspect (0 = all, default: 500)")
    ap.add_argument(
        "--min-et",
        type=float,
        default=20.0,
        help="crystal E_T threshold in GeV (20.0 = photon default; ~50.0 = jet-scale central PbPb)",
    )
    ap.add_argument("--swiss-threshold", type=float, default=0.95, help="max allowed 1 - E_4/E_1 (default: 0.95)")
    ap.add_argument("--time-threshold", type=float, default=4.0, help="max allowed |time| in ns (default: 4.0)")
    args = ap.parse_args(argv)

    if args.self_test:
        return run_self_test()
    if not args.input:
        ap.print_help(sys.stderr)
        return 2
    return inspect_reco_file(args.input, args.max_events, args.min_et, args.swiss_threshold, args.time_threshold)


if __name__ == "__main__":
    sys.exit(main())
