"""compare_outputs.py — Comprehensive statistical comparison of UTESpac Python vs MATLAB outputs.

After running the MATLAB pipeline (produces .mat files) and the Python pipeline (produces .pkl
files) for any site under UTESpac_MATLAB/, run this script to get an expansive field-by-field,
column-by-column statistical report.

Usage
-----
    python3 compare_outputs.py                         # all sites, GPF only (default)
    python3 compare_outputs.py --site siteFire1        # one site, GPF only
    python3 compare_outputs.py --site siteFire1 --lpf  # LPF only
    python3 compare_outputs.py --avg                   # 30-min averaged only, GPF
    python3 compare_outputs.py --raw                   # 20 Hz raw only
    python3 compare_outputs.py --csv results.csv       # save stats table to CSV
    python3 compare_outputs.py --tol-rel 0.05          # change PASS threshold (default 0.01)

Output
------
For each matched (mat, pkl) pair the script prints a table with one row per numeric column:

  Status  Field            Column                  N_valid  Bias        RMSE        max_abs
          p99_abs     max_rel     mean_rel    R2      NaN_mismatch

Status codes:
  OK   — max_rel < tol_rel  (default 1 %)
  WARN — 1 % ≤ max_rel < 10 %
  BIG  — max_rel ≥ 10 %      (or absolute error > tol_abs when ref ≈ 0)
"""

import os
import sys
import io
import argparse
import pickle
import warnings
import csv as _csv
from datetime import datetime
import numpy as np
import scipy.io as sio

# ── paths ────────────────────────────────────────────────────────────────────

ROOT_PY    = os.path.dirname(os.path.abspath(__file__))
MATLAB_DIR = os.path.join(ROOT_PY, "UTESpac_MATLAB")


# ── file discovery ───────────────────────────────────────────────────────────

def discover_pairs(matlab_dir, site_filter=None, lpf_only=False, gpf_only=False,
                   avg_only=False, raw_only=False):
    """Return list of dicts describing matched (pkl, mat) file pairs.

    .mat files are looked up in  matlab_dir/site/output/.
    .pkl files are looked up in  matlab_dir/site/output/  first, then in
    ROOT_PY/site/output/ (where the Python pipeline writes when rootFolder
    points at UTESpac_Python/ rather than UTESpac_MATLAB/).

    Each dict has keys: site, pf_mode, output_type, pkl_path, mat_path.
    """
    pairs = []

    # Collect candidate site names from both search roots
    candidate_sites = set()
    for search_root in (matlab_dir, ROOT_PY):
        try:
            candidate_sites.update(
                d for d in os.listdir(search_root)
                if d.startswith("site") and os.path.isdir(os.path.join(search_root, d))
            )
        except FileNotFoundError:
            pass

    if not candidate_sites:
        print(f"ERROR: no site* directories found under {matlab_dir} or {ROOT_PY}")
        return pairs

    for site in sorted(candidate_sites):
        if site_filter and site != site_filter:
            continue

        mat_dir = os.path.join(matlab_dir, site, "output")
        if not os.path.isdir(mat_dir):
            continue

        # .mat files always come from the MATLAB output folder
        mat_files = os.listdir(mat_dir)
        mat_stems = {os.path.splitext(f)[0]: f for f in mat_files if f.endswith(".mat")}

        # .pkl files: prefer matlab_dir copy, fall back to ROOT_PY copy
        pkl_stems: dict = {}
        for pkl_root in (matlab_dir, ROOT_PY):
            pkl_dir = os.path.join(pkl_root, site, "output")
            if not os.path.isdir(pkl_dir):
                continue
            for f in os.listdir(pkl_dir):
                if f.endswith(".pkl"):
                    stem = os.path.splitext(f)[0]
                    if stem not in pkl_stems:   # matlab_dir takes priority
                        pkl_stems[stem] = os.path.join(pkl_dir, f)

        for stem in sorted(set(pkl_stems) & set(mat_stems)):
            # Parse stem: e.g. Fire1_30minAvg_LPF_LinDet_2025_06_08
            #                or Fire1_raw_GPF_LinDet_2025_06_08
            parts = stem.split("_")
            pf_mode    = next((p for p in parts if p in ("LPF", "GPF")), "unknown")
            output_type = "raw" if "raw" in parts else "avg"

            if lpf_only and pf_mode != "LPF":
                continue
            if gpf_only and pf_mode != "GPF":
                continue
            if avg_only and output_type != "avg":
                continue
            if raw_only and output_type != "raw":
                continue

            pairs.append({
                "site":        site,
                "pf_mode":     pf_mode,
                "output_type": output_type,
                "stem":        stem,
                "pkl_path":    pkl_stems[stem],
                "mat_path":    os.path.join(mat_dir, mat_stems[stem]),
            })
    return pairs


# ── loaders ──────────────────────────────────────────────────────────────────

def load_pkl(path):
    with open(path, "rb") as fh:
        return pickle.load(fh)


def load_mat_avg(path):
    """Load averaged MATLAB output struct → {field: ndarray}."""
    mat = sio.loadmat(path, squeeze_me=True)
    ref = mat["output"]
    names = ref.dtype.names
    values = ref.item()
    result = {}
    for name, val in zip(names, values):
        try:
            arr = np.asarray(val, dtype=float)
            if arr.ndim == 1:
                arr = arr.reshape(-1, 1)
            result[name] = arr
        except (TypeError, ValueError):
            pass
    return result


def load_mat_raw(path):
    """Load raw MATLAB output struct → {field: ndarray}."""
    mat = sio.loadmat(path, squeeze_me=True)
    ref = mat["rawFlux"]
    names = ref.dtype.names
    values = ref.item()
    result = {}
    for name, val in zip(names, values):
        try:
            arr = np.asarray(val, dtype=float)
            if arr.ndim == 1:
                arr = arr.reshape(-1, 1)
            result[name] = arr
        except (TypeError, ValueError):
            pass
    return result


# ── header lookup ────────────────────────────────────────────────────────────

def get_header(pkl, field):
    """Return a flat column-label list for a field from the pkl dict.

    Handles both flat lists ['a','b',...] and nested [['a','b',...]] storage.
    """
    for suffix in (f"{field}Header", f"{field}header", f"{field}Header".lower()):
        if suffix in pkl:
            h = pkl[suffix]
            # unwrap one level of nesting when stored as [['col1', 'col2', ...]]
            if h and isinstance(h[0], list):
                h = h[0]
            return h
    return None


# ── per-column statistics ─────────────────────────────────────────────────────

def column_stats(py_col, ref_col, tol_rel=0.01, tol_abs=1e-9):
    """Compute comprehensive statistics comparing two 1-D arrays.

    Returns a dict; all float values are NaN when no valid pairs exist.
    """
    py  = np.asarray(py_col,  dtype=float).ravel()
    ref = np.asarray(ref_col, dtype=float).ravel()
    n   = min(len(py), len(ref))
    py, ref = py[:n], ref[:n]

    py_nan  = np.isnan(py)
    ref_nan = np.isnan(ref)
    nan_mismatch = int(np.sum(py_nan != ref_nan))
    mask = ~(py_nan | ref_nan)
    n_valid = int(mask.sum())

    null = float("nan")
    if n_valid == 0:
        return {
            "N_total": n, "N_valid": 0, "N_nan_mismatch": nan_mismatch,
            "bias": null, "RMSE": null,
            "max_abs": null, "mean_abs": null, "p95_abs": null, "p99_abs": null,
            "max_rel": null, "mean_rel": null, "R2": null,
            "status": "SKIP",
        }

    diff  = py[mask] - ref[mask]
    adiff = np.abs(diff)
    ref_m = ref[mask]
    denom = np.abs(ref_m).copy()
    denom[denom < tol_abs] = tol_abs
    rel   = adiff / denom

    bias     = float(diff.mean())
    rmse     = float(np.sqrt((diff**2).mean()))
    max_abs  = float(adiff.max())
    mean_abs = float(adiff.mean())
    p95_abs  = float(np.percentile(adiff, 95))
    p99_abs  = float(np.percentile(adiff, 99))
    max_rel  = float(rel.max())
    mean_rel = float(rel.mean())

    # R² (correlation coefficient squared) — meaningful only when ref has variance
    ref_std = float(ref_m.std())
    if ref_std > 1e-14 and n_valid > 2:
        r = float(np.corrcoef(py[mask], ref_m)[0, 1])
        r2 = r**2 if np.isfinite(r) else null
    else:
        r2 = null

    if max_rel < tol_rel:
        status = "OK  "
    elif max_rel < 0.10:
        status = "WARN"
    else:
        status = "BIG "

    return {
        "N_total": n, "N_valid": n_valid, "N_nan_mismatch": nan_mismatch,
        "bias": bias, "RMSE": rmse,
        "max_abs": max_abs, "mean_abs": mean_abs,
        "p95_abs": p95_abs, "p99_abs": p99_abs,
        "max_rel": max_rel, "mean_rel": mean_rel,
        "R2": r2,
        "status": status,
    }


# ── field comparison ──────────────────────────────────────────────────────────

def compare_field(field, py_arr, ref_arr, headers, tol_rel, tol_abs):
    """Compare two 2-D arrays column by column. Returns list of row dicts."""
    rows = []
    py  = np.asarray(py_arr,  dtype=float)
    ref = np.asarray(ref_arr, dtype=float)
    if py.ndim == 1:  py  = py.reshape(-1, 1)
    if ref.ndim == 1: ref = ref.reshape(-1, 1)

    nrow = min(py.shape[0], ref.shape[0])
    ncol = min(py.shape[1], ref.shape[1])
    py  = py[:nrow, :ncol]
    ref = ref[:nrow, :ncol]

    for c in range(ncol):
        col_label = str(headers[c]) if (headers and c < len(headers)) else f"col{c}"
        stats = column_stats(py[:, c], ref[:, c], tol_rel=tol_rel, tol_abs=tol_abs)
        rows.append({"field": field, "column": col_label, **stats})
    return rows


# ── averaged output comparison ────────────────────────────────────────────────

# Fields to skip (non-numeric metadata)
SKIP_FIELDS = {"tableNames", "warnings", "dataInfo"}

# Fields whose timestamp column (col 0) should be skipped
HAS_TIMESTAMP_COL0 = {
    "spdAndDir", "rotatedSonic", "PFSonic",
    "H", "Hlat", "tau", "tke", "sigma", "R", "L",
    "eta", "delta_flux_ctrb", "delta_time_ctrb", "turbtr", "epsilon",
    "skew", "H_SNSP", "Flux_lat", "LHflux", "CO2flux", "derivedT",
}


def compare_avg(pkl, mat_ref, tol_rel, tol_abs):
    """Compare averaged output dict (pkl) vs MATLAB struct (mat_ref). Returns stat rows."""
    rows = []
    all_fields = set(pkl) | set(mat_ref)

    for field in sorted(all_fields):
        if "Header" in field or "header" in field:
            continue
        if field in SKIP_FIELDS:
            continue
        if field not in pkl:
            rows.append({"field": field, "column": "(all)", "status": "MISS_PY",
                         **{k: float("nan") for k in
                            ["N_total","N_valid","N_nan_mismatch","bias","RMSE",
                             "max_abs","mean_abs","p95_abs","p99_abs","max_rel","mean_rel","R2"]}})
            continue
        if field not in mat_ref:
            rows.append({"field": field, "column": "(all)", "status": "MISS_REF",
                         **{k: float("nan") for k in
                            ["N_total","N_valid","N_nan_mismatch","bias","RMSE",
                             "max_abs","mean_abs","p95_abs","p99_abs","max_rel","mean_rel","R2"]}})
            continue

        try:
            py_arr  = np.asarray(pkl[field],     dtype=float)
            ref_arr = np.asarray(mat_ref[field], dtype=float)
        except (TypeError, ValueError):
            continue

        if py_arr.ndim == 1:  py_arr  = py_arr.reshape(-1, 1)
        if ref_arr.ndim == 1: ref_arr = ref_arr.reshape(-1, 1)

        headers = get_header(pkl, field)

        # Skip timestamp column 0 for fields that carry it
        col_start = 1 if field in HAS_TIMESTAMP_COL0 else 0
        if headers and col_start > 0:
            headers = headers[col_start:]

        rows.extend(compare_field(
            field,
            py_arr[:, col_start:],
            ref_arr[:, col_start:],
            headers, tol_rel, tol_abs,
        ))
    return rows


# ── raw output comparison ─────────────────────────────────────────────────────

# Map pkl raw keys → MATLAB rawFlux field names (when they differ)
RAW_KEY_MAP = {
    "uPF":       "uPF",
    "vPF":       "vPF",
    "wPF":       "wPF",
    "u_tilt":    "u_tilt",
    "v_tilt":    "v_tilt",
    "w_tilt":    "w_tilt",
    "sonTs":     "sonTs",
    "Theta_v_son": "Theta_v_son",
    "P":         "P",
    "rhov":      "rhov",
    "rhovPrime": "rhovPrime",
    "rhoCO2":    "rhoCO2",
    "rhoCO2Prime": "rhoCO2Prime",
}

# For 20 Hz raw, use every Nth sample to keep memory manageable
RAW_STRIDE = 100


def compare_raw(pkl, mat_ref, tol_rel, tol_abs):
    """Compare raw 20 Hz output (pkl) vs MATLAB rawFlux struct (mat_ref). Returns stat rows."""
    rows = []
    all_keys = set(RAW_KEY_MAP) | set(mat_ref)

    for py_key in sorted(all_keys):
        mat_key = RAW_KEY_MAP.get(py_key, py_key)
        py_present  = py_key  in pkl
        ref_present = mat_key in mat_ref

        if not py_present and not ref_present:
            continue
        if not py_present:
            rows.append({"field": py_key, "column": "(all)", "status": "MISS_PY",
                         **{k: float("nan") for k in
                            ["N_total","N_valid","N_nan_mismatch","bias","RMSE",
                             "max_abs","mean_abs","p95_abs","p99_abs","max_rel","mean_rel","R2"]}})
            continue
        if not ref_present:
            rows.append({"field": py_key, "column": "(all)", "status": "MISS_REF",
                         **{k: float("nan") for k in
                            ["N_total","N_valid","N_nan_mismatch","bias","RMSE",
                             "max_abs","mean_abs","p95_abs","p99_abs","max_rel","mean_rel","R2"]}})
            continue

        try:
            py_arr  = np.asarray(pkl[py_key],      dtype=float)
            ref_arr = np.asarray(mat_ref[mat_key], dtype=float)
        except (TypeError, ValueError):
            continue

        if py_arr.ndim  == 0: py_arr  = py_arr.reshape(1, 1)
        if ref_arr.ndim == 0: ref_arr = ref_arr.reshape(1, 1)
        if py_arr.ndim  == 1: py_arr  = py_arr.reshape(-1, 1)
        if ref_arr.ndim == 1: ref_arr = ref_arr.reshape(-1, 1)

        # Subsample to keep comparison fast
        py_arr  = py_arr [::RAW_STRIDE]
        ref_arr = ref_arr[::RAW_STRIDE]

        ncol = min(py_arr.shape[1], ref_arr.shape[1])
        for c in range(ncol):
            col_label = f"col{c}" if ncol > 1 else py_key
            stats = column_stats(py_arr[:, c], ref_arr[:, c], tol_rel=tol_rel, tol_abs=tol_abs)
            rows.append({"field": py_key, "column": col_label, **stats})

    return rows


# ── printing ──────────────────────────────────────────────────────────────────

STAT_KEYS  = ["bias", "RMSE", "max_abs", "p99_abs", "max_rel", "mean_rel", "R2"]
STAT_WIDTHS = [12, 12, 12, 12, 10, 10, 8]


def _fmt(val, width):
    if not isinstance(val, float) or not np.isfinite(val):
        return f"{'nan':>{width}}"
    if abs(val) == 0:
        return f"{'0':>{width}}"
    if abs(val) >= 1000 or (abs(val) < 0.001 and val != 0):
        return f"{val:>{width}.3e}"
    return f"{val:>{width}.6f}"


def print_table(rows, site, stem, tol_rel):
    col_w = 30
    field_w = 22
    n_w = 8
    nm_w = 7

    sep = "-" * (6 + field_w + col_w + n_w + nm_w + sum(STAT_WIDTHS) + len(STAT_WIDTHS) + 4)
    header_line = (
        f"  {'':6s}{'Field':<{field_w}}{'Column':<{col_w}}"
        f"{'N_valid':>{n_w}} {'NaN_mm':>{nm_w}}"
    )
    for key, w in zip(STAT_KEYS, STAT_WIDTHS):
        header_line += f"  {key:>{w}}"

    print(f"\n{'='*len(sep)}")
    print(f"  Site: {site}   File: {stem}   tol_rel={tol_rel:.1%}")
    print(f"{'='*len(sep)}")
    print(header_line)
    print(sep)

    for row in rows:
        status  = row.get("status", "")
        field   = str(row.get("field",  ""))[:field_w]
        column  = str(row.get("column", ""))[:col_w]
        n_valid = row.get("N_valid", float("nan"))
        nan_mm  = row.get("N_nan_mismatch", float("nan"))

        n_str  = f"{n_valid:>{n_w}d}" if isinstance(n_valid, int) else f"{'nan':>{n_w}}"
        nm_str = f"{nan_mm:>{nm_w}d}" if isinstance(nan_mm, int) else f"{'nan':>{nm_w}}"

        line = f"  {status:<6s}{field:<{field_w}}{column:<{col_w}}{n_str} {nm_str}"
        for key, w in zip(STAT_KEYS, STAT_WIDTHS):
            line += f"  {_fmt(row.get(key, float('nan')), w)}"
        print(line)

    print(sep)
    ok   = sum(1 for r in rows if r.get("status","").strip() == "OK")
    warn = sum(1 for r in rows if r.get("status","").strip() == "WARN")
    big  = sum(1 for r in rows if r.get("status","").strip() == "BIG")
    skip = sum(1 for r in rows if r.get("status","").strip() in ("SKIP","MISS_PY","MISS_REF"))
    print(f"  Summary: {ok} OK  {warn} WARN  {big} BIG  {skip} SKIP/MISSING\n")


# ── CSV export ────────────────────────────────────────────────────────────────

ALL_COLS = ["site", "pf_mode", "output_type", "stem", "field", "column", "status",
            "N_total", "N_valid", "N_nan_mismatch",
            "bias", "RMSE", "max_abs", "mean_abs", "p95_abs", "p99_abs",
            "max_rel", "mean_rel", "R2"]


def save_csv(all_rows, path):
    with open(path, "w", newline="") as fh:
        writer = _csv.DictWriter(fh, fieldnames=ALL_COLS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(all_rows)
    print(f"\nStats saved to: {path}")


# ── main ──────────────────────────────────────────────────────────────────────

LOG_DIR = os.path.join(ROOT_PY, "comparison_logs")


class _Tee:
    """Write to both stdout and a file simultaneously."""
    def __init__(self, fh):
        self._fh = fh
        self._stdout = sys.stdout
    def write(self, s):
        self._stdout.write(s)
        self._fh.write(s)
    def flush(self):
        self._stdout.flush()
        self._fh.flush()


def main():
    parser = argparse.ArgumentParser(
        description="Comprehensive statistical comparison of UTESpac Python vs MATLAB outputs.")
    parser.add_argument("--site",     default=None,  help="Limit to one site folder name")
    parser.add_argument("--lpf",      action="store_true", help="LPF files only")
    parser.add_argument("--gpf",      action="store_true", help="GPF files only")
    parser.add_argument("--avg",      action="store_true", help="30-min averaged files only")
    parser.add_argument("--raw",      action="store_true", help="20 Hz raw files only")
    parser.add_argument("--csv",      default=None,  metavar="FILE", help="Save stats to CSV")
    parser.add_argument("--tol-rel",  default=0.01,  type=float,
                        help="Relative error threshold for OK/WARN boundary (default 0.01)")
    parser.add_argument("--tol-abs",  default=1e-9,  type=float,
                        help="Absolute floor used when ref ≈ 0 (default 1e-9)")
    args = parser.parse_args()

    # Auto-generate log file name: comparison_logs/<site>_<YYYYMMDD_HHMMSS>.txt
    os.makedirs(LOG_DIR, exist_ok=True)
    site_tag  = args.site if args.site else "all"
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path  = os.path.join(LOG_DIR, f"{site_tag}_{timestamp}.txt")
    csv_path  = args.csv or os.path.join(LOG_DIR, f"{site_tag}_{timestamp}.csv")

    with open(log_path, "w") as log_fh:
        sys.stdout = _Tee(log_fh)

        # Default to GPF when neither --lpf nor --gpf is given (GPF is the mode
        # used for analyses and AmeriFlux submission).
        gpf_only = args.gpf or (not args.lpf and not args.gpf)
        pairs = discover_pairs(
            MATLAB_DIR,
            site_filter=args.site,
            lpf_only=args.lpf,
            gpf_only=gpf_only,
            avg_only=args.avg,
            raw_only=args.raw,
        )

        if not pairs:
            print("No matched (pkl, mat) pairs found. Check that both pipelines have been run.")
            sys.stdout = sys.stdout._stdout
            sys.exit(1)

        print(f"Found {len(pairs)} matched file pair(s).")
        all_csv_rows = []

        for pair in pairs:
            site   = pair["site"]
            stem   = pair["stem"]
            pf     = pair["pf_mode"]
            otype  = pair["output_type"]

            print(f"\nLoading: {stem}")
            try:
                pkl = load_pkl(pair["pkl_path"])
            except Exception as exc:
                print(f"  ERROR loading pkl: {exc}")
                continue
            try:
                if otype == "raw":
                    mat_ref = load_mat_raw(pair["mat_path"])
                else:
                    mat_ref = load_mat_avg(pair["mat_path"])
            except Exception as exc:
                print(f"  ERROR loading mat: {exc}")
                continue

            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                if otype == "raw":
                    rows = compare_raw(pkl, mat_ref, args.tol_rel, args.tol_abs)
                else:
                    rows = compare_avg(pkl, mat_ref, args.tol_rel, args.tol_abs)

            for r in rows:
                r.update({"site": site, "pf_mode": pf, "output_type": otype, "stem": stem})
            all_csv_rows.extend(rows)

            print_table(rows, site, stem, args.tol_rel)

        save_csv(all_csv_rows, csv_path)
        print(f"\nLog saved to:  {log_path}")
        print(f"CSV saved to:  {csv_path}")

        sys.stdout = sys.stdout._stdout


if __name__ == "__main__":
    main()
