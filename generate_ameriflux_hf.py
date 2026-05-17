"""generate_ameriflux_hf.py – AmeriFlux high-frequency CSV export from UTESpac raw pkl output.

Format references:
  https://ameriflux.lbl.gov/data/how-to-upload-data/uploading-high-frequency-data/
  https://ameriflux.lbl.gov/half-hourly-hourly-data-upload-format/

Output structure:
  One CSV per 30-min period: <SITE_ID>_HF_<YYYYMMDDHHMM>_<YYYYMMDDHHMM>.csv
  All CSVs for each site type are collected into a zip archive at the end.
  Individual CSVs are removed after zipping.

Variables per sonic height (V = 1 = lowest height, AmeriFlux positional convention):
  TIMESTAMP         – YYYYMMDDHHMMSS.cc  (local standard time; cc = centiseconds)
  U_1_{V}_1         – planar-fit u wind component          [m/s]
  V_1_{V}_1         – planar-fit v wind component          [m/s]
  W_1_{V}_1         – planar-fit w wind component          [m/s]
  T_SONIC_1_{V}_1   – sonic temperature                    [°C]
  H2O_1_{V}_1       – water vapour molar density           [mmol/m³]
                         (rhov [g/m³] / 18.015 × 1000)
  CO2_1_{V}_1       – CO2 molar density                    [µmol/m³]
                         (rhoCO2 [mg/m³] / 44.010 × 1000)
  PA_1_1_1          – atmospheric pressure (one barometer) [kPa]

Unit note — Python vs MATLAB raw pkl convention:
  rhov:    g/m³   in both Python and MATLAB
  rhoCO2:  mg/m³  in Python pkl  (MATLAB stores in kg/m³; see fluxes.py comment)

Missing data: -9999  (NaN and ±Inf replaced before writing)
"""

import os
import re
import glob
import pickle
import zipfile
from collections import defaultdict

import numpy as np
import pandas as pd

# ── configuration ─────────────────────────────────────────────────────────────

SITE_ID = "US-xFM"    # replace with official AmeriFlux site ID when registered

ROOT_PY = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(ROOT_PY, "ameriflux_hf_output")

PF_TYPE = "GPF"       # "LPF" or "GPF"
MISSING = -9999.0

MW_H2O          = 18.015    # g/mol
MW_CO2          = 44.010    # g/mol
MATLAB_EPOCH    = 719529.0  # MATLAB days to Unix epoch (1970-01-01)

# ── helpers ────────────────────────────────────────────────────────────────────

def matlab_to_unix(t_arr):
    """MATLAB serial date float array → Unix seconds."""
    return (np.asarray(t_arr) - MATLAB_EPOCH) * 86400.0


def make_timestamps(t_arr):
    """Vectorised MATLAB serial → 'YYYYMMDDHHMMSS.cc' pandas Series."""
    dti = pd.to_datetime(matlab_to_unix(t_arr), unit="s")
    cs  = (dti.microsecond // 10000).astype(str).str.zfill(2)
    return dti.strftime("%Y%m%d%H%M%S") + "." + cs


def file_ts(t_val, offset_s=0.0):
    """Single MATLAB serial → YYYYMMDDHHMM string for file naming."""
    unix = matlab_to_unix(t_val) + offset_s
    return pd.Timestamp(unix, unit="s").strftime("%Y%m%d%H%M")


def safe(arr):
    """Replace NaN / ±Inf with MISSING (-9999)."""
    out = np.where(np.isfinite(arr), arr, MISSING)
    return out


# ── main ───────────────────────────────────────────────────────────────────────

os.makedirs(OUT_DIR, exist_ok=True)

raw_pkls = sorted(glob.glob(
    os.path.join(ROOT_PY, "site*", "output", f"*_raw_{PF_TYPE}_*.pkl")))

if not raw_pkls:
    print(f"No raw {PF_TYPE} pkl files found under {ROOT_PY}/site*/output/")
    raise SystemExit(1)

print(f"Found {len(raw_pkls)} raw {PF_TYPE} pkl file(s).")

# site-type → list of written CSV paths (for zip grouping)
site_csvs = defaultdict(list)

for pkl_path in raw_pkls:
    # Determine site type prefix (e.g. "siteIRGA" or "siteGill")
    site_dir  = os.path.basename(os.path.dirname(os.path.dirname(pkl_path)))
    site_type = re.match(r"(site[A-Za-z]+)\d", site_dir)
    site_type = site_type.group(1) if site_type else site_dir

    print(f"\nProcessing {os.path.basename(pkl_path)}  [{site_dir}] …")

    with open(pkl_path, "rb") as fh:
        raw = pickle.load(fh)

    t     = raw["t"]           # (N,)        MATLAB serial dates
    z     = raw["z"]           # (n_sonics,) heights [m]
    uPF   = raw["uPF"]         # (N, n_sonics) m/s
    vPF   = raw["vPF"]
    wPF   = raw["wPF"]
    Ts    = raw["sonTs"]       # (N, n_sonics) °C
    P_raw = raw["P"]           # (N, 2): col0 = timestamp, col1 = pressure [kPa]
    rhov  = raw["rhov"]        # (N, n_sonics) g/m³
    rCO2  = raw["rhoCO2"]      # (N, n_sonics) mg/m³  (Python convention)

    n_sonics = len(z)
    N        = len(t)

    # Infer sampling frequency
    dt_s = float(np.nanmedian(np.diff(t[:min(2000, N)]))) * 86400.0
    hz   = int(round(1.0 / dt_s))
    block = hz * 30 * 60       # samples per 30-min
    n_blocks = N // block

    print(f"  {hz} Hz  ·  {n_sonics} sonic(s) at {list(z)} m  ·  {n_blocks} × 30-min files")

    # V-index: V=1=lowest height
    sorted_z = sorted(z)
    v_of     = {zi: sorted_z.index(zi) + 1 for zi in z}

    # Convert all data to output units (vectorised over full 48-h array)
    u_out   = safe(uPF)
    v_out   = safe(vPF)
    w_out   = safe(wPF)
    ts_out  = safe(Ts)
    pa_out  = safe(P_raw[:, 1])                      # kPa
    h2o_out = safe(rhov)  / MW_H2O * 1e3             # g/m³ → mmol/m³
    co2_out = safe(rCO2)  / MW_CO2 * 1e3             # mg/m³ → µmol/m³

    # Build column name list (order: variables grouped by height, then PA at end)
    col_names = ["TIMESTAMP"]
    for zi in z:
        v = v_of[zi]
        col_names += [
            f"U_1_{v}_1", f"V_1_{v}_1", f"W_1_{v}_1",
            f"T_SONIC_1_{v}_1",
            f"H2O_1_{v}_1",
            f"CO2_1_{v}_1",
        ]
    col_names.append("PA_1_1_1")

    for bi in range(n_blocks):
        s0, s1 = bi * block, (bi + 1) * block
        t_blk  = t[s0:s1]

        t_start = file_ts(t_blk[0])
        t_end   = file_ts(t_blk[-1], offset_s=1.0 / hz)

        # Build DataFrame
        data = {"TIMESTAMP": make_timestamps(t_blk)}
        for si, zi in enumerate(z):
            v = v_of[zi]
            data[f"U_1_{v}_1"]       = u_out[s0:s1, si]
            data[f"V_1_{v}_1"]       = v_out[s0:s1, si]
            data[f"W_1_{v}_1"]       = w_out[s0:s1, si]
            data[f"T_SONIC_1_{v}_1"] = ts_out[s0:s1, si]
            data[f"H2O_1_{v}_1"]     = h2o_out[s0:s1, si]
            data[f"CO2_1_{v}_1"]     = co2_out[s0:s1, si]
        data["PA_1_1_1"] = pa_out[s0:s1]

        df = pd.DataFrame(data, columns=col_names)

        fname    = f"{SITE_ID}_HF_{t_start}_{t_end}.csv"
        out_path = os.path.join(OUT_DIR, fname)
        df.to_csv(out_path, index=False, float_format="%.4f", lineterminator="\n")

        site_csvs[site_type].append(out_path)

    print(f"  Wrote {n_blocks} CSV files.")

# ── zip by site type ───────────────────────────────────────────────────────────
print("\nZipping …")
for site_type, csv_list in site_csvs.items():
    zip_name = f"{SITE_ID}_{site_type}_HF.zip"
    zip_path = os.path.join(OUT_DIR, zip_name)
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for csv_path in sorted(csv_list):
            zf.write(csv_path, os.path.basename(csv_path))
    for csv_path in csv_list:
        os.remove(csv_path)
    print(f"  → {zip_name}  ({len(csv_list)} files)")

print("\nDone.")
