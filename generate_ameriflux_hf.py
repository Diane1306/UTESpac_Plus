"""generate_ameriflux_hf.py – AmeriFlux high-frequency CSV export from UTESpac raw pkl output.

Format references:
  https://ameriflux.lbl.gov/data/how-to-upload-data/uploading-high-frequency-data/
  https://ameriflux.lbl.gov/half-hourly-hourly-data-upload-format/

Output structure:
  One CSV per 30-min period: <SITE_ID>_HF_<YYYYMMDDHHMM>_<YYYYMMDDHHMM>.csv
  Flat zip archive (no subdirectories) per site type — AmeriFlux HF standard.
  Individual CSVs are removed after zipping.

Variable convention (EddyPro / AmeriFlux standard; V = 1 = lowest height):
  TIMESTAMP         – YYYYMMDDHHMMSS.cc  (local standard time, no DST; cc = centiseconds)
  U_1_{V}_1         – planar-fit u wind component           [m s⁻¹]
  V_1_{V}_1         – planar-fit v wind component           [m s⁻¹]  (~0 after PF rotation)
  W_1_{V}_1         – planar-fit w wind component           [m s⁻¹]  (~0 after PF rotation)
  T_SONIC_1_{V}_1   – sonic temperature                     [K]
                         (sonTs [°C] + 273.15 — EddyPro / AmeriFlux convention is Kelvin)
  H2O_1_{V}_1       – water vapour dry mole fraction        [mmol mol⁻¹]
                         converted from rhov [g m⁻³] via ideal gas law using PA and T_SONIC
  CO2_1_{V}_1       – CO2 dry mole fraction                 [µmol mol⁻¹]
                         converted from rhoCO2 [mg m⁻³] via ideal gas law
  PA_1_1_1          – atmospheric pressure (one barometer)  [Pa]
                         (P [kPa] × 1000 — AmeriFlux / EddyPro convention is Pa)

Mole-fraction conversion (ideal gas, dry-air basis):
  n_moist = PA [Pa] / (R × T [K])        total molar density [mol m⁻³]
  n_H2O   = rhov [g m⁻³] / 18.015       H2O molar density   [mol m⁻³]
  n_dry   = n_moist − n_H2O              dry-air molar density [mol m⁻³]
  H2O     = (n_H2O / n_dry) × 1000      [mmol mol⁻¹ dry air]
  CO2     = (rCO2 [mg m⁻³] / 44010) / n_dry × 1e6   [µmol mol⁻¹ dry air]

Unit notes:
  rhov:    g m⁻³  in both Python and MATLAB pkl
  rhoCO2:  mg m⁻³ in Python pkl  (MATLAB stores kg m⁻³; see fluxes.py comment)

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

R_GAS        = 8.314    # universal gas constant [J mol⁻¹ K⁻¹]
MW_H2O       = 18.015   # g mol⁻¹
MW_CO2       = 44.010   # g mol⁻¹
MATLAB_EPOCH = 719529.0 # MATLAB days to Unix epoch (1970-01-01)

# ── helpers ────────────────────────────────────────────────────────────────────

def matlab_to_unix(t_arr):
    return (np.asarray(t_arr) - MATLAB_EPOCH) * 86400.0


def make_timestamps(t_arr):
    """Vectorised MATLAB serial → 'YYYYMMDDHHMMSS.cc' pandas Series."""
    dti = pd.to_datetime(matlab_to_unix(t_arr), unit="s")
    cs  = (dti.microsecond // 10000).astype(str).str.zfill(2)
    return dti.strftime("%Y%m%d%H%M%S") + "." + cs


def file_ts(t_val, offset_s=0.0):
    """Single MATLAB serial → YYYYMMDDHHMM for file naming."""
    return pd.Timestamp(matlab_to_unix(t_val) + offset_s, unit="s").strftime("%Y%m%d%H%M")


def safe(arr):
    """Replace NaN / ±Inf with MISSING (-9999)."""
    return np.where(np.isfinite(arr), arr, MISSING)


def density_to_molefrac(rhov_gm3, rco2_mgm3, T_K, PA_Pa):
    """Convert gas densities to dry mole fractions via ideal gas law.

    Parameters
    ----------
    rhov_gm3   : ndarray  water vapour density  [g m⁻³]
    rco2_mgm3  : ndarray  CO2 density           [mg m⁻³]
    T_K        : ndarray  sonic temperature      [K]
    PA_Pa      : ndarray  air pressure           [Pa]

    Returns
    -------
    h2o_mmol   : ndarray  H2O dry mole fraction  [mmol mol⁻¹]
    co2_umol   : ndarray  CO2 dry mole fraction  [µmol mol⁻¹]
    """
    valid = np.isfinite(T_K) & (T_K > 200) & np.isfinite(PA_Pa) & (PA_Pa > 1e4)
    n_moist = np.full_like(T_K, np.nan)
    n_moist[valid] = PA_Pa[valid] / (R_GAS * T_K[valid])   # mol m⁻³

    n_H2O = rhov_gm3 / MW_H2O                              # mol m⁻³
    n_dry = n_moist - n_H2O
    # Guard against non-positive n_dry
    n_dry = np.where(n_dry > 0.1, n_dry, np.nan)

    h2o_mmol = safe(n_H2O / n_dry * 1e3)                   # mmol mol⁻¹
    co2_umol = safe((rco2_mgm3 / (MW_CO2 * 1e3)) / n_dry * 1e6)  # µmol mol⁻¹
    return h2o_mmol, co2_umol


# ── main ───────────────────────────────────────────────────────────────────────

os.makedirs(OUT_DIR, exist_ok=True)

raw_pkls = sorted(glob.glob(
    os.path.join(ROOT_PY, "site*", "output", f"*_raw_{PF_TYPE}_*.pkl")))

if not raw_pkls:
    print(f"No raw {PF_TYPE} pkl files found under {ROOT_PY}/site*/output/")
    raise SystemExit(1)

print(f"Found {len(raw_pkls)} raw {PF_TYPE} pkl file(s).")

site_csvs = defaultdict(list)

for pkl_path in raw_pkls:
    site_dir  = os.path.basename(os.path.dirname(os.path.dirname(pkl_path)))
    site_type = re.match(r"(site[A-Za-z]+)\d", site_dir)
    site_type = site_type.group(1) if site_type else site_dir

    print(f"\nProcessing {os.path.basename(pkl_path)}  [{site_dir}] …")

    with open(pkl_path, "rb") as fh:
        raw = pickle.load(fh)

    t     = raw["t"]
    z     = raw["z"]
    uPF   = raw["uPF"]
    vPF   = raw["vPF"]
    wPF   = raw["wPF"]
    Ts_C  = raw["sonTs"]            # °C
    P_kPa = raw["P"][:, 1]         # kPa
    rhov  = raw["rhov"]             # g m⁻³
    rCO2  = raw["rhoCO2"]           # mg m⁻³  (Python convention)

    n_sonics = len(z)
    N        = len(t)

    dt_s  = float(np.nanmedian(np.diff(t[:min(2000, N)]))) * 86400.0
    hz    = int(round(1.0 / dt_s))
    block = hz * 30 * 60
    n_blocks = N // block

    print(f"  {hz} Hz  ·  {n_sonics} sonic(s) at {list(z)} m  ·  {n_blocks} × 30-min files")

    sorted_z = sorted(z)
    v_of     = {zi: sorted_z.index(zi) + 1 for zi in z}

    # ── convert to AmeriFlux units (full 48-h array, vectorised) ─────────────
    T_K   = safe(Ts_C) + 273.15     # °C → K  (use first sonic column for conversion)
    PA_Pa = safe(P_kPa) * 1e3       # kPa → Pa

    # Pre-compute mole fractions per sonic using per-height T and shared PA
    h2o_out = np.full_like(rhov, MISSING)
    co2_out = np.full_like(rCO2, MISSING)
    for si in range(n_sonics):
        T_si = safe(Ts_C[:, si]) + 273.15
        h2o_out[:, si], co2_out[:, si] = density_to_molefrac(
            safe(rhov[:, si]), safe(rCO2[:, si]), T_si, PA_Pa)

    u_out  = safe(uPF)
    v_out  = safe(vPF)
    w_out  = safe(wPF)
    ts_out = safe(Ts_C) + 273.15    # K

    # ── build column names ────────────────────────────────────────────────────
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

    # ── write one CSV per 30-min block ────────────────────────────────────────
    for bi in range(n_blocks):
        s0, s1 = bi * block, (bi + 1) * block
        t_blk  = t[s0:s1]

        t_start = file_ts(t_blk[0])
        t_end   = file_ts(t_blk[-1], offset_s=1.0 / hz)

        data = {"TIMESTAMP": make_timestamps(t_blk)}
        for si, zi in enumerate(z):
            v = v_of[zi]
            data[f"U_1_{v}_1"]       = u_out[s0:s1, si]
            data[f"V_1_{v}_1"]       = v_out[s0:s1, si]
            data[f"W_1_{v}_1"]       = w_out[s0:s1, si]
            data[f"T_SONIC_1_{v}_1"] = ts_out[s0:s1, si]
            data[f"H2O_1_{v}_1"]     = h2o_out[s0:s1, si]
            data[f"CO2_1_{v}_1"]     = co2_out[s0:s1, si]
        data["PA_1_1_1"] = PA_Pa[s0:s1]

        df = pd.DataFrame(data, columns=col_names)

        fname    = f"{SITE_ID}_HF_{t_start}_{t_end}.csv"
        out_path = os.path.join(OUT_DIR, fname)
        df.to_csv(out_path, index=False, float_format="%.4f", lineterminator="\n")
        site_csvs[site_type].append(out_path)

    print(f"  Wrote {n_blocks} CSV files.")

# ── zip by site type (flat — no subdirectories) ───────────────────────────────
print("\nZipping …")
for site_type, csv_list in site_csvs.items():
    zip_name = f"{SITE_ID}_{site_type}_HF.zip"
    zip_path = os.path.join(OUT_DIR, zip_name)
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for csv_path in sorted(csv_list):
            zf.write(csv_path, os.path.basename(csv_path))  # flat — no subdirs
    for csv_path in csv_list:
        os.remove(csv_path)
    print(f"  → {zip_name}  ({len(csv_list)} files)")

print("\nDone.")
