"""generate_ameriflux.py – Create AmeriFlux BASE half-hourly CSV from UTESpac output.

Format requirements (ameriflux.lbl.gov/half-hourly-hourly-data-upload-format/):
  - ASCII CSV, comma-delimited, period as decimal separator
  - First two columns: TIMESTAMP_START, TIMESTAMP_END (YYYYMMDDHHMM, local standard time)
  - Single header row of variable names, no units row
  - Missing value: -9999
  - File name: <SITE_ID>_HH_<START>_<END>.csv

Usage::

    python generate_ameriflux.py
"""

import os
import glob
import pickle
import warnings
import numpy as np
import pandas as pd
from datetime import datetime, timedelta

# ── configuration ────────────────────────────────────────────────────────────

SITE_ID   = "US-xFM"   # replace with official AmeriFlux site ID when registered

_cwd = os.getcwd()
if "diane_wt" in _cwd:
    ROOT_PY   = "/Users/diane_wt/Library/CloudStorage/Box-Box/Diane/code/UTESpac_Python"
    ROOT_1MIN = ("/Users/diane_wt/Library/CloudStorage/Box-Box/Lab Library/"
                 "French Meadows/Summer2025/data/processed/FM_DOL_1min")
elif "admin-dianew68" in _cwd:
    ROOT_PY   = "C:/Users/admin-dianew68/Box/Diane/code/UTESpac_Python"
    ROOT_1MIN = ("C:/Users/admin-dianew68/Box/Lab Library/"
                 "French Meadows/Summer2025/data/processed/FM_DOL_1min")
else:
    ROOT_PY   = os.path.dirname(os.path.abspath(__file__))
    ROOT_1MIN = ""

OUT_DIR   = os.path.join(ROOT_PY, "ameriflux_output")

PF_TYPE   = "GPF"   # "LPF" or "GPF"
AVG_PER   = 30      # minutes
MISSING   = -9999.0
GAMMA     = 0.0098  # dry adiabatic lapse rate [K/m]
Z_REF     = 4.42    # reference height [m] used in UTESpac for both sites

# EC heights in ascending order (V=1 = lowest); matched to site output
IRGA_HEIGHTS = [4.42, 6.35, 13.94, 32.18]   # siteIRGA
GILL_HEIGHT  = 51.5                           # siteGill
EC_HEIGHTS   = IRGA_HEIGHTS + [GILL_HEIGHT]  # V=1..5

# T/RH and radiation heights (ascending); heights match column names in 1-min file
HMP_HEIGHTS  = [6.35, 15.0, 30.0, 51.5]     # V=1..4
HMP_COLS     = {6.35: ("Temp_6.35", "RH_6.35"),
                15.0: ("Temp_15",   "RH_15"),
                30.0: ("Temp_30",   "RH_30"),
                51.5: ("Temp_51.5", "RH_51.5")}
RAD_HEIGHTS  = [4.22, 51.5]                  # V=1..2
RAD_COLS     = {4.22:  ("SWIN_4.22",  "SWOUT_4.22",  "LWIN_4.22",  "LWOUT_4.22",  "Rn_4.22"),
                51.5:  ("SWIN_51.5",  "SWOUT_51.5",  "LWIN_51.5",  "LWOUT_51.5",  "Rn_51.5")}

# ── helpers ───────────────────────────────────────────────────────────────────

MATLAB_EPOCH = datetime(1970, 1, 1)
MATLAB_OFFSET = 719529  # datenum(1970-01-01)

def matlab_to_dt(serial):
    """Convert MATLAB serial date float → Python datetime, rounded to nearest minute."""
    raw = datetime.fromordinal(int(serial)) + timedelta(days=serial % 1) \
          - timedelta(days=366)
    # Round to nearest minute to remove MATLAB floating-point jitter (~3 µs)
    seconds = raw.second + raw.microsecond / 1e6
    if seconds >= 30:
        raw = raw.replace(second=0, microsecond=0) + timedelta(minutes=1)
    else:
        raw = raw.replace(second=0, microsecond=0)
    return raw

def dt_to_ameriflux(dt):
    """Format datetime as AmeriFlux YYYYMMDDHHMM string."""
    return dt.strftime("%Y%m%d%H%M")

def _hstr(h):
    """Format height as UTESpac uses: 4.42 → '4.42', 51.5 → '51.5'."""
    return f"{h:g}"

def get_col(data_arr, header_list, pattern):
    """Return column from data_arr whose header contains pattern, or NaN array."""
    for i, h in enumerate(header_list):
        if pattern in h:
            return data_arr[:, i].astype(float)
    return np.full(data_arr.shape[0], np.nan)

def load_pkl_site(site_folder, pf_type, avg_per=30):
    """Load and concatenate all matching 30-min avg pkl files for a site."""
    pattern = os.path.join(site_folder, "output",
                           f"*_30minAvg_{pf_type}_LinDet_*.pkl")
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No pkl files matching {pattern}")
    parts = []
    for fp in files:
        with open(fp, "rb") as f:
            parts.append(pickle.load(f))
    # Concatenate arrays within each key
    combined = {}
    ref = parts[0]
    for key in ref:
        arrays = [p[key] for p in parts if key in p]
        if isinstance(arrays[0], np.ndarray) and arrays[0].ndim >= 1:
            combined[key] = np.concatenate(arrays, axis=0)
        else:
            combined[key] = arrays[0]   # headers, scalars — take first
    return combined

def to_series(data_arr, header_list, pattern, t_index):
    """Extract a column as a pd.Series indexed by t_index."""
    col = get_col(data_arr, header_list, pattern)
    if len(col) == len(t_index):
        return pd.Series(col, index=t_index)
    return pd.Series(np.nan, index=t_index)

# ── load UTESpac pkl data ─────────────────────────────────────────────────────

print("Loading IRGA data…")
irga = load_pkl_site(os.path.join(ROOT_PY, "siteIRGA"), PF_TYPE)

print("Loading Gill data…")
gill = load_pkl_site(os.path.join(ROOT_PY, "siteGill"), PF_TYPE)

# ── build continuous 30-min timestamp grid ────────────────────────────────────

# Timestamps are MATLAB serial dates at end of period (stored in H[:,0])
irga_times = irga["H"][:, 0]
gill_times = gill["H"][:, 0]

# Convert to datetime
irga_dt = [matlab_to_dt(t) for t in irga_times]
gill_dt = [matlab_to_dt(t) for t in gill_times]

# Build full continuous grid covering both sites
t_start = min(irga_dt[0], gill_dt[0])
t_end   = max(irga_dt[-1], gill_dt[-1])

# Snap start to 30-min boundary
t_start_end = t_start   # first timestamp is already end-of-period
t_end_end   = t_end

# Create end-of-period index
freq  = f"{AVG_PER}min"
t_idx = pd.date_range(start=t_start_end, end=t_end_end, freq=freq)

# Create lookup series: matlab serial → position in each site's arrays
irga_idx = pd.Series(np.arange(len(irga_dt)),
                     index=pd.DatetimeIndex(irga_dt)).reindex(t_idx)
gill_idx = pd.Series(np.arange(len(gill_dt)),
                     index=pd.DatetimeIndex(gill_dt)).reindex(t_idx)

def site_col(site_data, header_key, hdr_key, pattern, t_idx_series):
    """Extract a column from site data, aligned to t_idx."""
    arr     = site_data[header_key]
    headers = site_data[hdr_key]
    col     = get_col(arr, headers, pattern)
    valid   = t_idx_series.dropna().astype(int)
    out     = np.full(len(t_idx_series), np.nan)
    for pos_t, pos_s in zip(valid.index, valid.values):
        loc = t_idx_series.index.get_loc(pos_t)
        out[loc] = col[pos_s] if pos_s < len(col) else np.nan
    return out

# ── pre-extract PA (kPa) for use in EC loop molar density calculations ───────

irga_hdr_raw  = irga.get("FMDOL_20HzHeader", [])
irga_table_pa = irga.get("FMDOL_20Hz")
names_row_pa  = irga_hdr_raw[0] if (irga_hdr_raw and isinstance(irga_hdr_raw[0], list)) \
                else irga_hdr_raw
pa_col_idx_ec = next((i for i, h in enumerate(names_row_pa)
                      if "Pressure" in str(h)), None)
pa_timeseries = np.full(len(t_idx), np.nan)
if pa_col_idx_ec is not None and irga_table_pa is not None:
    pa_raw_ec = irga_table_pa[:, pa_col_idx_ec].astype(float)
    med_ec = np.nanmedian(pa_raw_ec)
    if med_ec > 50000:
        pa_raw_ec = pa_raw_ec / 1000.0
    elif med_ec > 200:
        pa_raw_ec = pa_raw_ec / 10.0
    pa_ser_ec = pd.Series(pa_raw_ec, index=pd.DatetimeIndex(irga_dt))
    pa_timeseries = pa_ser_ec.reindex(t_idx).values.astype(float)

# ── build output DataFrame ────────────────────────────────────────────────────

df = pd.DataFrame(index=t_idx)

# Timestamps
df["TIMESTAMP_END"]   = [dt_to_ameriflux(t) for t in t_idx]
df["TIMESTAMP_START"] = [dt_to_ameriflux(t - timedelta(minutes=AVG_PER)) for t in t_idx]

# ── turbulent fluxes ──────────────────────────────────────────────────────────

for v_idx, height in enumerate(EC_HEIGHTS, start=1):
    hn  = _hstr(height)
    src = gill if height == GILL_HEIGHT else irga
    hi  = gill_idx if height == GILL_HEIGHT else irga_idx

    def _col(arr_key, hdr_key, pattern):
        return site_col(src, arr_key, hdr_key, pattern, hi)

    # ---- wind sector flag (1 = tower-disturbed, 0 = clean) ----
    wd_flag = _col("spdAndDir", "spdAndDirHeader", f"{hn}m flag")
    bad_sector = np.isfinite(wd_flag) & (wd_flag != 0)   # True where wind is bad

    # ---- density / thermodynamic terms ----
    # Use specificHum (level-specific) when available, else global H[:,1]/[:,2]
    spec_hum_h = f"{height} m: "   # prefix used in specificHum headers
    if "specificHum" in src:
        sh      = src["specificHum"]
        sh_hdr  = src["specificHumHeader"]
        # For 4.42 m use 6.35 m specificHum row
        lookup_h = 6.35 if height == 4.42 else height
        lhn      = f"{lookup_h} m: "
        rho_arr  = get_col(sh, sh_hdr, f"{lhn}rho_airmoistAvg(kg/m^3)")
        r_arr    = get_col(sh, sh_hdr, f"{lhn}rAvg(g/kg)")
        vt_arr   = get_col(sh, sh_hdr, f"{lhn}virtualThetaAvg(K)")
        # align to t_idx
        rho_col  = np.full(len(t_idx), np.nan)
        r_col    = np.full(len(t_idx), np.nan)
        vt_col   = np.full(len(t_idx), np.nan)
        for pos_t, pos_s in zip(hi.dropna().index, hi.dropna().astype(int).values):
            loc = t_idx.get_loc(pos_t)
            if pos_s < len(rho_arr):
                rho_col[loc] = rho_arr[pos_s]
                r_col[loc]   = r_arr[pos_s]
                vt_col[loc]  = vt_arr[pos_s]
        cp_col = 1004.67 * (1.0 + 0.84 * r_col / 1000.0)
        Lv_col = (2.501 - 0.00237 * (vt_col - 273.15)) * 1e6  # J/kg
    else:
        # Fallback: use global rho and cp stored in H matrix (cols 1 & 2)
        h_arr   = src["H"]
        rho_raw = h_arr[:, 1]
        cp_raw  = h_arr[:, 2]
        lv_raw  = _col("LHflux", "LHfluxHeader", f"{hn}m Lv(J/g)")
        rho_col = np.full(len(t_idx), np.nan)
        cp_col  = np.full(len(t_idx), np.nan)
        Lv_col  = np.full(len(t_idx), np.nan)
        for pos_t, pos_s in zip(hi.dropna().index, hi.dropna().astype(int).values):
            loc = t_idx.get_loc(pos_t)
            if pos_s < len(rho_raw):
                rho_col[loc] = rho_raw[pos_s]
                cp_col[loc]  = cp_raw[pos_s]
                Lv_col[loc]  = lv_raw[loc] * 1000.0  # J/g → J/kg
        vt_col = np.full(len(t_idx), np.nan)  # not available without specificHum
        r_col  = np.full(len(t_idx), np.nan)

    # ---- H: sensible heat flux [W m⁻²] ----
    thv_wPF  = _col("H",      "Hheader",      f"{hn}m son:Theta_v'wPF'")
    H_flux   = rho_col * cp_col * thv_wPF
    df[f"H_1_{v_idx}_1"] = H_flux

    # ---- LE: latent heat flux [W m⁻²] ----
    q_wpl    = _col("LHflux", "LHfluxHeader", f"{hn}m wPF'q_WPL'(m/s kg/m3)")
    LE_flux  = Lv_col * q_wpl
    df[f"LE_1_{v_idx}_1"] = LE_flux

    # ---- FC: CO2 flux [µmol m⁻² s⁻¹] ----
    # wPF'CO2_WPL' is in kg/(m²·s); ÷44.01 g/mol ×1000 g/kg ×1e6 µmol/mol → µmol/(m²·s)
    co2_wpl  = _col("CO2flux", "CO2fluxHeader", f"{hn}m: wPF'CO2_WPL'(m/s kg/m^3)")
    df[f"FC_1_{v_idx}_1"] = co2_wpl / 44.01 * 1000.0 * 1e6

    # ---- USTAR [m s⁻¹] ----
    tau_pf   = _col("tau", "tauHeader", f"{hn}m :sqrt(uPF'wPF'^2+vPF'wPF'^2)")
    df[f"USTAR_1_{v_idx}_1"] = np.sqrt(np.abs(tau_pf)) * np.sign(tau_pf)

    # ---- WS: mean streamwise wind speed [m s⁻¹] after PF + yaw rotation ----
    # rotatedSonic stores 30-min mean of (PF+yaw)-rotated u; mean(vPF)≈0 by yaw
    ws_pf = _col("rotatedSonic", "rotatedSonicHeader", f"{hn}m:u")
    df[f"WS_1_{v_idx}_1"] = ws_pf

    # ---- WD: wind direction [°] — geographic, uses unrotated frame ----
    wd = _col("spdAndDir", "spdAndDirHeader", f"{hn}m direction")
    df[f"WD_1_{v_idx}_1"] = wd

    # ---- MO_LENGTH: Obukhov length [m] ----
    L_col = _col("L", "Lheader", f"{hn}m L:")
    df[f"MO_LENGTH_1_{v_idx}_1"] = L_col

    # ---- TKE [m² s⁻²]: 0.5*(σ_uPF² + σ_vPF² + σ_wPF²) after PF rotation ----
    su = _col("sigma", "sigmaHeader", f"{hn}m :sigma_uPF")
    sv = _col("sigma", "sigmaHeader", f"{hn}m :sigma_vPF")
    sw = _col("sigma", "sigmaHeader", f"{hn}m :sigma_wPF")
    df[f"TKE_1_{v_idx}_1"] = 0.5 * (su**2 + sv**2 + sw**2)

    # ---- TAU: Reynolds stress [Pa = kg m⁻¹ s⁻²] ----
    df[f"TAU_1_{v_idx}_1"] = rho_col * tau_pf

    # ---- Velocity standard deviations [m s⁻¹] ----
    df[f"U_SIGMA_1_{v_idx}_1"] = su
    df[f"V_SIGMA_1_{v_idx}_1"] = sv
    df[f"W_SIGMA_1_{v_idx}_1"] = sw

    # ---- T_SONIC_SIGMA: sonic temperature std dev [°C] ----
    ts_sig = _col("sigma", "sigmaHeader", f"{hn}m :sigma_Tson")
    df[f"T_SONIC_SIGMA_1_{v_idx}_1"] = ts_sig

    # ---- ZL: stability parameter z/L [nondimensional] ----
    with np.errstate(divide="ignore", invalid="ignore"):
        zl = np.where(np.isfinite(L_col) & (L_col != 0.0),
                      height / L_col, np.nan)
    df[f"ZL_1_{v_idx}_1"] = zl

    # ---- FH2O: H2O flux [mmol H2O m⁻² s⁻¹] ----
    df[f"FH2O_1_{v_idx}_1"] = q_wpl * 1e6 / 18.015

    # ---- H2O: mean mole fraction [mmol H2O mol⁻¹ air] ----
    # mixing ratio r [g/kg] → mmol/mol; only valid after UTESpac re-run with specificHum fix
    df[f"H2O_1_{v_idx}_1"] = r_col * 28.97 / 18.015  # NaN when specificHum absent

    # ---- Molar air density for gas-phase sigma conversions ----
    PA_Pa   = pa_timeseries * 1000.0   # kPa → Pa
    with np.errstate(divide="ignore", invalid="ignore"):
        rho_mol = np.where(
            np.isfinite(PA_Pa) & np.isfinite(vt_col) & (vt_col > 200),
            PA_Pa / (8.314 * vt_col), np.nan)

    # ---- CO2_SIGMA [µmol CO2 mol⁻¹ air] ----
    # sigma_CO2 [mg/m³] ÷ 44010 [mg/mol] ÷ rho_mol [mol/m³] × 1e6
    sig_co2 = _col("sigma", "sigmaHeader", f"{hn}m :sigma_CO2")  # mg/m³
    with np.errstate(divide="ignore", invalid="ignore"):
        df[f"CO2_SIGMA_1_{v_idx}_1"] = sig_co2 / (44.01 * 1000.0) / rho_mol * 1e6

    # ---- H2O_SIGMA [mmol H2O mol⁻¹ air] ----
    # sigma_H2O [g/m³] ÷ 18.015 [g/mol] ÷ rho_mol [mol/m³] × 1e3
    sig_h2o = _col("sigma", "sigmaHeader", f"{hn}m :sigma_H2O")  # g/m³
    with np.errstate(divide="ignore", invalid="ignore"):
        df[f"H2O_SIGMA_1_{v_idx}_1"] = sig_h2o / 18.015 / rho_mol * 1e3

    # ---- T_SONIC: raw sonic temperature [°C] ----
    # theta_v_son = T_son + Gamma*(height - zRef), so invert to recover T_son
    theta_v = _col("derivedT", "derivedTheader", f"{height} m: theta_v_son")
    df[f"T_SONIC_1_{v_idx}_1"] = theta_v - GAMMA * (height - Z_REF)

    # ---- CO2: mean CO2 mole fraction [µmol mol⁻¹] (total, not dry) ----
    co2_ppm = _col("CO2flux", "CO2fluxHeader", f"{hn}m: CO2 (ppm)")
    df[f"CO2_1_{v_idx}_1"] = co2_ppm

    # ---- WD_FILTER: tower-obstruction wind direction flag ----
    # 0 = clean sector, 1 = wind from tower direction (disturbed)
    df[f"WD_FILTER_1_{v_idx}_1"] = wd_flag.astype(float)

    # ---- Mask wind-disturbed flux periods to NaN (→ -9999 in output) ----
    # AmeriFlux BASE convention: set knowingly bad data to missing rather than
    # relying on users to interpret a separate flag column.
    flux_cols = [
        f"H_1_{v_idx}_1",    f"LE_1_{v_idx}_1",  f"FC_1_{v_idx}_1",
        f"USTAR_1_{v_idx}_1", f"TAU_1_{v_idx}_1",
    ]
    for col in flux_cols:
        if col in df.columns:
            vals = df[col].values.astype(float)
            vals[bad_sector] = np.nan
            df[col] = vals

    # ---- SSITC quality flags (0=high, 1=moderate, 2=poor) ----
    for flag_var in ("TAU", "H", "LE", "FC"):
        flag_col = site_col(src, "fluxQC", "fluxQCHeader",
                            f"{hn}m:{flag_var}_SSITC_TEST", hi)
        df[f"{flag_var}_SSITC_TEST_1_{v_idx}_1"] = flag_col

# ── slow meteorology from 1-min data ─────────────────────────────────────────

print("Loading 1-min met data…")
min1_files = sorted(glob.glob(os.path.join(ROOT_1MIN, "FM_DOL_1min_*.txt")))
if min1_files:
    min1_frames = []
    for fp in min1_files:
        try:
            df_tmp = pd.read_csv(fp, header=0)
            min1_frames.append(df_tmp)
        except Exception:
            pass
    min1_all = pd.concat(min1_frames, ignore_index=True)

    # Build datetime index from year, day-of-year, HHMM, second columns
    def _1min_dt(row):
        hm  = int(row["HM"])
        h   = hm // 100
        m   = hm % 100
        try:
            return datetime(int(row["year"]), 1, 1) \
                   + timedelta(days=int(row["day"]) - 1,
                               hours=h, minutes=m,
                               seconds=int(row.get("second", 0)))
        except Exception:
            return pd.NaT

    min1_all["dt"] = min1_all.apply(_1min_dt, axis=1)
    min1_all = min1_all.dropna(subset=["dt"]).set_index("dt").sort_index()

    # Resample 1-min → 30-min mean; closed='right' label='right' means each
    # bin covers (T-30min, T] and is labelled at T — matching UTESpac timestamps.
    min1_30 = min1_all.resample(f"{AVG_PER}min", closed="right",
                                label="right").mean()
    # Reindex to t_idx with 1-minute tolerance to absorb any sub-minute offset
    min1_30.index = pd.DatetimeIndex(min1_30.index).round(f"{AVG_PER}min")

    def _met(col_name, idx=t_idx):
        if col_name not in min1_30.columns:
            warnings.warn(f"Column '{col_name}' not found in 1-min data.")
            return np.full(len(idx), np.nan)
        return min1_30[col_name].reindex(idx).values.astype(float)

    # ---- TA [°C] ----
    for v_idx, height in enumerate(HMP_HEIGHTS, start=1):
        t_col, rh_col = HMP_COLS[height]
        ta  = _met(t_col)
        rh  = _met(rh_col)
        df[f"TA_1_{v_idx}_1"]  = ta
        df[f"RH_1_{v_idx}_1"]  = rh
        # VPD [hPa] = es × (1 - RH/100)
        es  = 6.1078 * np.exp(17.27 * ta / (ta + 237.3))
        vpd = es * (1.0 - rh / 100.0)
        df[f"VPD_1_{v_idx}_1"] = np.where(np.isfinite(vpd), vpd, np.nan)

    # ---- PA from UTESpac averaged table [kPa] ----
    df["PA_1_1_1"] = pa_timeseries

    # ---- Radiation [W m⁻²] ----
    for v_idx, height in enumerate(RAD_HEIGHTS, start=1):
        sw_in, sw_out, lw_in, lw_out, rn = RAD_COLS[height]
        df[f"SW_IN_1_{v_idx}_1"]   = _met(sw_in)
        df[f"SW_OUT_1_{v_idx}_1"]  = _met(sw_out)
        df[f"LW_IN_1_{v_idx}_1"]   = _met(lw_in)
        df[f"LW_OUT_1_{v_idx}_1"]  = _met(lw_out)
        df[f"NETRAD_1_{v_idx}_1"]  = _met(rn)

    # ---- ALB: albedo [%] from 30-min SW (filter: SW_IN > 5 W/m²) ----
    for v_idx, height in enumerate(RAD_HEIGHTS, start=1):
        sw_in_v  = df[f"SW_IN_1_{v_idx}_1"].values.astype(float)
        sw_out_v = df[f"SW_OUT_1_{v_idx}_1"].values.astype(float)
        with np.errstate(divide="ignore", invalid="ignore"):
            alb = np.where(
                (sw_in_v != MISSING) & (sw_out_v != MISSING) & (sw_in_v > 5.0),
                100.0 * sw_out_v / sw_in_v, np.nan)
        df[f"ALB_1_{v_idx}_1"] = alb
else:
    warnings.warn("No FM_DOL_1min files found — met/radiation columns will be NaN.")

# ── replace NaN with -9999 and enforce column order ──────────────────────────

# Put timestamps first
ts_cols  = ["TIMESTAMP_START", "TIMESTAMP_END"]
data_cols = [c for c in df.columns if c not in ts_cols]
df = df[ts_cols + data_cols]

# Replace NaN with -9999
df = df.fillna(MISSING)

# Round numeric columns to 6 significant figures to avoid floating-point noise
for col in data_cols:
    vals = df[col].values.astype(float)
    mask = vals != MISSING
    if mask.any():
        vals[mask] = np.round(vals[mask], 6)
    df[col] = vals

# ── write CSV ────────────────────────────────────────────────────────────────

os.makedirs(OUT_DIR, exist_ok=True)

t_file_start = df["TIMESTAMP_START"].iloc[0]
t_file_end   = df["TIMESTAMP_END"].iloc[-1]
fname = f"{SITE_ID}_HH_{t_file_start}_{t_file_end}.csv"
out_path = os.path.join(OUT_DIR, fname)

df.to_csv(out_path, index=False)
print(f"\nWrote {len(df)} half-hourly rows × {len(df.columns)} columns")
print(f"Output: {out_path}")

# ── summary of coverage ───────────────────────────────────────────────────────

print("\nVariable coverage (non-missing rows / total):")
total = len(df)
for col in data_cols:
    n_valid = (df[col].values != MISSING).sum()
    if n_valid < total:
        pct = 100 * n_valid / total
        print(f"  {col:35s}  {n_valid:4d}/{total}  ({pct:.0f}%)")
print("  (columns with 100% coverage not shown)")
