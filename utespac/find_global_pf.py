"""findGlobalPF – compute or load multi-sector global planar-fit coefficients."""

import os
import pickle
import re
from datetime import datetime, timedelta
from typing import Dict, Optional

import numpy as np

from .get_data import get_data
from .pf_coefficients import pf_coefficients
from .stp_dn import stp_dn

_MATLAB_EPOCH = 719529  # MATLAB datenum for 1970-01-01


def _matlab_to_datetime(serial):
    """Convert a MATLAB serial date float to a Python datetime."""
    return datetime(1970, 1, 1) + timedelta(days=float(serial) - _MATLAB_EPOCH)


def _compute_direction(u, v, bearing, manufact):
    """Compute wind direction [0, 360) from u/v, matching MATLAB findGlobalPF.

    Manufacturer corrections (matching sonicRotation / MATLAB):
      1 (Campbell)  : u, v as-is
      0 (RMYoung)   : swap u↔v, negate new v  → uDir=v, vDir=-u  (wait, MATLAB: uDir=v, vDir=u*-1)
      2 (Gill)      : negate both              → uDir=-u, vDir=-v
    """
    if manufact == 0:       # RMYoung
        u_dir, v_dir = v, u * -1
    elif manufact == 2:     # Gill WindmasterPro
        u_dir, v_dir = -u, -v
    else:                   # Campbell (default)
        u_dir, v_dir = u, v
    return np.mod(np.arctan2(-v_dir, u_dir) * 180.0 / np.pi + bearing, 360.0)


def _show_selection_figure(z, site_name, direction, t_dn, u_dn, v_dn, w_dn, spd_dn):
    """Create the two-panel figure (histogram + time series) and return (fig, ax1, ax2)."""
    try:
        import matplotlib
        matplotlib.use("TkAgg")   # use a GUI backend; falls back gracefully
    except Exception:
        pass
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8))
    fig.suptitle(f"{site_name}  {z} m", fontsize=13)

    # ── top: direction histogram ──────────────────────────────────────────────
    valid = direction[~np.isnan(direction)]
    ax1.hist(valid, bins=100, range=(0, 360), color="steelblue", edgecolor="none")
    ax1.set_xlim(0, 360)
    ax1.set_xlabel("Direction (°)")
    ax1.set_ylabel(f"{round(1440 / (len(t_dn) / max(len(valid), 1)), 0):.0f}-min occurrences"
                   if len(t_dn) > 0 else "Occurrences")
    ax1.set_title("Wind direction — enter bin boundaries below, then press Enter with no value")

    # ── bottom: time series ───────────────────────────────────────────────────
    try:
        dates = [_matlab_to_datetime(d) for d in t_dn]
    except Exception:
        dates = list(range(len(t_dn)))

    ax2.plot(dates, u_dn,   "b.", markersize=3, label="u")
    ax2.plot(dates, v_dn,   "g.", markersize=3, label="v")
    ax2.plot(dates, w_dn,   "r.", markersize=3, label="w")
    ax2.plot(dates, spd_dn, "c.", markersize=3, label="spd")
    ax2.legend(loc="upper right", fontsize=8)
    ax2.set_title("Time series — enter date barriers after selecting direction bins")
    if isinstance(dates[0], datetime):
        ax2.xaxis.set_major_formatter(mdates.DateFormatter("%m/%d"))
        fig.autofmt_xdate(rotation=30)
    ax2.autoscale(axis="x", tight=True)

    plt.tight_layout()
    plt.show(block=False)
    plt.pause(0.1)
    return fig, ax1, ax2


def find_global_pf(info: Dict, template: Dict, sensor_info: Dict) -> Dict:
    """Compute (or load) global planar-fit coefficients for each sonic height.

    For each height the user is shown a wind-direction histogram and a velocity
    time-series figure (matching MATLAB findGlobalPF), then prompted to:
      1. Choose direction-bin boundaries (drawn live as vertical lines).
      2. Choose date-range boundaries (drawn live on the time series).

    The resulting ``PFinfo`` dict is saved to ``<siteFolder>/PFinfo.pkl``.
    """
    pf_path = os.path.join(info["rootFolder"], info["siteFolder"], "PFinfo.pkl")

    if not info["PF"]["recalculateGlobalCoefficients"] and os.path.isfile(pf_path):
        with open(pf_path, "rb") as fh:
            return pickle.load(fh)

    print("Finding Global Planar Fit Coefficients (b0, b1, b2)")

    # ── load LPF-averaged data ────────────────────────────────────────────────
    all_sites = sorted(d for d in os.listdir(info["rootFolder"]) if d.startswith("site"))
    site_num  = all_sites.index(info["siteFolder"]) + 1
    site_name = info["siteFolder"][4:]  # strip leading "site"

    data = get_data(
        info["rootFolder"],
        site=site_num,
        avg_per=info["avgPer"],
        qualifier="LPF",
        rows=0,
    )

    # ── determine sonic heights ───────────────────────────────────────────────
    hdr    = data.get("spdAndDirHeader", [])
    z_vals = []
    for h in hdr:
        nums = re.findall(r"[\d.]+", str(h))
        if nums:
            z_vals.append(float(nums[0]))
    z_vals = sorted(set(z_vals[1:]), reverse=True)  # skip timestamp at index 0

    table_names = data.get("tableNames", [])
    pf_info: Dict = {}

    for ii, z in enumerate(z_vals):
        skip = input(f"\nSkip {z} m? (1=process, 0=skip): ").strip()
        if skip == "0":
            if os.path.isfile(pf_path):
                with open(pf_path, "rb") as fh:
                    old = pickle.load(fh)
                cm_key = f"cm_{round(z * 100)}"
                if cm_key in old:
                    pf_info[cm_key] = old[cm_key]
                    continue
            print("  No existing PFinfo.pkl found — cannot skip this level.")
            skip = "1"

        varname_u = template["u"].replace("*", str(z))
        varname_v = template["v"].replace("*", str(z))
        varname_w = template["w"].replace("*", str(z))

        # ── find sonic table and columns ──────────────────────────────────────
        local_table = u_col = v_col = w_col = None
        for tname in table_names:
            hdr_key = f"{tname}Header"
            if hdr_key not in data:
                continue
            names = data[hdr_key][0] if isinstance(data[hdr_key][0], list) else data[hdr_key]
            if varname_u in names:
                local_table = tname
                u_col = names.index(varname_u)
                v_col = names.index(varname_v) if varname_v in names else None
                w_col = names.index(varname_w) if varname_w in names else None
                break

        if local_table is None or u_col is None:
            print(f"  Could not find sonic data at {z} m — skipping.")
            continue

        n_periods = data[local_table].shape[0]

        # ── quality flags (matching MATLAB findGlobalPF lines 69-89) ─────────
        def _flag_col(flag_key, col):
            arr = data.get(flag_key)
            if arr is None:
                return np.zeros(n_periods, dtype=bool)
            a = np.asarray(arr)
            return a[:, col].astype(bool) if a.ndim == 2 and col < a.shape[1] \
                   else np.zeros(n_periods, dtype=bool)

        nan_u   = _flag_col(f"{local_table}NanFlag",   u_col)
        nan_v   = _flag_col(f"{local_table}NanFlag",   v_col if v_col is not None else u_col)
        nan_w   = _flag_col(f"{local_table}NanFlag",   w_col if w_col is not None else u_col)
        spike_u = _flag_col(f"{local_table}SpikeFlag", u_col)
        spike_v = _flag_col(f"{local_table}SpikeFlag", v_col if v_col is not None else u_col)
        spike_w = _flag_col(f"{local_table}SpikeFlag", w_col if w_col is not None else u_col)

        wind_flag_col = next(
            (j for j, h in enumerate(hdr) if str(h).startswith(f"{z}m flag")), None)
        wind_flag = data["spdAndDir"][:, wind_flag_col].astype(bool) \
                    if wind_flag_col is not None else np.zeros(n_periods, dtype=bool)

        # sonic diagnostic flag
        varname_diag = template.get("sonDiagnostic", "diagnostic_*").replace("*", str(z))
        tbl_hdr = data.get(f"{local_table}Header", [])
        if tbl_hdr and isinstance(tbl_hdr[0], list):
            tbl_hdr = tbl_hdr[0]
        diag_col = tbl_hdr.index(varname_diag) if varname_diag in tbl_hdr else None
        if diag_col is not None:
            diag_vals = data[local_table][:, diag_col].copy()
            diag_vals = np.where(np.isnan(diag_vals), 0.0, diag_vals)
            diag_limit = info.get("diagnosticTest", {}).get("meanSonicDiagnosticLimit", 50)
            diag_vals[diag_vals < diag_limit] = 0.0
            diag_flag = diag_vals.astype(bool)
        else:
            diag_flag = np.zeros(n_periods, dtype=bool)

        total_flag = nan_u | spike_u | nan_v | spike_v | nan_w | spike_w \
                     | wind_flag | diag_flag

        # ── extract flagged arrays ────────────────────────────────────────────
        u_raw = data[local_table][:, u_col].copy()
        v_raw = data[local_table][:, v_col].copy() if v_col is not None else np.full(n_periods, np.nan)
        w_raw = data[local_table][:, w_col].copy() if w_col is not None else np.full(n_periods, np.nan)
        t_raw = data[local_table][:, 0].copy()

        u_raw[total_flag] = np.nan
        v_raw[total_flag] = np.nan
        w_raw[total_flag] = np.nan

        # ── wind speed column (for speed filter) ─────────────────────────────
        spd_col_idx = next((j for j, h in enumerate(hdr) if h == f"{z}m speed"), None)
        spd_raw = data["spdAndDir"][:, spd_col_idx].copy() if spd_col_idx is not None \
                  else np.full(n_periods, np.nan)
        spd_raw[total_flag] = np.nan

        # ── speed filter ──────────────────────────────────────────────────────
        max_w = info["PF"]["globalCalcMaxWind"]
        min_w = info["PF"]["globalCalcMinWind"]
        speed_mask = (spd_raw > max_w) | (spd_raw < min_w)
        u_raw[speed_mask] = v_raw[speed_mask] = w_raw[speed_mask] = np.nan

        # ── step down to PF averaging period ─────────────────────────────────
        factor = max(1, round(info["PF"]["avgPer"] / info["avgPer"]))
        u_dn   = stp_dn(u_raw.reshape(-1, 1), factor)[:, 0]
        v_dn   = stp_dn(v_raw.reshape(-1, 1), factor)[:, 0]
        w_dn   = stp_dn(w_raw.reshape(-1, 1), factor)[:, 0]
        t_dn   = stp_dn(t_raw.reshape(-1, 1), factor)[:, 0]
        spd_dn = stp_dn(spd_raw.reshape(-1, 1), factor)[:, 0]

        # ── compute direction from u/v with manufacturer correction ──────────
        # sensor_info["u"] columns: [table_idx, col_idx, height_m, bearing_deg, manufact_code]
        bearing  = float(sensor_info["u"][ii, 3]) if sensor_info["u"].shape[1] > 3 else 0.0
        manufact = int(sensor_info["u"][ii, 4])   if sensor_info["u"].shape[1] > 4 else 1
        direction = _compute_direction(u_dn, v_dn, bearing, manufact)

        # ── show figure ───────────────────────────────────────────────────────
        try:
            fig, ax1, ax2 = _show_selection_figure(
                z, site_name, direction, t_dn, u_dn, v_dn, w_dn, spd_dn)
            has_figure = True
        except Exception as exc:
            print(f"  (Could not show figure: {exc})")
            has_figure = False

        # ── direction bin selection ───────────────────────────────────────────
        print(f"\nSelect direction bin boundaries for {z} m.")
        print("  Enter one value per prompt, press Enter with no value when done.")
        bins = []
        y_lim = ax1.get_ylim() if has_figure else (0, 1)
        while True:
            val = input("  Bin boundary (° 0–360, or Enter to finish): ").strip()
            if not val:
                break
            try:
                b = float(val)
            except ValueError:
                continue
            if 0 <= b <= 360:
                bins.append(b)
                if has_figure:
                    import matplotlib.pyplot as plt
                    ax1.plot([b, b], [y_lim[0], y_lim[1]], "g--", linewidth=2)
                    plt.draw()
                    plt.pause(0.05)
        bins = sorted(bins)

        # ── date barrier selection ────────────────────────────────────────────
        all_days = np.unique(np.floor(t_dn[~np.isnan(t_dn)]))
        print(f"\nData spans {len(all_days)} day(s).")
        print("  Optionally split data by date for different PF calculations.")
        print("  Enter MATLAB serial dates (e.g. 739776), press Enter to skip.")
        date_barriers = []
        if has_figure:
            import matplotlib.pyplot as plt
            ax2_ylim = ax2.get_ylim()
        while True:
            val = input("  Date barrier (MATLAB serial date, or Enter to finish): ").strip()
            if not val:
                break
            try:
                db = float(val)
            except ValueError:
                continue
            if all_days[0] <= db <= all_days[-1]:
                date_barriers.append(db)
                if has_figure:
                    dt = _matlab_to_datetime(db)
                    ax2.plot([dt, dt], [ax2_ylim[0], ax2_ylim[1]], "g--", linewidth=2)
                    plt.draw()
                    plt.pause(0.05)
        date_barriers = sorted(date_barriers)

        if has_figure:
            import matplotlib.pyplot as plt
            plt.close(fig)

        # ── build date windows and compute PF coefficients ────────────────────
        if not date_barriers:
            date_bins = [(all_days[0], all_days[-1])]
        else:
            edges = [all_days[0]] + date_barriers + [all_days[-1]]
            date_bins = [(edges[k], edges[k + 1]) for k in range(len(edges) - 1)]

        cm_key = f"cm_{round(z * 100)}"
        pf_info[cm_key] = {}

        for d0, d1 in date_bins:
            key  = f"day_{int(d0)}to{int(d1)}"
            mask = (t_dn >= d0) & (t_dn <= d1 + 1)
            coef = pf_coefficients(
                np.column_stack([u_dn[mask], v_dn[mask], w_dn[mask], direction[mask]]),
                bins,
            )
            pf_info[cm_key][key] = coef

    # ── save ──────────────────────────────────────────────────────────────────
    with open(pf_path, "wb") as fh:
        pickle.dump(pf_info, fh)
    print(f"\nPFinfo saved to {pf_path}")

    # ── display summary ───────────────────────────────────────────────────────
    print("\nVerify Global Planar Fit Coefficients")
    for cm_key, date_dict in pf_info.items():
        print(f"\n  Sonic height: {cm_key}")
        for day_key, coef_dict in date_dict.items():
            d0 = int(day_key.split("_")[1].split("to")[0])
            d1 = int(day_key.split("to")[1])
            print(f"    {_matlab_to_datetime(d0).strftime('%Y-%m-%d')} → "
                  f"{_matlab_to_datetime(d1).strftime('%Y-%m-%d')}")
            for bin_key, coef in coef_dict.items():
                b0, b1, b2 = coef[0], coef[1], coef[2]
                pitch = np.degrees(np.arcsin(-b1 / np.sqrt(1 + b1**2)))
                roll  = np.degrees(np.arcsin( b2 / np.sqrt(1 + b2**2)))
                print(f"      {bin_key}: b0={b0:.3g}  b1={b1:.3g}  b2={b2:.3g}"
                      f"  pitch={pitch:.3g}°  roll={roll:.3g}°")

    ok = input("\nIs this correct? (1=yes, 0=no): ").strip()
    if ok != "1":
        raise RuntimeError("PF coefficients rejected by user. Check inputs and rerun.")

    ok2 = input("Okay to begin analysis? (1=yes, 0=no): ").strip()
    if ok2 != "1":
        raise RuntimeError("Program stopped by user.")

    return pf_info
