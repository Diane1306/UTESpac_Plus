"""calc_ssitc_flags – EddyPro-like Mauder & Foken (2004) SSITC quality flags.

Flags:  0 = high quality, 1 = moderate (usable), 2 = poor (discard / gap-fill)

Steady-state test: compares full-period covariance to mean of sub-period
covariances (default 5-min sub-periods within the 30-min averaging window).

ITC test: uses sigma_w / u* relative to similarity theory expectation.
"""

import numpy as np


def calc_ssitc_flags(
    wPF_P, uPF_P, vPF_P, ThvP,
    H2Op, rhov_ext,
    rho_CO2p, rhoc_ext,
    ustar, L_val,
    rot_flag, Ts_flag, h2o_flag, co2_flag,
    n_sub, height, d=0.0,
):
    """Return (tau_flag, H_flag, LE_flag, FC_flag) each 0/1/2 or NaN.

    Parameters
    ----------
    wPF_P, uPF_P, vPF_P : ndarray
        Detrended PF-rotated wind fluctuations for the averaging period.
    ThvP : ndarray
        Detrended virtual potential temperature fluctuations.
    H2Op : ndarray or None
        Detrended H2O fluctuations [g/m³]; None if no H2O sensor.
    rhov_ext : ndarray or None
        WPL H2O external fluctuation [kg/m³]; None if no H2O sensor.
    rho_CO2p : ndarray or None
        Detrended CO2 fluctuations [kg/m³]; None if no CO2 sensor.
    rhoc_ext : ndarray or None
        WPL CO2 external fluctuation [kg/m³]; None if no CO2 sensor.
    ustar : float
        Friction velocity [m/s] for this period.
    L_val : float
        Obukhov length [m] for this period.
    rot_flag, Ts_flag, h2o_flag, co2_flag : bool
        Period-level QC flags.
    n_sub : int
        Number of sub-periods (avgPer / SSITC_subAvgMin).
    height : float
        Sonic measurement height [m].
    d : float
        Displacement height [m] (default 0).
    """
    z_eff = max(height - d, 0.1)

    # ── ITC_w ─────────────────────────────────────────────────────────────────
    n_w = np.sum(~np.isnan(wPF_P))
    sigma_w = np.nanstd(wPF_P, ddof=1) if n_w > 1 else np.nan

    if (np.isnan(sigma_w) or np.isnan(ustar) or np.isnan(L_val)
            or not np.isfinite(L_val) or ustar <= 0):
        itc_dev = np.nan
    else:
        zeta = z_eff / L_val
        itc_measured = sigma_w / ustar
        itc_model = 1.25 * (1.0 - 3.0 * zeta) ** (1.0 / 3.0) if zeta < 0 else 1.25
        if np.isnan(itc_model) or itc_model <= 0:
            itc_dev = np.nan
        else:
            itc_dev = abs((itc_measured - itc_model) / itc_model) * 100.0

    # ── TAU ──────────────────────────────────────────────────────────────────
    ss_uw = _ss_dev(uPF_P, wPF_P, n_sub)
    ss_vw = _ss_dev(vPF_P, wPF_P, n_sub)
    devs  = [v for v in [ss_uw, ss_vw] if not np.isnan(v)]
    tau_ss_dev = max(devs) if devs else np.nan

    tau_flag = _combine(tau_ss_dev, itc_dev)
    if rot_flag:
        tau_flag = 2

    # ── H ────────────────────────────────────────────────────────────────────
    H_flag = _combine(_ss_dev(wPF_P, ThvP, n_sub), itc_dev)
    if rot_flag or Ts_flag:
        H_flag = 2

    # ── LE ───────────────────────────────────────────────────────────────────
    if H2Op is not None and rhov_ext is not None:
        H2O_WPL = H2Op + rhov_ext * 1e3   # g/m³
        LE_flag = _combine(_ss_dev(wPF_P, H2O_WPL, n_sub), itc_dev)
        if rot_flag or h2o_flag:
            LE_flag = 2
    else:
        LE_flag = np.nan

    # ── FC ───────────────────────────────────────────────────────────────────
    if rho_CO2p is not None and rhoc_ext is not None:
        CO2_WPL = rho_CO2p + rhoc_ext     # kg/m³
        FC_flag = _combine(_ss_dev(wPF_P, CO2_WPL, n_sub), itc_dev)
        if rot_flag or co2_flag:
            FC_flag = 2
    else:
        FC_flag = np.nan

    return tau_flag, H_flag, LE_flag, FC_flag


def _ss_dev(xP, yP, n_sub):
    """Steady-state deviation [%] between full-period and mean sub-period covariance."""
    if xP is None or yP is None:
        return np.nan
    valid = ~(np.isnan(xP) | np.isnan(yP))
    if valid.sum() < 10:
        return np.nan
    cov_full = np.nanmean(xP * yP)
    if np.isnan(cov_full) or abs(cov_full) < np.finfo(float).eps:
        return np.nan
    n   = len(xP)
    bps = np.round(np.linspace(0, n, n_sub + 1)).astype(int)
    cov_sub = np.array([
        np.nanmean(xP[bps[k]:bps[k + 1]] * yP[bps[k]:bps[k + 1]])
        for k in range(n_sub)
    ])
    return abs((np.nanmean(cov_sub) - cov_full) / cov_full) * 100.0


def _combine(ss_dev, itc_dev):
    """Combine SS and ITC deviations into 0/1/2 flag."""
    if np.isnan(ss_dev) or np.isnan(itc_dev):
        return 2
    if ss_dev < 30 and itc_dev < 30:
        return 0
    if ss_dev < 100 and itc_dev < 100:
        return 1
    return 2
