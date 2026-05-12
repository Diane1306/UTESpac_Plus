"""find_delta_time – time-fraction asymmetry (ejections minus sweeps)."""

import numpy as np


def find_delta_time(ww: np.ndarray, uu: np.ndarray) -> float:
    """Compute Δ_time = (ejection fraction – sweep fraction) of total time.

    Parameters
    ----------
    ww : 1-D ndarray   Vertical velocity perturbations w'.
    uu : 1-D ndarray   Scalar/velocity perturbations u' (or T', q', etc.).

    Returns
    -------
    delta_time : float
    """
    flux_all  = ww * uu
    flux_total = np.nansum(flux_all)
    if flux_total == 0:
        return np.nan
    downgradient = (flux_all * flux_total) > 0
    ejection  = downgradient & (ww > 0)
    sweep     = downgradient & (ww < 0)
    n = len(ww)
    return float(np.sum(ejection)) / n - float(np.sum(sweep)) / n
