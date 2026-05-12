"""find_delta_flux – flux contribution asymmetry (ejections minus sweeps)."""

import numpy as np


def find_delta_flux(ww: np.ndarray, uu: np.ndarray) -> float:
    """Compute Δ_flux = (ejection contribution – sweep contribution) / total flux.

    Parameters
    ----------
    ww : 1-D ndarray   Vertical velocity perturbations w'.
    uu : 1-D ndarray   Scalar/velocity perturbations u' (or T', q', etc.).

    Returns
    -------
    delta_flux : float
    """
    flux_all  = ww * uu
    flux_total = np.nansum(flux_all)
    if flux_total == 0:
        return np.nan
    downgradient = (flux_all * flux_total) > 0
    ejection  = downgradient & (ww > 0)
    sweep     = downgradient & (ww < 0)
    return (np.nansum(flux_all[ejection]) / flux_total
            - np.nansum(flux_all[sweep])  / flux_total)
