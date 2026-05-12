"""find_eta – transport efficiency η (total / downgradient flux ratio)."""

import numpy as np


def find_eta(ww: np.ndarray, uu: np.ndarray) -> float:
    """Compute η = total_flux / downgradient_flux.

    Parameters
    ----------
    ww : 1-D ndarray   w' perturbations.
    uu : 1-D ndarray   Scalar/velocity perturbations.

    Returns
    -------
    eta : float
    """
    flux_all   = ww * uu
    flux_total = np.nansum(flux_all)
    downgradient = (flux_all * flux_total) > 0
    dg_total = np.nansum(flux_all[downgradient])
    if dg_total == 0:
        return np.nan
    return flux_total / dg_total
