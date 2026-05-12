"""PF_coefficients – compute planar-fit (b0, b1, b2) for one or more direction bins."""

from typing import Dict, List
import numpy as np


def pf_coefficients(input_mat: np.ndarray, bins: List[float]) -> Dict[str, np.ndarray]:
    """Solve for planar-fit coefficients that force w̄ → 0.

    Parameters
    ----------
    input_mat : ndarray, shape (N, 4)
        Columns: [u, v, w, direction].
    bins : list of float
        Direction bin boundaries in degrees.  Empty list = single sector.

    Returns
    -------
    coef : dict
        Keys like ``'degrees_0_to_90'``; values are 3-element arrays [b0, b1, b2].
    """
    coef: Dict[str, np.ndarray] = {}
    num_bins = len(bins)
    if num_bins == 0:
        num_bins = 1

    for i in range(num_bins):
        if len(bins) == 0:
            matrix = input_mat
            key = "degrees_0_to_0"
        elif i < len(bins) - 1:
            mask = (input_mat[:, 3] > bins[i]) & (input_mat[:, 3] <= bins[i + 1])
            matrix = input_mat[mask, :]
            key = f"degrees_{bins[i]:g}_to_{bins[i+1]:g}"
        else:
            mask = (input_mat[:, 3] > bins[i]) | (input_mat[:, 3] <= bins[0])
            matrix = input_mat[mask, :]
            key = f"degrees_{bins[i]:g}_to_{bins[0]:g}"

        ubar = matrix[:, 0]
        vbar = matrix[:, 1]
        wbar = matrix[:, 2]

        # Remove NaNs
        valid = ~(np.isnan(ubar) | np.isnan(vbar) | np.isnan(wbar))
        ubar, vbar, wbar = ubar[valid], vbar[valid], wbar[valid]

        n = len(ubar)
        if n < 4:
            coef[key] = np.array([np.nan, np.nan, np.nan])
            continue

        su  = np.sum(ubar);  sv  = np.sum(vbar);  sw  = np.sum(wbar)
        suv = np.sum(ubar * vbar)
        suw = np.sum(ubar * wbar)
        svw = np.sum(vbar * wbar)
        su2 = np.sum(ubar ** 2)
        sv2 = np.sum(vbar ** 2)

        H = np.array([[n,   su,  sv ],
                      [su,  su2, suv],
                      [sv,  suv, sv2]], dtype=float)
        g = np.array([sw, suw, svw], dtype=float)

        try:
            coef[key] = np.linalg.solve(H, g)
        except np.linalg.LinAlgError:
            coef[key] = np.array([np.nan, np.nan, np.nan])

    return coef
