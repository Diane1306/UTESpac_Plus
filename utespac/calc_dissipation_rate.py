"""calc_dissipation_rate – TKE dissipation rate via structure function."""

import numpy as np


def calc_structure_function(u: np.ndarray) -> np.ndarray:
    """Second-order longitudinal structure function D_LL(r).

    Parameters
    ----------
    u : 1-D ndarray  velocity time series

    Returns
    -------
    sf : ndarray, length floor(N/2)
    """
    n = len(u)
    r_max = n // 2
    sf = np.zeros(r_max)
    for ri in range(1, r_max + 1):
        diff = (u[ri:] - u[:-ri]) ** 2
        sf[ri - 1] = np.nanmean(diff)
    return sf


def calc_dissipation_rate(u_prime: np.ndarray, u_mean: float, dt: float) -> float:
    """Estimate TKE dissipation rate ε from the structure function.

    Parameters
    ----------
    u_prime : 1-D ndarray   Velocity perturbation time series [m/s].
    u_mean  : float          Mean velocity to convert time lags to spatial lags [m/s].
    dt      : float          Time step [s].

    Returns
    -------
    epsilon : float   Dissipation rate [m²/s³].
    """
    u_prime = np.asarray(u_prime, dtype=float)
    n = len(u_prime)
    sf = calc_structure_function(u_prime)

    r_values = np.linspace(dt * u_mean, (n // 2) * dt * u_mean, n // 2)
    y2 = r_values ** (2.0 / 3.0)

    # Find the index where the 2/3 power law best matches the structure function
    n_fit = min(20, len(y2))
    with np.errstate(divide="ignore", invalid="ignore"):
        diff = np.abs(np.log(y2[:n_fit]) - np.log(np.where(sf[:n_fit] > 0, sf[:n_fit], np.nan)))
    I = int(np.nanargmin(diff))

    ratio = y2[I] / sf[I] if sf[I] > 0 else np.nan
    c_A2  = (y2 / ratio) / y2
    epsilon_arr = (c_A2 ** (3.0 / 2.0)) * 0.35
    return float(epsilon_arr[0]) if not np.isnan(epsilon_arr[0]) else np.nan
