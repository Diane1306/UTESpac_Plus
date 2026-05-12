"""RHtoSpecHum – convert relative humidity to specific humidity."""

import numpy as np


def rh_to_spec_hum(rh: np.ndarray, P: np.ndarray, T: np.ndarray) -> np.ndarray:
    """Convert relative humidity to specific humidity.

    Parameters
    ----------
    rh : array_like
        Relative humidity [%].
    P : array_like
        Pressure [kPa].
    T : array_like
        Temperature [K].

    Returns
    -------
    q : ndarray
        Specific humidity [kg/kg].
    """
    rh = np.asarray(rh, dtype=float)
    P  = np.asarray(P,  dtype=float)
    T  = np.asarray(T,  dtype=float)

    l_v = 2.5e6    # Latent heat of vaporization [J/kg]
    R_v = 461.5    # Gas constant for water vapour [J/(kg·K)]
    T_0 = 273.15   # Reference temperature [K]
    es0 = 6.11     # Saturated vapour pressure at T_0 [hPa]

    es = es0 * np.exp((l_v / R_v) * (1.0 / T_0 - 1.0 / T))  # [hPa]
    e  = rh / 100.0 * es                                      # actual vapour pressure [hPa]
    q  = 0.622 * (e / (P * 10.0))                             # [kg/kg]
    return q
