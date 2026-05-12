"""get_virtulPotTemp – virtual potential temperature and air densities from HMP data."""

import numpy as np


def get_virtual_pot_temp(altitude, level, Tair, RH):
    """Compute virtual potential temperature and air densities.

    Matches MATLAB get_virtulPotTemp: pressure is computed from the
    hypsometric formula using site elevation + sensor height rather than
    from a measured barometer, so the result is independent of barometer
    availability.

    Parameters
    ----------
    altitude : float
        Site elevation above sea level [m]  (info["siteElevation"]).
    level : float
        Sensor height above the reference level [m]  (sonHeight - zRef).
    Tair : array_like
        Air temperature [°C or K – auto-detected].
    RH : array_like
        Relative humidity [%].

    Returns
    -------
    virtual_theta  : ndarray  Virtual potential temperature [K]
    r              : ndarray  Mixing ratio [g/kg]
    rho_airmoist   : ndarray  Density of moist air [kg/m³]
    rho_airdry     : ndarray  Density of dry air [kg/m³]
    rho_H2O        : ndarray  Density of water vapour [kg/m³]
    """
    Tair = np.asarray(Tair, dtype=float)
    RH   = np.asarray(RH,   dtype=float)

    if np.nanmedian(Tair) < 200:
        Tair = Tair + 273.15          # → K

    # Pressure from hypsometric formula [Pa] — matches MATLAB
    P_air = 101325.0 * (1.0 - 2.25577e-5 * (altitude + level)) ** 5.25588

    # Saturated vapour pressure [Pa]
    e_sat = 1000.0 * np.exp(52.57633 - 6790.4985 / Tair - 5.02808 * np.log(Tair))
    e_air = e_sat * (RH / 100.0)      # actual vapour pressure [Pa]

    r = 621.97 * e_air / (P_air - e_air)  # mixing ratio [g/kg]

    Gamma = 0.0098   # dry adiabatic lapse rate [K/m]
    theta = Tair + Gamma * level

    virtual_theta = theta * (1.0 + 0.61 * r / 1000.0)  # [K]

    Rd = 287.058   # [J/(kg·K)]
    Rv = 461.495   # [J/(kg·K)]
    rho_airdry   = (P_air - e_air) / (Rd * Tair)
    rho_H2O      = e_air / (Rv * Tair)
    rho_airmoist = rho_airdry + rho_H2O

    return virtual_theta, r, rho_airmoist, rho_airdry, rho_H2O
