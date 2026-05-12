"""calc_SNSP_angle – slope-normal / slope-parallel angle decomposition."""

import numpy as np


def calc_snsp_angle(phi: float, alpha: float):
    """Decompose wind direction angles in the slope coordinate system.

    Parameters
    ----------
    phi   : float   Meteorological wind direction from north [degrees].
    alpha : float   Slope angle [degrees].

    Returns
    -------
    a1 : float  Angle between u_slope-normal and true vertical [degrees].
    a2 : float  Angle between v_slope-normal and true vertical [degrees].
    """
    kesi = phi - 30.0   # French Meadows downslope direction is 30° from north
    a1 = np.degrees(np.arcsin(np.cos(np.radians(kesi))      * np.sin(np.radians(alpha))))
    a2 = np.degrees(np.arcsin(np.cos(np.radians(kesi - 90)) * np.sin(np.radians(alpha))))
    return a1, a2
