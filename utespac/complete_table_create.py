"""completeTableCreate – create a NaN-filled table with a complete timestamp grid."""

import numpy as np
from .campbell_date import serial_date_to_campbell_date, datetime_to_matlab_datenum
from datetime import datetime


def complete_table_create(
    begin_serial_day: float,
    end_serial_day: float,
    sample_frequency: float,
    expected_columns: int,
) -> np.ndarray:
    """Create a NaN-padded table with a uniform timestamp grid.

    The timestamp grid runs from ``begin_serial_day + dt`` to
    ``end_serial_day + 1`` (inclusive), matching MATLAB behaviour.

    Parameters
    ----------
    begin_serial_day : float
        MATLAB serial date of the first day (integer part).
    end_serial_day : float
        MATLAB serial date of the last day (integer part).
    sample_frequency : float
        Sampling frequency in Hz.
    expected_columns : int
        Number of columns in the output table.  The first 4 columns hold
        the Campbell date vector [year, DOY, HHMM, seconds].

    Returns
    -------
    filled_table : ndarray, shape (N, expected_columns)
    """
    dt_days = 1.0 / sample_frequency / 86400.0
    t = np.arange(begin_serial_day + dt_days, end_serial_day + 1.0 + dt_days / 2.0, dt_days)

    filled = np.full((len(t), expected_columns), np.nan)
    camp_dates = serial_date_to_campbell_date(t)
    filled[:, :4] = camp_dates
    return filled
