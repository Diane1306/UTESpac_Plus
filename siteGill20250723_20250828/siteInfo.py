"""Site configuration for siteGill (converted from siteInfo.m)."""

sonicOrientation     = [36]  # UU1
sonicManufact        = [2]
tower                = 210  # FM
siteElevation        = 1980  # (m)
angle                = 8.2  # mean of 20 m radius buffer
tableNames           = ["FMDOL_10Hz", "FMDOL_1min"]  # modified by Diane
tableScanFrequency   = [10, 1/60]  # [Hz]
tableNumberOfColumns = [13, 12]  # modified by Diane (21 = +Temp_6.35, RH_6.35)

# Single sonic at 51.5 m; reference T/RH/P to lowest sonic on full tower (4.42 m).
useTrefHMP    = True
avgSlowFreq   = 1
shiftzRef     = True
zRefLowestSon = 4.42
ascending     = False  # heights stored high → low in input tables

# SSITC quality flag settings
SSITC_subAvgMin   = 5     # sub-period length [min] for steady-state test
displacementHeight = 0    # zero-plane displacement height [m]
canopyHeight      = 19.3  # canopy height [m] for in-canopy ITC model
useCanopyITC      = True  # use Rannik et al. canopy σ_w/u* for z ≤ canopyHeight
