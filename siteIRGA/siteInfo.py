"""Site configuration for siteIRGA (converted from siteInfo.m)."""

sonicOrientation     = [243, 128, 134, 139]  # UU1
sonicManufact        = [1, 1, 1, 1]  # 3 IRGASONs
tower                = 210  # FM
siteElevation        = 1980  # (m)
angle                = 8.2  # mean of 20 m radius buffer
tableNames           = ["FMDOL_20Hz"]  # modified by Diane
tableScanFrequency   = [20]  # [Hz]
tableNumberOfColumns = [54]  # modified by Diane

useTrefHMP    = True
avgSlowFreq   = 1
shiftzRef     = False
zRefLowestSon = 4.42
ascending     = False  # heights stored high → low in input tables

# SSITC quality flag settings
SSITC_subAvgMin   = 5     # sub-period length [min] for steady-state test
displacementHeight = 0    # zero-plane displacement height [m]
canopyHeight      = 19.3  # canopy height [m] for in-canopy ITC model
useCanopyITC      = True  # use Rannik et al. canopy σ_w/u* for z ≤ canopyHeight
