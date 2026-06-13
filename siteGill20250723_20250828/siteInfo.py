"""Site configuration for siteGill (converted from siteInfo.m)."""

sonicOrientation     = [36]  # degrees from true north
sonicManufact        = [2]   # 1 Gill WindmasterPro: 0 as RMYoung, 1 as IRGASON, 2 as Gill WindmasterPro
tower                = 210   # Tower orientation relative to Sonic
siteElevation        = 1980  # (m)
angle                = 8.2  # mean of 20 m radius buffer
tableNames           = ["FMDOL_10Hz", "FMDOL_1min"]  # fast and slow table names
tableScanFrequency   = [10, 1/60]  # [Hz]
tableNumberOfColumns = [13, 12]  # Note that the number of columns in the output structure will be 3 less than the number in the .csv file.  This is because the 4 column date vector is replaced with a single-column serial time.

# Single sonic at 51.5 m; reference T/RH/P to lowest sonic on full tower (4.42 m).
useTrefHMP    = True   # Use slow-response sensor (e.g. HMP, EE181) if exits
avgSlowFreq   = 1      # min; used to compute EC pressure at HMP freq
shiftzRef     = True   # True if needs to prelocate reference height
zRefLowestSon = 4.42   # Reference height
ascending     = False  # heights stored high → low in input tables

# SSITC quality flag settings
SSITC_subAvgMin   = 5     # sub-period length [min] for steady-state test
displacementHeight = 0    # zero-plane displacement height [m]
canopyHeight      = 19.3  # canopy height [m] for in-canopy ITC model
useCanopyITC      = True  # use Rannik et al. canopy σ_w/u* for z ≤ canopyHeight
