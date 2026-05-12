"""siteInfo – template for site-specific configuration (French Meadows / UU1).

Copy this file to ``<siteFolder>/siteInfo.py`` and adjust the values.
The main script executes it with ``exec()`` and reads the variables listed below.
"""

# Sonic orientations (azimuth from north, degrees).
# Order: tables sorted alphabetically, then columns in ascending height order.
sonicOrientation = [247, 120, 119]  # UU1

# Sonic manufacturer code per sonic:
#   1 = Campbell (CSAT3 / IRGASON)
#   0 = RMYoung  (v_RMY = u_Campbell !)
#   2 = Gill WindmasterPro
sonicManufact = [1, 1, 1]

# Orientation of the tower body relative to true north [degrees].
tower = 210

# Site elevation above sea level [m].
siteElevation = 1980

# Slope angle [degrees] – used for SNSP heat flux decomposition.
angle = 0.0

# Expected table names (must match the CSV filename table identifier).
tableNames = ["FM_20Hz_FluxTower"]

# Sampling frequency for each table [Hz].
tableScanFrequency = [20]

# Number of columns expected in each raw CSV table (before date compression).
tableNumberOfColumns = [19]
