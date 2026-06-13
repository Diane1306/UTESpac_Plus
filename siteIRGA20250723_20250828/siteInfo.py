"""Site configuration for siteIRGA (converted from siteInfo.m)."""

sonicOrientation     = [243, 128, 134, 139]  # degrees from true north
sonicManufact        = [1, 1, 1, 1]  		 # 4 IRGASONs: 0 as RMYoung, 1 as IRGASON, 2 as Gill WindmasterPro
tower                = 210  			     # Tower orientation relative to Sonic
siteElevation        = 1980  				 # (m)
angle                = 8.2         			 # mean of 20 m radius buffer
tableNames           = ["FMDOL_20Hz", "FMDOL_1min"]  # fast and slow table names
tableScanFrequency   = [20, 1/60]  # [Hz]
tableNumberOfColumns = [48, 10]    # Note that the number of columns in the output structure will be 3 less than the number in the .csv file.  This is because the 4 column date vector is replaced with a single-column serial time.

# NOTE — T/RH height labelling in FMDOL_1min_header:
#   The LICOR daqm HMP sensors are physically at 15 m and 30 m, but the
#   column names in the header use the co-located sonic heights (13.94 m and
#   32.18 m) so that findInstruments can pair T/RH with the correct sonic.
#   shiftsSonHeight/shiftsHMPHeight tells fluxes.py to use the true HMP
#   heights for the altitude correction inside get_virtual_pot_temp.

useTrefHMP      = True 			   # Use slow-response sensor (e.g. HMP, EE181) if exits
avgSlowFreq     = 1                # min; used to compute EC pressure at HMP freq
shiftsSonHeight = [13.94, 32.18]   # sonic heights that have paired HMP sensors; Can skip this and the following line if no need to shift
shiftsHMPHeight = [15,    30   ]   # corresponding physical HMP heights [m]
shiftzRef     = False              # True if needs to prelocate reference height
zRefLowestSon = 4.42               # Reference height
ascending     = False  			   # heights stored high → low in input tables

# SSITC quality flag settings
SSITC_subAvgMin   = 5     # sub-period length [min] for steady-state test
displacementHeight = 0    # zero-plane displacement height [m]
canopyHeight      = 19.3  # canopy height [m] for in-canopy ITC model
useCanopyITC      = True  # use Rannik et al. canopy σ_w/u* for z ≤ canopyHeight
