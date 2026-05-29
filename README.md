# UTESpac Python

Python port of UTESpac (Utah Turbulence in Environmental Studies Process and Analysis Code).  
Original MATLAB code by Derek Jensen & Eric Pardyjak, modified by Diane Wang.

## Installation

```bash
pip3 install numpy scipy matplotlib
pip3 install netCDF4   # optional — only needed when saveNetCDF=True
```

## Running the pipeline

`raw_processing` provides sample code for generating formatted input for UTESpac package. Site folders includes sample `siteinfo` files. After generating formatted input and creating corresponding site folder, edit the `info` dict at the top of `utespac_main.py` to set paths, averaging period, QC settings, and output options, then run:

```bash
python3 utespac_main.py
```
Note that `python utespac_main.py --globalCalculation local` needs to be run first before running `python utespac_main.py` (default GPF).

---

## Utility scripts

### Convert siteInfo.m → siteInfo.py

Converts a MATLAB `siteInfo.m` site configuration file into the Python equivalent that
`find_files()` expects.  Run from the `UTESpac_Python/` directory.

```bash
# Convert an explicit .m file
python3 convert_siteinfo.py UTESpac_MATLAB/siteGill/siteInfo.m

```
The script handles scalar values, numeric arrays (`[1 2 3]` or `[1, 2, 3]`), string cell
arrays (`{'name'}`), and preserves inline `%` comments as `#` comments.

---

### Generate AmeriFlux BASE data

`generate_ameriflux.py` loads GPF `.pkl` output files from `siteIRGA` and `siteGill`,
aggregates slow meteorology and radiation from the 1-min data files, and writes a
half-hourly AmeriFlux BASE CSV.  Edit the path constants at the top of the script before
running.

```bash
python3 generate_ameriflux.py
```

The output CSV is written to `ameriflux_output/` and follows the
[AmeriFlux BASE format](https://ameriflux.lbl.gov/half-hourly-hourly-data-upload-format/):
comma-delimited, `TIMESTAMP_START` / `TIMESTAMP_END` in `YYYYMMDDHHMM` local standard time,
missing values as `-9999`.

Variables included:

| Group | Variables |
|---|---|
| Turbulent fluxes | `H`, `LE`, `FC`, `TAU`, `USTAR` (×5 EC heights) |
| Wind | `WS`, `WD` (×5 heights) |
| Stability | `MO_LENGTH`, `ZL`, `TKE` (×5 heights) |
| Sonic temperature | `T_SONIC`, `T_SONIC_SIGMA` (×5 heights) |
| Gas scalars | `CO2`, `CO2_SIGMA`, `H2O`, `H2O_SIGMA`, `FH2O` (×5 heights) |
| Velocity variances | `U_SIGMA`, `V_SIGMA`, `W_SIGMA` (×5 heights) |
| Wind direction QC | `WD_FILTER` (×5 heights): 0=clean sector, 1=tower-disturbed |
| Quality flags | `TAU_SSITC_TEST`, `H_SSITC_TEST`, `LE_SSITC_TEST`, `FC_SSITC_TEST` (×5 heights) |
| Slow met | `TA`, `RH`, `VPD` (×4 HMP heights), `PA` |
| Radiation | `SW_IN`, `SW_OUT`, `LW_IN`, `LW_OUT`, `NETRAD`, `ALB` (×2 rad heights) |

**Note on SSITC flags:** SSITC flags are interpreted as diagnostic indicators of nonstationarity
and similarity-theory departure, rather than as direct indicators of instrument failure or
unusable observations over forested complex terrain.
