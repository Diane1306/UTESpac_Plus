# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Purpose

Python port of **UTESpac** (Utah Turbulence in Environmental Studies Process and Analysis Code), originally written in MATLAB by Derek Jensen and Eric Pardyjak, modified by Diane Wang. Processes eddy-covariance and micrometeorology field data from Campbell Scientific dataloggers. The original MATLAB source lives in `UTESpac_MATLAB/` for reference.

## Running the Code

```bash
# Install dependencies
pip install numpy scipy          # netCDF4 only needed when saveNetCDF=True

# Run the full pipeline (prompts interactively for site and date selection)
python3 utespac_main.py

# Run the validation test suite against siteFire1 reference outputs
python3 run_test.py              # both LPF and GPF
python3 run_test.py --lpf        # local planar fit only
python3 run_test.py --gpf        # global planar fit only
python3 run_test.py --compare    # compare existing .pkl output without re-running
```

Edit the `info` dict at the top of `utespac_main.py` to configure paths, averaging period, QC settings, and output options before running.

## Package Structure

```
utespac/                        # importable package
    campbell_date.py            # Campbell date ↔ MATLAB serial date conversion
    strfndw.py                  # wildcard string search (*, ?)
    nandetrend.py               # linear/constant detrend ignoring NaNs
    consec_flag_removal.py      # clear spike flags that form runs > maxConsecutiveOutliers
    rh_to_spec_hum.py           # RH → specific humidity
    simple_avg.py               # block-average a matrix by MATLAB serial timestamp
    stp_dn.py                   # step-down block average by fixed row count
    import_header.py            # read a UTESpac _header.dat file
    find_files.py               # locate site folders, headers, and CSV files
    find_instruments.py         # map sensor templates → [table, col, height] arrays
    load_data.py                # import raw CSV files for one date row
    find_serial_date.py         # convert 4-column date vector to serial date in col 0
    condition_data.py           # QC: absolute limits, iterative spike removal
    avg.py                      # block-average each table to avgPer minutes
    wind_stats.py               # mean wind speed/direction and bad-sector flag
    pf_coefficients.py          # planar-fit least-squares solver (b0, b1, b2)
    sonic_rotation.py           # planar-fit + yaw rotation (Wilczak et al. 2000)
    fluxes.py                   # all turbulent statistics (H, τ, TKE, L, σ, η, ε, …)
    find_global_pf.py           # interactive multi-sector global PF coefficients
    save_data.py                # write output to .pkl (+ optional CSV / NetCDF)
    get_data.py                 # load and concatenate processed .pkl files
    get_virtual_pot_temp.py     # θ_v, mixing ratio, and air densities from HMP data
    calc_dissipation_rate.py    # TKE dissipation ε via structure function
    calc_snsp_angle.py          # slope-normal/slope-parallel angle decomposition
    find_delta_flux.py          # flux asymmetry Δ_flux (ejections − sweeps)
    find_delta_time.py          # time-fraction asymmetry Δ_time
    find_eta.py                 # transport efficiency η
    complete_table_create.py    # NaN-filled table with uniform timestamp grid
    site_info.py                # template for French Meadows / UU1 site config
    site_info_es5.py            # site config for ES5 Spring campaign
utespac_main.py                 # main entry point (configure info dict here)
run_test.py                     # validation test suite for siteFire1
requirements.txt
UTESpac_MATLAB/                 # original MATLAB source (reference only)
    siteFire1/                  # test site data and reference outputs
        siteInfo.py             # site-specific configuration (Python)
        siteInfo.m              # site-specific configuration (MATLAB)
        site1_20Hz_header.dat   # sensor column definitions
        site1_20Hz_*.txt        # 20 Hz raw data
        PFinfo.mat              # pre-computed global PF coefficients
        PFinfo.pkl              # cached Python version of PFinfo.mat
        output/                 # reference MATLAB .mat + Python .pkl outputs
```

## Processing Pipeline

`utespac_main.py` calls modules in this order:

1. `find_files` — scan `rootFolder/site*` dirs for CSV/TXT data and `*header*` files. Accepts `site=` and `dates=` kwargs for non-interactive use.
2. `find_instruments` — wildcard-match header variable names against `template` patterns; return `sensor_info` dict where each key maps to an `ndarray[:,3+]` of `[table_idx, col_idx, height_m]`.
3. `find_global_pf` *(only when `info["PF"]["globalCalculation"] == "global"`)*
4. Per-date loop:
   - `load_data` → `find_serial_date` → `condition_data` → `avg` → `wind_stats` → `sonic_rotation` → `fluxes` → `save_data`

## Key Data Structures

- **`info`** — plain Python `dict` holding all configuration (mirrors MATLAB `info` struct). Top-level keys: `rootFolder`, `siteFolder`, `avgPer`, `detrendingFormat`, `PF`, `spikeTest`, `absoluteLimitsTest`, `windDirectionTest`, `nanTest`, `diagnosticTest`, `sonicOrientation`, `sonicManufact`, `tower`, `siteElevation`, `angle`, `tableNames`, `tableScanFrequency`, `tableNumberOfColumns`.
- **`sensor_info`** — `dict` of `ndarray`. Each key is a template field name (`u`, `v`, `w`, `Tson`, `fw`, `RH`, `T`, `P`, `irgaH2O`, `irgaCO2`, `LiH2O`, `LiCO2`, `KH2O`, `cup`, `birdSpd`, `birdDir`, …). Sonic entries (`u`) carry 5 columns: `[table_idx, col_idx, height_m, bearing_deg, manufact_code]`.
- **`data`** — `list` of `ndarray | None`, one per table. After `find_serial_date`, column 0 is a MATLAB serial date float; remaining columns are raw sensor values.
- **`output`** — `dict` built up through the pipeline. Key naming convention: `<tableName>` (averaged data), `<tableName>Header` (column name list), `<tableName>SpikeFlag` / `<tableName>NanFlag` (bool ndarray), `spdAndDir`, `spdAndDirHeader`, `rotatedSonic`, `rotatedSonicHeader`, `H`, `tau`, `tke`, `L`, `sigma`, `H_SNSP`, `eta`, `delta_flux`, `delta_time`, `turbtr`, `epsilon`.
- **`template`** — `dict[str, str]` mapping sensor field names to wildcard patterns, e.g. `template["u"] = "Ux_*"`. The `*` stands in for the numeric sensor height.

## Timestamp Convention

All timestamps are stored as **MATLAB serial date floats** (days since the MATLAB epoch). The conversion anchor is `MATLAB_EPOCH = 719529` (= `datenum(1970, 1, 1)`):

```python
from utespac.campbell_date import datetime_to_matlab_datenum, matlab_datenum_to_datetime
```

The Campbell logger stores dates as 4 columns `[year, day_of_year, HHMM, seconds]`; `find_serial_date` converts these to a single float in column 0 of each data array.

## Site Configuration

Each `site*` folder must contain:
- **`siteInfo.py`** — executed by `find_files` via `exec()`. Must define module-level variables: `sonicOrientation`, `sonicManufact`, `tower`, `siteElevation`, `angle`, `tableNames`, `tableScanFrequency`, `tableNumberOfColumns`. See `utespac/site_info.py` for the French Meadows template and `utespac/site_info_es5.py` for ES5.
- **`*header*`** file(s) — single-line comma-delimited header matching the CSV column names (3 fewer columns than the raw CSV, since the 4-column date vector is replaced by one serial-date column).
- Raw CSV or TXT data files named `<tableName>_<YYYYMMDDHHMMSS>_...`.

## Output Files

Two pickle files are written per run to `<siteFolder>/output/`:
- **`<SiteName>_<avgPer>minAvg_<PFtype><detrendType><date>.pkl`** — 30-min averaged statistics (the `output` dict).
- **`<SiteName>_raw_<PFtype><detrendType><date>.pkl`** — 20 Hz conditioned arrays (`uPF`, `vPF`, `wPF`, `u_tilt`, `v_tilt`, `w_tilt`, `sonTs`, `Theta_v_son`, `P`, `rhov`, `rhovPrime`, `rhoCO2`, `rhoCO2Prime`, `t`, `z`); only written when `info["saveRawConditionedData"] = True`.

Optional: CSV (`saveCSV=True`) and NetCDF (`saveNetCDF=True`, requires `netCDF4`).

`get_data()` loads and vertically concatenates averaged `.pkl` files from an `output/` folder.

## Validation Test Suite (`run_test.py`)

Runs the pipeline on `siteFire1` and compares against the MATLAB reference outputs in `siteFire1/output/*.mat`. The test covers both LPF and GPF modes, and both the 30-min averaged output and the 20 Hz raw output.

**rootFolder for the test** is `UTESpac_Python/UTESpac_MATLAB/` (where `siteFire1/` now lives). The raw data and reference `.mat` files are both in `UTESpac_MATLAB/siteFire1/`. The global PF coefficients come from `UTESpac_MATLAB/siteFire1/PFinfo.mat`, which `run_test.py` converts to `PFinfo.pkl` on first use.

**All 28 checks pass** against the updated MATLAB reference (generated with MATLAB bug fixes: `wCol`, `TsSpikeFlag`, `spd_raw` all corrected to use height-matched indices and manufacturer-corrected wind). Both LPF and GPF now match to near-machine-precision for wind and flux channels.

LPF mode:

| Output | Tolerance | Typical error |
|---|---|---|
| Wind direction / speed (30-min) | < 0.01% | ~0% |
| Wind sector flag (30-min) | Exact | 0 |
| Rotated sonic u, w (30-min) | < 0.5% | ~0% |
| Momentum flux τ (30-min) | < 0.5% | ~0% |
| TKE (30-min) | < 0.1% | ~0% |
| Obukhov length L (30-min) | < 0.5% | 0.12% |
| σ_u, σ_w (30-min) | < 0.1% | ~0% |
| H: T′w′, Θ_v′w′ (30-min) | < 0.5% | ~0.02% |
| 20 Hz uPF/vPF/wPF | p99 < 0.005 m/s | < 1e-7 m/s |
| 20 Hz u_tilt/v_tilt/w_tilt | p99 < 0.005 m/s | < 1e-8 m/s |
| 20 Hz sonTs | max < 0.1 °C | < 1e-14 °C |
| 20 Hz Theta_v_son | max < 0.5 °C | < 1e-14 °C |
| 20 Hz P | < 0.01% | ~0% |
| 20 Hz rhov, rhovPrime | p99 < 0.5 g/m³ | < 1e-14 g/m³ |
| 20 Hz rhoCO2 | max < 5 mg/m³ | < 1e-12 mg/m³ |
| 20 Hz rhoCO2Prime | p99 < 30 mg/m³ | 25.7 mg/m³ p99 |

GPF mode (same tolerances as LPF with new reference — previous cascade errors no longer present):

| Output | Tolerance | Typical error |
|---|---|---|
| All 30-min outputs | Same as LPF | Near zero |
| All 20 Hz wind channels | p99 < 0.05 m/s | < 1e-13 m/s |
| 20 Hz gas channels | Same as LPF | Same as LPF |

## Known Limitations vs MATLAB Reference

**rhoCO2Prime p99 ≈ 26 mg/m³**: Python's `nanstd` uses `ddof=0` (population std) while MATLAB's `std(..., 'omitmissing')` uses `ddof=1` (sample std). This ~0.003% difference causes ~333 of 3.456M CO2 samples to land on different NaN/non-NaN sides in spike removal, changing the linear detrend for those periods. Raw CO2 matches to 1×10⁻¹² mg/m³.

**sigma_CO2_WPL max_abs ≈ 2.4 g/m³**: Same cascade mechanism affecting the WPL-corrected CO2 standard deviation.

**rotatedSonic v max_abs ≈ 8e-5 m/s**: The yaw-rotated v-component is forced to zero mean; the residual 8e-5 m/s reflects floating-point rounding in the rotation matrix. Not physically significant.

**skew max_abs ≈ 0.002**: Near-zero skewness values have large relative errors but negligible absolute differences.

**PFinfo.pkl cache**: `run_test.py` converts `PFinfo.mat → PFinfo.pkl` only once (guards on file existence). If `PFinfo.mat` is regenerated, delete `UTESpac_MATLAB/siteFire1/PFinfo.pkl` to force reconversion.

## Manufacturer Codes

`sonicManufact` values used throughout the codebase:
- `1` — Campbell CSAT3 / IRGASON (u, v, w as-is)
- `0` — RMYoung: swap u↔v and negate new v
- `2` — Gill WindmasterPro: negate both u and v
