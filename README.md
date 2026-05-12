# UTESpac Python

Python port of UTESpac (Utah Turbulence in Environmental Studies Process and Analysis Code).  
Original MATLAB code by Derek Jensen & Eric Pardyjak, modified by Diane Wang.

## Installation

```bash
pip3 install numpy scipy matplotlib
pip3 install netCDF4   # optional — only needed when saveNetCDF=True
```

## Running the pipeline

Edit the `info` dict at the top of `utespac_main.py` to set paths, averaging period, QC
settings, and output options, then run:

```bash
python3 utespac_main.py
```

## Validation test suite

Runs the pipeline on `siteFire1` and compares output against the MATLAB reference files in
`UTESpac_MATLAB/siteFire1/output/`.

```bash
python3 run_test.py              # LPF + GPF
python3 run_test.py --lpf        # local planar fit only
python3 run_test.py --gpf        # global planar fit only
python3 run_test.py --compare    # compare existing .pkl output without re-running
```

---

## Utility scripts

### Convert siteInfo.m → siteInfo.py

Converts a MATLAB `siteInfo.m` site configuration file into the Python equivalent that
`find_files()` expects.  Run from the `UTESpac_Python/` directory.

```bash
# Convert all site*/siteInfo.m files found under UTESpac_MATLAB/
python3 convert_siteinfo.py

# Convert a single site folder
python3 convert_siteinfo.py UTESpac_MATLAB/siteFire1/

# Convert an explicit .m file
python3 convert_siteinfo.py UTESpac_MATLAB/siteGill/siteInfo.m

# Preview the output without writing any files
python3 convert_siteinfo.py --dry-run UTESpac_MATLAB/siteFire1/

# Overwrite existing siteInfo.py without prompting
python3 convert_siteinfo.py --force
```

The script handles scalar values, numeric arrays (`[1 2 3]` or `[1, 2, 3]`), string cell
arrays (`{'name'}`), and preserves inline `%` comments as `#` comments.

---

### Compare Python vs MATLAB outputs

After running both the MATLAB pipeline (producing `.mat` files) and the Python pipeline
(producing `.pkl` files) for any site under `UTESpac_MATLAB/`, run this script to get a
field-by-field, column-by-column statistical comparison.  Run from the `UTESpac_Python/`
directory.

```bash
# Compare all sites, all modes (LPF + GPF, averaged + raw)
python3 compare_outputs.py

# Limit to one site
python3 compare_outputs.py --site siteFire1

# LPF or GPF only
python3 compare_outputs.py --site siteFire1 --lpf
python3 compare_outputs.py --site siteFire1 --gpf

# 30-min averaged output only
python3 compare_outputs.py --avg

# 20 Hz raw output only
python3 compare_outputs.py --raw

# Save the full stats table to a CSV file
python3 compare_outputs.py --csv results.csv

# Change the OK/WARN threshold (default is 1 %)
python3 compare_outputs.py --tol-rel 0.05
```

Each row in the output table covers one numeric column and reports:

| Statistic | Meaning |
|---|---|
| `N_valid` | Non-NaN pairs compared |
| `NaN_mm` | Positions where one side is NaN and the other is not |
| `bias` | mean(Python − MATLAB) |
| `RMSE` | Root-mean-square error |
| `max_abs` / `mean_abs` | Absolute error extremes |
| `p95_abs` / `p99_abs` | Percentile absolute errors |
| `max_rel` / `mean_rel` | Relative errors |
| `R²` | Squared correlation |

Status flags: `OK` (max\_rel < tol\_rel), `WARN` (1–10 %), `BIG` (> 10 %),
`MISS_PY` / `MISS_REF` (field absent from one side), `SKIP` (all NaN).
