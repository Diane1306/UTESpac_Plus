"""getData – load and concatenate processed output files from a site folder."""

import os
import glob
import pickle
from typing import Dict, Optional
import numpy as np


def get_data(
    root_folder: str,
    site: str = None,
    avg_per: int = None,
    qualifier: str = None,
    rows=None,
) -> Dict:
    """Load and vertically concatenate processed UTESpac output files.

    Parameters
    ----------
    root_folder : str
        Root directory containing ``site*`` sub-directories.
    site : str or int or None
        Site name (without ``'site'`` prefix) or 1-based integer index.
        If None the user is prompted to select.
    avg_per : int or None
        Averaging period in minutes (used to filter file names).
    qualifier : str or None
        Substring that must appear in output file names (e.g. ``'LPF'``).
    rows : int, list, or 0
        0 → all rows; integer → first N rows; list → specific row indices.

    Returns
    -------
    output_struct : dict
        Concatenated structure with the same keys as individual output files.
    """
    sites = sorted(
        d for d in os.listdir(root_folder)
        if os.path.isdir(os.path.join(root_folder, d)) and d.startswith("site")
    )

    if site is None:
        for i, s in enumerate(sites):
            print(f"  {i + 1}. {s}")
        choice = int(input("Please indicate site of interest: ")) - 1
        site_dir = sites[choice]
    elif isinstance(site, int):
        site_dir = sites[site - 1]
    elif site.startswith("site"):
        site_dir = site
    else:
        site_dir = f"site{site}"

    site_path = os.path.join(root_folder, site_dir, "output")
    if not os.path.isdir(site_path):
        raise FileNotFoundError(f"Output folder not found: {site_path}")

    # Build glob pattern
    pattern_parts = ["*"]
    if avg_per is not None:
        pattern_parts.append(f"_{avg_per}*")
    if qualifier is not None:
        pattern_parts.append(f"*{qualifier}*")
    pattern_parts.append("*.pkl")
    pattern = os.path.join(site_path, "".join(pattern_parts))

    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No output files found matching {pattern}")

    # Select rows
    if rows is None or rows == 0:
        selected_files = files
    elif isinstance(rows, int):
        selected_files = files[:rows]
    else:
        selected_files = [files[r - 1] for r in rows if 0 < r <= len(files)]

    # Load and concatenate
    output_struct: Dict = {}
    for fpath in selected_files:
        with open(fpath, "rb") as fh:
            d = pickle.load(fh)
        for key, val in d.items():
            if key not in output_struct:
                output_struct[key] = val
            elif isinstance(val, np.ndarray) and isinstance(output_struct[key], np.ndarray):
                if val.ndim == output_struct[key].ndim == 2 and val.shape[1] == output_struct[key].shape[1]:
                    output_struct[key] = np.vstack([output_struct[key], val])
                elif val.ndim == output_struct[key].ndim == 1:
                    output_struct[key] = np.concatenate([output_struct[key], val])

    return output_struct
