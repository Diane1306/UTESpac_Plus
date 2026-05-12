"""windStats – compute mean wind speed, direction, and direction QC flag."""

from typing import Dict, List
import numpy as np


def wind_stats(output: Dict, sensor_info: Dict, table_names: List[str], info: Dict) -> Dict:
    """Compute wind speed and direction for each sonic and flag bad-sector data.

    Populates ``output['spdAndDir']`` and ``output['spdAndDirHeader']``.
    Also initialises ``output['warnings']``.
    """
    output.setdefault("warnings", [])

    if "u" not in sensor_info:
        return output

    num_sonics   = sensor_info["u"].shape[0]
    tower_bearing = float(info.get("tower", 0))
    wind_env      = float(info["windDirectionTest"]["envelopeSize"])

    for ii in range(num_sonics):
        try:
            tbl_idx  = int(sensor_info["u"][ii, 0])
            bearing  = float(sensor_info["u"][ii, 3])
            height   = float(sensor_info["u"][ii, 2])
            manufact = int(sensor_info["u"][ii, 4]) if sensor_info["u"].shape[1] > 4 else 1

            u_col = int(sensor_info["u"][sensor_info["u"][:, 2] == height, 1][0])
            v_col = int(sensor_info["v"][sensor_info["v"][:, 2] == height, 1][0])
            tname = table_names[tbl_idx]

            t = output[tname][:, 0]
            u_raw = output[tname][:, u_col]
            v_raw = output[tname][:, v_col]

            # Manufacturer-specific component swap
            if manufact == 1:
                u_dir, v_dir = u_raw, v_raw
            elif manufact == 0:   # RMYoung: v = Campbell u
                u_dir, v_dir = v_raw, u_raw * -1.0
            else:                 # Gill WindmasterPro
                u_dir, v_dir = u_raw * -1.0, v_raw * -1.0

            direction = (np.degrees(np.arctan2(-v_dir, u_dir)) + bearing) % 360.0
            speed     = np.hypot(u_dir, v_dir)

            # Wind direction quality flag
            total_bearing = (tower_bearing + bearing) % 360.0
            flag = np.zeros(len(direction))

            if total_bearing + wind_env > 360.0:
                min_a = total_bearing - wind_env
                max_a = (total_bearing + wind_env) - 360.0
                flag[(direction > min_a) | (direction < max_a)] = 1
            elif total_bearing - wind_env < 0.0:
                min_a = 360.0 + total_bearing - wind_env
                max_a = total_bearing + wind_env
                flag[(direction > min_a) | (direction < max_a)] = 1
            else:
                min_a = total_bearing - wind_env
                max_a = total_bearing + wind_env
                flag[(direction > min_a) & (direction < max_a)] = 1

            # Build spdAndDir matrix (timestamps in col 0, then [dir, spd, flag] * nSonics)
            n = len(t)
            if "spdAndDir" not in output:
                output["spdAndDir"]       = np.full((n, 1 + num_sonics * 3), np.nan)
                output["spdAndDirHeader"] = ["timeStamp"] + [""] * (num_sonics * 3)
                output["spdAndDir"][:, 0] = t

            c0 = 1 + ii * 3
            output["spdAndDir"][:n, c0]     = direction
            output["spdAndDir"][:n, c0 + 1] = speed
            output["spdAndDir"][:n, c0 + 2] = flag
            output["spdAndDirHeader"][c0]     = f"{height}m direction"
            output["spdAndDirHeader"][c0 + 1] = f"{height}m speed"
            output["spdAndDirHeader"][c0 + 2] = (
                f"{height}m flag {min_a:.3g}<dir<{max_a:.3g}"
            )

        except Exception as exc:
            msg = f"Unable to find wind stats at {height}m: {exc}"
            import warnings; warnings.warn(msg)
            output["warnings"].append(msg)

    return output
