"""importHeader – read a UTESpac header .dat file."""

import re
from typing import List


def import_header(filename: str) -> List[List]:
    """Read a single-line, comma-delimited header file.

    Returns a list of two sublists:
        result[0] : list of variable name strings (row 0 in MATLAB)
        result[1] : list of sensor heights (float) extracted from each name,
                    or None when no numeric height is found.
    """
    with open(filename, "r") as fh:
        line = fh.readline()

    # Strip surrounding quotes and whitespace from each field
    fields = [f.strip().strip('"').strip("'") for f in line.split(",")]
    names = fields

    heights = []
    for name in names:
        # Find all digit/decimal runs; the last consecutive chunk is the height
        digit_positions = [i for i, c in enumerate(name) if c.isdigit()]
        decimal_positions = [i for i, c in enumerate(name) if c == "."]
        all_pos = sorted(set(digit_positions + decimal_positions))

        if not all_pos:
            heights.append(None)
            continue

        # Find the last chunk of consecutive positions
        chunks = []
        chunk = [all_pos[0]]
        for pos in all_pos[1:]:
            if pos == chunk[-1] + 1:
                chunk.append(pos)
            else:
                chunks.append(chunk)
                chunk = [pos]
        chunks.append(chunk)
        last_chunk = chunks[-1]
        height_str = "".join(name[p] for p in last_chunk)
        try:
            heights.append(float(height_str))
        except ValueError:
            heights.append(None)

    return [names, heights]
