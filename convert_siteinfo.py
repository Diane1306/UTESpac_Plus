"""convert_siteinfo.py — Convert a MATLAB siteInfo.m into a Python siteInfo.py.

Reads every  info.fieldName = value; % comment  assignment from the .m file
and writes module-level variable assignments compatible with find_files()'s exec().

Usage
-----
    python3 convert_siteinfo.py siteFire1/              # folder containing siteInfo.m
    python3 convert_siteinfo.py siteFire1/siteInfo.m    # explicit file path
    python3 convert_siteinfo.py                         # converts all site*/siteInfo.m
                                                        #   found under UTESpac_MATLAB/
    python3 convert_siteinfo.py --dry-run siteFire1/    # print output without writing

The output file is written next to the input as siteInfo.py.  If siteInfo.py already
exists the script asks for confirmation before overwriting (skipped with --force).
"""

import os
import re
import sys
import argparse

# ── MATLAB value parsers ──────────────────────────────────────────────────────

def _to_number(s):
    """Return int or float for a scalar numeric string."""
    s = s.strip()
    try:
        v = float(s)
        return int(v) if v == int(v) else v
    except ValueError:
        return s  # leave as-is if unparseable


def _parse_numeric_array(inner):
    """Parse the inside of [ ... ] — comma- or space-separated numbers."""
    tokens = re.split(r"[,\s]+", inner.strip())
    tokens = [t for t in tokens if t]
    values = [_to_number(t) for t in tokens]
    return values


def _parse_cell_array(inner):
    """Parse the inside of { ... } — quoted strings."""
    return re.findall(r"'([^']*)'", inner)


def _parse_matlab_value(raw):
    """Convert a MATLAB RHS expression to a Python object.

    Handles:
      scalar           210  /  8.2
      numeric array    [35]  /  [1, 1, 1]  /  [247 120 119]
      cell string arr  {'site1_20Hz'}  /  {'a', 'b'}
    """
    s = raw.strip()

    # cell array of strings: {'a', 'b'}
    m = re.fullmatch(r"\{(.+)\}", s, re.DOTALL)
    if m:
        strings = _parse_cell_array(m.group(1))
        return strings

    # numeric array: [1 2 3] or [1, 2, 3]
    m = re.fullmatch(r"\[(.+)\]", s, re.DOTALL)
    if m:
        return _parse_numeric_array(m.group(1))

    # bare scalar
    return _to_number(s)


# ── Python value formatter ────────────────────────────────────────────────────

def _format_python_value(val):
    """Format a Python object as source code."""
    if isinstance(val, list):
        if all(isinstance(v, str) for v in val):
            items = ", ".join(f'"{v}"' for v in val)
        else:
            items = ", ".join(
                str(v) if isinstance(v, int) else repr(v)
                for v in val
            )
        return f"[{items}]"
    if isinstance(val, float):
        # keep trailing .0 so it's obviously a float
        return repr(val) if val != int(val) else f"{int(val)}.0"
    return repr(val) if isinstance(val, str) else str(val)


# ── parser ────────────────────────────────────────────────────────────────────

# Matches:  info.fieldName = <value> ;   optionally followed by  % comment
_ASSIGN_RE = re.compile(
    r"^\s*info\.(\w+)\s*=\s*(.+?)\s*;?\s*(?:%(.*))?$"
)

# Field display order (fields not in this list are appended in parse order)
_FIELD_ORDER = [
    "sonicOrientation",
    "sonicManufact",
    "tower",
    "siteElevation",
    "angle",
    "tableNames",
    "tableScanFrequency",
    "tableNumberOfColumns",
]


def parse_m_file(path):
    """Parse a siteInfo.m file.

    Returns a list of (field_name, python_value, inline_comment) tuples,
    in the order they appear in the file.
    """
    entries = {}  # field → (value, comment) — last assignment wins
    order   = []  # preserve first-seen order

    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            m = _ASSIGN_RE.match(line)
            if not m:
                continue
            field   = m.group(1)
            raw_val = m.group(2).strip()
            comment = (m.group(3) or "").strip()

            try:
                value = _parse_matlab_value(raw_val)
            except Exception as exc:
                print(f"  WARNING: could not parse '{field}' ({raw_val!r}): {exc}")
                continue

            if field not in entries:
                order.append(field)
            entries[field] = (value, comment)

    result = []
    # emit in preferred order, then anything extra
    seen = set()
    for f in _FIELD_ORDER + order:
        if f in entries and f not in seen:
            val, cmt = entries[f]
            result.append((f, val, cmt))
            seen.add(f)
    return result


# ── formatter ─────────────────────────────────────────────────────────────────

def format_py_file(entries, m_path):
    """Render the Python source string for the converted siteInfo."""
    site_name = os.path.basename(os.path.dirname(os.path.abspath(m_path)))
    lines = [f'"""Site configuration for {site_name} (converted from siteInfo.m)."""', ""]

    # measure alignment: longest field name
    max_len = max((len(f) for f, _, _ in entries), default=10)
    pad = max_len + 1  # +1 for at least one space

    for field, value, comment in entries:
        py_val = _format_python_value(value)
        assignment = f"{field:<{pad}}= {py_val}"
        if comment:
            # clean up trailing whitespace from MATLAB comment
            assignment += f"  # {comment.strip()}"
        lines.append(assignment)

    lines.append("")  # trailing newline
    return "\n".join(lines)


# ── file I/O ──────────────────────────────────────────────────────────────────

def convert_file(m_path, dry_run=False, force=False):
    """Convert one siteInfo.m → siteInfo.py."""
    m_path  = os.path.abspath(m_path)
    py_path = os.path.join(os.path.dirname(m_path), "siteInfo.py")

    print(f"\nConverting: {m_path}")
    entries = parse_m_file(m_path)
    if not entries:
        print("  No info.* assignments found — skipping.")
        return False

    source = format_py_file(entries, m_path)

    if dry_run:
        print(f"  [dry-run] would write {py_path}")
        print("  " + "\n  ".join(source.splitlines()))
        return True

    if os.path.isfile(py_path) and not force:
        ans = input(f"  {py_path} already exists. Overwrite? [y/N] ").strip().lower()
        if ans != "y":
            print("  Skipped.")
            return False

    with open(py_path, "w", encoding="utf-8") as fh:
        fh.write(source)
    print(f"  Written: {py_path}")
    return True


# ── discovery ─────────────────────────────────────────────────────────────────

def find_m_files(root):
    """Return all siteInfo.m paths under site* subdirectories of root."""
    found = []
    try:
        entries = sorted(os.listdir(root))
    except FileNotFoundError:
        return found
    for name in entries:
        if not name.startswith("site"):
            continue
        candidate = os.path.join(root, name, "siteInfo.m")
        if os.path.isfile(candidate):
            found.append(candidate)
    return found


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Convert siteInfo.m to siteInfo.py.")
    parser.add_argument(
        "paths", nargs="*",
        help="Site folder(s) or siteInfo.m file(s). "
             "Omit to convert all site*/siteInfo.m under UTESpac_MATLAB/.",
    )
    parser.add_argument("--dry-run", action="store_true",
                        help="Print output without writing files.")
    parser.add_argument("--force", action="store_true",
                        help="Overwrite existing siteInfo.py without prompting.")
    args = parser.parse_args()

    script_dir  = os.path.dirname(os.path.abspath(__file__))
    matlab_dir  = os.path.join(script_dir, "UTESpac_MATLAB")

    if args.paths:
        m_files = []
        for p in args.paths:
            p = os.path.abspath(p)
            if os.path.isdir(p):
                candidate = os.path.join(p, "siteInfo.m")
                if os.path.isfile(candidate):
                    m_files.append(candidate)
                else:
                    print(f"WARNING: no siteInfo.m found in {p}")
            elif os.path.isfile(p):
                m_files.append(p)
            else:
                print(f"WARNING: path not found: {p}")
    else:
        m_files = find_m_files(matlab_dir)
        if not m_files:
            print(f"No site*/siteInfo.m files found under {matlab_dir}")
            sys.exit(1)
        print(f"Found {len(m_files)} siteInfo.m file(s) under {matlab_dir}")

    converted = sum(
        convert_file(p, dry_run=args.dry_run, force=args.force)
        for p in m_files
    )
    print(f"\nDone. {converted}/{len(m_files)} file(s) converted.")


if __name__ == "__main__":
    main()
