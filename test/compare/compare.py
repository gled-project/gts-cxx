#!/usr/bin/env python3
"""Compares the summaries of gts-ref (C GTS) and gts-cxx (C++) case by case.

usage: compare.py REFDIR CXXDIR [RTOL] [--classify DIR]

For every case run by both programs:
- unsupported in gts-cxx: skipped and counted;
- the C library crashed: reported, not compared;
- otherwise the statuses must agree (ok or error), and every summary line
  present on both sides must match: integers exactly, floating-point values
  to the relative tolerance RTOL (default 1e-9) plus an absolute 1e-12.
Lines present on one side only are listed once per key, not as failures.
The "command", "status", "exit" and "glib_*" lines are not compared.

Cases listed in ties.txt, next to this script, are reported but not counted
as differences: there GTS breaks ties by memory address, and compare.sh
checks their surfaces by distance instead.

A boolean case is degenerate when gts-cxx reports orientation tests that
were exactly zero ("sos_ties"), so that Simulation of Simplicity decided
them. GTS broke those ties by memory address, which depends on everything
the program allocated before, and its results change with it; gts-cxx
breaks them by vertex index. A degenerate case that differs is not counted
as a difference. Instead, when gts-cxx succeeds, its result must be closed,
manifold and orientable. The same check applies where the C library
crashed. --classify DIR takes the sos_ties lines from the summaries in DIR,
for comparisons where CXXDIR holds summaries the C library made of the C++
surfaces.

Exits with 1 if any other case differs, or a degenerate result is not valid.
"""

import math
import pathlib
import sys

SKIP = {"command", "status", "exit", "glib_critical", "glib_warning"}


def parse(path):
    lines = {}
    for line in path.read_text().splitlines():
        key, _, rest = line.partition(" ")
        lines[key] = rest
    return lines


def values(text):
    out = []
    for tok in text.split():
        name, eq, v = tok.partition("=")
        out.append((name if eq else None, v if eq else tok))
    return out


def close(a, b, rtol):
    try:
        ia, ib = int(a), int(b)
        return ia == ib
    except ValueError:
        pass
    try:
        fa, fb = float(a), float(b)
    except ValueError:
        return a == b
    if math.isnan(fa) or math.isnan(fb):
        return math.isnan(fa) and math.isnan(fb)
    return abs(fa - fb) <= rtol * max(abs(fa), abs(fb)) + 1e-12


def ties():
    path = pathlib.Path(__file__).with_name("ties.txt")
    names = set()
    for line in path.read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            names.add(line.split()[0])
    return names


def valid(c):
    return all(c.get(k) == "1" for k in ("closed", "manifold", "orientable"))


def main():
    args = sys.argv[1:]
    classify_dir = None
    if "--classify" in args:
        i = args.index("--classify")
        classify_dir = pathlib.Path(args[i + 1])
        del args[i:i + 2]
    ref_dir, cxx_dir = pathlib.Path(args[0]), pathlib.Path(args[1])
    rtol = float(args[2]) if len(args) > 2 else 1e-9
    classify_dir = classify_dir or cxx_dir

    tied = ties()
    n_same = n_diff = n_tie = n_unsupported = n_ref_crash = 0
    n_degenerate = n_invalid = 0
    one_sided = {}
    for cxx in sorted(cxx_dir.glob("*.txt")):
        ref = ref_dir / cxx.name
        if not ref.exists():
            continue
        r, c = parse(ref), parse(cxx)
        name = cxx.stem
        cls = classify_dir / cxx.name
        sos = parse(cls).get("sos_ties", "0") if cls.exists() else "0"
        degenerate = sos != "0"
        cxx_ok = c.get("status", "").startswith("ok")
        if c.get("status", "").startswith("unsupported"):
            n_unsupported += 1
            continue
        # Summaries made outside run-cases.sh have no exit line.
        if r.get("exit", "0") not in ("0", "1"):
            n_ref_crash += 1
            note = ""
            if cxx_ok and not valid(c):
                n_invalid += 1
                note = "  INVALID: not closed, manifold and orientable"
            print(f"C crashed ({r.get('exit')}), C++: status {c.get('status', '?')}  {name}{note}")
            continue

        diffs = []
        rs, cs = r.get("status", "?").split()[0], c.get("status", "?").split()[0]
        if rs != cs:
            diffs.append(f"status: C {r.get('status')} | C++ {c.get('status')}")
        for key in r.keys() | c.keys():
            if key in SKIP:
                continue
            if key not in r or key not in c:
                side = "C" if key in r else "C++"
                one_sided.setdefault((key, side), 0)
                one_sided[(key, side)] += 1
                continue
            vr, vc = values(r[key]), values(c[key])
            if len(vr) != len(vc) or any(nr != nc or not close(xr, xc, rtol)
                                         for (nr, xr), (nc, xc) in zip(vr, vc)):
                diffs.append(f"{key}: C {r[key]} | C++ {c[key]}")
        if diffs and degenerate:
            n_degenerate += 1
            if cxx_ok and not valid(c):
                n_invalid += 1
                print(f"INVALID {name}: degenerate ({sos} ties), C++ result not closed, manifold and orientable")
            else:
                print(f"DEGEN {name}: {sos} ties, {len(diffs)} lines differ, C {rs}, C++ {cs}")
        elif diffs and name in tied and rs == cs:
            n_tie += 1
            print(f"TIE  {name}: {len(diffs)} lines differ; compared by distance in compare.sh")
        elif diffs:
            n_diff += 1
            print(f"DIFF {name}")
            for d in sorted(diffs):
                print(f"    {d}")
        else:
            n_same += 1

    for (key, side), n in sorted(one_sided.items()):
        print(f"note: '{key}' only from {side}, in {n} cases")
    print(f"same {n_same}, different {n_diff}, tied {n_tie}, degenerate {n_degenerate} "
          f"(invalid {n_invalid}), unsupported in C++ {n_unsupported}, "
          f"C crashed {n_ref_crash}  (rtol {rtol})")
    return 1 if n_diff or n_invalid else 0


if __name__ == "__main__":
    sys.exit(main())
