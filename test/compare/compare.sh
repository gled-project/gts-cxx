#!/bin/bash
# Runs all cases through the C reference and the C++ library and compares the
# summaries (compare.py). Then reads every surface the C++ library wrote back
# through the C library, and compares those summaries with the ones of the
# surfaces the C library wrote itself: this checks the .gts writer. The cases
# in ties.txt are compared by the distance between the two surfaces, and
# those in fixed.txt by the check given there.
#
# usage: test/compare/compare.sh [OUTDIR]     (default test/compare/out)

set -u
cd "$(dirname "$0")/../.." || exit 1
OUT=${1:-test/compare/out}

test/compare/run-cases.sh test/compare/gts-ref "$OUT/ref" 2>/dev/null
test/compare/run-cases.sh test/compare/gts-cxx "$OUT/cxx" 2>/dev/null

echo "=== C and C++ summaries"
python3 test/compare/compare.py "$OUT/ref" "$OUT/cxx"
rc=$?

echo "=== cases with tied costs, by distance between the C and the C++ surface"
rc3=0
while read -r name tol; do
  case "$name" in ''|'#'*) continue;; esac
  a=$OUT/ref/$name.gts b=$OUT/cxx/$name.gts
  if [ ! -f "$a" ] || [ ! -f "$b" ]; then echo "skipped $name (no surface from both)"; continue; fi
  d=$(setarch -R test/compare/gts-ref distance "$a" "$b" 2>/dev/null |
      sed -n 's/^distance_[ab]_to_[ab] .*max=\([^ ]*\).*/\1/p' | sort -g | tail -1)
  if python3 -c "import sys; sys.exit(0 if float('$d') <= float('$tol') else 1)"; then
    echo "ok   $name: max distance $d <= $tol"
  else
    echo "FAIL $name: max distance $d > $tol"; rc3=1
  fi
done < test/compare/ties.txt

echo "=== surfaces written by C and by C++, both read by C"
mkdir -p "$OUT/rt-ref" "$OUT/rt-cxx"
rm -f "$OUT"/rt-ref/*.txt "$OUT"/rt-cxx/*.txt
for f in "$OUT"/cxx/*.gts; do
  b=$(basename "$f" .gts)
  [ -f "$OUT/ref/$b.gts" ] || continue
  setarch -R test/compare/gts-ref stats "$OUT/ref/$b.gts" > "$OUT/rt-ref/$b.txt" 2>/dev/null
  setarch -R test/compare/gts-ref stats "$f" > "$OUT/rt-cxx/$b.txt" 2>/dev/null
done
python3 test/compare/compare.py "$OUT/rt-ref" "$OUT/rt-cxx" --classify "$OUT/cxx"
rc2=$?

exit $(( rc || rc2 || rc3 ))
