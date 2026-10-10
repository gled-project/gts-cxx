#!/bin/bash
# Runs every case through one program: gts-ref, which uses the C GTS library,
# or gts-cxx, which uses the C++ one. For each case NAME it writes, in the
# output directory:
#   NAME.txt   the summary printed by gts-ref, then "exit N"
#   NAME.err   anything gts-ref printed to stderr
#   NAME.gts   the resulting surface, when the command produces one
# Also runs every surface in test/data/gts-boolean through stats, and every
# pair in test/data/gts-boolean/tests through union, intersection and
# difference.
#
# GTS keeps faces and vertices in hash tables keyed by pointer, so its results
# depend on memory addresses. Coarsening chooses differently between equal
# costs, and some degenerate booleans abort in one run and succeed in the
# next. Each case therefore runs under "setarch -R", which turns off address
# space randomisation, so that the reference repeats exactly from run to run.
# The output directory in the "command" line is replaced by @OUT@.
#
# In cases.txt, @DATA@ stands for test/data, @OUT@ for this program's output
# directory and @REF@ for gts-ref's, so that both programs can read a
# surface written by the C library.
#
# A program that does not support a case yet exits with 3, which is counted
# as unsupported.
#
# usage: test/compare/run-cases.sh PROGRAM OUTDIR

set -u
cd "$(dirname "$0")/../.."

PROG=${1:?usage: run-cases.sh PROGRAM OUTDIR}
OUT=${2:?usage: run-cases.sh PROGRAM OUTDIR}
DATA=test/data
REF=$(dirname "$OUT")/ref
TIMEOUT=${TIMEOUT:-120}

mkdir -p "$OUT"
rm -f "$OUT"/*.txt "$OUT"/*.err "$OUT"/*.gts

n_ok=0; n_error=0; n_crash=0; n_unsupported=0
crashed=()

run_case()
{
  local name=$1; shift
  # An input made by an earlier case that this program does not support.
  local x
  for x in "$@"; do
    if [[ $x == "$OUT"/* && ! -f $x ]]; then
      printf 'command %s\nstatus unsupported input %s missing\nexit 3\n' "$*" "${x#"$OUT"/}" > "$OUT/$name.txt"
      sed -i "1s|$OUT/|@OUT@/|g" "$OUT/$name.txt"
      n_unsupported=$((n_unsupported + 1))
      return
    fi
  done
  setarch -R timeout "$TIMEOUT" "$PROG" -o "$OUT/$name.gts" "$@" > "$OUT/$name.txt" 2> "$OUT/$name.err"
  local rc=$?
  echo "exit $rc" >> "$OUT/$name.txt"
  sed -i "1s|$OUT/|@OUT@/|g" "$OUT/$name.txt"
  [ -s "$OUT/$name.gts" ] || rm -f "$OUT/$name.gts"
  [ -s "$OUT/$name.err" ] || rm -f "$OUT/$name.err"
  if   [ $rc -eq 0 ]; then n_ok=$((n_ok + 1))
  elif [ $rc -eq 1 ]; then n_error=$((n_error + 1))
  elif [ $rc -eq 3 ]; then n_unsupported=$((n_unsupported + 1))
  else n_crash=$((n_crash + 1)); crashed+=("$name (exit $rc)")
  fi
}

# Listed cases
while read -r name rest; do
  case "$name" in ''|'#'*) continue;; esac
  rest=${rest//@DATA@/$DATA}
  rest=${rest//@OUT@/$OUT}
  rest=${rest//@REF@/$REF}
  # shellcheck disable=SC2086
  run_case "$name" $rest
done < test/compare/cases.txt

# Every boolean test surface through stats
for f in "$DATA"/gts-boolean/*; do
  b=$(basename "$f")
  case "$b" in tests|README.md) continue;; esac
  run_case "stats-${b%.gts}" stats "$f"
done

# Every pair in the boolean test list, through three operations
while read -r a b _; do
  case "$a" in ''|'#'*) continue;; esac
  [ -f "$DATA/gts-boolean/$a" ] && [ -f "$DATA/gts-boolean/$b" ] || continue
  for op in union inter diff; do
    run_case "bool-$op-${a%.gts}-${b%.gts}" bool $op "$DATA/gts-boolean/$a" "$DATA/gts-boolean/$b"
  done
done < "$DATA/gts-boolean/tests"

echo "$(basename "$PROG"): ok $n_ok, error $n_error, unsupported $n_unsupported, crash or timeout $n_crash  (in $OUT)"
for c in "${crashed[@]}"; do echo "  crashed: $c"; done
