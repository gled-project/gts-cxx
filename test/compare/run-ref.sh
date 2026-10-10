#!/bin/bash
# Runs every case through the C reference; see run-cases.sh.
# usage: test/compare/run-ref.sh [OUTDIR]     (default test/compare/out/ref)
cd "$(dirname "$0")/../.." || exit 1
exec test/compare/run-cases.sh test/compare/gts-ref "${1:-test/compare/out/ref}"
