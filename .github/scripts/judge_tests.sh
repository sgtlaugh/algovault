#!/bin/bash
#
# Runs judge tests locally: each file in judge_tests solves a Library Checker or Aizu problem with a library file,
# and competitive-verifier runs it on the judge's own test data (Library Checker data is generated locally)
#
# Usage: .github/scripts/judge_tests.sh [files or directories under judge_tests ...]   (default: all of judge_tests)
# Needs: pip install competitive-verifier==4.2.2   (or COMPETITIVE_VERIFIER=/path/to/competitive-verifier)
#        new judge test files must be git added first, the resolver only sees tracked files
#
set -uo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
CV="${COMPETITIVE_VERIFIER:-competitive-verifier}"
command -v "$CV" > /dev/null || { echo "competitive-verifier not found: pip install competitive-verifier==4.2.2"; exit 2; }

cd "$ROOT" || exit 1
OUT=.competitive-verifier/local
mkdir -p "$OUT"
targets=("$@")
[[ ${#targets[@]} -eq 0 ]] && targets=(judge_tests)

"$CV" oj-resolve --include "${targets[@]}" --config .competitive-verifier/config.toml > "$OUT/verify_files.json" || exit 1
"$CV" verify --verify-json "$OUT/verify_files.json" --output "$OUT/result.json" --check-error
