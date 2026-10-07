#!/bin/bash
#
# self (default): compiles and runs every C/C++ self-test in code_library under AddressSanitizer and
#                 UndefinedBehaviorSanitizer, then runs every Python file
# stress:         does the same for every test in stress_tests, then reports library files without a stress test
#
# Reports all failures instead of stopping at the first one.
#
# Usage: .github/scripts/run_tests.sh [self|stress]       (override compilers with CC=... CXX=...)
#        STRESS_SEED / STRESS_SCALE are passed through to the stress tests, see stress_tests/common.h
#
set -uo pipefail

MODE="${1:-self}"
[[ "$MODE" == "self" || "$MODE" == "stress" ]] || { echo "usage: $0 [self|stress]"; exit 2; }

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
CC="${CC:-gcc}"
CXX="${CXX:-g++}"
PYTHON="${PYTHON:-python3}"
TEST_TIMEOUT="${TEST_TIMEOUT:-300}"

# Sanitizers abort on the first error so UB fails the test, the glibc++ assertions catch out of range container access
SAN_FLAGS=(-O1 -g -fno-omit-frame-pointer "-fsanitize=address,undefined" -fno-sanitize-recover=all -Wall -Wextra)
CFLAGS=(-std=c11 "${SAN_FLAGS[@]}")
CXXFLAGS=(-std=c++17 "${SAN_FLAGS[@]}" -D_GLIBCXX_ASSERTIONS)
export UBSAN_OPTIONS="print_stacktrace=1"

# Self-tests that cannot run in CI, with the reason
SKIP=(
    "code_library/hacking/anti_double_hash.cpp"  # the self-test runs a full collision search, ~150 s
)

# Library files that intentionally have no stress test yet, with the reason
STRESS_SKIP=(
    "code_library/hacking/anti_double_hash.cpp"  # a collision search takes ~150 s, too slow to verify
    "code_library/2SAT_tarjan.cpp"               # lexicographic 2SAT variant, pending a rename and doc review
)

# Contest judges usually give 8 MB of stack, tests must pass with it
ulimit -s 8192

BUILD="$(mktemp -d)"
trap 'rm -rf "$BUILD"' EXIT

passed=0 failed=0 skipped=0 warned=0
failures=()

contains(){  # value, array elements...
    local value="$1" item
    shift
    for item in "$@"; do [[ "$item" == "$value" ]] && return 0; done
    return 1
}

report_failure(){  # file, reason, log
    failures+=("$1 ($2)")
    failed=$((failed + 1))
    echo "FAIL  $1  ($2)"
    if [[ -n "${GITHUB_ACTIONS:-}" ]]; then
        echo "::group::log for $1"; tail -n 60 "$3"; echo "::endgroup::"
        echo "::error file=$1::$2"
    else
        tail -n 30 "$3" | sed 's/^/      /'
    fi
}

run_one(){  # repo relative path
    local rel="$1" src="$ROOT/$1" exe log start secs rc link=()
    exe="$BUILD/$(echo "$rel" | tr '/' '_').out"
    log="$exe.log"

    if contains "$rel" "${SKIP[@]}"; then
        echo "SKIP  $rel"
        skipped=$((skipped + 1))
        return
    fi

    start=$(date +%s)
    if [[ "$rel" == *.py ]]; then
        timeout "$TEST_TIMEOUT" "$PYTHON" "$src" < /dev/null > "$log" 2>&1
        rc=$?
    else
        # A test can ask for extra libraries with a line like: // LINK: -lgmpxx -lgmp
        read -ra link <<< "$(sed -n 's|^// LINK: *||p' "$src")"
        local compile=("$CXX" "${CXXFLAGS[@]}" -o "$exe" "$src" "${link[@]}")
        [[ "$rel" == *.c ]] && compile=("$CC" "${CFLAGS[@]}" -o "$exe" "$src" -lm "${link[@]}")
        if ! "${compile[@]}" > "$log" 2>&1; then
            report_failure "$rel" "compile error" "$log"
            return
        fi
        grep -q "warning:" "$log" && warned=$((warned + 1))
        timeout "$TEST_TIMEOUT" "$exe" < /dev/null >> "$log" 2>&1
        rc=$?
    fi
    secs=$(( $(date +%s) - start ))

    if [[ $rc -eq 124 ]]; then
        report_failure "$rel" "timed out after ${TEST_TIMEOUT}s" "$log"
    elif [[ $rc -ne 0 ]]; then
        report_failure "$rel" "exit code $rc" "$log"
    else
        echo "ok    $rel  (${secs}s)"
        passed=$((passed + 1))
    fi
}

check_coverage(){
    local lib base missing=()
    while IFS= read -r lib; do
        contains "$lib" "${STRESS_SKIP[@]}" && continue
        base="stress_tests/${lib#code_library/}"
        base="${base%.*}"
        [[ -e "$base.cpp" || -e "$base.c" || -e "$base.py" ]] || missing+=("$lib")
    done < <(find code_library -type f \( -name '*.c' -o -name '*.cpp' -o -name '*.py' \) | sort)

    echo "stress test coverage: ${#missing[@]} library files without a stress test"
    [[ ${#missing[@]} -eq 0 ]] && return
    printf '  missing: %s\n' "${missing[@]}"
    failed=$((failed + ${#missing[@]}))
    failures+=("${missing[@]/%/ (no stress test)}")
}

echo "$("$CXX" --version | head -1) | $("$CC" --version | head -1) | $("$PYTHON" --version)"

cd "$ROOT" || exit 1
TESTS_DIR="code_library"
[[ "$MODE" == "stress" ]] && TESTS_DIR="stress_tests"

while IFS= read -r file; do
    run_one "$file"
done < <(find "$TESTS_DIR" -type f \( -name '*.c' -o -name '*.cpp' -o -name '*.py' \) ! -name 'stress.py' | sort)

echo
[[ "$MODE" == "stress" ]] && check_coverage
echo "passed: $passed  failed: $failed  skipped: $skipped  (files with compiler warnings: $warned)"
if [[ $failed -gt 0 ]]; then
    printf '  failed: %s\n' "${failures[@]}"
    exit 1
fi
