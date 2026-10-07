#!/bin/bash
#
# Compiles and runs every C/C++ self-test under AddressSanitizer and UndefinedBehaviorSanitizer,
# then runs every Python file. Reports all failures instead of stopping at the first one.
#
# Usage: .github/scripts/run_tests.sh            (override compilers with CC=... CXX=...)
#
set -uo pipefail

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

# Files whose self-test cannot run in CI, with the reason
SKIP=(
    "code_library/hacking/anti_double_hash.cpp"  # the self-test runs a full collision search, ~150 s
)

# Contest judges usually give 8 MB of stack, tests must pass with it
ulimit -s 8192

BUILD="$(mktemp -d)"
trap 'rm -rf "$BUILD"' EXIT

passed=0 failed=0 skipped=0 warned=0
failures=()

is_skipped(){
    local f
    for f in "${SKIP[@]}"; do [[ "$1" == "$f" ]] && return 0; done
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
    local rel="$1" src="$ROOT/$1" exe log start secs rc
    exe="$BUILD/$(echo "$rel" | tr '/' '_').out"
    log="$exe.log"

    if is_skipped "$rel"; then
        echo "SKIP  $rel"
        skipped=$((skipped + 1))
        return
    fi

    start=$(date +%s)
    if [[ "$rel" == *.py ]]; then
        timeout "$TEST_TIMEOUT" "$PYTHON" "$src" < /dev/null > "$log" 2>&1
        rc=$?
    else
        local compile=("$CXX" "${CXXFLAGS[@]}" -o "$exe" "$src")
        [[ "$rel" == *.c ]] && compile=("$CC" "${CFLAGS[@]}" -o "$exe" "$src" -lm)
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

echo "$("$CXX" --version | head -1) | $("$CC" --version | head -1) | $("$PYTHON" --version)"

cd "$ROOT" || exit 1
while IFS= read -r file; do
    run_one "$file"
done < <(find code_library -type f \( -name '*.c' -o -name '*.cpp' -o -name '*.py' \) | sort)

echo
echo "passed: $passed  failed: $failed  skipped: $skipped  (files with compiler warnings: $warned)"
if [[ $failed -gt 0 ]]; then
    printf '  failed: %s\n' "${failures[@]}"
    exit 1
fi
