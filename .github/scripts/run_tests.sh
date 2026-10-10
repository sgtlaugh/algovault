#!/bin/bash
#
# self (default): compiles and runs every C++ self-test in code_library under AddressSanitizer and
#                 UndefinedBehaviorSanitizer, then runs every Python file
# stress:         does the same for every test in stress_tests, then reports library files without a stress test
# both modes also report library files missing from the README index, README ✔ marks that disagree with judge_tests,
# templates whose header does not state their complexity, and whitespace or style errors in tracked files
#
# Reports all failures instead of stopping at the first one.
#
# Usage: .github/scripts/run_tests.sh [self|stress]       (override the compiler with CXX=...)
#        STRESS_SEED / STRESS_SCALE are passed through to the stress tests, see stress_tests/common.h
#        JOBS=n             run n tests at a time (default: all cores)
#        BUILD=fast         -O2 without sanitizers, for high-iteration runs where wrong answers, not UB, are the target
#        CHANGED_SINCE=ref  only run tests whose library file, stress test or a library file that stress test includes
#                           changed since ref, everything when a shared file changed or ref is unknown
#        TEST_BUDGET=s      fail tests that pass but take longer than s seconds, so slow tests cannot pile up
#
set -uo pipefail

MODE="${1:-self}"
[[ "$MODE" == "self" || "$MODE" == "stress" ]] || { echo "usage: $0 [self|stress]"; exit 2; }

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
CXX="${CXX:-g++}"
PYTHON="${PYTHON:-python3}"
TEST_TIMEOUT="${TEST_TIMEOUT:-300}"
TEST_BUDGET="${TEST_BUDGET:-0}"
JOBS="${JOBS:-$(nproc)}"
BUILD_MODE="${BUILD:-sanitize}"
[[ "$BUILD_MODE" == "sanitize" || "$BUILD_MODE" == "fast" ]] || { echo "BUILD must be sanitize or fast"; exit 2; }

# Sanitizers abort on the first error so UB fails the test, the glibc++ assertions catch out of range container access
# _FORTIFY_SOURCE is pinned because Ubuntu's gcc enables it by default and other builds don't,
# and it adds warnings (unused scanf results) that -Werror turns into CI-only failures
# -U first: Ubuntu defines it built in, redefining to another value is itself an error under -Werror
SAN_FLAGS=(-O1 -g -fno-omit-frame-pointer "-fsanitize=address,undefined" -fno-sanitize-recover=all -Wall -Wextra -Werror -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=3)
CXXFLAGS=(-std=c++17 "${SAN_FLAGS[@]}" -D_GLIBCXX_ASSERTIONS)
[[ "$BUILD_MODE" == "fast" ]] && CXXFLAGS=(-std=c++17 -O2 -Wall -Wextra -Werror -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=3)
export UBSAN_OPTIONS="print_stacktrace=1"

# Self-tests that cannot run in CI, with the reason
SKIP=(
    "code_library/hacking/anti_double_hash.cpp"  # the self-test runs a full collision search, ~150 s
)

# Library files that intentionally have no stress test yet, with the reason
STRESS_SKIP=(
    "code_library/hacking/anti_double_hash.cpp"  # a collision search takes ~150 s, too slow to verify
)

# Files quoting outside code verbatim, kept diffable against their source
STYLE_SKIP=(
    "stress_tests/graphs/edge_coloring.cpp"  # KACTL EdgeColoring.h as a second reference
    "stress_tests/graphs/stable_marriage.cpp"  # the original Library implementation as a second reference
)

# Changing any of these can affect every test, so CHANGED_SINCE then runs everything
SHARED=(
    "stress_tests/common.h"
    "stress_tests/python/stress.py"
    ".github/scripts/run_tests.sh"
    ".github/workflows/ci.yml"
)

# Contest judges usually give 8 MB of stack, tests must pass with it
ulimit -s 8192

if [[ -z "${RUN_TESTS_WORKER:-}" ]]; then
    BUILD_DIR="$(mktemp -d)"
    trap 'rm -rf "$BUILD_DIR"' EXIT
    export BUILD_DIR
fi

passed=0 failed=0 skipped=0
failures=()

contains(){  # value, array elements...
    local value="$1" item
    shift
    for item in "$@"; do [[ "$item" == "$value" ]] && return 0; done
    return 1
}

report_failure(){  # file, reason, log
    printf 'FAIL\t%s\t%s\n' "$1" "$2" >> "$BUILD_DIR/results"
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
    exe="$BUILD_DIR/$(echo "$rel" | tr '/' '_').out"
    log="$exe.log"

    if contains "$rel" "${SKIP[@]}"; then
        echo "SKIP  $rel"
        printf 'SKIP\t%s\t\n' "$rel" >> "$BUILD_DIR/results"
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
        if ! "${compile[@]}" > "$log" 2>&1; then
            report_failure "$rel" "compile error" "$log"
            return
        fi
        timeout "$TEST_TIMEOUT" "$exe" < /dev/null >> "$log" 2>&1
        rc=$?
    fi
    secs=$(( $(date +%s) - start ))

    if [[ $rc -eq 124 ]]; then
        report_failure "$rel" "timed out after ${TEST_TIMEOUT}s" "$log"
    elif [[ $rc -ne 0 ]]; then
        report_failure "$rel" "exit code $rc" "$log"
    elif [[ $TEST_BUDGET -gt 0 && $secs -gt $TEST_BUDGET ]]; then
        report_failure "$rel" "passed but took ${secs}s, over the ${TEST_BUDGET}s budget" "$log"
    else
        echo "ok    $rel  (${secs}s)"
        printf 'OK\t%s\t\n' "$rel" >> "$BUILD_DIR/results"
    fi
}

check_coverage(){
    local lib base missing=()
    while IFS= read -r lib; do
        contains "$lib" "${STRESS_SKIP[@]}" && continue
        base="stress_tests/${lib#code_library/}"
        base="${base%.*}"
        [[ -e "$base.cpp" || -e "$base.py" ]] || missing+=("$lib")
    done < <(find code_library -type f \( -name '*.cpp' -o -name '*.py' \) | sort)

    echo "stress test coverage: ${#missing[@]} library files without a stress test"
    [[ ${#missing[@]} -eq 0 ]] && return
    printf '  missing: %s\n' "${missing[@]}"
    failed=$((failed + ${#missing[@]}))
    failures+=("${missing[@]/%/ (no stress test)}")
}

check_readme(){
    local lib missing=()
    while IFS= read -r lib; do
        grep -qF "($lib)" README.md || missing+=("$lib")
    done < <(find code_library -type f \( -name '*.cpp' -o -name '*.py' \) | sort)

    echo "readme index: ${#missing[@]} library files not linked from README.md"
    [[ ${#missing[@]} -eq 0 ]] && return
    printf '  missing: %s\n' "${missing[@]}"
    failed=$((failed + ${#missing[@]}))
    failures+=("${missing[@]/%/ (not in README index)}")
}

check_judge_marks(){  # a template carries a README ✔ exactly when a judge test's first library include is that template
    local file lib line target tested=() bad=()
    while IFS= read -r file; do
        lib="$(grep -m1 -oP '#include "(\.\./)+\Kcode_library/[^"]+' "$file")"
        [[ -n "$lib" ]] && tested+=("$lib")
    done < <(git ls-files 'judge_tests/*.cpp')

    while IFS= read -r lib; do
        line="$(grep -F "($lib)" README.md | head -1)"
        target="$(grep -oP '\[✔\]\(\K[^)]+' <<< "$line")"
        if contains "$lib" "${tested[@]}"; then
            [[ -n "$target" && -f "$target" ]] || bad+=("$lib")
        elif [[ -n "$target" ]]; then
            bad+=("$lib")
        fi
    done < <(find code_library -type f \( -name '*.cpp' -o -name '*.py' \) | sort)

    echo "judge marks: ${#bad[@]} index entries whose ✔ disagrees with judge_tests"
    [[ ${#bad[@]} -eq 0 ]] && return
    printf '  bad: %s\n' "${bad[@]}"
    failed=$((failed + ${#bad[@]}))
    failures+=("${bad[@]/%/ (judge mark)}")
}

check_headers(){  # every template opens with a header that states its complexity
    local file bad=()
    while IFS= read -r file; do
        if [[ "$file" == *.py ]]; then
            python3 -c 'import ast, sys; doc = ast.get_docstring(ast.parse(open(sys.argv[1]).read())) or ""; sys.exit("Complexity" not in doc)' "$file" || bad+=("$file")
        else
            awk 'NR == 1 && !/^\/\*\*\*/ {exit 1} /^\*\*\*\// {exit !found} /^ \* Complexity/ {found = 1} END {if (!found) exit 1}' "$file" || bad+=("$file")
        fi
    done < <(git ls-files 'code_library/*.cpp' 'code_library/*.py')

    echo "headers: ${#bad[@]} templates without a header stating their complexity"
    [[ ${#bad[@]} -eq 0 ]] && return
    printf '  bad: %s\n' "${bad[@]}"
    failed=$((failed + ${#bad[@]}))
    failures+=("${bad[@]/%/ (no complexity in header)}")
}

check_whitespace(){  # tabs, trailing whitespace, CRLF or a missing final newline
    local file bad=()
    while IFS= read -r file; do
        if grep -qP '\t|[ \t]$|\r' "$file" || [[ -s "$file" && -n "$(tail -c1 "$file")" ]]; then
            bad+=("$file")
        fi
    done < <(git ls-files '*.cpp' '*.h' '*.py' '*.md' '*.sh' '*.yml')

    echo "whitespace: ${#bad[@]} tracked files with tabs, trailing whitespace, CRLF or no final newline"
    [[ ${#bad[@]} -eq 0 ]] && return
    printf '  bad: %s\n' "${bad[@]}"
    failed=$((failed + ${#bad[@]}))
    failures+=("${bad[@]/%/ (whitespace)}")
}

check_style(){  # ){ and struct X{ open a body, ) {} is empty, keywords take a space, template arguments do not
    local line bad=()
    while IFS= read -r line; do
        contains "${line%%:*}" "${STYLE_SKIP[@]}" || bad+=("$line")
    done < <(git ls-files '*.cpp' '*.h' | xargs grep -nP '\) \{(?!\})|\b(struct|class|namespace)\s+\w+(\s*:\s*[^{;]+)? \{(?!\})|\)\{\}|\b(for|if|while|switch)\(|\b(?!template\b)\w+ <(?=[\w:]+[\s\w:,<>*&]*>)')

    echo "style: ${#bad[@]} lines off the brace or spacing convention"
    [[ ${#bad[@]} -eq 0 ]] && return
    printf '  bad: %s\n' "${bad[@]}"
    failed=$((failed + ${#bad[@]}))
    failures+=("${bad[@]/%/ (style)}")
}

select_tests(){  # prints the tests to run, one per line
    local all changed file lib stress deps
    all="$(find "$TESTS_DIR" -type f \( -name '*.cpp' -o -name '*.py' \) ! -name 'stress.py' | sort)"
    if [[ -z "${CHANGED_SINCE:-}" ]] || ! git cat-file -e "${CHANGED_SINCE}^{commit}" 2> /dev/null; then
        [[ -n "${CHANGED_SINCE:-}" ]] && echo "CHANGED_SINCE=$CHANGED_SINCE is not a known commit, running everything" >&2
        echo "$all"
        return
    fi

    changed="$(git diff --name-only "$CHANGED_SINCE" HEAD)"
    for file in "${SHARED[@]}"; do
        grep -qxF "$file" <<< "$changed" && { echo "$file changed, running everything" >&2; echo "$all"; return; }
    done

    while IFS= read -r file; do
        [[ -z "$file" ]] && continue
        lib="code_library/${file#*/}" stress="stress_tests/${file#*/}"
        deps="$lib"$'\n'"$stress"
        # A stress test may build on other library files, e.g. the hld test runs against segment_tree.cpp
        [[ "$TESTS_DIR" == "stress_tests" && -f "$stress" ]] && deps+=$'\n'"$(grep -oP '#include "(\.\./)+\Kcode_library/[^"]+' "$stress")"
        if grep -qxFf <(grep . <<< "$deps") <<< "$changed"; then echo "$file"; fi
    done <<< "$all"
}

if [[ -n "${RUN_TESTS_WORKER:-}" ]]; then  # one test, run by the parallel driver below
    cd "$ROOT" || exit 1
    run_one "$2"
    exit 0
fi

echo "$("$CXX" --version | head -1) | $("$PYTHON" --version) | build: $BUILD_MODE | jobs: $JOBS"

cd "$ROOT" || exit 1
TESTS_DIR="code_library"
[[ "$MODE" == "stress" ]] && TESTS_DIR="stress_tests"

: > "$BUILD_DIR/results"
tests="$(select_tests)"
echo "running $(grep -c . <<< "$tests") of $(find "$TESTS_DIR" -type f \( -name '*.cpp' -o -name '*.py' \) ! -name 'stress.py' | wc -l) tests"
grep . <<< "$tests" | RUN_TESTS_WORKER=1 xargs -P "$JOBS" -I{} bash "$0" "$MODE" {}

passed=$(grep -c '^OK' "$BUILD_DIR/results")
skipped=$(grep -c '^SKIP' "$BUILD_DIR/results")
failed=$(grep -c '^FAIL' "$BUILD_DIR/results")
failures=()
while IFS=$'\t' read -r _ file reason; do failures+=("$file ($reason)"); done < <(grep '^FAIL' "$BUILD_DIR/results")
selected=$(grep -c . <<< "$tests")
if [[ $((passed + skipped + failed)) -ne $selected ]]; then  # a worker died without reporting, never pass silently
    failures+=("$((selected - passed - skipped - failed)) of $selected selected tests reported no result")
    failed=$((failed + 1))
fi

echo
[[ "$MODE" == "stress" ]] && check_coverage
check_readme
check_judge_marks
check_headers
check_whitespace
check_style
echo "passed: $passed  failed: $failed  skipped: $skipped"
if [[ $failed -gt 0 ]]; then
    printf '  failed: %s\n' "${failures[@]}"
    exit 1
fi
