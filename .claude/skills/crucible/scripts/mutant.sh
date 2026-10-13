#!/usr/bin/env bash
# Plants one bug in a scratch copy of a library file and confirms that its self-test or stress test fails.
# Run from the repo root: bash .claude/skills/crucible/scripts/mutant.sh <path under code_library> '<sed expression>'
# Exit 0: the bug was caught. Exit 1: it survived. Exit 3: the mutant is invalid (no change, or it does not compile).
set -uo pipefail
REL="$1" EXPR="$2"
[[ -f "code_library/$REL" ]] || { echo "no such file: code_library/$REL (run from the repo root)"; exit 2; }

# The CI flags, read from run_tests.sh so this never drifts from CI
eval "$(grep -E '^SAN_FLAGS=' .github/scripts/run_tests.sh)"
export LD_LIBRARY_PATH="$(dirname "$(g++ -print-file-name=libstdc++.so)")"
ulimit -s 8192

M=$(mktemp -d /tmp/crucible_mutant.XXXX)
cp -r stress_tests "$M/" && mkdir -p "$M/code_library" && cp -r code_library/. "$M/code_library/"
sed "$EXPR" "code_library/$REL" > "$M/code_library/$REL"
if cmp -s "code_library/$REL" "$M/code_library/$REL"; then
    echo "INVALID MUTANT: the sed expression changed nothing"; exit 3
fi
echo "planted: $(diff "code_library/$REL" "$M/code_library/$REL" | grep '^>' | head -2 | sed 's/^> *//')"

cd "$M" || exit 2
base="stress_tests/${REL%.*}"
caught=""

run(){  # label, then the command to run
    local label="$1"; shift
    timeout 300 "$@" < /dev/null > "$M/$label.out" 2>&1
    local rc=$?
    if [[ $rc -ne 0 ]]; then
        caught+="$label "
        echo "$label failed (exit $rc): $(grep -m1 -E 'Assertion|AssertionError|runtime error|ERROR' "$M/$label.out" | cut -c1-200)"
    fi
}

if [[ "$REL" == *.py ]]; then
    run self python3 "code_library/$REL"
    run stress python3 "$base.py"
else
    compiler=(g++ -std=c++17); [[ "$REL" == *.c ]] && compiler=(gcc -std=c11)
    "${compiler[@]}" "${SAN_FLAGS[@]}" $([[ "$REL" == *.c ]] || echo -D_GLIBCXX_ASSERTIONS) "code_library/$REL" -o "$M/self" -lm > "$M/self.log" 2>&1 \
        || { echo "INVALID MUTANT: the library does not compile"; tail -5 "$M/self.log"; exit 3; }
    run self "$M/self"

    stress="$base.cpp"
    link=$(sed -n 's|^// LINK: *||p' "$stress")
    g++ -std=c++17 "${SAN_FLAGS[@]}" -D_GLIBCXX_ASSERTIONS "$stress" $link -o "$M/stress" > "$M/stress.log" 2>&1 \
        || { echo "INVALID MUTANT: the stress test does not compile"; tail -5 "$M/stress.log"; exit 3; }
    run stress "$M/stress"
fi

if [[ -z "$caught" ]]; then echo "MUTANT SURVIVED: neither the self-test nor the stress test failed"; exit 1; fi
echo "mutant caught by: $caught"
