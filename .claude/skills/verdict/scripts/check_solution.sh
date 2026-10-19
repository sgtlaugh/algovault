#!/bin/bash
# check_solution.sh <solution> <template> [<template> ...]
# Checks a practice solution before it is submitted: the template copies come first and unchanged (POJ, UVa and HDU may adapt
# them for C++11), the code after them has no comments and no warnings, and the file builds for its judge
set -uo pipefail
solution=$1; shift
bad=0
cxx="${CXX:-g++}"

prefix=""
for template in "$@"; do
    prefix+="$(awk '/^int main\(\)\{/{exit} {print}' "$template")"$'\n'
done
copy_lines=$(grep -c '' <<< "$prefix")
copy_lines=$((copy_lines - 1))
old_judge=0
[[ "$(basename "$solution")" =~ ^(poj|uva|hdu)_ ]] && old_judge=1

if [[ "$(head -n "$copy_lines" "$solution")"$'\n' != "$prefix" ]]; then
    changed=$(diff <(printf '%s' "$prefix") <(head -n "$copy_lines" "$solution") | grep -c '^[<>]')
    if [[ $old_judge -eq 1 ]]; then
        echo "note: template copy differs in $changed lines, allowed only for C++11"
    else
        echo "FAIL: template copy differs from code_library in $changed lines"
        bad=1
    fi
fi

comments=$(awk -v n="$copy_lines" 'NR > n && /\/\/|\/\*/ {print NR": "$0}' "$solution")
[[ -n "$comments" ]] && { echo "FAIL: comments after the template copy:"; echo "$comments"; bad=1; }

warnings=$("$cxx" -std=c++17 -Wall -Wextra -fsyntax-only "$solution" 2>&1 | grep -P "^\Q$solution\E:\d+:\d+: warning" |
    awk -F: -v n="$copy_lines" '$2 > n')
[[ -n "$warnings" ]] && { echo "FAIL: warnings after the template copy:"; echo "$warnings"; bad=1; }

if [[ $old_judge -eq 1 ]]; then
    # GCC 5 in C++11 mode on the judge; -O0 linking catches static constexpr members used without a definition
    flags=(-std=c++11 -Werror=c++14-extensions -Werror=c++17-extensions -O0)
else
    flags=(-std=c++17 -O2)
fi
mkdir -p /tmp/verdict
out=$("$cxx" "${flags[@]}" "$solution" -o /tmp/verdict/solution 2>&1) || { echo "FAIL: does not build with ${flags[*]}:"; grep -m5 -E 'error|undefined' <<< "$out"; bad=1; }

# GCC 14 flags miss GCC 5 library differences, like tuple's explicit constructor or glibc's minor() macro
if [[ $old_judge -eq 1 ]] && docker image inspect gcc:5 > /dev/null 2>&1; then
    dir="$(cd "$(dirname "$solution")" && pwd)"
    out=$(docker run --rm -v "$dir":/src:ro gcc:5 g++ -std=c++11 -O0 "/src/$(basename "$solution")" -o /dev/null 2>&1) ||
        { echo "FAIL: does not build on GCC 5:"; grep -m5 -E 'error|undefined' <<< "$out"; bad=1; }
elif [[ $old_judge -eq 1 ]]; then
    echo "note: docker pull gcc:5 checks against the judge's real compiler"
fi

[[ $bad -eq 0 ]] && echo "OK: built as /tmp/verdict/solution with ${flags[*]}"
exit $bad
