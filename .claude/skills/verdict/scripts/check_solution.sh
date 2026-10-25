#!/bin/bash
# check_solution.sh <solution> <template> [<template> ...]
# Checks a practice solution before it is submitted: the template copies come first and unchanged (POJ, UVa and HDU may adapt
# them for C++11), fast_io.cpp is copied last, the code after them has no comments, no stdio or iostream and no warnings,
# and the file builds for its judge
set -uo pipefail
solution=$1; shift
bad=0
cxx="${CXX:-g++}"

# The trailing x keeps each copy's closing blank line, which command substitution would strip and so misjoin two copies
prefix=""
for template in "$@"; do
    part=$(awk '/^int main\(\)\{/{exit} {print}' "$template"; printf x)
    prefix+="${part%x}"
done
copy_lines=$(printf '%s' "$prefix" | grep -c '')
old_judge=0
[[ "$(basename "$solution")" =~ ^(poj|uva|hdu)_ ]] && old_judge=1

if ! cmp -s <(head -n "$copy_lines" "$solution") <(printf '%s' "$prefix"); then
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

# All input and output go through fio: scanf, printf, cin or cout next to its buffers would reorder the stream
[[ " $* " == *"misc/fast_io.cpp "* ]] || { echo "FAIL: copy code_library/misc/fast_io.cpp last, after the other templates"; bad=1; }
stdio=$(grep -nP '\b(scanf|printf|puts|gets|getchar|putchar|cin|cout)\b' "$solution" | awk -F: -v n="$copy_lines" '$1 > n')
[[ -n "$stdio" ]] && { echo "FAIL: stdio or iostream after the template copy, use fio read and write:"; echo "$stdio"; bad=1; }

# Layout rules from STYLE.md that a script can check without judgment; the solution part starts after the copies
layout=$(python3 -I - "$solution" "$copy_lines" "$old_judge" <<'PY'
import re, sys
path, start, old_judge = sys.argv[1], int(sys.argv[2]), sys.argv[3] == '1'
lines = open(path).read().split('\n')
types = r'(const |static |unsigned |signed )*(int|long|short|double|float|char|bool|string|vector|array|pair|map|set|' \
        r'unordered_map|unordered_set|queue|deque|stack|priority_queue|bitset|uint32_t|uint64_t|size_t)\b'
decl = re.compile(r'^    ' + types + r'[^(]*;$')
fail, note = [], []

for i in range(start, len(lines)):
    no, line = i + 1, lines[i]
    if line == '' and i > start and lines[i - 1] == '':
        fail.append(f'{no}: two blank lines in a row')
    m = re.match(r'#define (\w+)\(', line)
    if m and m.group(1) not in ('clr', 'dbg'):
        fail.append(f'{no}: macro {m.group(1)}() takes arguments, write a function')
    if re.search(r"(\b[a-z_]\w*(\[[^]]*\])?\s*[-+]\s*|[<>]=?\s*)(48|57|65|90|97|122)\b", line):
        note.append(f'{no}: an ASCII code where a character literal would read better: {line.strip()}')

    if re.match(r'^\S.*\)\{$', line) and not line.startswith(('struct', 'namespace')):
        end = i + 1
        while end < len(lines) and lines[end] != '}':
            end += 1
        is_main = line.startswith('int main')
        # Blank lines split steps, so a short helper (one declaration, one loop, the return) needs none
        if not is_main and end - i - 1 < 8:
            continue
        j = i + 1
        while j < end and decl.match(lines[j]):
            j += 1
        if j > i + 1 and lines[j] != '':
            fail.append(f'{j + 1}: no blank line after the declarations at the top of the function')
        if not is_main and re.match(r'^    return \S', lines[end - 1]) and lines[end - 2] == '    }':
            fail.append(f'{end}: no blank line before the final return after a loop')

if old_judge:
    text = '\n'.join(lines[start:])
    if re.search(r'\bread\(t\);', text) and not re.search(r'while \(read\(', text):
        note.append('the test count is read once, UVa, POJ and HDU files often repeat the block: while (read(t)){ while (t--){')

for f in fail:
    print('FAIL ' + f)
for n in note:
    print('note ' + n)
PY
)
if [[ -n "$layout" ]]; then
    grep '^FAIL' <<< "$layout" | sed 's/^FAIL /FAIL: layout, line /' && bad=1
    grep '^note' <<< "$layout" | sed 's/^note /note: /'
fi

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
# The -O0 build only proves linking: the binary that gets timed is rebuilt at -O2, an -O0 one runs about 10x slower
if [[ $old_judge -eq 1 && $bad -eq 0 ]]; then
    flags=(-std=c++11 -O2)
    out=$("$cxx" "${flags[@]}" "$solution" -o /tmp/verdict/solution 2>&1) || { echo "FAIL: does not build with ${flags[*]}:"; grep -m5 -E 'error|undefined' <<< "$out"; bad=1; }
fi

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
