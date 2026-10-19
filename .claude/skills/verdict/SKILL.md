---
name: verdict
description: >
  Solve an algovault practice problem on its template in the user's own contest style, check it, hand it over for
  submission, and once the user reports it accepted, commit it under practice/ with the ✅ link. POJ, UVa and HDU
  solutions get a C++11 copy of the template.
trigger: /verdict
allowed-tools: Read, Write, Edit, Glob, Grep, WebFetch, Bash(git grep:*), Bash(git diff:*), Bash(git log:*), Bash(git status:*), Bash(git show:*), Bash(ls:*), Bash(g++:*), Bash(python3:*), Bash(bash .claude/skills/verdict/scripts/check_solution.sh:*), Bash(python3 .claude/skills/verdict/scripts/commit_one.py:*), Bash(env CHANGED_SINCE=:*), Bash(mkdir:*), Bash(cmp:*), Bash(diff:*), Bash(setsid:*), Bash(xclip:*), Bash(/tmp/verdict/:*)
argument-hint: "<problem URL or statement> [<topic>/<template file>] [--commit]"
---

# /verdict - an accepted practice solution in the user's style

Read the practice bullets of `CONTRIBUTING.md` (`## Adding a template`) and all of `STYLE.md` next to this file on every run. Scratch files go in `/tmp/verdict/`. Never look up the user's old solutions: solve the problem fresh.

## Steps

1. **Problem.** Fetch the statement; when the judge blocks fetching, ask the user to paste it instead of guessing. Find its entry in `code_library/<topic>/README.md`: the `#### ` section header and the exact entry line. A problem missing from the README needs the user's approval before it is added.
2. **Template.** Copy the section's template from `code_library/` byte for byte up to its `int main(){`, never from `practice/`. A solution on two templates copies both, in the order they are used.
3. **C++11 for POJ, UVa and HDU.** These judges compile C++11 with GCC 5. Change the copy only where `check_solution.sh` reports an error, with the smallest fix from the table below, and keep every other line identical.
4. **Solve.** Correctness first, then speed: pick the algorithm that fits the limits with room, and use the template's API as its header documents it. Write everything after the copies to `STYLE.md`: no comments, the user's habits, blank lines between steps.
5. **Check:**
   ```bash
   bash .claude/skills/verdict/scripts/check_solution.sh <solution> <template> [<template> ...]
   ```
   It must print `OK`. Then run the binary it built, `/tmp/verdict/solution`:
   - every sample, compared with `diff` against the expected output, exactly
   - one random input at the maximum limits, timed, at most half the time limit
   - a brute force on small random inputs, compared with `diff`, unless `main` only reads the input, makes the template's calls and prints their results
6. **Hand over.** Copy the file to the clipboard, then confirm it in a separate call, since a plain `xclip` dies with its shell:
   ```bash
   setsid -f xclip -selection clipboard -i <solution> < /dev/null > /dev/null 2>&1
   xclip -selection clipboard -o | cmp - <solution>
   ```
   Wait for the user's verdict. Wrong answer or time limit: fix, check again and hand over again.
7. **Commit**, only after the user reports it accepted and said commit, or passed `--commit`. Save it as `practice/<topic>/<judge>_<code>.cpp` (naming per CONTRIBUTING, `_<template>` suffix when the problem has solutions on several templates), run `env CHANGED_SINCE=HEAD JOBS=4 bash .github/scripts/run_tests.sh self` until it exits 0, then:
   ```bash
   python3 .claude/skills/verdict/scripts/commit_one.py <readme> '<section header>' '<entry line>' <solution> '<message>'
   ```
   It stages only the solution and its README line, so another session's edits stay out. Message `feat: add <judge> <code> <title> on the <template>`, at most 72 characters, no body. No AI attribution, never `--no-verify`, never amend, never push.

## C++11 fixes for POJ, UVa and HDU

| Error | Smallest fix |
|---|---|
| `structured bindings only available with -std=c++17` | `.first` / `.second`, or `tie(a, b) = f()`; a binding that writes back becomes a reference to the pair |
| `wrong number of template arguments`, `missing template arguments` | Spell them: `SegmentTree<int> st(a)` |
| `undefined reference to X::NAME` | Define it after the struct: `constexpr long long X::NAME;`, or `template <typename T> constexpr int X<T>::NAME;` |
| `inline variables are only available` | Drop `inline` and define it once after the struct as above |
| `optional does not name a type` | Delete the function when the solution never calls it, else return `pair<bool, T>` |
| `gcd` or `lcm` not declared | `__gcd(a, b)`, `a / __gcd(a, b) * b`, with `abs` on any argument that can be negative: `__gcd(-6, 4)` is -2 where `gcd` gives 2 |
| `is_same_v`, `void_t`, `if constexpr` | Keep only the branch the solution's types take, or `is_same<A, B>::value` |
| `auto` in a lambda parameter or as a return type | The real type; a lambda that recursed through an `auto&& self` parameter becomes a `function<...>` calling itself by name |
| `converting to tuple from initializer list would use explicit constructor` (GCC 5) | `make_tuple(...)` instead of `{...}` |
| `constexpr` function body not a single return, call to non-constexpr | Drop that `constexpr` unless a compile-time constant needs it |
| `clamp` not declared | `max(lo, min(x, hi))` |
| `could not convert {...}` for a struct with default member values | Add a constructor taking the fields |

## Output

The report MUST contain, in this order:

1. **Problem:** judge, code, title, and the README section and entry.
2. **Approach:** the algorithm in two or three plain sentences and why it fits the limits.
3. **C++11 changes:** each changed line of the template copy with its table row, or "none".
4. **Checks:** the `check_solution.sh` line, each sample's result, and the max input's time.
5. **Status:** "in clipboard, waiting for the verdict", or the commit sha and subject.
