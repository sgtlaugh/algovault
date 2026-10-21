---
name: crucible
description: >
  Add or rework an algovault template: port an algorithm from a file of the user's own code, reimplement one from a
  web page, or refactor an existing code_library file to be simpler, easier to use or faster. Every result ships with a
  brute-force stress test proven to catch a planted bug, docs from /lucid, and a green CI suite.
trigger: /crucible
allowed-tools: Read, Write, Edit, Glob, Grep, WebFetch, Skill, Bash(git grep:*), Bash(git diff:*), Bash(git log:*), Bash(git status:*), Bash(git rev-parse:*), Bash(git stash create:*), Bash(git add:*), Bash(git commit:*), Bash(ls:*), Bash(g++:*), Bash(python3:*), Bash(env CHANGED_SINCE=:*), Bash(bash .github/scripts/run_tests.sh:*), Bash(bash .github/scripts/judge_tests.sh:*), Bash(bash .claude/skills/crucible/scripts/mutant.sh:*), Bash(env LD_LIBRARY_PATH=:*), Bash(mkdir:*), Bash(cp:*), Bash(/tmp/crucible_bench/:*)
argument-hint: "<file path | URL | code_library file> [--name <file name>] [--commit]"
---

# /crucible - port, reimplement or refactor a template

The repo's rules live in `CONTRIBUTING.md`, read it in full first on every run. Documentation is delegated to `/lucid`. This skill owns what is specific to algorithms: the stress test, the planted bug, the benchmark and the judge test.

**Why the stress test is the core.** Of the 10 files ported in June 2026, 4 shipped silent wrong answers that an October 2026 audit found (sparse 2D Fenwick key collisions above the documented range, an off-by-one at `y2 = n`, a range-set recipe ignoring 0, a suffix array sentinel clash). Self-tests check what the author expected; a brute force checks what is true.

## Mode, from the argument

| Argument | Mode | Result |
|---|---|---|
| A file or directory outside `code_library/` | **port** | a new template from the user's code |
| A URL | **web** | a new template reimplemented from the page |
| A file in `code_library/` | **refactor** | the same template, simpler, easier to use or faster |

One template per run. A directory in port mode means one algorithm spread over its files.

Run shell commands under bash, edit files with Edit or Write (never `sed -i`), and keep scratch work such as benchmarks in `/tmp/crucible_bench/`.

## Steps

1. **Understand.**
   - Read the source in full, two recent `code_library/` files of the same topic for style, and `.github/scripts/run_tests.sh` for the CI flags.
   - port and web: `ls -R code_library/` and `git grep -il <algorithm name>`; if the algorithm already exists, stop and report it, offering refactor mode instead.
   - web: find the page's licence. Code under a permissive licence (CC0, MIT, Apache, public domain) may be adapted, credited in the header. Code under CC BY-SA, GPL or no stated licence is never copied: reimplement it from the idea and cite the page. A run without a licence verdict in the report is not done.
   - refactor: read the file's stress test, judge tests and callers (`git grep -n <struct or function> -- code_library stress_tests judge_tests`), and benchmark the current version at -O2 before changing anything, on its largest documented input, or else the stress test's largest input or the judge problem's limit.
   - Ambiguous or suspicious source (an off-by-one in a known algorithm, a complexity claim the code does not meet): stop and ask. Porting a suspected bug silently is not allowed.
2. **Write the code** to the CONTRIBUTING rules, sections `## API and code` and `## Style`: a struct with per-instance state, 0-based positions with inclusive ranges `[l, r]`, no mutable globals and no generic global names, `T x(args)`. For each item of the edge-case list below that applies, the code handles it or the header states the limit.
   - port and web: the source's API is not a constraint. Design the call a contestant would want, then fit the algorithm to it.
   - refactor: an API change must break old calls at compile time, never change their meaning silently; update every caller. Renaming a template parameter or a local that callers never spell needs no caller change. A new requirement on a type (such as default constructible) counts as an API change and goes in the header's notes.
3. **Simplicity gate**, before any testing. Write the common case as the 1 to 3 lines a contestant would type, then answer each question with yes, or change the code until it is yes:
   - Does the common case need no configuration and no parameter the code could do without?
   - Does every public name say what the call does, with arguments in the order the sentence reads?
   - Does it reuse the names its family already uses (`find`, `same`, `size`, `combine`, `query`, `update`)?
   - Is customizing a matter of passing a lambda, functor or `Ops` struct, never editing the template?
   - Does every call return its answer rather than fill an out-parameter?
   - Is everything public, with no `private:` section and no `class`?
   - Are loops plain and variables named for what they hold, with no abstraction used once and no language machinery that does not prevent a silent wrong answer?
   - Does the struct read top down (types and constants, fields, constructors, the API in header order, helpers in call order), with named options as a nested enum rather than int constants?
   - Would deleting each comment lose something the code does not already show? A comment restating the next line goes, one naming what a variable holds becomes a rename, and debugging history, failing inputs and judge anecdotes go in the commit or the stress test. Answer it from the comment table in the report, never from memory.
   - Do blank lines split each function body into its steps?
   A no that stays a no is allowed only for a higher priority (correct, fast, easy to call, easy to read inside), with a one-line reason in the report, such as a measured speed cost.
   The questions mirror CONTRIBUTING.md `## API and code` and `## Style`, so a change there updates them in the same commit.
4. **Stress test**, at `stress_tests/<same path>` (`stress_tests/python/` for Python), mandatory:
   - an independent reference: a brute force, and also the original source when it compiles (port mode)
   - exhaustive on small inputs, random on large ones, every boundary the header claims (max n, max value, extreme moduli)
   - randomness through `stress::rand_int` and loop counts through `stress::scaled` (`rng` and `scaled` in Python)
   - refactor: keep the existing test and extend it for every new behaviour
   - no practical reference: stop and ask; `STRESS_SKIP` needs the user's approval and a reason
   - a test that expects an abort must check that the library's own assert fired (its stderr names the library file and `Assertion`): CI's `-D_GLIBCXX_ASSERTIONS` also aborts on an out-of-range access, so a bare SIGABRT check passes even with the assert deleted
5. **Prove the stress test can fail.** Plant one realistic bug (an off-by-one, a dropped edge case, a swapped operand) for every logic change, and at least one even when no logic changed:
   ```bash
   bash .claude/skills/crucible/scripts/mutant.sh <path under code_library> '<sed expression>'
   ```
   It must print `mutant caught by: ...`. A mutant that survives means the test is too weak: strengthen it and plant again. `INVALID MUTANT` does not count.
6. **Performance.**
   - refactor: benchmark old against new at -O2 on the same input. More than 5% slower is reverted, unless the user asked for the trade. A kept slower choice gets a `///` comment with the numbers. A reverted attempt is reported with its numbers.
   - port: compare with the source when it compiles; a port slower than its source reports both timings, and the header notes state the trade.
7. **Judge test**, when a Library Checker or Aizu problem fits: `judge_tests/<topic>/<problem>.cpp` per CONTRIBUTING, plus the ✔ on the README index entry; one ✔ linking any one of the file's judge tests is enough. refactor: run every judge test that includes a file whose code changed. The first run prepares the judge data and can take minutes.
8. **Docs.** Give `/lucid` the state after the code change as its base, so its scope check covers only the docs: the commit of the code change with `--commit`, otherwise a snapshot from `git stash create`, which records the working tree without changing any file. Run `/lucid` on the file (header, `main()`, README index entry and practice section heading). A new template also gets its README index entry, alphabetical, variants nested under their base. Practice problems are optional; any you add must be verified per CONTRIBUTING.
9. **Verify**, from the repo root:
   ```bash
   env CHANGED_SINCE=<commit before the run> JOBS=4 bash .github/scripts/run_tests.sh self
   env CHANGED_SINCE=<commit before the run> JOBS=4 bash .github/scripts/run_tests.sh stress
   env LD_LIBRARY_PATH="$(dirname "$(g++ -print-file-name=libstdc++.so)")" bash .github/scripts/judge_tests.sh <the judge tests from step 7>
   ```
   When `competitive-verifier` is not on the PATH, add `COMPETITIVE_VERIFIER=<path to the executable>` just before `bash`. When it is not installed at all, stop before committing and ask the user to install it (`pip install competitive-verifier==4.2.2`): an unrun judge test is not a pass, and indirect checks such as comparing disassembly do not replace it.
   Done only when both suites exit 0 with the file's tests among those passed (`passed: 0` tested nothing and does not count) and every judge test succeeds.
10. **Commit**, only with `--commit`: one commit per library file together with its stress test, judge test and README lines. Conventional subject (`feat:` new template, `refactor:`, `perf:`, `fix:`), at most 72 characters, no body. No AI attribution, never `--no-verify`, never amend, never push.

## Edge cases to check

The bug classes the October 2026 audit found in this repo:
- empty input, `n <= 1`, duplicates, collinear or degenerate geometry
- a value of 0 colliding with a sentinel; negative inputs and `INT_MIN`
- `mod == 1`, inputs `>= mod` or negative not reduced
- overflow: `a * b % m` with 64-bit operands needs `(unsigned __int128)a * b % m` and a "Requires __int128" note; `int` parameters truncating `long long`
- a template type deduced from a literal, like `f(n, 0, cost)` making `T = int` while the cost is `long long`
- precision: `acosl(-1.0L)` in long double code, relative eps
- recursion as deep as `n` must fit the 8 MB stack CI uses: the stress test's largest `n` passes under it

## Output

The report MUST contain, in this order:

1. **Mode and source:** the file or URL; for web mode the licence and whether the code was adapted or reimplemented.
2. **API:** the common case as written in step 3, the simplicity gate answers (any no with its reason), and for refactor the API before and after with every caller updated.
3. **Comments:** a table of every `///` line outside the header and `main()`, each with `why:` and the information it adds that the code does not show, marked kept, added or removed. A comment with no such information is removed, not listed as kept.
4. **Header claims:** each range and complexity claim with the stress case or self-test that checks it.
5. **Stress test:** the reference used and what it covers.
6. **Planted bugs:** each sed expression with its `mutant caught by:` line, and any mutant that first survived with how the test was strengthened.
7. **Performance:** the benchmark numbers, including reverted attempts, or "not measured" with the reason (port with no compilable source).
8. **Judge tests:** each one with its result, or "none fits".
9. **Docs:** sections 1 to 4 of the `/lucid` report (README entry, claims, main(), proposed); its verification and commit are items 10 and 11 here.
10. **Verification:** the last 3 lines of both suite runs.
11. **Commits:** sha and subject for each, or "not committed".
