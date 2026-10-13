---
name: crucible
description: >
  Add or rework an algovault template: port an algorithm from a file of the user's own code, reimplement one from a
  web page, or refactor an existing code_library file to be simpler, easier to use or faster. Every result ships with a
  brute-force stress test proven to catch a planted bug, docs from /lucid, and a green CI suite.
trigger: /crucible
allowed-tools: Read, Write, Edit, Glob, Grep, WebFetch, Skill, Bash(git grep:*), Bash(git diff:*), Bash(git log:*), Bash(git status:*), Bash(git rev-parse:*), Bash(git add:*), Bash(git commit:*), Bash(ls:*), Bash(g++:*), Bash(python3:*), Bash(env CHANGED_SINCE=:*), Bash(bash .github/scripts/run_tests.sh:*), Bash(bash .github/scripts/judge_tests.sh:*), Bash(bash .claude/skills/crucible/scripts/mutant.sh:*)
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

## Steps

1. **Understand.**
   - Read the source in full, two recent `code_library/` files of the same topic for style, and `.github/scripts/run_tests.sh` for the CI flags.
   - port and web: `ls -R code_library/` and `git grep -il <algorithm name>`; if the algorithm already exists, stop and report it, offering refactor mode instead.
   - web: find the page's licence. Code under a permissive licence (CC0, MIT, Apache, public domain) may be adapted, credited in the header. Code under CC BY-SA, GPL or no stated licence is never copied: reimplement it from the idea and cite the page. A run without a licence verdict in the report is not done.
   - refactor: read the file's stress test, judge tests and callers (`git grep -n <struct or function> -- code_library stress_tests judge_tests`), and benchmark the current version at -O2 on its largest documented input before changing anything.
   - Ambiguous or suspicious source (an off-by-one in a known algorithm, a complexity claim the code does not meet): stop and ask. Porting a suspected bug silently is not allowed.
2. **Write the code** to the CONTRIBUTING style: a struct with per-instance state, 0-based positions with inclusive ranges `[l, r]`, no mutable globals and no generic global names, `T x(args)`. For each item of the edge-case list below that applies, the code handles it or the header states the limit.
   - refactor: an API change must break old calls at compile time, never change their meaning silently; update every caller. Keep the library's simple patterns (a functor template with a default, an Ops struct).
3. **Stress test**, at `stress_tests/<same path>` (`stress_tests/python/` for Python), mandatory:
   - an independent reference: a brute force, and also the original source when it compiles (port mode)
   - exhaustive on small inputs, random on large ones, every boundary the header claims (max n, max value, extreme moduli)
   - randomness through `stress::rand_int` and loop counts through `stress::scaled` (`rng` and `scaled` in Python)
   - refactor: keep the existing test and extend it for every new behaviour
   - no practical reference: stop and ask; `STRESS_SKIP` needs the user's approval and a reason
4. **Prove the stress test can fail.** For every logic change, plant one realistic bug (an off-by-one, a dropped edge case, a swapped operand):
   ```bash
   bash .claude/skills/crucible/scripts/mutant.sh <path under code_library> '<sed expression>'
   ```
   It must print `mutant caught`. A mutant that survives means the test is too weak: strengthen it and plant again. `INVALID MUTANT` does not count.
5. **Performance.**
   - refactor: benchmark old against new at -O2 on the same input. More than 5% slower is reverted, unless the user asked for the trade. A kept slower choice gets a `///` comment with the numbers.
   - port: compare with the source when it compiles; a slower port needs a measured reason.
6. **Judge test**, when a Library Checker or Aizu problem fits: `judge_tests/<topic>/<problem>.cpp` per CONTRIBUTING, plus the ✔ on the README index entry. refactor: run every judge test that includes a file whose API changed.
7. **Docs.** Run `/lucid` on the file (header, `main()`, README index entry and practice section heading). A new template also gets its README index entry, alphabetical, variants nested under their base. Practice problems are optional; any you add must be verified per CONTRIBUTING.
8. **Verify**, from the repo root:
   ```bash
   env CHANGED_SINCE=<commit before the run> JOBS=4 bash .github/scripts/run_tests.sh self
   env CHANGED_SINCE=<commit before the run> JOBS=4 bash .github/scripts/run_tests.sh stress
   bash .github/scripts/judge_tests.sh <the judge tests from step 6>
   ```
   Done only when both suites exit 0 with the file's tests among those passed (`passed: 0` tested nothing and does not count) and every judge test succeeds.
9. **Commit**, only with `--commit`: one commit per library file together with its stress test, judge test and README lines. Conventional subject (`feat:` new template, `refactor:`, `perf:`, `fix:`), at most 72 characters, no body. No AI attribution, never `--no-verify`, never amend, never push.

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
2. **API:** for a new template, its public calls; for refactor, before and after with every caller updated.
3. **Header claims:** each range and complexity claim with the stress case or self-test that checks it.
4. **Stress test:** the reference used and what it covers.
5. **Planted bugs:** each sed expression with the `mutant caught` line.
6. **Performance:** the benchmark numbers, or "not measured" with the reason (port with no compilable source).
7. **Judge tests:** each one with its result, or "none fits".
8. **Docs:** the `/lucid` report.
9. **Verification:** the last 3 lines of both suite runs.
10. **Commits:** sha and subject for each, or "not committed".
