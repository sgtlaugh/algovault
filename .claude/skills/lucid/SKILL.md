---
name: lucid
description: >
  Rewrite the documentation of algovault templates to the CONTRIBUTING.md rules - the header comment,
  the usage example in main() and the README index entry - so a reader knows what a template does and
  how to call it within a minute. Docs only: never changes library behaviour or public names.
trigger: /lucid
allowed-tools: Read, Edit, Write, Glob, Grep, Bash(git grep:*), Bash(git diff:*), Bash(git log:*), Bash(git status:*), Bash(git rev-parse:*), Bash(git add:*), Bash(git commit:*), Bash(env CHANGED_SINCE=:*), Bash(bash .github/scripts/run_tests.sh:*)
argument-hint: "<library file or directory> [--commit], e.g. 'code_library/data_structures/hashmap.cpp --commit'"
---

# /lucid - clear headers and usage examples

The rules live in the `## Documentation` section of `CONTRIBUTING.md`. Read it first on every run; this skill is the workflow that applies it, not a second copy of the rules.

## Scope

- **Edits only:** the `/*** ... ***/` header (the module docstring for Python), the body of `main()`, and the file's README index entry and practice section heading.
- **Never:** library code outside the header and `main()`, stress tests, judge tests, public names. A rename that would help is reported under "Proposed, not applied", never made: renames break pasted code and callers.
- **Commits only with `--commit`.** Without it, leave the edits in the working tree and say "not committed".
- A directory argument means every library file in it, one file at a time, each finished (verified, and committed with `--commit`) before the next.

## Steps, per file

1. **Read** the whole file, its stress test (`stress_tests/<same path>`), its README index entry, and its callers: `git grep -n <struct or function name> -- code_library stress_tests judge_tests`. Record `base=$(git rev-parse HEAD)` before the first edit.
2. **Header.** Rewrite it to the CONTRIBUTING rules. Every statement in it is a claim to check against the code:
   - complexity from the loops and recursion, memory from the allocations
   - each API line from the actual signature, its return value and its asserts
   - value ranges from the types and asserts
   A claim you cannot confirm from the code is removed, not guessed.
3. **main().** Rewrite it to the CONTRIBUTING rules: the first block is the worked example, one call per assert, a blank line before each new object. Each further block shows one more use of the API, like a second constructor or a mode; a block that only checks an edge case or scale the stress test covers is removed. Every expected value is worked out by hand from the input; a value taken from running the program is forbidden.
4. **README.** Make the index entry agree with the header's second line: the same meaning, often shorter. If the template's display name changed, update its practice section heading in `code_library/<topic>/README.md` to match.
5. **Verify.** Run both, from the repo root:
   ```bash
   env CHANGED_SINCE=$base JOBS=4 bash .github/scripts/run_tests.sh self
   env CHANGED_SINCE=$base JOBS=4 bash .github/scripts/run_tests.sh stress
   ```
   Then check the scope: `git diff $base --stat` lists only the library file and READMEs, and every changed hunk of the library file lies inside the header or `main()`.
   A run that reports `passed: 0` tested nothing and does not count.
   The file is not done until both suites exit 0 with the file's tests among those passed, and the scope check holds. A failure is fixed in the header or `main()`; if the fix needs a library change, stop and report it.
6. **Commit** (only with `--commit`): one commit per library file, together with its README lines. Subject `docs: <template> header and usage follow CONTRIBUTING` or similar, conventional, at most 72 characters, no body. No AI attribution, never `--no-verify`, never amend, never push.

## Checklist before finishing a file

Every item is a yes, or the file is not done. The list mirrors CONTRIBUTING.md, so a change to its Documentation section updates this list in the same commit.

- [ ] Line 1 is the name, line 2 one plain sentence; blank lines separate name, complexity, API and notes
- [ ] Complexity is one line for one cost, otherwise one bullet per operation or struct, time first, memory only when it matters
- [ ] One API line per public call, constructor first, each one sentence: `call(args): what it returns or does, preconditions`
- [ ] No sentence is split across lines and no line holds two calls
- [ ] The header has no usage example; a snippet appears only where calling the template needs code
- [ ] `main()` has no chained `assert(a && b)` over several calls, no `auto x = T(...)`, no timing, no printing
- [ ] A blank line separates blocks and precedes each new object in `main()`
- [ ] Both suite runs exited 0, and the diff touches only the header, `main()` and README lines

## Output, per file

The report MUST contain, in this order:

1. **File** and its README entry, before and after.
2. **Claims:** a table of every header line, each with the code that confirms it (`file:line` or the signature). A line without evidence is not allowed in the header.
3. **main():** each kept block with the API use it shows, and each removed block with the stress test check that covers it.
4. **Proposed, not applied:** renames, missing derivation comments and library defects noticed, or "none".
5. **Verification:** the last 3 lines of both suite runs and the `git diff --stat` output.
6. **Commit:** the sha and subject, or "not committed".
