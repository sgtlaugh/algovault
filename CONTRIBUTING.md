# Contributing

Bug reports, fixes and new templates are welcome as issues or pull requests.

## Adding a template

- One standalone C++17 file under `code_library/<topic>/` that compiles on its own with GCC
- A header and a usage example in `main()` as described in [Documentation](#documentation), with edge cases in the stress test
- A brute-force stress test at the same path under `stress_tests/`, against an independent reference
- The test includes `../common.h` and the library with `#define main library_main`
- Randomness through `stress::rand_int`, loops sized by `stress::scaled`, a few seconds by default
- An entry in its README index section, alphabetically, variants nested under their base
- Verified practice problems in `code_library/<topic>/README.md`, as `[<Judge> <code> - <title>](<url>)` in natural order,
  each on a live judge or with official data and a judge, like the Code Jam archive
- An accepted solution built on the template goes in `practice/<topic>/<judge>_<code>.cpp` as submitted, linked from its
  practice problem as ` [✅](../../practice/<topic>/<file> "Accepted solution")`; templates are copied from `code_library/`, never from `practice/`
- A problem solved on several templates gets one file per template, `<judge>_<code>_<template>.cpp`, each linked from its section
- Python templates go in `code_library/python/`, tests in `stress_tests/python/` via `stress.py`
- A file CI truly cannot stress test goes in `STRESS_SKIP` in `run_tests.sh`, with the reason
- A judge test under `judge_tests/<topic>/<problem>.cpp` when a Library Checker or Aizu problem fits, starting with
  `// competitive-verifier: PROBLEM <url>` and `// competitive-verifier: TLE <the problem's time limit>`, plus a ✔ on the index entry

## Documentation

Every template reads the same way, so a reader knows what it does and how to call it within a minute.
The `/lucid` skill in `.claude/skills/lucid/` applies these rules to existing files.

### Header

```cpp
/***
 *
 * Disjoint Sparse Table
 * O(1) range queries on a static array for any associative operation
 *
 * Complexity:
 *   - build: O(n log n) time and memory
 *   - query: O(1)
 *
 * DisjointSparseTable<T> table(a): range sums of a
 * DisjointSparseTable table(a, op): range folds with op, which must be associative but need not commute
 * table.query(l, r): op over a[l..r], 0 <= l <= r < n
 *
 * Use sparse_table.cpp for min, max or gcd, and a prefix array for plain sums
 *
***/
```

- The name on the first line, one plain sentence on what it is for on the second
- Blank lines between the blocks: name, complexity, API, notes
- Complexity on one line when there is one cost, `Complexity: O(α(n)) amortized per operation`
- Otherwise one bullet per operation, or per struct when a file has several, time first
- Memory only when it is not the obvious O(n) or its real size matters, like `about 205 MB by default`
- The API lists the constructor first, then one line per public call: `call(args): what it returns or does, preconditions`
- One sentence per line, a sentence never split across lines, never two calls on one line
- No usage example, `main()` is the example; a short snippet only when calling it needs code, like a functor
- Notes on value ranges, traps and which template to use instead, one short line each
- Plain words: no filler, no marketing, no lists of applications

### main()

```cpp
int main(){
    /// Sum and product over the same array
    vector<int> a = {5, 6, 1, 13, 7};
    DisjointSparseTable<int> sum(a);
    assert(sum.query(1, 3) == 20);
    assert(sum.query(4, 4) == 7);

    DisjointSparseTable prod(a, [](int x, int y){ return x * y % 1000; });
    assert(prod.query(0, 3) == 390);

    return 0;
}
```

- The first block is the worked example of the main API, under a `///` line saying what it shows
- One call per assert, never `assert(a && b)` over several calls
- `T x(args)`, not `auto x = T(args)`
- A blank line between blocks, and before each new object inside a block
- A `///` comment only where the expected value is not obvious
- A few blocks only: edge cases and scale belong in the stress test, no timing and no printing
- Expected values worked out by hand, never copied from the program's output

### README

- An index entry is `[Name](path) - short description`, matching the header's second line, variants nested under their base
- Practice section headings match the index names

## Style

- Comments explain why, never what
- `){` and `struct X{` open a body, `) {}` is empty, `for (` takes a space, `vector<int>` does not
- A struct when there is state, a namespace only for stateless helpers like `fft`
- Recursion must fit an 8 MB stack, which CI enforces as most judges do
- A template that needs another one carries a verbatim copy between `/// BEGIN COPY <name> from <file>` and
  `/// END COPY <name>`, matching `/// BEGIN SHARED <name>` ... `/// END SHARED <name>` in the original, which CI checks

## Checks

CI fails on a failing test, a compiler warning, a missing stress test or index entry, or a whitespace or style error.

```bash
bash .github/scripts/run_tests.sh self      # every self-test
bash .github/scripts/run_tests.sh stress    # every stress test
bash .github/scripts/judge_tests.sh         # every judge test, or name files; needs pip install competitive-verifier==4.2.2
```

`run_tests.sh` builds with `-std=c++17 -Wall -Wextra -Werror` under AddressSanitizer and UndefinedBehaviorSanitizer and takes:

- `JOBS=n` tests at a time, all cores by default
- `CHANGED_SINCE=<commit>` only the tests of files changed since that commit, which is what pushes run
- `BUILD=fast` for `-O2` without sanitizers, with `STRESS_SCALE=20` the nightly high-iteration run
- `STRESS_SEED` and `STRESS_SCALE` to reproduce a stress failure, every stress test prints its seed when it starts

Judge tests compile with `-O2` and run on the judge's own data; new judge test files must be `git add`ed first.

## Complexity notation

- Inline code: `O(n log n)`, `O(n^(2/3))`, `O(m sqrt(n))`
- A space for products, not `*`
- `^` for powers, `sqrt(x)` for roots
- Lowercase `n`, `m` for input sizes, naming any other variable
- In the README index, only where it helps choose a template

## The Zen Of Contributing
Inspired from [The Zen Of Python](https://www.python.org/dev/peps/pep-0020/#id2)

```python
Beautiful is better than ugly

Simple is way better than complex

Consistency matters

Flat is preferred over nested

Typing is better than incomprehensible macros

Correct is better than fast

Because readability counts

But not as much as speed

Untested is broken, and a brute force is the best reviewer

Warnings are errors

Measure, don't guess

Document the contract, not the code

A black box needs its limits written down

Comments say why, the code already says what

If the implementation is hard to explain, it's a bad idea

If the implementation is easy to explain, it may be a good idea

Structs are one honking great idea, let's do more of those!
```
