# Contributing

Bug reports, fixes and new templates are welcome as issues or pull requests.

## Adding a template

The `/crucible` skill in `.claude/skills/crucible/` walks through all of this, for a new template or a refactor of an existing one.

- One standalone C++17 file under `code_library/<topic>/` that compiles on its own with GCC
- A header and a usage example in `main()` as described in [Documentation](#documentation), with edge cases in the stress test
- A brute-force stress test at the same path under `stress_tests/`, against an independent reference
- The test includes `../common.h` and the library with `#define main library_main`
- Randomness through `stress::rand_int`, loops sized by `stress::scaled`, a few seconds by default
- An entry in its README index section, alphabetically, variants nested under their base
- Verified practice problems in `code_library/<topic>/README.md`, as `[<Judge> <code> - <title>](<url>)` in natural order,
  each on a live judge or with official data and a judge, like the Code Jam archive
- A problem goes under the template at the core of its solution, not under setup steps like compression or sorting
- An accepted solution built on the template goes in `practice/<topic>/<judge>_<code>.cpp` (`.py` for Python templates) as submitted, linked from its
  practice problem as ` [✅](../../practice/<topic>/<file> "Accepted solution")`; templates are copied from `code_library/`, never from `practice/`
- POJ, UVa and HDU solutions compile as C++11, which CI checks: their template copy changes only where C++11 needs it, as the
  `/verdict` skill lists
- A problem solved on several templates gets one file per template, `<judge>_<code>_<template>.cpp`, each linked from its section
- `<judge>` and `<code>` are the practice entry's in lowercase with underscores for spaces, `codeforces_25d.cpp`, `codeforces_gym_100570d.cpp`
- Python templates go in `code_library/python/`, tests in `stress_tests/python/` via `stress.py`
- A file CI truly cannot stress test goes in `STRESS_SKIP` in `run_tests.sh`, with the reason
- A judge test under `judge_tests/<topic>/<problem>.cpp` when a Library Checker or Aizu problem fits, starting with
  `// competitive-verifier: PROBLEM <url>` and `// competitive-verifier: TLE <the problem's time limit>`, plus a ✔ on the index entry

## API and code

A template is used under contest pressure by someone who did not write it, so the simplest correct design wins over the cleverest one.
The priorities, in order: correct, fast, easy to call, easy to read inside; a rule below gives way only to a higher priority, with a comment saying so.

### The API

- The common case is one or two lines with no configuration: `DSU dsu(n)`, `SegmentTree<long long> st(n)`, `DisjointSparseTable<int> sum(a)`
- Ask for nothing the code can do without; `DisjointSparseTable` and `SegmentTree` need no identity because their folds start from the data
- Names say what a call does and read as a sentence: `dsu.relate(a, b, d)` records `x[b] = x[a] + d`, `same(a, b)`, `size(a)`
- One family, one vocabulary: `find`, `same`, `size` in every DSU, `combine` for every fold functor
- Customize by passing a lambda, a functor or an `Ops` struct, never by editing the template
- Return the answer instead of filling an out-parameter: `mo()` returns its answers, `exact_cover()` an `optional`
- 0-based positions and inclusive ranges `[l, r]` everywhere
- Everything is public, no `private:` and no `class`: contestants read and extend a template's state, and the header says which calls are the API
- An API change breaks old calls at compile time; a call that still compiles never changes meaning

### The code

- Plain loops and named variables: `v, va, vb, vab` with a derivation comment, not `t0..t7` or `x, y, z, c, d`
- Standard terms over clever ones: `__lg` over `31 - __builtin_clz`, a `while` loop over a comma chain in a `for` header
- No abstraction for one use, and no language machinery (template template parameters, SFINAE traits, CRTP, macros) unless
  it turns a silent wrong answer into a compile error, with a comment saying so
- The one non-obvious line gets a one-line derivation; everything else should not need a comment
- Blank lines split a function body into its steps, like setup, the main loop and the answer
- A trick for speed stays only with a measured reason in a comment
- A struct reads top down: types and constants, every field, constructors, the API in header order, then helpers in the
  order they are called; CI fails a field or type declared after a method
- Related fields or constants share a line, `int n, m;` or `vector<int> parent, depth;`, unrelated ones do not
- A set of named options is a nested enum, `enum Cmp{LESSEQ, EQUAL, GREATEQ};`, not int constants, so passing an option
  from the wrong set or a bare number fails to compile
- Loop variables live in their loop; short conventional names (`i, j, u, v, cur, res`) suit locals, while fields and API
  names say what they hold

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
- The header states what a caller must pass and gets back; how the code works goes in a comment at that code
- Plain words: no filler, no marketing, no lists of applications
- A term outside standard contest vocabulary is defined in plain words where it first appears, and math is one sentence, never pseudocode

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

- An index entry is `[Name](path) - description`, the description one short phrase on what the template is for, at most 50 characters
- Details live in the header, not the index: no feature lists, no implementation names, a complexity only when it sets the template apart
- A variant is nested under its base and says only how it differs, like `every version stays queryable`
- Practice section headings match the index names

## Style

- Comments explain why, never what: a comment whose deletion loses nothing the code already shows is deleted, and one naming
  what a variable holds becomes a better name
- A comment is one line at the code it explains, with no debugging history, failing inputs or judge anecdotes: those go in
  the commit message or the stress test
- A header note states the condition instead of hedging: `beyond double precision`, not `may get either answer`
- `){` and `struct X{` open a body, `) {}` is empty, `for (` takes a space, `vector<int>` does not
- A struct when there is state, a namespace only for stateless helpers like `fft`
- Recursion must fit an 8 MB stack, which CI enforces as most judges do
- A template that needs another one carries a verbatim copy between `/// BEGIN COPY <name> from <file>` and
  `/// END COPY <name>`, matching `/// BEGIN SHARED <name>` ... `/// END SHARED <name>` in the original, which CI checks

## Checks

CI fails on a failing test, a compiler warning, a missing stress test or index entry, an index description over 50 characters, a ✔ that disagrees with `judge_tests/`, a ✅ that disagrees with `practice/`, a header without a complexity, a copied region that drifted, a struct member declared after a method, or a whitespace or style error.

```bash
bash .github/scripts/run_tests.sh self      # every self-test
bash .github/scripts/run_tests.sh stress    # every stress test
bash .github/scripts/judge_tests.sh         # every judge test, or name files; needs pip install competitive-verifier==4.2.2
```

`run_tests.sh` builds with `-std=c++17 -Wall -Wextra -Werror` under AddressSanitizer and UndefinedBehaviorSanitizer and takes:

- `JOBS=n` tests at a time, all cores by default
- `CHANGED_SINCE=<commit>` only the tests of files changed since that commit, uncommitted edits included, which is what pushes run
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
