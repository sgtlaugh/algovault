# Contributing

Bug reports, fixes and new templates are welcome as issues or pull requests.

## The Zen Of Contributing

Inspired by [The Zen Of Python](https://www.python.org/dev/peps/pep-0020/#id2)

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

## Where to start

| Adding | Read | Skill |
|---|---|---|
| A template, or reworking one | [Templates](#templates) to [Judge tests](#judge-tests), [API and code](#api-and-code) | `/crucible` |
| Documentation for a template | [Documentation](#documentation) | `/lucid` |
| A practice problem or solution | [Practice problems](#practice-problems) | `/verdict` |

The skills live in `.claude/skills/` and apply these rules step by step.

## Templates

Each template is one C++17 file under `code_library/<topic>/` that compiles on its own with GCC. It comes with:

- **Header and `main()`:** as in [Documentation](#documentation)
- **Stress test:** as in [Stress tests](#stress-tests)
- **Index entry:** in its README index section, alphabetical, variants nested under their base
- **Judge test:** when one fits, as in [Judge tests](#judge-tests)

Also:

- **Python:** templates go in `code_library/python/`, their tests in `stress_tests/python/` using `stress.py`
- **Recursion:** must fit an 8 MB stack, as on most judges; CI runs with that limit
- **Shared code:** a template built on another carries a copy of it, as in [Shared code](#shared-code)

### Shared code

Every template must compile on its own, so it can be pasted into a solution as one file. A template built on another, like offline dynamic connectivity on the rollback DSU, therefore cannot `#include` it: it carries a copy instead.

1. In the original, wrap the code others need in `/// BEGIN SHARED <name>` and `/// END SHARED <name>`
2. In the template that needs it, paste that code between `/// BEGIN COPY <name> from <path>` and `/// END COPY <name>`
3. CI fails when a copy differs from its original, so a fix to shared code goes into every copy:
   `git grep "BEGIN COPY <name>"` lists them

```cpp
// rollback_dsu.cpp, the original
/// BEGIN SHARED rollback_dsu
struct RollbackDSU{ ... };
/// END SHARED rollback_dsu

// offline_dynamic_connectivity.cpp, the copy
/// BEGIN COPY rollback_dsu from code_library/data_structures/rollback_dsu.cpp
struct RollbackDSU{ ... };
/// END COPY rollback_dsu
```

## Stress tests

Every template has a stress test at the same path under `stress_tests/`. It checks the template against an independent brute force, a few seconds by default.

```cpp
#include "../common.h"

#define main library_main
#include "../../code_library/dp/convex_hull_trick.cpp"
#undef main
```

- **Randomness:** values from `stress::rand_int`, loop counts from `stress::scaled`
- **Exceptions:** a file CI truly cannot stress test goes in `STRESS_SKIP` in `run_tests.sh`, with the reason

## Judge tests

When a Library Checker or Aizu problem fits, add `judge_tests/<topic>/<problem>.cpp`, starting as below, and a ✔ on the template's index entry:

```cpp
// competitive-verifier: PROBLEM <url>
// competitive-verifier: TLE <the problem's time limit>
```

## Practice problems

Each template's section in `code_library/<topic>/README.md` lists problems it solves, with a ✅ linking an accepted solution:

```markdown
- [Codeforces 932F - Escape Through Leaf](https://codeforces.com/problemset/problem/932/F)
- [Codeforces 1083E - The Fair Nut and Rectangles](https://codeforces.com/problemset/problem/1083/E)
- [CSES 2084 - Monster Game I](https://cses.fi/problemset/task/2084) [✅](../../practice/data_structures/cses_2084.cpp "Accepted solution")
```

**Problems**

- **Format:** `[<Judge> <code> - <title>](<url>)`
- **Order:** case-insensitive, numbers by value: `Codeforces 932F` before `Codeforces 1083E`
- **Verified:** each one can be submitted, to a live judge or to official data with a checker, like the Code Jam archive
- **Placement:** under the template at the core of the solution, not under setup steps like compression or sorting
- **Mirrors:** a problem on several judges is listed once, preferring UVa, then SPOJ, then Codeforces
- **Library Checker:** never a practice problem; one the template solves becomes its judge test instead

**Solutions**

- **File:** `practice/<topic>/<judge>_<code>.cpp`, `.py` for Python templates, exactly as submitted
- **Name:** the entry's judge and code in lowercase, spaces as underscores: `codeforces_25d.cpp`, `codeforces_gym_100570d.cpp`
- **Several templates:** one file per template, `<judge>_<code>_<template>.cpp`, each linked from its own section
- **Link:** ` [✅](../../practice/<topic>/<file> "Accepted solution")` after the entry
- **Template copy:** taken from `code_library/`, never from another practice file
- **POJ, UVa and HDU:** these judges compile C++11, and so does CI; change the template copy only where C++11 needs it, as `/verdict` lists

## API and code

A template is used under contest pressure by someone who did not write it, so the simplest correct design wins.

**Priorities, in order:** correct, fast, easy to call, easy to read inside. A rule below gives way only to a higher priority, with a comment saying so.

### The API

- **Common case:** one or two lines with no configuration: `DSU dsu(n)`, `SegmentTree<long long> st(n)`, `DisjointSparseTable<int> sum(a)`
- **Parameters:** ask for nothing the code can do without; `SegmentTree` takes no identity, since its folds start from the data
- **Names:** say what a call does and read as a sentence: `dsu.relate(a, b, d)` records `x[b] = x[a] + d`
- **Vocabulary:** one family, one set of names: `find`, `same`, `size` in every DSU, `combine` in every fold functor
- **Customizing:** by passing a lambda, a functor or an `Ops` struct, never by editing the template
- **Results:** returned, never written to an out-parameter: `mo()` returns its answers, `exact_cover()` an `optional`
- **Indexing:** 0-based positions, inclusive ranges `[l, r]`
- **Visibility:** everything public, no `private:` and no `class`, because contestants read and extend a template's state; the header lists the API
- **Changes:** an API change breaks old calls at compile time; a call that still compiles never changes meaning

### The code

- **Plain:** plain loops, and variables named for what they hold: `v, va, vb, vab`, not `t0..t7`
- **Steps:** blank lines split a function body into its steps, like setup, the main loop and the answer
- **No machinery:** no abstraction used once, and no template template parameters, SFINAE, CRTP or macros, unless they turn a silent wrong answer into a compile error, with a comment saying so
- **Speed tricks:** stay only with a measured reason in a comment
- **Structs and namespaces:** a struct when there is state, a namespace only for stateless helpers like `fft`
- **Options:** a set of named options is a nested enum, `enum Cmp{LESSEQ, EQUAL, GREATEQ};`, so a wrong option or a bare number fails to compile
- **Fields:** related ones share a line, `int n, m;` or `vector<int> parent, depth;`, unrelated ones do not
- **Short names:** `i, j, u, v, cur, res` suit locals, and loop variables live in their loop; fields and API names say what they hold

**Struct order**, which CI checks:

1. Types and constants
2. Fields
3. Constructors
4. The API, in header order
5. Helpers, in the order they are called

Free functions follow header order too, except that a helper comes before the function that calls it.

## Documentation

Every template reads the same way, so a reader knows what it does and how to call it within a minute.

### Header

```cpp
/***
 *
 * Partially Persistent DSU
 * Union-find that remembers its past: were a and b in the same set after the first t connect calls?
 *
 * Complexity: O(log n) per operation
 *
 * PartiallyPersistentDSU dsu(n): nodes 0..n, each in its own set
 * dsu.connect(a, b): merge the sets of a and b, false if they already were one set
 * dsu.find(a, t): representative of a's set after the first t connect calls
 * dsu.same(a, b, t): true if a and b were in the same set after the first t connect calls
 * dsu.size(a, t): number of nodes in a's set after the first t connect calls
 * dsu.now: number of connect calls so far, failed ones included
 *
 * Notes:
 *   - t = 0 is the state before any call, any t >= now is the current state
 *   - The first call that joined a and b is the smallest t with same(a, b, t), found by binary search
 *
***/
```

**Blocks**, separated by blank lines:

1. **Name**, then one plain sentence on what the template is for
2. **Complexity**, time first: one line for a single cost, otherwise one bullet per operation, or per struct when a file has several. Memory only when it is not the obvious O(n) or its size matters, like `about 205 MB`
3. **API**, one line per public call: `call(args): what it returns or does, preconditions`
4. **Notes**, one bullet each: value ranges, traps, and which template to use instead

**API order**, as in the example:

1. Constructors, simplest first: `PartiallyPersistentDSU dsu(n)`
2. Calls that change the structure, alphabetical: `connect`
3. Calls that only read it, alphabetical: `find`, `same`, `size`
4. Public fields, in the order a reader needs them: `now`

A call goes by its main effect: `connect` returns a bool but changes the structure, and a DSU's `find` compresses paths but only reads.

**API lines**

- **What, not how:** say what the call returns or does, and nothing the signature or an assert already shows
- **Self-contained:** define each call on its own: `a becomes the f with divisor_zeta(f) = a`, not `undoes divisor_zeta`
- **Math names:** name the object as the math does, so calls read as the formula: `Lagrange f(terms, mod)`, then `f(x)`

**Every line**

- **One sentence:** never split over two lines, never two calls on one line
- **120 characters at most:** a long line is cut, never wrapped
- **New facts only:** delete a line that repeats another, merge two that only make sense together
- **The contract:** what a caller passes and gets back; how the code works goes in a comment at that code
- **Plain words:** no filler, no marketing, no lists of applications; define a term outside standard contest vocabulary where it first appears, and write math as one sentence, never pseudocode
- **No usage example:** `main()` is the example; a short snippet only when calling the template needs code, like a functor

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

- **First block:** the worked example of the main API, under a `///` line saying what it shows
- **Asserts:** one call each, never `assert(a && b)` over several calls
- **Objects:** `T x(args)`, never `auto x = T(args)`
- **Blank lines:** between blocks, and before each new object inside a block
- **Comments:** only where an expected value is not obvious
- **Expected values:** worked out by hand, never copied from the program's output
- **Size:** a few blocks; edge cases and scale belong in the stress test, and nothing is timed or printed

### README

- **Index entry:** `[Name](path) - description`, the description one phrase on what the template is for, 50 characters at most
- **Details:** live in the header, so no feature lists and no implementation names; a complexity only when it sets the template apart
- **Variants:** nested under their base, saying only how they differ, like `every version stays queryable`
- **Practice headings:** match the index names

## Style

**Comments**

- **Why, never what:** delete a comment whose deletion loses nothing; replace one naming what a variable holds with a better name
- **Invariants:** comment a rule a caller or editor could break unknowingly, like a field that must stay false for zero
- **Commit message instead:** derivations, proofs, debugging history, measurements, failing inputs and judge anecdotes
- **One line:** at most 120 characters, at the code it explains, never two `///` lines in a row
- **Tests too:** a test comment names the case it covers and the bug it would catch

**Write this, not that**

| Write | Not |
|---|---|
| `f(x){` and `struct X{` | `f(x) {` and `struct X {` |
| `X() {}` | `X(){}` |
| `for (`, `if (`, `while (` | `for(`, `if(`, `while(` |
| `vector<int>` | `vector <int>` |
| `__lg(x)` | `31 - __builtin_clz(x)` |
| a `while` loop | a comma chain in a `for` header |
| `beyond double precision` | `may get either answer` |

## Complexity notation

- **Format:** inline code, `^` for powers, `sqrt(x)` for roots: `O(n log n)`, `O(n^(2/3))`, `O(m sqrt(n))`
- **Products:** a space, not `*`: `O(n m)`
- **Variables:** lowercase `n`, `m` for input sizes, any other letter defined beside it: `O(n C) for capacity C`

## Checks

```bash
bash .github/scripts/run_tests.sh self      # every self-test
bash .github/scripts/run_tests.sh stress    # every stress test
bash .github/scripts/judge_tests.sh         # every judge test, or name files; needs pip install competitive-verifier==4.2.2
```

Both `run_tests.sh` modes build with `-std=c++17 -Wall -Wextra -Werror` under AddressSanitizer and UndefinedBehaviorSanitizer. CI fails on:

- **Tests:** a failing test, a compiler warning, a library file without a stress test, a `fast_io.cpp` that does not build as C++20
- **Index:** a library file missing from the index, a description over 50 characters, a ✔ that disagrees with `judge_tests/`
- **Practice:** a ✅ that disagrees with `practice/`, copies of a problem that differ in title or URL, a POJ, UVa or HDU solution that does not build as C++11
- **Headers:** no complexity, a `*` in an `O(...)` product, a line over 120 characters
- **Code:** a copied region that drifted, a field or type after a method, API calls out of header order, two `///` lines in a row, a whitespace or style error

| Variable | Effect |
|---|---|
| `JOBS=n` | runs n tests at a time, all cores by default |
| `CHANGED_SINCE=<commit>` | tests only files changed since that commit, uncommitted edits included, as pushes do |
| `BUILD=fast` | builds with `-O2` and no sanitizers; with `STRESS_SCALE=20` it is the nightly run |
| `STRESS_SEED`, `STRESS_SCALE` | reproduce a stress failure from the seed each stress test prints |

Judge tests compile with `-O2` and run on the judge's own data. New judge test files must be `git add`ed first.
