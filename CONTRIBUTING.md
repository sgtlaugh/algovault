# Contributing

Bug reports, fixes and new templates are welcome as issues or pull requests.

## Adding or changing a template

- **One standalone file.** A C++17 file under `code_library/<topic>/` that compiles on its own with GCC, so it can be copied into a solution as is
- **A header that documents it.** What it computes, its complexity, its API and its limits, so the file works as a black box
- **Self-tests in `main()`.** Assert-based checks of known answers
- **A brute-force stress test.** Under `stress_tests/` at the same relative path, comparing against an independent reference on random inputs
  - Include `../common.h` and the library with `#define main library_main`
  - Draw randomness through `stress::rand_int` and scale loops with `stress::scaled`, so `STRESS_SEED` and `STRESS_SCALE` reproduce and lengthen runs
  - Keep the default run to a few seconds under sanitizers
- **An entry in the README index.** In its section, alphabetically
- **Comments explain why.** Constraints, tradeoffs and pitfalls, never what the code already says

## Running the checks

CI fails on any test failure, any compiler warning, a library file without a stress test, or one missing from the README index. Run the same checks locally:

```bash
bash .github/scripts/run_tests.sh self      # every self-test
bash .github/scripts/run_tests.sh stress    # every stress test
```

Both compile with `-std=c++17` under AddressSanitizer and UndefinedBehaviorSanitizer with `-Wall -Wextra -Werror`.

## Complexity notation

In headers and the README index, write complexities as inline code: `O(n log n)`, `O(n^(2/3))`, `O(m sqrt(n))`, `O(2^n n^2)`.
Use a space for products instead of `*`, `^` for powers, `sqrt(x)` for roots, and lowercase `n`, `m` for input sizes, naming any other variable such as `k` terminals or `W` capacity.
In the index, give a complexity only where it helps choose a template.

## The Zen Of Contributing
Inspired from [The Zen Of Python](https://www.python.org/dev/peps/pep-0020/#id2)

```python
Beautiful is better than ugly

Simple is way better than complex

Consistency matters

Flat is preferred over nested

Typing is better than incomprehensible macros

Because readability counts

But not as much as speed

Some documentation is better than no documentation

No documentation is better than extensive documentation

But not as important as ease of reusing as a black box

Four spaces are better than tabs

Tabs are better than no spaces

If the implementation is hard to explain, it's a bad idea

If the implementation is easy to explain, it may be a good idea

Structs are one honking great idea, let's do more of those!
```
