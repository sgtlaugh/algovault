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
- **Python templates** go under `code_library/python/`, with tests under `stress_tests/python/` that start with `from stress import rng, scaled`
- **An 8 MB stack.** CI runs tests with `ulimit -s 8192`, as many judges do, so deep recursion must fit or be iterative
- **Exceptions are explicit.** A file that truly cannot be stress tested in CI goes in `STRESS_SKIP` in `.github/scripts/run_tests.sh`, with the reason
- **An entry in the README index.** In its section, alphabetically
- **Comments explain why.** Constraints, tradeoffs and pitfalls, never what the code already says
- **A struct when there is state**, such as a graph, tree or table; a namespace only for stateless functions that share helpers, like `fft` or `ntt`
- **Four spaces, no tabs**

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
