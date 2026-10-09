# Contributing

Bug reports, fixes and new templates are welcome as issues or pull requests.

## Adding a template

- One standalone C++17 file under `code_library/<topic>/` that compiles on its own with GCC
- A header stating what it computes, its complexity, its API and its limits
- Assert-based self-tests in `main()`
- A brute-force stress test at the same path under `stress_tests/`, against an independent reference
- The test includes `../common.h` and the library with `#define main library_main`
- Randomness through `stress::rand_int`, loops sized by `stress::scaled`, a few seconds by default
- An entry in its README index section, alphabetically, variants nested under their base
- Python templates go in `code_library/python/`, tests in `stress_tests/python/` via `stress.py`
- A file CI truly cannot stress test goes in `STRESS_SKIP` in `run_tests.sh`, with the reason

## Style

- Comments explain why, never what
- `){` and `struct X{` open a body, `) {}` is empty, `for (` takes a space, `vector<int>` does not
- A struct when there is state, a namespace only for stateless helpers like `fft`
- Recursion must fit an 8 MB stack, which CI enforces as most judges do

## Checks

CI fails on a failing test, a compiler warning, a missing stress test or index entry, or a whitespace or style error.

```bash
bash .github/scripts/run_tests.sh self      # every self-test
bash .github/scripts/run_tests.sh stress    # every stress test
```

Both build with `-std=c++17 -Wall -Wextra -Werror` under AddressSanitizer and UndefinedBehaviorSanitizer.

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
