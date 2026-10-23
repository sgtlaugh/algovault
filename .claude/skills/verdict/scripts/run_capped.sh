#!/bin/bash
# run_capped.sh <command> [args...]
# Runs a solution, brute force or generator with 4 GiB of address space and 120 s of CPU, so a runaway allocation dies
# alone with bad_alloc instead of pushing the whole desktop into swap. Not for sanitizer builds, which reserve far more
# virtual memory than they use.
exec prlimit --as=$((4 << 30)) --cpu=120 -- "$@"
