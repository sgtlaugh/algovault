"""
Shared helpers for the Python stress tests, the Python counterpart of stress_tests/common.h

STRESS_SEED fixes the random seed, the default is fixed so CI runs are reproducible
STRESS_SCALE multiplies iteration counts for longer runs, the default is 1
Importing this module also makes code_library/python importable
"""

import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'code_library', 'python'))

SEED = int(os.environ.get('STRESS_SEED') or 20260105)
SCALE = max(1, int(os.environ.get('STRESS_SCALE') or 1))

rng = random.Random(SEED)
print(f'STRESS_SEED={SEED}', file=sys.stderr)


def scaled(iterations):
    return iterations * SCALE
