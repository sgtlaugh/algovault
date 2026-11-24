"""check_order.py: prints entries that break natural order, practice entries within a section and README index siblings

An index entry is compared with the previous entry at its own indent in the same list, so variants sort under their base.
"""
import glob
import re

from practice_toc import natural


def misordered(path, item):
    """Each line matching item that sorts before the previous sibling at its indent in the same block."""
    previous = {}
    for number, line in enumerate(open(path, encoding='utf-8'), 1):
        m = re.match(item, line)
        if not m:
            if line.startswith('#') or not line.strip():
                previous = {}
            continue
        indent = len(m.group(1))
        previous = {k: v for k, v in previous.items() if k <= indent}
        if indent in previous and natural(line) < natural(previous[indent]):
            print(f'{path}:{number}: {line.strip()[:80]}')
        previous[indent] = line


for path in sorted(glob.glob('code_library/*/README.md')):
    misordered(path, r'()- \[.+? - ')
misordered('README.md', r'( *)- \[[^\]]+\]\(code_library/')
