"""practice_audit.py <topic>: signals for spotting low-value practice problems in code_library/<topic>/README.md

Per problem, one row per README section:
  loc   - lines of the solution after the template copies, trivial lines like braces excluded
  sim   - highest token similarity (0 to 1) to another solution in the same section, and which one
  pop   - UVa distinct accepted users from uHunt, or Codeforces rating/solved count; high means easy
  flags - thin (loc <= 8), similar (sim >= 0.5), easy (UVa pop > 3000 or Codeforces rating <= 1400),
          judge test (the same URL is a judge test), and other sections listing the problem
The flags only point at candidates: a judgment of whether the problem teaches a new use of the template decides.
Popularity is fetched once and cached in /tmp/practice_audit_cache.json.
"""
import collections
import glob
import json
import os
import re
import sys
import urllib.request

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../..')
CACHE = '/tmp/practice_audit_cache.json'
TRIVIAL = {'', '{', '}', 'return 0;', 'using namespace fio;', 'int main(){'}
KEYWORDS = set('''int long double char bool void auto const static return if else for while do break continue struct
using namespace true false string vector pair size push_back read write main unsigned'''.split())
ENTRY = re.compile(r'^- \[(.+?) (\S+) - (.+?)\]\((.+?)\)(?: \[✔\]\(\.\./\.\./(practice/\S+?) )?')


def body(path):
    """The solution after its template copies, which end with the fast_io namespace."""
    text = open(path).read()
    m = re.search(r'^namespace fio\b.*?^}$', text, re.S | re.M)
    return text[m.end():] if m else text[text.rfind('int main(){'):]


def shingles(code):
    """Token 6-grams with identifiers and numbers erased, so renamed variables still match."""
    toks = re.findall(r'[A-Za-z_]\w*|\d+|\S', re.sub(r'"[^"]*"', '"s"', code))
    toks = ['N' if t.isdigit() else t if t in KEYWORDS or not t[0].isalpha() else 'ID' for t in toks]
    return {tuple(toks[i:i + 6]) for i in range(len(toks) - 5)}


def fetch(cache, url):
    if url not in cache:
        with urllib.request.urlopen(url, timeout=20) as r:
            cache[url] = json.loads(r.read())
        with open(CACHE, 'w') as f:
            json.dump(cache, f)
    return cache[url]


def popularity(cache, judge, code):
    """(text, easy) for UVa and Codeforces, None for judges without public statistics."""
    try:
        if judge == 'UVa':
            dacu = fetch(cache, f'https://uhunt.onlinejudge.org/api/p/num/{code}').get('dacu')
            return str(dacu), dacu is not None and dacu > 3000
        if judge == 'Codeforces':
            cf = fetch(cache, 'https://codeforces.com/api/problemset.problems')['result']
            m = re.match(r'(\d+)(\w+)', code)
            key = (int(m.group(1)), m.group(2))
            stats = cf['problemStatistics']
            rating = next((p.get('rating') for p in cf['problems'] if (p['contestId'], p['index']) == key), None)
            solved = next((s['solvedCount'] for s in stats if (s['contestId'], s['index']) == key), None)
            return f'{rating}r/{solved}', rating is not None and rating <= 1400
    except (OSError, ValueError, KeyError) as e:
        return f'error {type(e).__name__}', False
    return None


def main():
    topic = sys.argv[1]
    cache = json.load(open(CACHE)) if os.path.exists(CACHE) else {}
    judge_urls = set()
    for path in glob.glob(f'{ROOT}/judge_tests/**/*.cpp', recursive=True):
        judge_urls |= set(re.findall(r'PROBLEM (\S+)', open(path).read()))

    sections, current = collections.defaultdict(list), None
    for line in open(f'{ROOT}/code_library/{topic}/README.md').read().split('\n'):
        m = re.match(r'^#### \[(.+?)\]', line)
        if m:
            current = m.group(1)
        elif current and ENTRY.match(line):
            sections[current].append(ENTRY.match(line).groups())
    listed_in = collections.defaultdict(list)
    for name, entries in sections.items():
        for judge, code, _, _, _ in entries:
            listed_in[f'{judge} {code}'].append(name)

    for name, entries in sections.items():
        rows = []
        for judge, code, _, url, solution in entries:
            row = {'id': f'{judge} {code}', 'url': url, 'pop': popularity(cache, judge, code), 'solution': solution}
            if solution:
                code_text = body(f'{ROOT}/{solution}')
                row['loc'] = sum(1 for line in code_text.split('\n') if line.strip() not in TRIVIAL)
                row['shingles'] = shingles(code_text)
            rows.append(row)

        print(f'\n### {name} ({len(rows)} problems)')
        print(f"{'problem':34} {'loc':>4} {'sim':>5} {'pop':>12}  flags")
        for row in rows:
            flags = [f'also in {other}' for other in listed_in[row['id']] if other != name]
            if row['url'] in judge_urls:
                flags.append('judge test')
            pop, easy = row['pop'] or ('-', False)
            if easy:
                flags.append('easy')
            if not row['solution']:
                print(f"{row['id']:34} {'-':>4} {'-':>5} {pop:>12}  {'; '.join(['unsolved'] + flags)}")
                continue

            sim, closest = 0.0, ''
            for other in rows:
                if other is not row and other.get('shingles') and row['shingles']:
                    s = len(row['shingles'] & other['shingles']) / len(row['shingles'] | other['shingles'])
                    if s > sim:
                        sim, closest = s, other['id']
            if row['loc'] <= 8:
                flags.append('thin')
            if sim >= 0.5:
                flags.append(f'similar to {closest}')
            print(f"{row['id']:34} {row['loc']:>4} {sim:>5.2f} {pop:>12}  {'; '.join(flags)}")


if __name__ == '__main__':
    main()
