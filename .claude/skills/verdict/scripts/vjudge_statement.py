#!/usr/bin/env python3
"""vjudge_statement.py <problem URL or vjudge id like SPOJ-QTREE>: a problem's limits and statement as plain text.

Reads it through vjudge.net, which serves the Codeforces and SPOJ statements their own sites hide behind Cloudflare.
Credentials: ~/.config/vjudge/credentials, username then password on two lines, mode 600, never inside a git work tree.
The password is never printed; the session cookie is kept in ~/.config/vjudge/cookies, mode 600.
Statements are cached in ~/.cache/algovault/statements, so parallel agents fetch each problem once.
"""
import html
import html.parser
import http.cookiejar
import json
import os
import re
import subprocess
import sys

import requests

BASE = 'https://vjudge.net'
DIR = os.path.expanduser('~/.config/vjudge')
CACHE = os.path.expanduser('~/.cache/algovault/statements')

# A URL names its problem exactly, so no README title or judge display name is ever parsed
VJUDGE_IDS = [
    (r'codeforces\.com/(?:problemset/problem|contest)/(\d+)/(?:problem/)?(\w+)', 'CodeForces-{0}{1}'),
    (r'codeforces\.com/gym/(\d+)/problem/(\w+)', 'Gym-{0}{1}'),
    (r'spoj\.com/problems/(\w+)', 'SPOJ-{0}'),
    (r'onlinejudge\.org/external/\d+/(\d+)\.pdf', 'UVA-{0}'),
    (r'onlinejudge\.org/.*[?&]problem=(\d+)', 'UVA-{0}'),
    (r'cses\.fi/problemset/task/(\d+)', 'CSES-{0}'),
    (r'atcoder\.jp/contests/\w+/tasks/(\w+)', 'AtCoder-{0}'),
    (r'acm\.timus\.ru/problem\.aspx\?(?:space=1&)?num=(\d+)', 'URAL-{0}'),
    (r'onlinejudge\.u-aizu\.ac\.jp/(?:problems|courses/[\w/]+/description\.jsp\?id=)/?(\w+)', 'Aizu-{0}'),
    (r'judge\.u-aizu\.ac\.jp/onlinejudge/description\.jsp\?id=(\w+)', 'Aizu-{0}'),
    (r'acm\.hdu\.edu\.cn/showproblem\.php\?pid=(\d+)', 'HDU-{0}'),
    (r'poj\.org/problem\?id=(\d+)', 'POJ-{0}'),
    (r'codechef\.com/(?:problems|[A-Z0-9]+/problems)/(\w+)', 'CodeChef-{0}'),
    (r'open\.kattis\.com/problems/(\w+)', 'Kattis-{0}'),
    (r'hackerrank\.com/(?:contests/[\w-]+/)?challenges/([\w-]+)', 'HackerRank-{0}'),
    (r'usaco\.org/index\.php\?page=viewproblem2&cpid=(\d+)', 'USACO-{0}'),
    (r'oj\.uz/problem/view/(\w+)', 'OJUZ-{0}'),
    (r'dmoj\.ca/problem/(\w+)', 'DMOJ-{0}'),
    (r'loj\.ac/p/(\d+)', 'LibreOJ-{0}'),
    (r'qoj\.ac/problem/(\d+)', 'QOJ-{0}'),
    (r'csacademy\.com/contest/archive/task/(\w+)', 'CSAcademy-{0}'),
    (r'toph\.co/p/([\w-]+)', 'Toph-{0}'),
    (r'yukicoder\.me/problems/no/(\d+)', 'Yukicoder-{0}'),
    (r'uoj\.ac/problem/(\d+)', 'UniversalOJ-{0}'),
    (r'luogu\.com\.cn/problem/(\w+)', '洛谷-{0}'),
    (r'judge\.yosupo\.jp/problem/(\w+)', 'Yosupo-{0}'),
    (r'eolymp\.com/\w+/problems/(\d+)', 'EOlymp-{0}'),
]


class Text(html.parser.HTMLParser):
    """Plain text from statement HTML: <pre> blocks verbatim, since samples must match byte for byte."""

    BLOCK = {'p', 'div', 'br', 'h1', 'h2', 'h3', 'h4', 'h5', 'tr', 'ul', 'ol', 'table', 'section'}

    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.out, self.pre = [], 0

    def handle_starttag(self, tag, attrs):
        if tag == 'pre':
            self.pre += 1
            self.out.append('\n```\n')
        elif self.pre:
            if tag == 'br':
                self.out.append('\n')
        elif tag == 'li':
            self.out.append('\n- ')
        elif tag == 'sup':
            self.out.append('^')
        elif tag == 'img':
            self.out.append('[image: ' + dict(attrs).get('src', '') + ']')
        elif tag in self.BLOCK:
            self.out.append('\n')

    def handle_endtag(self, tag):
        if tag == 'pre':
            self.pre -= 1
            if not ''.join(self.out[-1:]).endswith('\n'):
                self.out.append('\n')
            self.out.append('```\n')
        elif self.pre:
            # Codeforces puts each sample line in its own div inside the pre
            if tag == 'div' and not ''.join(self.out[-1:]).endswith('\n'):
                self.out.append('\n')
        elif tag in self.BLOCK or tag in ('td', 'th'):
            self.out.append('\n' if tag in self.BLOCK else ' ')

    def handle_data(self, data):
        if self.pre:
            data = data.replace('\r\n', '\n')
            if self.out and self.out[-1].endswith('```\n'):
                data = data.lstrip('\n')
        self.out.append(data if self.pre else re.sub(r'\s+', ' ', data))


def outside_git(path):
    """Refuses a secrets path inside any git work tree, so credentials and cookies can never be committed."""
    d = os.path.dirname(os.path.realpath(path))
    while d != os.path.dirname(d):
        if os.path.exists(os.path.join(d, '.git')):
            sys.exit(f'{path} is inside the git work tree {d}, keep vjudge secrets in ~/.config/vjudge')
        d = os.path.dirname(d)
    return path


def vjudge_id(arg):
    if re.fullmatch(r'[^/\s]+-[^/\s]+', arg) and '://' not in arg:
        return arg
    for pattern, fmt in VJUDGE_IDS:
        m = re.search(pattern, arg)
        if m:
            return fmt.format(*m.groups())
    sys.exit(f'no vjudge mapping for {arg}, fetch the judge page instead')


def credentials():
    path = outside_git(os.path.join(DIR, 'credentials'))
    if not os.path.exists(path):
        sys.exit(f'{path} is missing, CONTRIBUTING.md "Judge credentials" shows how to create it')
    if os.stat(path).st_mode & 0o077:
        sys.exit(f'{path} must be readable by you only: chmod 600 {path}')
    with open(path) as f:
        user, password = f.read().split('\n')[:2]
    return user.strip(), password


def session():
    s = requests.Session()
    s.headers['User-Agent'] = 'Mozilla/5.0 (X11; Linux x86_64) algovault-statement'
    jar = http.cookiejar.LWPCookieJar(outside_git(os.path.join(DIR, 'cookies')))
    if os.path.exists(jar.filename):
        jar.load(ignore_discard=True)
    s.cookies = jar
    return s


def login(s):
    user, password = credentials()
    s.get(BASE + '/', timeout=30)
    s.post(BASE + '/user/login', data={'username': user, 'password': password},
           headers={'X-Requested-With': 'XMLHttpRequest', 'Referer': BASE + '/'}, timeout=30)
    s.cookies.save(ignore_discard=True)
    os.chmod(s.cookies.filename, 0o600)


def problem_page(s, pid):
    """The problem's dataJson, logging in again when the session cookie has expired."""
    for attempt in range(2):
        r = s.get(f'{BASE}/problem/{pid}', timeout=30, allow_redirects=False)
        if r.status_code == 200:
            m = re.search(r'name="dataJson">(.*?)</textarea>', r.text, re.S)
            if not m:
                sys.exit(f'{pid}: no problem data on the page, the site layout changed')
            return json.loads(html.unescape(m.group(1)))
        if r.status_code == 404:
            sys.exit(f'{pid}: not found on vjudge')
        if attempt == 0:
            login(s)
    sys.exit('vjudge login failed, check ~/.config/vjudge/credentials')


def pdf_text(s, path):
    """UVa and oj.uz statements are PDFs on vjudge's CDN, read through pdftotext from a pipe so no file is written."""
    cdn = s.get(BASE + '/util/cfg', timeout=30).json()['cdnBaseUrl']
    pdf = s.get(cdn + path, timeout=60).content
    text = subprocess.run(['pdftotext', '-layout', '-', '-'], input=pdf, capture_output=True, check=True).stdout
    return text.decode()


def statement(pid):
    s = session()
    data = problem_page(s, pid)
    brief = next((b for b in data['descBriefs'] if b.get('official')), data['descBriefs'][0])
    page = s.get(f"{BASE}/problem/description/{brief['key']}?{brief['version']}", timeout=30).text
    desc = json.loads(html.unescape(re.search(r'data-json-container"[^>]*>(.*?)</textarea>', page, re.S).group(1)))

    lines = [pid]
    for p in data.get('properties', []):
        if not p.get('spoiler'):
            lines.append(f"{p['title']}: {re.sub(r'<[^>]+>', '', p['content'])}")
    for section in desc['sections']:
        value = section['value']
        if section.get('title'):
            lines.append('\n## ' + section['title'])
        if value.get('format') == 'HTML':
            parser = Text()
            parser.feed(value['content'])
            body = ''.join(parser.out)
        else:
            body = value['content']
        body = re.sub(r'\[pdf:CDN_BASE_URL/([^\]]+)\]', lambda m: pdf_text(s, m.group(1)), body)
        body = re.sub(r'(?m)^[ \t]+$', '', body)
        lines.append(re.sub(r'\n{3,}', '\n\n', body).strip())
    return '\n'.join(lines) + '\n'


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    pid = vjudge_id(sys.argv[1])
    cached = os.path.join(CACHE, pid.replace('/', '_') + '.txt')
    if not os.path.exists(cached):
        text = statement(pid)
        os.makedirs(CACHE, exist_ok=True)
        with open(cached, 'w') as f:
            f.write(text)
    with open(cached) as f:
        sys.stdout.write(f.read())


if __name__ == '__main__':
    main()
