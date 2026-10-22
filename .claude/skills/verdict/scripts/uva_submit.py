#!/usr/bin/env python3
"""UVa Online Judge client: login, submit, verdict.

uva_submit.py login                      check that the stored credentials log in
uva_submit.py submit <problem> <file>    submit as C++11 and wait for the verdict
uva_submit.py last [count]               latest verdicts from uHunt

Credentials: ~/.config/uva/credentials, username then password on two lines, mode 600, never inside a git work tree.
The password is never printed; the session cookie is kept in ~/.config/uva/cookies, mode 600.
Submits only when run with submit: /verdict calls it after the user says submit, never on its own.
"""
import html.parser
import http.cookiejar
import os
import re
import sys
import time
import urllib.parse

import requests

BASE = 'https://onlinejudge.org'
UHUNT = 'https://uhunt.onlinejudge.org/api'
DIR = os.path.expanduser('~/.config/uva')
LANGUAGES = {'c': 1, 'java': 2, 'c++': 3, 'pascal': 4, 'c++11': 5, 'python': 6}
VERDICTS = {10: 'Submission error', 15: 'Cannot be judged', 20: 'In queue', 30: 'Compile error',
            35: 'Restricted function', 40: 'Runtime error', 45: 'Output limit', 50: 'Time limit',
            60: 'Memory limit', 70: 'Wrong answer', 80: 'Presentation error', 90: 'Accepted'}


class LoginForm(html.parser.HTMLParser):
    """Collects the action and every input of the form with id mod_loginform."""

    def __init__(self):
        super().__init__()
        self.action, self.fields, self.inside = None, {}, False

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if tag == 'form' and a.get('id') == 'mod_loginform':
            self.inside, self.action = True, a.get('action')
        elif tag == 'input' and self.inside and a.get('name') and a.get('type') != 'checkbox':
            self.fields[a['name']] = a.get('value', '')

    def handle_endtag(self, tag):
        if tag == 'form':
            self.inside = False


def outside_git(path):
    """Refuses a secrets path inside any git work tree, so credentials and cookies can never be committed."""
    d = os.path.dirname(os.path.realpath(path))
    while d != os.path.dirname(d):
        if os.path.exists(os.path.join(d, '.git')):
            sys.exit(f'{path} is inside the git work tree {d}, keep UVa secrets in ~/.config/uva')
        d = os.path.dirname(d)
    return path


def credentials():
    path = outside_git(os.path.join(DIR, 'credentials'))
    if os.stat(path).st_mode & 0o077:
        sys.exit(f'{path} must be readable by you only: chmod 600 {path}')
    with open(path) as f:
        user, password = f.read().split('\n')[:2]
    return user.strip(), password


def session():
    s = requests.Session()
    s.headers['User-Agent'] = 'Mozilla/5.0 (X11; Linux x86_64) uva-cli'
    jar = http.cookiejar.LWPCookieJar(outside_git(os.path.join(DIR, 'cookies')))
    if os.path.exists(jar.filename):
        jar.load(ignore_discard=True)
    s.cookies = jar
    return s


def save(s):
    s.cookies.save(ignore_discard=True)
    os.chmod(s.cookies.filename, 0o600)


def logged_in(page):
    return 'mod_loginform' not in page and re.search(r'task=logout|Logout', page) is not None


def login(s):
    page = s.get(BASE + '/', timeout=30).text
    if '<title>Online Judge - Offline</title>' in page:
        sys.exit('UVa is offline, try again later')
    if logged_in(page):
        return True
    form = LoginForm()
    form.feed(page)
    if not form.action:
        sys.exit('login form not found on the home page, the site layout changed')
    user, password = credentials()
    form.fields.update(username=user, passwd=password, remember='yes')
    s.post(form.action, data=form.fields, timeout=30)
    ok = logged_in(s.get(BASE + '/', timeout=30).text)
    if ok:
        save(s)
    return ok


def uid():
    user, _ = credentials()
    return requests.get(f'{UHUNT}/uname2uid/{user}', timeout=30).text.strip()


def last(count=5):
    data = requests.get(f'{UHUNT}/subs-user-last/{uid()}/{count}', timeout=30).json()
    problems = {}
    for sub in data['subs']:
        sid, pid, verdict, runtime = sub[0], sub[1], sub[2], sub[3]
        if pid not in problems:
            problems[pid] = requests.get(f'{UHUNT}/p/id/{pid}', timeout=30).json()
        p = problems[pid]
        print(f"{sid}  UVa {p['num']}  {VERDICTS.get(verdict, verdict)}  {runtime} ms")
    return data['subs']


def submit(s, problem, path, language='c++11'):
    if not login(s):
        sys.exit('login failed, check ~/.config/uva/credentials')
    with open(path) as f:
        code = f.read()
    r = s.post(BASE + '/index.php?option=com_onlinejudge&Itemid=25&page=save_submission',
               data={'problemid': '', 'category': '', 'localid': str(problem),
                     'language': str(LANGUAGES[language]), 'code': code}, timeout=60)
    m = re.search(r'Submission received with ID (\d+)', urllib.parse.unquote_plus(r.url) + ' ' + r.text)
    if not m:
        msg = re.search(r'mosmsg=([^&]+)', r.url)
        sys.exit('submission not confirmed: ' + (urllib.parse.unquote_plus(msg.group(1)) if msg else r.url))
    sid = int(m.group(1))
    print(f'submitted UVa {problem} as {language}, id {sid}')
    for _ in range(60):
        time.sleep(5)
        for sub in requests.get(f'{UHUNT}/subs-user-last/{uid()}/10', timeout=30).json()['subs']:
            if sub[0] == sid and sub[2] not in (0, 20):
                print(f'{VERDICTS.get(sub[2], sub[2])}  {sub[3]} ms')
                return sub[2]
    print('no verdict after 5 minutes, check later with: uva_submit.py last')


def main():
    if len(sys.argv) < 2 or sys.argv[1] not in ('login', 'submit', 'last'):
        sys.exit(__doc__)
    s = session()
    if sys.argv[1] == 'login':
        print('login ok' if login(s) else 'login failed')
    elif sys.argv[1] == 'submit':
        sys.exit(0 if submit(s, sys.argv[2], sys.argv[3]) == 90 else 1)
    else:
        last(int(sys.argv[2]) if len(sys.argv) > 2 else 5)


if __name__ == '__main__':
    main()
