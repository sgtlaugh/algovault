import io
import os
import subprocess
import sys

from stress import rng, scaled

from fast_io import Scanner

SEPARATORS = [' ', '  ', '\t', '\n', '\r\n', '\n\n', ' \n ', '　', '\x0b']
TOKEN_CHARS = 'az09-+.é\U0001f9a6'


def random_text(count):
    tokens = [''.join(rng.choice(TOKEN_CHARS) for _ in range(rng.randint(1, 8))) for _ in range(count)]
    text = rng.choice(['', '\n', '  ']) + ''.join(t + rng.choice(SEPARATORS) for t in tokens)
    return text if rng.randint(0, 1) else text.rstrip()  # often no whitespace after the last token


def read_all(scanner):
    tokens = []
    while True:
        tokens.append(scanner.next_token())
        if tokens[-1] is None:
            return tokens


def main():
    for it in range(scaled(600)):
        text = random_text(rng.randint(0, 30 if it % 10 else 3000))
        sys.stdin = io.StringIO(text)
        scanner = Scanner()
        assert read_all(scanner) == text.split() + [None]
        try:
            scanner.next_token()
            assert False, 'the token stream must end after None'
        except StopIteration:
            pass

    # The real stdin of a fresh interpreter, with universal newlines and utf-8 decoding
    text = random_text(20000)
    script = 'import sys; sys.path.insert(0, sys.argv[1]); from fast_io import Scanner; s = Scanner(); t = []\n' \
             'while True:\n    t.append(s.next_token())\n    if t[-1] is None: break\nprint(repr(t))'
    library = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'code_library', 'python')
    out = subprocess.run([sys.executable, '-c', script, library], input=text.encode(), capture_output=True, check=True,
                         env=dict(os.environ, PYTHONIOENCODING='utf-8', PYTHONUTF8='1'))
    assert out.stdout.decode().strip() == repr(text.split() + [None])


if __name__ == '__main__':
    main()
