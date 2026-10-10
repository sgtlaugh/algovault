"""
Fast input scanner

Scanner().next_token(): the next whitespace separated token from sys.stdin, then None at the end of input

Complexity: O(total input length) over all calls, memory holds one line and its tokens at a time
"""

import io
import sys


class Scanner():
    def __init__(self):
        self.tokens = self.get_tokens()

    def next_token(self):
        return next(self.tokens)

    def get_tokens(self):
        for line in sys.stdin:
            for token in line.split():
                yield token
        yield None


def main():
    # Typical usage reads from stdin, or from a file after sys.stdin = open('input.txt'):
    #     sc = Scanner()
    #     t = int(sc.next_token())

    sys.stdin = io.StringIO("3\nwombat  42\n\n   quokka\n")
    sc = Scanner()
    assert [sc.next_token() for _ in range(5)] == ['3', 'wombat', '42', 'quokka', None]


if __name__ == '__main__':
    main()
