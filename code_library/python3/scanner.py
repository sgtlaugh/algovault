import sys


class Scanner():
    def __init__(self):
        self.tokens = []
        self.index = -1

        for line in sys.stdin:
            self.tokens.extend(line.split())

    def next_token(self):
        self.index += 1
        return None if self.index == len(self.tokens) \
            else self.tokens[self.index]


def main():
    # uncomment if need to read input from file, otherwise input ends with EOF from stdin
    # sys.stdin = open('input.txt')

    sc = Scanner()
    t = int(sc.next_token())

    for _ in range(t):
        line = sc.next_token()
        print(line)


if __name__ == '__main__':
    main()
