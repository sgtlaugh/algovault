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
    # uncomment if need to read input from file, otherwise input ends with EOF from stdin
    # sys.stdin = open('input.txt')

    sc = Scanner()
    t = int(sc.next_token())

    for _ in range(t):
        line = sc.next_token()
        print(line)


if __name__ == '__main__':
    main()
