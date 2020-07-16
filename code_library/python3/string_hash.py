class StringHash:
    def __init__(self, text, base=997, mod=10 ** 18 + 3):
        n = len(text)

        self.text = text
        self.base = base
        self.mod = mod
        self.power = [1] * (n + 2)
        self.prefix_hash = [0] * (n + 2)

        for i in range(1, n + 1):
            self.power[i] = self.power[i - 1] * self.base % self.mod
            self.prefix_hash[i] = (self.prefix_hash[i - 1] * self.base + ord(self.text[i - 1]) + 13) % self.mod

    def get_hash(self, l, r):
        h = (self.prefix_hash[r + 1] - self.power[r - l + 1] * self.prefix_hash[l]) % self.mod
        return h + self.mod if h < 0 else h


def main():
    sh = StringHash('abracadabra')

    assert sh.get_hash(0, 0) == sh.get_hash(3, 3)
    assert sh.get_hash(0, 3) == sh.get_hash(7, 10)
    assert sh.get_hash(0, 4) != sh.get_hash(5, 9)


if __name__ == '__main__':
    main()
