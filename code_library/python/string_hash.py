BASE = 997
MOD = 10**18 + 3

class StringHash:
    _powers = [1]  # shared by all instances, grown to the longest text seen

    def __init__(self, text):
        powers = StringHash._powers
        while len(powers) <= len(text):
            powers.append(powers[-1] * BASE % MOD)

        self.prefix_hash = [0] * (len(text) + 2)
        for i in range(1, len(text) + 1):
            self.prefix_hash[i] = (self.prefix_hash[i - 1] * BASE + ord(text[i - 1]) + 13) % MOD

    def get_hash(self, l, r):
        return (self.prefix_hash[r + 1] - StringHash._powers[r - l + 1] * self.prefix_hash[l]) % MOD


def main():
    sh = StringHash('abracadabra')

    assert sh.get_hash(0, 0) == sh.get_hash(3, 3)
    assert sh.get_hash(0, 3) == sh.get_hash(7, 10)
    assert sh.get_hash(0, 4) != sh.get_hash(5, 9)


if __name__ == '__main__':
    main()
