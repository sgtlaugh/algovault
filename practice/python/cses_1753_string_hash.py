"""
String Hash

Polynomial rolling hash, get_hash(l, r) is the hash of text[l..r] inclusive

Complexity: O(n) to build, O(1) per get_hash, n + 1 prefix hashes per instance plus powers shared up to the longest text
"""


class StringHash:
    # Class attributes so contest code that assigns a global MOD or BASE cannot corrupt live instances
    _BASE = 997
    _MOD = 10**18 + 3
    _powers = [1]  # shared by all instances, grown to the longest text seen

    def __init__(self, text):
        base, mod, powers = StringHash._BASE, StringHash._MOD, StringHash._powers
        while len(powers) <= len(text):
            powers.append(powers[-1] * base % mod)

        h, prefix_hash = 0, [0]
        for c in text:
            h = (h * base + ord(c) + 13) % mod
            prefix_hash.append(h)
        self.prefix_hash = prefix_hash

    def get_hash(self, l, r):
        return (self.prefix_hash[r + 1] - StringHash._powers[r - l + 1] * self.prefix_hash[l]) % StringHash._MOD


def main():
    import sys

    text, pattern = sys.stdin.read().split()
    n, m = len(text), len(pattern)
    target = StringHash(pattern).get_hash(0, m - 1)
    sh = StringHash(text)
    print(sum(1 for i in range(n - m + 1) if sh.get_hash(i, i + m - 1) == target))


if __name__ == '__main__':
    main()
