from stress import rng, scaled

from string_hash import StringHash


def main():
    """Hash equality against substring equality, a collision has probability ~1e-18 per pair and should never appear"""
    for it in range(scaled(300)):
        alphabet = 'ab'[:rng.randint(1, 2)] if it % 3 else 'ab\u0161\U0001f600'  # U+0161 is 'a' + 256, any byte truncation collides
        text = ''.join(rng.choice(alphabet) for _ in range(rng.randint(1, 40 if it % 10 else 3000)))
        h = StringHash(text)
        for _ in range(100):
            length = rng.randint(1, len(text))
            l1, l2 = rng.randint(0, len(text) - length), rng.randint(0, len(text) - length)
            equal = text[l1:l1 + length] == text[l2:l2 + length]
            assert (h.get_hash(l1, l1 + length - 1) == h.get_hash(l2, l2 + length - 1)) == equal

    # Longer than the old fixed power table
    text = 'ab' * 75000 + 'a'
    h = StringHash(text)
    assert h.get_hash(0, len(text) - 3) == h.get_hash(2, len(text) - 1)
    assert h.get_hash(0, len(text) - 2) != h.get_hash(1, len(text) - 1)


if __name__ == '__main__':
    main()
