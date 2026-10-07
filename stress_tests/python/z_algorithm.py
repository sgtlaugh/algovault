from stress import rng, scaled

from z_algorithm import z_function


def brute(text):
    z = []
    for i in range(len(text)):
        k = 0
        while i + k < len(text) and text[k] == text[i + k]:
            k += 1
        z.append(k)
    return z


def main():
    for it in range(scaled(3000)):
        alphabet = 'ab'[:rng.randint(1, 2)] if it % 4 else 'abcz'
        text = ''.join(rng.choice(alphabet) for _ in range(rng.randint(0, 30 if it % 10 else 300)))
        if it % 5 == 0:
            text = (text[:rng.randint(1, 5)] * 100)[:len(text)]  # periodic, long z boxes
        assert z_function(text) == brute(text), text


if __name__ == '__main__':
    main()
