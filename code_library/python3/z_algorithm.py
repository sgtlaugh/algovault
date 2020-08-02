def z_function(text):
    """
    Z Algorithm in O(n)
    :param text: text string to process
    :return: the z_array, where z_array[i] = length of the longest common prefix of text[i:] and text
    """

    n = len(text)
    z_array = [0 for _ in range(n)]
    l = r = 0

    for i in range(1, n):
        z = z_array[i - l]
        if i + z >= r:
            z = max(r - i, 0)
            while i + z < n and text[z] == text[i + z]:
                z += 1

            l, r = i, i + z

        z_array[i] = z

    z_array[0] = n
    return z_array


def main():
    assert z_function('a') == [1]
    assert z_function('aaaa') == [4, 3, 2, 1]
    assert z_function('aaabaab') == [7, 2, 1, 0, 2, 1, 0]
    assert z_function('abababab') == [8, 0, 6, 0, 4, 0, 2, 0]
    assert z_function('abracadabra') == [11, 0, 0, 1, 0, 1, 0, 4, 0, 0, 1]


if __name__ == '__main__':
    main()
