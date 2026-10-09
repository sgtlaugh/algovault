import random
from functools import lru_cache

from stress import SEED, rng, scaled

from graph_matching import maximum_matching


def brute(n, edges):
    adj = [0] * n
    for u, v in edges:
        adj[u] |= 1 << v
        adj[v] |= 1 << u

    @lru_cache(maxsize=None)
    def best(used):
        """Decide the lowest unmatched node: leave it single or pair it with a free neighbour"""
        if used == (1 << n) - 1:
            return 0
        x = (~used & (used + 1)).bit_length() - 1
        res = best(used | 1 << x)
        for y in range(n):
            if adj[x] >> y & 1 and not used >> y & 1:
                res = max(res, 1 + best(used | 1 << x | 1 << y))
        return res

    return best(0)


def main():
    """The rank test errs with probability at most n / mod per call, ~1e-8 here, so any mismatch is a bug"""
    random.seed(SEED)  # the library draws its Tutte matrix from the global generator
    for _ in range(scaled(300)):
        n = rng.randint(1, 12)
        density = rng.random()
        edges = [(u, v) if rng.random() < 0.5 else (v, u) for u in range(n) for v in range(u + 1, n) if rng.random() < density]
        if rng.random() < 0.5:
            # Disjoint odd cycles, each leaves a node unmatched, plus a few chords that may or may not help
            order = rng.sample(range(n), n)
            edges, start = [], 0
            while start + 3 <= n:
                size = rng.choice([3, 3, 5])
                cycle = order[start:start + size]
                edges += [(cycle[i], cycle[(i + 1) % len(cycle)]) for i in range(len(cycle))]
                start += size
            edges += [tuple(rng.sample(range(n), 2)) for _ in range(rng.randint(0, 2)) if n >= 2]
            edges = list({(min(e), max(e)) for e in edges})

        rng.shuffle(edges)
        nodes = max((max(e) for e in edges), default=-1) + 1  # the library sizes the graph by the largest label
        assert maximum_matching(edges) == brute(nodes, edges)


if __name__ == '__main__':
    main()
