from stress import rng, scaled

from dijkstra import Graph, dijkstra

INF = float('inf')


def floyd_warshall(nodes, arcs):
    dist = {u: {v: (0 if u == v else INF) for v in nodes} for u in nodes}
    for u, v, w in arcs:
        dist[u][v] = min(dist[u][v], w)

    for k in nodes:
        for i in nodes:
            for j in nodes:
                dist[i][j] = min(dist[i][j], dist[i][k] + dist[k][j])

    return dist


def main():
    for it in range(scaled(1000)):
        n = rng.randint(1, 9)
        names = [f'stop-{i}' for i in range(n)] if it % 2 else list(range(n))
        graph, arcs = Graph(), []
        for _ in range(rng.randint(0, 3 * n)):
            u, v, w = rng.choice(names), rng.choice(names), rng.choice([0, rng.randint(0, 20), rng.randint(0, 10**12)])
            directed = rng.random() < 0.5
            graph.add_edge(u, v, w, directed)
            arcs.append((u, v, w))
            if not directed:
                arcs.append((v, u, w))

        dist = floyd_warshall(names, arcs)
        for s in names:
            for t in names:
                expected = dist[s][t]
                assert dijkstra(graph, s, t) == (None if expected == INF else expected)
        assert dijkstra(graph, names[0], 'nowhere') is None

    try:
        Graph().add_edge(0, 1, -1)
        assert False, 'a negative cost must be rejected'
    except ValueError:
        pass


if __name__ == '__main__':
    main()
