"""
Dijkstra

Graph().add_edge(u, v, cost, directed=False), costs must be non-negative
dijkstra(graph, source, destination): shortest distance, or None when unreachable

Complexity: O((V + E) log E), a lazy-deletion heap holds one entry per edge relaxation
"""

from collections import defaultdict
from heapq import heappop, heappush


class Graph:
    def __init__(self):
        self.adj_list = defaultdict(list)

    def add_edge(self, from_node, to_node, cost, directed=False):
        """
        Adds an edge in the graph between from_node and to_node

        Parameters
        ----------
        from_node : Outward edge from this node
        to_node : Incoming edge to this node
        cost: Cost of traversing this edge
        directed: True if the edge is uni-directional, False if its bi-directional

        """

        if cost < 0: raise ValueError("Costs must be non-negative for dijkstra to work")

        self.adj_list[from_node].append((to_node, cost))
        if not directed:
            self.add_edge(to_node, from_node, cost, True)


def dijkstra(graph, source, destination):
    """
    Returns the shortest distance in the graph between source and destination nodes
    Returns None if source and destination are not connected or does not exist

    Parameters
    ----------
    graph : Graph object containing the nodes and edges
    source: Source node
    destination: Destination node

    """

    distance = {source: 0}
    queue = [(0, source)]

    while queue:
        (d, u) = heappop(queue)
        if d > distance[u]: continue
        if u == destination: return d

        for v, w in graph.adj_list[u]:
            if (v not in distance) or ((d + w) < distance[v]):
                distance[v] = d + w
                heappush(queue, (distance[v], v))


def main():
    import sys

    data = sys.stdin.buffer.read().split()
    pos = 0
    t = int(data[pos])
    pos += 1
    out = []
    for case in range(1, t + 1):
        n, m, s, e = map(int, data[pos:pos + 4])
        pos += 4
        graph = Graph()
        for _ in range(m):
            u, v, w = int(data[pos]), int(data[pos + 1]), int(data[pos + 2])
            pos += 3
            graph.add_edge(u, v, w)

        res = dijkstra(graph, s, e)
        out.append('Case #%d: %s' % (case, 'unreachable' if res is None else res))
    print('\n'.join(out))


if __name__ == '__main__':
    main()
