from collections import defaultdict
from heapq import heapify, heappop, heappush


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

    Complexity
    ----------
    O((E + V) log(V)) roughly
    E : Number of edges in the graph
    V : Number of nodes in the graph

    Parameters
    ----------
    graph : Graph object containing the nodes and edges
    source: Source node
    destination: Destination node

    """

    distance = {source: 0}
    queue = [(distance[source], source)]
    heapify(queue)

    while queue:
        (d, u) = heappop(queue)
        if u == destination: return d

        for v, w in graph.adj_list[u]:
            if (v not in distance) or ((d + w) < distance[v]):
                distance[v] = d + w
                heappush(queue, (distance[v], v))


if __name__ == "__main__":
    graph = Graph()
    graph.add_edge("Redwood City", "Emerald Hills", 75)
    graph.add_edge("Redwood City", "Menlo Park", 100)
    graph.add_edge("Menlo Park", "Stanford", 50)
    graph.add_edge("Stanford", "Portola Valley", 100)
    graph.add_edge("Emerald Hills", "Portola Valley", 80)
    graph.add_edge("Emerald Hills", "Woodside", 50)
    graph.add_edge("Menlo Park", "Woodside", 200)

    assert dijkstra(graph, "Stanford", "Woodside") == 230
    assert dijkstra(graph, "Menlo Park", "Sunnyvale") is None
