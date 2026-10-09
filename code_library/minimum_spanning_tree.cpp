/***
 *
 * Minimum Spanning Tree (Kruskal)
 * Minimum spanning forest of an undirected weighted graph
 *
 * Complexity: O(m log m)
 *
 * minimum_spanning_tree(n, edges): edges[i] = {u, v, w}, 0-based vertices, long long weights (negative allowed)
 *     the total weight of the forest must fit in long long
 * Returns {total weight, indices of the chosen edges, connected}
 *     connected is true when the forest is one tree spanning all n vertices
 *     for a disconnected graph the result is a minimum spanning forest, one tree per component
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct MSTResult{
    long long weight;
    vector<int> chosen;
    bool connected;
};

MSTResult minimum_spanning_tree(int n, const vector<array<long long, 3>>& edges){
    vector<int> order(edges.size()), parent(n), rank(n, 0);
    iota(order.begin(), order.end(), 0);
    iota(parent.begin(), parent.end(), 0);
    stable_sort(order.begin(), order.end(), [&](int a, int b){ return edges[a][2] < edges[b][2]; });

    auto find = [&](int x){
        while (x != parent[x]) x = parent[x] = parent[parent[x]];
        return x;
    };

    MSTResult res{0, {}, false};
    for (int i : order){
        int a = find(edges[i][0]), b = find(edges[i][1]);
        if (a == b) continue;
        if (rank[a] < rank[b]) swap(a, b);
        parent[b] = a;
        if (rank[a] == rank[b]) rank[a]++;
        res.weight += edges[i][2];
        res.chosen.push_back(i);
    }
    res.connected = (int)res.chosen.size() == max(0, n - 1);
    return res;
}

int main(){
    /***
     *   0 --1-- 1
     *   |  \    |
     *   4   3   2
     *   |    \  |
     *   3 --5-- 2
    ***/
    vector<array<long long, 3>> edges = {{0, 1, 1}, {1, 2, 2}, {0, 2, 3}, {0, 3, 4}, {3, 2, 5}};
    MSTResult r = minimum_spanning_tree(4, edges);
    assert(r.weight == 7 && r.connected);
    assert((r.chosen == vector<int>{0, 1, 3}));

    MSTResult forest = minimum_spanning_tree(5, {{0, 1, -4}, {2, 3, 6}, {3, 2, 1}});
    assert(forest.weight == -3 && !forest.connected && forest.chosen.size() == 2);

    assert(minimum_spanning_tree(1, {}).connected && minimum_spanning_tree(1, {}).weight == 0);
    assert(minimum_spanning_tree(0, {}).connected);
    assert(!minimum_spanning_tree(2, {}).connected);
    return 0;
}
