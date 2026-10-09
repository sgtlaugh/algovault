/***
 *
 * K Shortest Walks (Eppstein)
 * Lengths of the k shortest walks from s to t in a directed graph, vertices and edges may repeat
 *
 * Complexity: O((n + m) log m + k log k) time, O((n + m) log m + k) memory
 *
 * KShortestWalks g(n); g.add_edge(u, v, w): directed edge, 0-based, w >= 0
 * g.solve(s, t, k): the k smallest walk lengths in non-decreasing order, counted with multiplicity
 *     fewer than k values when fewer walks exist (empty when t is unreachable from s)
 *     when s == t the empty walk counts, so the first value is 0
 *     parallel edges give distinct walks
 *     the largest reported length plus n * max(w) must fit in long long: unreported candidates reach that far
 *
 * Dijkstra from t on the reversed graph gives a shortest path tree, every other edge u -> v is a sidetrack
 * costing w + d[v] - d[u] >= 0, and a walk is the tree path plus a sequence of sidetracks
 * h[u] is a persistent leftist heap of the sidetracks leaving the tree path from u, sharing h[parent]
 * The answers are popped best-first from a heap of (length, heap node): from a node, either replace its
 * sidetrack with a child in the same heap, or append the cheapest sidetrack after it
 *
 * Example:
 *   KShortestWalks g(3);
 *   g.add_edge(0, 1, 1), g.add_edge(1, 2, 1), g.add_edge(0, 2, 3), g.add_edge(2, 0, 1);
 *   g.solve(0, 2, 4);  // {2, 3, 5, 6}
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct KShortestWalks{
    static constexpr long long INF = numeric_limits<long long>::max();

    struct Edge{
        int u, v;
        long long w;
    };

    struct Node{
        long long key;
        int to, left, right, rank;
    };

    int n;
    vector<Edge> edges;
    vector<Node> pool;

    KShortestWalks(int n) : n(n) {}

    void add_edge(int u, int v, long long w){
        edges.push_back({u, v, w});
    }

    vector<long long> solve(int s, int t, int k){
        vector<vector<int>> rev(n), out(n);
        for (int id = 0; id < (int)edges.size(); id++){
            rev[edges[id].v].push_back(id);
            out[edges[id].u].push_back(id);
        }

        vector<long long> d(n, INF);
        vector<int> tree_edge(n, -1), order;
        priority_queue<pair<long long, int>, vector<pair<long long, int>>, greater<>> heap;
        d[t] = 0, heap.push({0, t});
        while (!heap.empty()){
            auto [du, u] = heap.top();
            heap.pop();
            if (du != d[u]) continue;
            order.push_back(u);
            for (int id : rev[u]){
                auto [x, v, w] = edges[id];
                if (du + w < d[x]) d[x] = du + w, tree_edge[x] = id, heap.push({d[x], x});
            }
        }

        if (d[s] == INF || k <= 0) return {};

        /// Dijkstra pops a vertex after its tree parent, so h[parent] is complete when h[u] needs it
        pool.clear();
        vector<int> h(n, -1);
        for (int u : order){
            vector<pair<long long, int>> sidetracks;
            for (int id : out[u]){
                auto [x, v, w] = edges[id];
                if (id != tree_edge[u] && d[v] != INF) sidetracks.push_back({w + d[v] - d[u], v});
            }
            sort(sidetracks.rbegin(), sidetracks.rend());

            /// A sorted chain hung on left children is already a valid leftist heap: every right child is empty
            int own = -1;
            for (auto [key, v] : sidetracks){
                pool.push_back({key, v, own, -1, 1});
                own = pool.size() - 1;
            }
            h[u] = meld(own, tree_edge[u] == -1 ? -1 : h[edges[tree_edge[u]].v]);
        }

        vector<long long> result = {d[s]};
        priority_queue<pair<long long, int>, vector<pair<long long, int>>, greater<>> best;
        if (h[s] != -1) best.push({d[s] + pool[h[s]].key, h[s]});
        while ((int)result.size() < k && !best.empty()){
            auto [len, x] = best.top();
            best.pop();
            result.push_back(len);

            const Node& node = pool[x];
            if (node.left != -1) best.push({len - node.key + pool[node.left].key, node.left});
            if (node.right != -1) best.push({len - node.key + pool[node.right].key, node.right});
            if (h[node.to] != -1) best.push({len + pool[h[node.to]].key, h[node.to]});
        }
        return result;
    }

private:
    int rank(int x) const{
        return x == -1 ? 0 : pool[x].rank;
    }

    /// Copies only the nodes on the merged right spines, so recursion depth and new nodes are O(log size)
    int meld(int a, int b){
        if (a == -1) return b;
        if (b == -1) return a;
        if (pool[b].key < pool[a].key) swap(a, b);

        const Node copy = pool[a];
        int c = pool.size();
        pool.push_back(copy);
        int right = meld(copy.right, b);

        Node& node = pool[c];
        node.right = right;
        if (rank(node.left) < rank(node.right)) swap(node.left, node.right);
        node.rank = rank(node.right) + 1;
        return c;
    }
};

int main(){
    KShortestWalks g(3);
    g.add_edge(0, 1, 1), g.add_edge(1, 2, 1), g.add_edge(0, 2, 3), g.add_edge(2, 0, 1);
    assert((g.solve(0, 2, 7) == vector<long long>{2, 3, 5, 6, 6, 7, 8}));
    assert((g.solve(0, 2, 1) == vector<long long>{2}));
    assert((g.solve(0, 2, 0) == vector<long long>{}));
    assert((g.solve(1, 0, 3) == vector<long long>{2, 5, 6}));

    KShortestWalks chain(3);
    chain.add_edge(0, 1, 4);
    assert((chain.solve(0, 2, 5) == vector<long long>{}));
    assert((chain.solve(1, 0, 5) == vector<long long>{}));
    assert((chain.solve(2, 2, 5) == vector<long long>{0}));

    KShortestWalks loop(1);
    loop.add_edge(0, 0, 2);
    assert((loop.solve(0, 0, 4) == vector<long long>{0, 2, 4, 6}));

    KShortestWalks parallel(2);
    parallel.add_edge(0, 1, 5), parallel.add_edge(0, 1, 5), parallel.add_edge(0, 1, 2);
    assert((parallel.solve(0, 1, 4) == vector<long long>{2, 5, 5}));

    KShortestWalks zero(3);
    zero.add_edge(0, 1, 0), zero.add_edge(1, 0, 0), zero.add_edge(1, 2, 4);
    assert((zero.solve(0, 2, 3) == vector<long long>{4, 4, 4}));
    return 0;
}
