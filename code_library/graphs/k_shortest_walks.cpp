/***
 *
 * K Shortest Walks (Eppstein)
 * Lengths of the k shortest walks from s to t in a directed graph, vertices and edges may repeat
 *
 * Complexity: O((n + m) log m + k log k) time, O((n + m) log m + k) memory
 *
 * KShortestWalks g(n): vertices 0 to n - 1, no edges
 * g.add_edge(u, v, w): directed edge, w >= 0
 * g.solve(s, t, k): the k smallest walk lengths in non-decreasing order, fewer if fewer walks exist
 *
 * Notes:
 *   - Lengths count with multiplicity and parallel edges give distinct walks
 *   - When s == t the empty walk counts, so the first value is 0
 *   - The largest reported length plus n max(w) must fit in long long
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// BEGIN COPY LeftistHeap from code_library/data_structures/leftist_heap.cpp
template <typename T, typename Compare = less<T>>
struct LeftistHeap{
    static constexpr int EMPTY = -1;

    struct Node{
        T val;
        int left, right, rank;
    };

    vector<Node> nodes;
    Compare cmp;

    /// Copies only the nodes on the merged right spines, whose length is bounded by the ranks
    int meld(int a, int b){
        if (a == EMPTY) return b;
        if (b == EMPTY) return a;
        if (cmp(nodes[b].val, nodes[a].val)) swap(a, b);

        Node copy = nodes[a];
        nodes.push_back(copy);
        int c = (int)nodes.size() - 1;

        int right = meld(copy.right, b);  /// the call can reallocate nodes, assign after it
        nodes[c].right = right;
        if (rank(nodes[c].left) < rank(right)) swap(nodes[c].left, nodes[c].right);
        nodes[c].rank = rank(nodes[c].right) + 1;
        return c;
    }

    int pop(int h){
        assert(h != EMPTY);
        return meld(nodes[h].left, nodes[h].right);
    }

    int push(int h, const T& x){
        nodes.push_back({x, EMPTY, EMPTY, 1});
        return meld(h, (int)nodes.size() - 1);
    }

    T top(int h) const{
        assert(h != EMPTY);
        return nodes[h].val;
    }

    int rank(int h) const{
        return h == EMPTY ? 0 : nodes[h].rank;
    }
};
/// END COPY LeftistHeap

struct KShortestWalks{
    static constexpr long long INF = numeric_limits<long long>::max();

    struct Edge{
        int u, v;
        long long w;
    };

    int n;
    vector<Edge> edges;
    LeftistHeap<pair<long long, int>> pool;

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
        pool.nodes.clear();
        vector<int> h(n, pool.EMPTY);
        for (int u : order){
            vector<pair<long long, int>> sidetracks;
            for (int id : out[u]){
                auto [x, v, w] = edges[id];
                if (id != tree_edge[u] && d[v] != INF) sidetracks.push_back({w + d[v] - d[u], v});
            }
            sort(sidetracks.rbegin(), sidetracks.rend());

            /// Sorted decreasing, each push of a smaller value becomes the root above the old heap, copying O(1) nodes
            int own = pool.EMPTY;
            for (auto sidetrack : sidetracks) own = pool.push(own, sidetrack);
            h[u] = pool.meld(own, tree_edge[u] == -1 ? pool.EMPTY : h[edges[tree_edge[u]].v]);
        }

        vector<long long> result = {d[s]};
        priority_queue<pair<long long, int>, vector<pair<long long, int>>, greater<>> best;
        if (h[s] != pool.EMPTY) best.push({d[s] + pool.top(h[s]).first, h[s]});
        while ((int)result.size() < k && !best.empty()){
            auto [len, x] = best.top();
            best.pop();
            result.push_back(len);

            const auto& node = pool.nodes[x];
            auto [key, to] = node.val;
            if (node.left != pool.EMPTY) best.push({len - key + pool.top(node.left).first, node.left});
            if (node.right != pool.EMPTY) best.push({len - key + pool.top(node.right).first, node.right});
            if (h[to] != pool.EMPTY) best.push({len + pool.top(h[to]).first, h[to]});
        }

        return result;
    }
};

int main(){
    /// 0 -> 1 -> 2 weighs 2, the direct 0 -> 2 weighs 3, and 2 -> 0 of weight 1 lets walks loop
    KShortestWalks g(3);
    g.add_edge(0, 1, 1), g.add_edge(1, 2, 1), g.add_edge(0, 2, 3), g.add_edge(2, 0, 1);
    assert((g.solve(0, 2, 4) == vector<long long>{2, 3, 5, 6}));  /// 0 1 2, 0 2, 0 1 2 0 1 2, then 2 walks of length 6
    assert((g.solve(1, 0, 3) == vector<long long>{2, 5, 6}));      /// walks may repeat vertices and edges

    KShortestWalks chain(3);
    chain.add_edge(0, 1, 4);
    assert(chain.solve(0, 2, 5).empty());                          /// 2 is unreachable
    assert((chain.solve(2, 2, 5) == vector<long long>{0}));        /// only the empty walk
    return 0;
}
