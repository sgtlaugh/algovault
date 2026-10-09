#include "../common.h"

#define main library_main
#include "../../code_library/graphs/directed_mst.cpp"
#undef main

/// Every choice of one incoming edge per non-root node, kept when all parent chains reach the root
pair<bool, long long> brute(int n, int root, const vector<Edge>& edges){
    vector<vector<int>> incoming(n);
    for (auto& e : edges) if (e.u != e.v) incoming[e.v].push_back(&e - &edges[0]);

    /// rest[x] is the cheapest any completion of nodes x..n-1 can cost, which keeps pruning valid with negative weights
    vector<long long> rest(n + 1, 0);
    for (int x = n - 1; x >= 0; x--){
        long long cheapest = 0;
        for (int id : incoming[x]) cheapest = min(cheapest, edges[id].w);
        rest[x] = rest[x + 1] + (x == root ? 0 : cheapest);
    }

    bool found = false;
    long long best = 0;
    vector<int> parent(n, -1);
    function<void(int, long long)> choose = [&](int x, long long total){
        if (found && total + rest[x] >= best) return;
        if (x == n){
            for (int s = 0; s < n; s++){
                int y = s;
                for (int steps = 0; y != root && steps <= n; steps++) y = parent[y];
                if (y != root) return;
            }
            found = true, best = total;
            return;
        }
        if (x == root) return choose(x + 1, total);
        for (int id : incoming[x]){
            parent[x] = edges[id].u;
            choose(x + 1, total + edges[id].w);
        }
    };

    choose(0, 0);
    return {found, best};
}

/// Checks res.parent is an arborescence from root made of given edges and returns its cheapest weight, O(n + m)
long long tree_weight(int n, int root, const vector<Edge>& edges, const Arborescence& res){
    assert((int)res.parent.size() == n && res.parent[root] == -1);

    const long long none = LLONG_MAX;
    vector<long long> cheapest(n, none);
    for (auto& e : edges){
        if (e.v != root && e.u == res.parent[e.v]) cheapest[e.v] = min(cheapest[e.v], e.w);
    }

    vector<int> state(n, 0), chain;
    state[root] = 2;
    for (int s = 0; s < n; s++){
        int y = s;
        while (state[y] == 0) state[y] = 1, chain.push_back(y), y = res.parent[y];
        assert(state[y] == 2);
        for (int x : chain) state[x] = 2;
        chain.clear();
    }

    long long total = 0;
    for (int v = 0; v < n; v++){
        if (v == root) continue;
        assert(cheapest[v] != none);
        total += cheapest[v];
    }
    return total;
}

/// Cheap 2-cycles i <-> i + 1 behind expensive root edges, each contraction forms a new 2-cycle with the next node,
/// so round-based Chu-Liu/Edmonds needs about n / 2 full passes (11 s at -O2 for n = 2e4) while O(m log m) is instant;
/// the 5 s bound leaves headroom for slow CI runners and still separates the two
void nested_two_cycles(){
    int n = 20000;
    DirectedMST g(n);
    for (int i = 1; i < n; i++) g.add_edge(0, i, 1000000 + i);
    for (int i = 1; i + 1 < n; i++) g.add_edge(i + 1, i, 1), g.add_edge(i, i + 1, 2);

    auto start = chrono::steady_clock::now();
    Arborescence res = g.solve(0);
    assert(chrono::steady_clock::now() - start < chrono::seconds(5));

    /// Any single root edge 0 -> j completes the tree with j - 1 edges of weight 1 and n - 1 - j of weight 2
    assert(res.exists && res.weight == 1000000 + 2LL * n - 3);
    assert(tree_weight(n, 0, g.edges, res) == res.weight);
}

/// Larger random graphs have no brute force, but the parents must still form a tree that costs exactly the reported weight
void large_random(){
    for (long long it = 0; it < stress::scaled(60); it++){
        int n = stress::rand_int(1, 3000), root = stress::rand_int(0, n - 1), m = stress::rand_int(0, 4 * n);
        long long max_w = it % 2 ? 5 : 1000000000;
        DirectedMST g(n);
        for (int x = 0; x < n; x++){
            if (x != root && it % 3) g.add_edge(stress::rand_int(0, n - 1), x, stress::rand_int(-max_w, max_w));
        }
        for (int i = 0; i < m; i++) g.add_edge(stress::rand_int(0, n - 1), stress::rand_int(0, n - 1), stress::rand_int(-max_w, max_w));

        Arborescence res = g.solve(root);
        if (res.exists) assert(tree_weight(n, root, g.edges, res) == res.weight);
        else assert(res.weight == 0 && res.parent.empty());
    }
}

int main(){
    nested_two_cycles();
    large_random();

    for (long long it = 0; it < stress::scaled(40000); it++){
        int n = stress::rand_int(1, 7), root = stress::rand_int(0, n - 1), mode = it % 4;
        /// mode 3 sits on the header's n * max|w| <= 1e18 boundary
        long long max_w = mode == 3 ? 1000000000000000000LL / n : mode == 2 ? 1000000000 : 10;
        long long min_w = mode == 0 ? 0 : -max_w;
        vector<Edge> edges(stress::rand_int(0, 15));
        for (auto& e : edges) e = Edge(stress::rand_int(0, n - 1), stress::rand_int(0, n - 1), stress::rand_int(min_w, max_w));
        if (it % 2){
            /// Plant an arborescence so most instances have an answer and cycles must be contracted
            for (int x = 0; x < n; x++){
                if (x != root) edges.push_back(Edge(x ? stress::rand_int(0, x - 1) : root, x, stress::rand_int(min_w, max_w)));
            }
            shuffle(edges.begin(), edges.end(), stress::rng());
        }

        DirectedMST g(n);
        for (auto& e : edges) g.add_edge(e.u, e.v, e.w);
        Arborescence res = g.solve(root);
        auto [found, best] = brute(n, root, edges);
        assert(res.exists == found);
        if (found) assert(res.weight == best && tree_weight(n, root, edges, res) == best);
        else assert(res.weight == 0 && res.parent.empty());
    }
    return 0;
}
