#include "common.h"

#define main library_main
#include "../code_library/directed_mst.cpp"
#undef main

/// Every choice of one incoming edge per non-root node, kept when all parent chains reach the root
long long brute(int n, int root, const vector<Edge>& edges){
    vector<vector<int>> incoming(n);
    for (auto& e : edges) if (e.u != e.v) incoming[e.v].push_back(&e - &edges[0]);

    long long best = -1;
    vector<int> parent(n, -1);
    function<void(int, long long)> choose = [&](int x, long long total){
        if (best != -1 && total >= best) return;
        if (x == n){
            for (int s = 0; s < n; s++){
                int y = s;
                for (int steps = 0; y != root && steps <= n; steps++) y = parent[y];
                if (y != root) return;
            }
            best = total;
            return;
        }
        if (x == root) return choose(x + 1, total);
        for (int id : incoming[x]){
            parent[x] = edges[id].u;
            choose(x + 1, total + edges[id].w);
        }
    };
    choose(0, 0);
    return best;
}

int main(){
    for (long long it = 0; it < stress::scaled(40000); it++){
        int n = stress::rand_int(1, 7), root = stress::rand_int(0, n - 1), max_w = it % 4 ? 10 : 1000000000;
        vector<Edge> edges(stress::rand_int(0, 15));
        for (auto& e : edges) e = Edge(stress::rand_int(0, n - 1), stress::rand_int(0, n - 1), stress::rand_int(0, max_w));
        if (it % 2){
            /// Plant an arborescence so most instances have an answer and cycles must be contracted
            for (int x = 0; x < n; x++){
                if (x != root) edges.push_back(Edge(x ? stress::rand_int(0, x - 1) : root, x, stress::rand_int(0, max_w)));
            }
            shuffle(edges.begin(), edges.end(), stress::rng());
        }
        assert(directed_mst(n, root, edges) == brute(n, root, edges));
    }
    return 0;
}
