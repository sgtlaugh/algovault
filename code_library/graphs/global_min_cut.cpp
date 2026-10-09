/***
 *
 * Global minimum cut with the Stoer-Wagner algorithm
 * Minimum total weight of edges to remove so that an undirected graph splits into two non-empty parts
 * 0 based indexing for nodes, so nodes are numbered from 0 to n-1
 *
 * Complexity: O(n^3) time, O(n^2) memory for the adjacency matrix
 *
 * add_edge(u, v, w) adds an undirected edge, parallel edges add up and self-loops are ignored
 * Weights must be non-negative and the sum of all weights must fit in a long long
 * min_cut() needs n >= 2 and returns {weight, side}: side holds the nodes of one part, the rest form the other
 * A disconnected graph has cut weight 0 and side is a union of components
 * min_cut() works on a copy of the matrix, so edges can be added and the cut queried again
 *
 * Example:
 *   GlobalMinCut g(4);
 *   g.add_edge(0, 1, 5), g.add_edge(1, 2, 1), g.add_edge(2, 3, 5), g.add_edge(3, 0, 1);
 *   auto [weight, side] = g.min_cut();  // weight = 2, side = {0, 1} or {2, 3}
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct GlobalMinCut{
    int n;
    vector <vector<long long>> adj;

    GlobalMinCut(int n) : n(n), adj(n, vector<long long>(n, 0)) {}

    void add_edge(int u, int v, long long w){
        assert(w >= 0);
        if (u == v) return;
        adj[u][v] += w, adj[v][u] += w;
    }

    /// Each phase orders the live nodes by maximum adjacency, cuts the last one off, then merges it into the one before it
    pair<long long, vector<int>> min_cut() const{
        assert(n >= 2);
        auto mat = adj;
        vector <vector<int>> group(n);
        vector <bool> merged(n, false);
        for (int i = 0; i < n; i++) group[i] = {i};
        long long best = LLONG_MAX;
        vector <int> side;

        for (int phase = 1; phase < n; phase++){
            vector <long long> key(n, 0);
            vector <bool> added(n, false);
            int prev = -1, last = -1;

            for (int it = 0; it <= n - phase; it++){
                int sel = -1;
                for (int i = 0; i < n; i++){
                    if (!merged[i] && !added[i] && (sel == -1 || key[i] > key[sel])) sel = i;
                }
                if (it == n - phase && key[sel] < best) best = key[sel], side = group[sel];
                added[sel] = true, prev = last, last = sel;
                for (int i = 0; i < n; i++) key[i] += mat[sel][i];
            }

            group[prev].insert(group[prev].end(), group[last].begin(), group[last].end());
            for (int i = 0; i < n; i++) mat[prev][i] += mat[last][i], mat[i][prev] = mat[prev][i];
            mat[prev][prev] = 0, merged[last] = true;
        }

        sort(side.begin(), side.end());
        return {best, side};
    }
};

int main(){
    auto is_side = [](const vector<int>& side, const vector<int>& expected, int n){
        vector <int> complement;
        for (int i = 0; i < n; i++){
            if (!binary_search(expected.begin(), expected.end(), i)) complement.push_back(i);
        }
        return side == expected || side == complement;
    };

    {
        GlobalMinCut g(8);
        vector <array<int, 3>> edges = {{0, 1, 2}, {0, 4, 3}, {1, 2, 3}, {1, 4, 2}, {1, 5, 2}, {2, 3, 4}, {2, 6, 2}, {3, 6, 2}, {3, 7, 2}, {4, 5, 3}, {5, 6, 1}, {6, 7, 3}};
        for (auto [u, v, w]: edges) g.add_edge(u, v, w);
        auto [weight, side] = g.min_cut();
        assert(weight == 4 && is_side(side, {2, 3, 6, 7}, 8));
    }

    {
        GlobalMinCut g(4);
        g.add_edge(0, 1, 5), g.add_edge(1, 2, 1), g.add_edge(2, 3, 5), g.add_edge(3, 0, 1);
        auto [weight, side] = g.min_cut();
        assert(weight == 2 && is_side(side, {0, 1}, 4));
    }

    {
        GlobalMinCut g(2);
        assert(g.min_cut().first == 0);
        g.add_edge(0, 1, 2), g.add_edge(1, 0, 3), g.add_edge(1, 1, 100);
        auto [weight, side] = g.min_cut();
        assert(weight == 5 && is_side(side, {0}, 2));
    }

    {
        GlobalMinCut g(5);
        g.add_edge(0, 1, 7), g.add_edge(1, 2, 7), g.add_edge(3, 4, 9);
        auto [weight, side] = g.min_cut();
        assert(weight == 0 && is_side(side, {3, 4}, 5));
    }

    {
        GlobalMinCut g(3);
        g.add_edge(0, 1, 4000000000000000000LL), g.add_edge(1, 2, 2000000000000000000LL), g.add_edge(0, 2, 1000000000000000000LL);
        auto [weight, side] = g.min_cut();
        assert(weight == 3000000000000000000LL && is_side(side, {2}, 3));
        g.add_edge(0, 2, 2200000000000000000LL);
        tie(weight, side) = g.min_cut();
        assert(weight == 5200000000000000000LL && is_side(side, {2}, 3));
    }

    return 0;
}
