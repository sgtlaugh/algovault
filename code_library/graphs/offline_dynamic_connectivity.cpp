/***
 *
 * Offline Dynamic Connectivity
 * Answers connectivity and component count queries on an undirected graph under edge insertions and deletions
 *
 * Complexity: O(n + (m + q) log(m + q) log n) time, O(n + q + m log q) memory, for m edge operations and q queries
 *
 * Every edge copy is alive for a contiguous range of queries, so it is inserted into the O(log q) segment tree
 * nodes covering that range; a DFS over the tree unites a node's edges on entry and rolls them back on exit,
 * using a DSU with union by size and no path compression so each union is undone in O(1)
 *
 * Nodes are numbered from 0 to n - 1, parallel edges and self loops are allowed
 * remove_edge deletes one live copy of (u, v) in either orientation, the copy must exist
 * Queries return their index, solve() returns the answers in query order:
 *   query_connected(u, v) -> 1 if u and v are connected, else 0
 *   query_components()    -> number of connected components
 * solve() does not consume the operations, more can be added and solve() called again
 *
 * Example:
 *   DynamicConnectivity dc(3);
 *   dc.add_edge(0, 1);
 *   int a = dc.query_connected(0, 1), b = dc.query_components();
 *   dc.remove_edge(1, 0);
 *   int c = dc.query_components();
 *   auto res = dc.solve();  // res[a] = 1, res[b] = 2, res[c] = 3
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct DynamicConnectivity{
    int n, components;
    vector<array<int, 2>> queries; /// (u, v), or (-1, -1) for a component count
    vector<array<int, 4>> intervals; /// edge (u, v) alive for queries [l, r)
    map<pair<int, int>, vector<int>> open; /// first query index of each live copy of an edge
    vector<vector<pair<int, int>>> tree;
    vector<int> parent, set_size, history;

    DynamicConnectivity(int n): n(n) {}

    void add_edge(int u, int v){
        if (u > v) swap(u, v);
        open[{u, v}].push_back(queries.size());
    }

    int query_components(){
        queries.push_back({-1, -1});
        return queries.size() - 1;
    }

    int query_connected(int u, int v){
        queries.push_back({u, v});
        return queries.size() - 1;
    }

    void remove_edge(int u, int v){
        if (u > v) swap(u, v);
        auto it = open.find({u, v});
        assert(it != open.end());

        int l = it->second.back();
        it->second.pop_back();
        if (it->second.empty()) open.erase(it);
        if (l < (int)queries.size()) intervals.push_back({u, v, l, (int)queries.size()});
    }

    vector<int> solve(){
        int q = queries.size();
        vector<int> res(q);
        if (q == 0) return res;

        tree.assign(4 * q, {});
        for (auto [u, v, l, r]: intervals) insert(1, 0, q - 1, l, r - 1, {u, v});
        for (auto& [edge, starts]: open){
            for (int l: starts){
                if (l < q) insert(1, 0, q - 1, l, q - 1, edge);
            }
        }

        components = n;
        parent.resize(n);
        iota(parent.begin(), parent.end(), 0);
        set_size.assign(n, 1);
        history.clear();
        dfs(1, 0, q - 1, res);

        tree.clear();
        return res;
    }

private:
    void dfs(int node, int lo, int hi, vector<int>& res){
        int saved = history.size();
        for (auto [u, v]: tree[node]) unite(u, v);

        if (lo == hi){
            auto [u, v] = queries[lo];
            res[lo] = u == -1 ? components : find_root(u) == find_root(v);
        }
        else{
            int mid = (lo + hi) / 2;
            dfs(2 * node, lo, mid, res);
            dfs(2 * node + 1, mid + 1, hi, res);
        }

        while ((int)history.size() > saved){
            int c = history.back();
            history.pop_back();
            set_size[parent[c]] -= set_size[c];
            parent[c] = c;
            components++;
        }
    }

    int find_root(int u){
        while (parent[u] != u) u = parent[u];
        return u;
    }

    void insert(int node, int lo, int hi, int l, int r, pair<int, int> edge){
        if (r < lo || hi < l) return;
        if (l <= lo && hi <= r){
            tree[node].push_back(edge);
            return;
        }

        int mid = (lo + hi) / 2;
        insert(2 * node, lo, mid, l, r, edge);
        insert(2 * node + 1, mid + 1, hi, l, r, edge);
    }

    void unite(int u, int v){
        u = find_root(u), v = find_root(v);
        if (u == v) return;
        if (set_size[u] > set_size[v]) swap(u, v);

        parent[u] = v;
        set_size[v] += set_size[u];
        history.push_back(u);
        components--;
    }
};

int main(){
    DynamicConnectivity dc(4);
    int q0 = dc.query_components();
    dc.add_edge(0, 1);
    int q1 = dc.query_connected(1, 0);
    dc.add_edge(2, 1);
    int q2 = dc.query_components(), q3 = dc.query_connected(0, 2), q4 = dc.query_connected(0, 3);
    dc.remove_edge(1, 0);
    int q5 = dc.query_connected(0, 2), q6 = dc.query_components();
    dc.add_edge(0, 1);
    dc.add_edge(1, 0);
    dc.remove_edge(0, 1);
    int q7 = dc.query_connected(0, 2), q8 = dc.query_components();
    dc.add_edge(3, 3);
    dc.remove_edge(1, 2);
    dc.add_edge(2, 3);
    dc.remove_edge(3, 2);
    int q9 = dc.query_components(), q10 = dc.query_connected(3, 3), q11 = dc.query_connected(1, 2);

    auto res = dc.solve();
    assert((res == vector<int>{4, 1, 2, 1, 0, 0, 3, 1, 2, 3, 1, 0}));
    assert(q0 == 0 && q1 == 1 && q2 == 2 && q3 == 3 && q4 == 4 && q5 == 5);
    assert(q6 == 6 && q7 == 7 && q8 == 8 && q9 == 9 && q10 == 10 && q11 == 11);

    dc.add_edge(2, 3);
    dc.query_connected(1, 3);
    dc.query_components();
    assert((dc.solve() == vector<int>{4, 1, 2, 1, 0, 0, 3, 1, 2, 3, 1, 0, 0, 2}));

    DynamicConnectivity empty(0);
    assert(empty.solve().empty());

    DynamicConnectivity single(1);
    single.add_edge(0, 0);
    single.query_connected(0, 0);
    single.query_components();
    assert((single.solve() == vector<int>{1, 1}));

    DynamicConnectivity no_queries(3);
    no_queries.add_edge(0, 1);
    no_queries.remove_edge(0, 1);
    assert(no_queries.solve().empty());

    return 0;
}
