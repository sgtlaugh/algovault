/***
 *
 * Offline Dynamic Connectivity
 * Answers connectivity and component count queries on an undirected graph under edge insertions and deletions
 *
 * Complexity: O(n + (m + q) log(m + q) log n) time, O(n + q + m log q) memory, for m edge operations and q queries
 *
 * DynamicConnectivity dc(n): an empty graph on vertices 0 to n - 1, parallel edges and self loops allowed
 * dc.add_edge(u, v): adds a copy of the edge (u, v)
 * dc.query_components(): records a query for the number of connected components, returns its index
 * dc.query_connected(u, v): records a query for 1 if u and v are connected, else 0, returns its index
 * dc.remove_edge(u, v): deletes one live copy of (u, v) in either orientation, the copy must exist
 * dc.solve(): the answers in query order, indexed by what the query calls returned
 *
 * Notes:
 *   - solve() keeps the operations, so more can be added and solve() called again
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// BEGIN COPY rollback_dsu from code_library/data_structures/rollback_dsu.cpp
struct RollbackDSU{
    vector<int> sz, parent, history;

    RollbackDSU(int n) : sz(n + 1, 1), parent(n + 1){
        iota(parent.begin(), parent.end(), 0);
    }

    bool connect(int a, int b){
        a = find(a), b = find(b);
        if (a == b) return false;

        if (sz[a] > sz[b]) swap(a, b);
        parent[a] = b, sz[b] += sz[a];
        history.push_back(a);
        return true;
    }

    void rollback(int t){
        assert(0 <= t && t <= time());
        while ((int)history.size() > t){
            int a = history.back();
            history.pop_back();
            sz[parent[a]] -= sz[a], parent[a] = a;
        }
    }

    int find(int i){
        while (i != parent[i]) i = parent[i];
        return i;
    }

    bool same(int a, int b){
        return find(a) == find(b);
    }

    int size(int i){
        return sz[find(i)];
    }

    int time(){
        return history.size();
    }
};
/// END COPY rollback_dsu

struct DynamicConnectivity{
    int n;
    vector<array<int, 2>> queries; /// (u, v), or (-1, -1) for a component count
    vector<array<int, 4>> intervals; /// edge (u, v) alive for queries [l, r)
    map<pair<int, int>, vector<int>> open; /// first query index of each live copy of an edge
    vector<vector<pair<int, int>>> tree;
    RollbackDSU dsu;

    DynamicConnectivity(int n): n(n), dsu(n) {}

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

        dfs(1, 0, q - 1, res);

        tree.clear();
        return res;
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

    void dfs(int node, int lo, int hi, vector<int>& res){
        int saved = dsu.time();
        for (auto [u, v]: tree[node]) dsu.connect(u, v);

        if (lo == hi){
            auto [u, v] = queries[lo];
            res[lo] = u == -1 ? n - dsu.time() : dsu.same(u, v);
        }
        else{
            int mid = (lo + hi) / 2;
            dfs(2 * node, lo, mid, res);
            dfs(2 * node + 1, mid + 1, hi, res);
        }

        dsu.rollback(saved);
    }
};

int main(){
    DynamicConnectivity dc(4);
    int before = dc.query_components();
    dc.add_edge(0, 1), dc.add_edge(2, 1);
    int linked = dc.query_connected(0, 2), apart = dc.query_connected(0, 3);
    dc.remove_edge(1, 0);                           /// either orientation removes the edge
    int cut = dc.query_connected(0, 2), after = dc.query_components();

    /// Parallel copies: removing one leaves the other live
    dc.add_edge(0, 1), dc.add_edge(1, 0);
    dc.remove_edge(0, 1);
    int doubled = dc.query_connected(0, 2);

    auto res = dc.solve();
    assert(res[before] == 4);
    assert(res[linked] == 1);
    assert(res[apart] == 0);
    assert(res[cut] == 0);
    assert(res[after] == 3);                        /// {0}, {1, 2}, {3}
    assert(res[doubled] == 1);
    return 0;
}
