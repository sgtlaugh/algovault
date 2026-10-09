#include "../common.h"

#define main library_main
#include "../../code_library/graphs/offline_dynamic_connectivity.cpp"
#undef main

/// BFS component labels over the live edge copies
vector<int> bfs_labels(int n, const vector<pair<int, int>>& live){
    vector<vector<int>> adj(n);
    for (auto [u, v]: live){
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    vector<int> label(n, -1);
    for (int s = 0; s < n; s++){
        if (label[s] != -1) continue;
        queue<int> bfs;
        bfs.push(s), label[s] = s;
        while (!bfs.empty()){
            int u = bfs.front();
            bfs.pop();
            for (int v: adj[u]){
                if (label[v] == -1) label[v] = s, bfs.push(v);
            }
        }
    }
    return label;
}

/// Random operations, every query whose index check_every divides is compared against a BFS recount
/// On reuse runs solve() is also called halfway, so the final solve() checks it consumed nothing
void run(int n, int ops, int check_every, bool reuse){
    DynamicConnectivity dc(n);
    vector<pair<int, int>> live;
    vector<int> expected;
    auto check = [&](){
        auto res = dc.solve();
        assert(res.size() == expected.size());
        for (size_t i = 0; i < res.size(); i++) assert(expected[i] == -1 || res[i] == expected[i]);
    };

    for (int i = 0; i < ops; i++){
        if (reuse && i == ops / 2) check();

        int type = n == 0 ? 4 : stress::rand_int(0, 5);

        if (type <= 1){
            /// Mostly short hops so long paths and cycles form, sometimes any pair including self loops
            int u = stress::rand_int(0, n - 1);
            int v = stress::rand_int(0, 2) ? (u + stress::rand_int(1, 3)) % n : stress::rand_int(0, n - 1);
            dc.add_edge(u, v);
            live.push_back({u, v});
        }
        else if (type == 2 && !live.empty()){
            int k = stress::rand_int(0, live.size() - 1);
            auto [u, v] = live[k];
            swap(live[k], live.back());
            live.pop_back();
            if (stress::rand_int(0, 1)) swap(u, v);
            dc.remove_edge(u, v);
        }
        else{
            bool counting = type == 4;
            int u = counting ? -1 : stress::rand_int(0, n - 1), v = counting ? -1 : stress::rand_int(0, n - 1);
            int idx = counting ? dc.query_components() : dc.query_connected(u, v);
            assert(idx == (int)expected.size());

            if (idx % check_every){
                expected.push_back(-1);
                continue;
            }
            auto label = bfs_labels(n, live);
            if (!counting) expected.push_back(label[u] == label[v]);
            else{
                int roots = 0;
                for (int x = 0; x < n; x++) roots += label[x] == x;
                expected.push_back(roots);
            }
        }
    }

    check();
}

int main(){
    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = it < 200 ? it % 7 : stress::rand_int(1, it % 20 ? 10 : 120);
        run(n, stress::rand_int(0, it % 20 ? 40 : 600), 1, it % 3 == 0);
    }

    /// A path linked in order hangs each old root under the new node without union by size, making every find O(n)
    int len = 1 << 16;
    DynamicConnectivity path(len);
    for (int i = 0; i + 1 < len; i++) path.add_edge(i, i + 1);
    for (int i = 0; i < len; i++) path.query_connected(0, len - 1);

    auto start = chrono::steady_clock::now();
    assert(path.solve() == vector<int>(len, 1));
    assert(chrono::steady_clock::now() - start < chrono::seconds(2));

    /// Two permanent sets of size half + 1 and half are joined in every other leaf and rolled back each time
    /// A rollback that leaves the larger root's size inflated grows it by half per leaf, past INT_MAX here
    /// Within int range stale sizes still double per level of depth, so only the overflow is observable
    int half = 1 << 16;
    DynamicConnectivity sets(2 * half + 1);
    for (int i = 1; i <= half; i++) sets.add_edge(0, i), sets.add_edge(half + 1, half + i);
    vector<int> joined;
    for (int i = 0; i < half; i++){
        sets.add_edge(0, half + 1);
        sets.query_connected(1, 2 * half);
        sets.remove_edge(half + 1, 0);
        sets.query_components();
        joined.insert(joined.end(), {1, 2});
    }
    assert(sets.solve() == joined);

    for (int n: {1000, 100000}){
        auto start = chrono::steady_clock::now();
        run(n, 100000, n / 20, n == 1000);
        assert(chrono::steady_clock::now() - start < chrono::seconds(10));
    }

    return 0;
}
