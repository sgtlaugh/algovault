#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/partially_persistent_dsu.cpp"
#undef main

/// Applies one union to a label array in O(n)
void relabel(vector<int>& label, int a, int b){
    int from = label[a], to = label[b];
    if (from != to) for (auto& x : label) if (x == from) x = to;
}

/// Component labels recomputed from scratch over the edge list
vector<int> naive_labels(int n, const vector<pair<int, int>>& edges){
    vector<int> label(n + 1);
    iota(label.begin(), label.end(), 0);
    for (auto [a, b] : edges) relabel(label, a, b);
    return label;
}

/// Component sizes indexed by label
vector<int> naive_sizes(const vector<int>& label){
    vector<int> cnt(label.size());
    for (int x : label) cnt[x]++;
    return cnt;
}

/// Builds the whole edge history, then checks every time t against labels rebuilt from the first t edges
int main(){
    for (int n : {1, 2, 1000, 100000}){
        PartiallyPersistentDSU chain(n);
        for (int i = 0; i < n; i++) assert(chain.connect(i, i + 1));
        assert(chain.size(0, n) == n + 1 && chain.size(n, n / 2) == 1);
        assert(chain.size(0, n / 2) == n / 2 + 1 && !chain.same(0, n, n - 1));
    }

    for (long long it = 0; it < stress::scaled(1500); it++){
        int n = stress::rand_int(0, it % 10 ? 12 : 150), m = stress::rand_int(0, 3 * n + 3);
        PartiallyPersistentDSU dsu(n);
        vector<pair<int, int>> edges;
        vector<int> label = naive_labels(n, {});

        for (int k = 0; k < m; k++){
            int a = stress::rand_int(0, n), b = stress::rand_int(0, n);
            assert(dsu.connect(a, b) == (label[a] != label[b]));
            edges.push_back({a, b});
            relabel(label, a, b);
        }
        assert(dsu.now == m);

        vector<int> times(m + 3);
        iota(times.begin(), times.end(), 0);
        times.push_back(INT_MAX);

        label = naive_labels(n, {});
        for (int t : times){
            if (t >= 1 && t <= m) relabel(label, edges[t - 1].first, edges[t - 1].second);

            vector<int> cnt = naive_sizes(label);
            for (int i = 0; i <= n; i++){
                int j = stress::rand_int(0, n), root = dsu.find(i, t);
                assert(label[root] == label[i] && dsu.find(root, t) == root);
                assert(dsu.same(i, j, t) == (label[i] == label[j]));
                assert(dsu.size(i, t) == cnt[label[i]]);
            }
        }

        for (int i = 0; i <= n; i++){
            int depth = 0, root = i;
            while (dsu.joined[root] <= m) root = dsu.parent[root], depth++;
            assert(depth < 31 && (1 << depth) <= dsu.size(i, m));
        }
    }

    return 0;
}
