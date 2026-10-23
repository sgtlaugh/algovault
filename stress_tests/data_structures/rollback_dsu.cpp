#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/rollback_dsu.cpp"
#undef main

/// Union by size keeps every node within log2(component size) parent hops of its root
template <typename Dsu>
void check_depths(const Dsu& dsu){
    for (int i = 0; i < (int)dsu.parent.size(); i++){
        int depth = 0, root = i;
        while (dsu.parent[root] != root) root = dsu.parent[root], depth++;
        assert(depth < 31 && (1 << depth) <= dsu.sz[root]);
    }
}

void relabel(vector<int>& label, int a, int b){
    int from = label[a], to = label[b];
    if (from != to) for (auto& x : label) if (x == from) x = to;
}

vector<int> naive_labels(int n, const vector<pair<int, int>>& edges){
    vector<int> label(n + 1);
    iota(label.begin(), label.end(), 0);
    for (auto [a, b] : edges) relabel(label, a, b);
    return label;
}

vector<int> sizes_by_label(const vector<int>& label){
    vector<int> cnt(label.size());
    for (int x : label) cnt[x]++;
    return cnt;
}

/// Connects, checkpoints and rollbacks to random checkpoints, against labels recomputed from the surviving unions
int main(){
    /// A checkpoint in the future would otherwise roll back nothing and return silently
    stress::expect_assert("data_structures/rollback_dsu.cpp", []{
        RollbackDSU dsu(3);
        dsu.connect(0, 1);
        dsu.rollback(2);
    });
    stress::expect_assert("data_structures/rollback_dsu.cpp", []{ RollbackDSU(3).rollback(-1); });

    for (int n : {1, 2, 1000, 100000}){
        RollbackDSU chain(n);
        for (int i = 0; i < n; i++) assert(chain.connect(i, i + 1));
        check_depths(chain);
        assert(chain.size(n) == n + 1);

        chain.rollback(n / 2);
        assert(chain.time() == n / 2 && chain.size(0) == n / 2 + 1 && chain.size(n) == 1);
    }

    for (long long it = 0; it < stress::scaled(1500); it++){
        int n = stress::rand_int(0, it % 10 ? 12 : 120);
        RollbackDSU dsu(n);
        vector<pair<int, int>> unions;
        vector<int> checkpoints = {0}, label = naive_labels(n, unions);

        for (int q = stress::rand_int(0, 6 * n + 5); q; q--){
            int type = stress::rand_int(0, 5);
            if (type <= 2){
                int a = stress::rand_int(0, n), b = stress::rand_int(0, n);
                bool merged = label[a] != label[b];
                assert(dsu.connect(a, b) == merged);
                if (merged) unions.push_back({a, b}), relabel(label, a, b);
            }
            else if (type == 3) checkpoints.push_back(dsu.time());
            else{
                int t = checkpoints[stress::rand_int(0, (int)checkpoints.size() - 1)];
                if (t > (int)unions.size()) continue;
                dsu.rollback(t);
                unions.resize(t);
                label = naive_labels(n, unions);
            }

            vector<int> cnt = sizes_by_label(label);
            assert(dsu.time() == (int)unions.size());
            for (int i = 0; i <= n; i++){
                int j = stress::rand_int(0, n);
                assert(dsu.same(i, j) == (label[i] == label[j]));
                assert(dsu.size(i) == cnt[label[i]]);
            }
        }
        check_depths(dsu);
    }

    return 0;
}
