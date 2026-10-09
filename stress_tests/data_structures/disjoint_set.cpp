#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/disjoint_set.cpp"
#undef main

/// Union by size keeps every node within log2(component size) parent hops of its root
template <typename Dsu>
void check_depths(const Dsu& dsu){
    for (int i = 0; i < (int)dsu.parent.size(); i++){
        int depth = 0, root = i;
        while (dsu.parent[root] != root) root = dsu.parent[root], depth++;
        assert(depth < 31 && (1 << depth) <= dsu.counter[root]);
    }
}

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

/// Random unions and queries against a naive label array, merged by relabeling
void stress_dsu(){
    for (int n : {1, 2, 1000, 100000}){
        DSU chain(n);
        for (int i = 0; i < n; i++) assert(chain.connect(i, i + 1));  /// links the grown root under a fresh node each time
        check_depths(chain);
        assert(chain.component_size(0) == n + 1 && !chain.connect(0, n));
    }

    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(0, it % 10 ? 15 : 400);
        DSU dsu(n);
        vector<int> label = naive_labels(n, {});

        for (int q = stress::rand_int(0, 4 * n + 5); q; q--){
            int a = stress::rand_int(0, n), b = stress::rand_int(0, n);
            if (stress::rand_int(0, 2) == 0){
                assert(dsu.connect(a, b) == (label[a] != label[b]));
                relabel(label, a, b);
            }
            assert(dsu.is_connected(a, b) == (label[a] == label[b]));
            assert(dsu.component_size(a) == count(label.begin(), label.end(), label[a]));
        }
        check_depths(dsu);
    }
}

/// Connects, checkpoints and rollbacks to random checkpoints, against labels recomputed from the surviving unions
void stress_rollback(){
    for (int n : {1, 2, 1000, 100000}){
        RollbackDSU chain(n);
        for (int i = 0; i < n; i++) assert(chain.connect(i, i + 1));
        check_depths(chain);
        assert(chain.component_size(n) == n + 1);

        chain.rollback(n / 2);
        assert(chain.time() == n / 2 && chain.component_size(0) == n / 2 + 1 && chain.component_size(n) == 1);
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

            vector<int> cnt = naive_sizes(label);
            assert(dsu.time() == (int)unions.size());
            for (int i = 0; i <= n; i++){
                int j = stress::rand_int(0, n);
                assert(dsu.is_connected(i, j) == (label[i] == label[j]));
                assert(dsu.component_size(i) == cnt[label[i]]);
            }
        }
        check_depths(dsu);
    }
}

/// Relations generated from a hidden assignment, some perturbed into contradictions once the pair is connected
void stress_weighted(){
    for (int n : {1000, 100000}){
        WeightedDSU chain(n);
        for (int i = 0; i < n; i++) assert(chain.connect(i, i + 1, 1));
        check_depths(chain);
        assert(chain.diff(0, n) == n && chain.diff(n, 0) == -n && !chain.connect(n, 0, n));
    }

    /// Pairwise merges of equal blocks never compress a path, so node i sits popcount(i) hops deep
    /// The first diff from the deepest node sums 17 potentials near 1e18, which overflows unless reduced per hop
    {
        const int n = (1 << 17) - 1;
        const long long mod = 1e18;
        vector<long long> x(n + 1);
        for (auto& v : x) v = stress::rand_int(0, mod - 1);

        WeightedDSU binomial(n, mod);
        for (int step = 1; step <= n; step *= 2){
            for (int i = 0; i + step <= n; i += 2 * step) assert(binomial.connect(i, i + step, x[i + step] - x[i]));
        }
        check_depths(binomial);

        for (int i = n; i >= 0; i--) assert(binomial.diff(i, 0) == ((x[0] - x[i]) % mod + mod) % mod);
    }

    const long long lo = LLONG_MIN / 2, hi = LLONG_MAX / 2;
    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(0, it % 10 ? 15 : 300);
        long long mod = vector<long long>{0, 0, 1, 2, 3, stress::rand_int(4, 1e9), (long long)1e18}[it % 7];
        bool extreme = mod == 0 && it % 2;

        vector<long long> x(n + 1);
        for (auto& v : x) v = mod ? stress::rand_int(0, mod - 1) : extreme ? stress::rand_int(lo, hi) : stress::rand_int(-1e18, 1e18);
        auto truth = [&](int a, int b){
            if (!mod) return x[b] - x[a];
            return ((x[b] - x[a]) % mod + mod) % mod;
        };

        WeightedDSU dsu(n, mod);
        vector<int> label = naive_labels(n, {});

        for (int q = stress::rand_int(0, 4 * n + 5); q; q--){
            int a = stress::rand_int(0, n), b = stress::rand_int(0, n);
            long long d = x[b] - x[a];
            bool lie = label[a] == label[b] && mod != 1 && stress::rand_int(0, 2) == 0;

            if (lie && mod) d += stress::rand_int(1, mod - 1);
            if (lie && !mod) d += d >= 0 ? -stress::rand_int(1, 5) : stress::rand_int(1, 5);
            if (mod) d += mod * stress::rand_int(-3, 3);

            assert(dsu.connect(a, b, d) == !lie);
            relabel(label, a, b);

            int i = stress::rand_int(0, n), j = stress::rand_int(0, n);
            assert(dsu.is_connected(i, j) == (label[i] == label[j]));
            assert(dsu.component_size(i) == count(label.begin(), label.end(), label[i]));
            if (label[i] == label[j]) assert(dsu.diff(i, j) == truth(i, j));
        }

        for (int i = 0; i <= n; i++){
            for (int j = 0; j <= n; j++) if (label[i] == label[j]) assert(dsu.diff(i, j) == truth(i, j));
        }
        check_depths(dsu);
    }
}

/// Builds the whole edge history, then checks every time t against labels rebuilt from the first t edges
void stress_persistent(){
    for (int n : {1, 2, 1000, 100000}){
        PartiallyPersistentDSU chain(n);
        for (int i = 0; i < n; i++) assert(chain.connect(i, i + 1));
        assert(chain.component_size(0, n) == n + 1 && chain.component_size(n, n / 2) == 1);
        assert(chain.component_size(0, n / 2) == n / 2 + 1 && !chain.is_connected(0, n, n - 1));
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
                int j = stress::rand_int(0, n), root = dsu.find_root(i, t);
                assert(label[root] == label[i] && dsu.find_root(root, t) == root);
                assert(dsu.is_connected(i, j, t) == (label[i] == label[j]));
                assert(dsu.component_size(i, t) == cnt[label[i]]);
            }
        }

        for (int i = 0; i <= n; i++){
            int depth = 0, root = i;
            while (dsu.joined[root] <= m) root = dsu.parent[root], depth++;
            assert(depth < 31 && (1 << depth) <= dsu.component_size(i, m));
        }
    }
}

int main(){
    stress_dsu();
    stress_rollback();
    stress_weighted();
    stress_persistent();

    return 0;
}
