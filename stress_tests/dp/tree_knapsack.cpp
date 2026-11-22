#include "../common.h"

#define main library_main
#include "../../code_library/dp/tree_knapsack.cpp"
#undef main

/// Every subset, kept when closed upward or connected (members minus one equals the member-parent edges)
pair<vector<long long>, vector<long long>> brute(const vector<int>& parent, const vector<int>& weight,
                                                 const vector<long long>& value, int W){
    int n = parent.size();
    vector<long long> closed(W + 1, 0), connected(W + 1, 0);
    for (int mask = 1; mask < (1 << n); mask++){
        long long w = 0, v = 0;
        int members = 0, edges = 0;
        bool is_closed = true;
        for (int i = 0; i < n; i++){
            if (!(mask >> i & 1)) continue;
            w += weight[i], v += value[i], members++;
            if (parent[i] != -1 && (mask >> parent[i] & 1)) edges++;
            else if (parent[i] != -1) is_closed = false;
        }
        if (w > W) continue;

        for (int c = w; c <= W; c++){
            if (is_closed) closed[c] = max(closed[c], v);
            if (edges == members - 1) connected[c] = max(connected[c], v);
        }
    }

    return {closed, connected};
}

/// Children merged pairwise in O(W^2), the classic tree DP, rooted at each vertex for connected sets
vector<long long> merge_dp(const vector<vector<int>>& adj, const vector<int>& weight, const vector<long long>& value,
                           int W, int v, int from){
    const long long NONE = LLONG_MIN / 4;
    vector<long long> f(W + 1, NONE);
    if (weight[v] <= W) f[weight[v]] = value[v];
    for (int u : adj[v]){
        if (u == from) continue;
        vector<long long> g = merge_dp(adj, weight, value, W, u, v), h = f;
        for (int a = 0; a <= W; a++){
            if (f[a] == NONE) continue;
            for (int b = 0; a + b <= W; b++){
                if (g[b] != NONE) h[a + b] = max(h[a + b], f[a] + g[b]);
            }
        }
        f = h;
    }

    return f;
}

vector<long long> prefix_max(vector<long long> f){
    f[0] = max(f[0], 0LL);
    for (int c = 1; c < (int)f.size(); c++) f[c] = max(f[c], f[c - 1]);
    return f;
}

vector<int> random_forest(int n, int shape){
    vector<int> parent(n, -1);
    for (int v = 1; v < n; v++){
        if (shape == 0) parent[v] = stress::rand_int(0, 5) ? stress::rand_int(0, v - 1) : -1;
        else if (shape == 1) parent[v] = v - 1;
        else if (shape == 2) parent[v] = 0;
        else parent[v] = stress::rand_int(max(0, v - 2), v - 1);
    }

    /// Relabel so parents are not always smaller than their children
    vector<int> label(n);
    iota(label.begin(), label.end(), 0);
    shuffle(label.begin(), label.end(), stress::rng());
    vector<int> res(n);
    for (int v = 0; v < n; v++) res[label[v]] = parent[v] == -1 ? -1 : label[parent[v]];
    return res;
}

void check_small(int n, int shape, int W, int max_weight, long long max_value){
    vector<int> parent = random_forest(n, shape), weight(n);
    vector<long long> value(n);
    for (int v = 0; v < n; v++){
        weight[v] = stress::rand_int(0, max_weight);
        value[v] = stress::rand_int(-max_value / 3, max_value);
    }

    auto [closed, connected] = brute(parent, weight, value, W);
    TreeKnapsack tk(parent, weight, value, W);
    assert(tk.ancestor_closed() == closed);
    assert(tk.connected() == connected);
}

void check_medium(int n, int W){
    vector<int> parent = random_forest(n, stress::rand_int(0, 3)), weight(n);
    vector<long long> value(n);
    for (int v = 0; v < n; v++){
        weight[v] = stress::rand_int(0, 6);
        value[v] = stress::rand_int(-20, 100);
    }

    vector<vector<int>> adj(n);
    for (int v = 0; v < n; v++){
        if (parent[v] != -1) adj[v].push_back(parent[v]), adj[parent[v]].push_back(v);
    }

    /// A virtual root n of weight 0 ties the forest together for the closed upward reference
    vector<vector<int>> rooted = adj;
    rooted.push_back({});
    vector<int> weight_root = weight;
    vector<long long> value_root = value;
    weight_root.push_back(0), value_root.push_back(0);
    for (int v = 0; v < n; v++){
        if (parent[v] == -1) rooted[n].push_back(v), rooted[v].push_back(n);
    }

    vector<long long> connected(W + 1, 0);
    for (int r = 0; r < n; r++){
        vector<long long> f = prefix_max(merge_dp(adj, weight, value, W, r, -1));
        for (int c = 0; c <= W; c++) connected[c] = max(connected[c], f[c]);
    }

    TreeKnapsack tk(parent, weight, value, W);
    assert(tk.ancestor_closed() == prefix_max(merge_dp(rooted, weight_root, value_root, W, n, -1)));
    assert(tk.connected() == connected);
}

int main(){
    TreeKnapsack empty({}, {}, {}, 3);
    assert((empty.ancestor_closed() == vector<long long>{0, 0, 0, 0}));
    assert((empty.connected() == vector<long long>{0, 0, 0, 0}));

    /// Small forests of every shape against all subsets, with zero weights, W = 0 and negative values
    for (int it = 0; it < stress::scaled(1500); it++){
        int n = stress::rand_int(1, 11), max_weight = stress::rand_int(0, 4);
        check_small(n, it % 4, stress::rand_int(0, 12), max_weight, it % 5 ? 20 : (long long)1e15);
    }

    for (int it = 0; it < stress::scaled(30); it++) check_medium(stress::rand_int(20, 70), stress::rand_int(0, 40));

    /// A path of 2e5 vertices of weight 1 and value 1: closed upward and connected both take best[c] = c
    int n = 200000, W = 5;
    vector<int> parent(n), weight(n, 1);
    for (int v = 0; v < n; v++) parent[v] = v - 1;
    TreeKnapsack path(parent, weight, vector<long long>(n, 1), W);
    vector<long long> diagonal = {0, 1, 2, 3, 4, 5};
    assert(path.ancestor_closed() == diagonal);
    assert(path.connected() == diagonal);

    /// Parent cycles are unreachable from every root, so both calls must reject them
    stress::expect_assert("dp/tree_knapsack.cpp", []{
        TreeKnapsack bad({1, 0}, {1, 1}, {1, 1}, 2);
        bad.ancestor_closed();
    });
    stress::expect_assert("dp/tree_knapsack.cpp", []{
        TreeKnapsack bad({-1, 2, 1}, {1, 1, 1}, {1, 1, 1}, 2);
        bad.connected();
    });
    stress::expect_assert("dp/tree_knapsack.cpp", []{ TreeKnapsack bad({-1}, {-1}, {1}, 2); });

    return 0;
}
