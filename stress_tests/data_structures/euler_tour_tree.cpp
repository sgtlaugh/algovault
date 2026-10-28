#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/euler_tour_tree.cpp"
#undef main

struct NaiveForest{
    int n;
    vector<set<int>> adj;
    vector<long long> value;

    NaiveForest(const vector<long long>& value) : n(value.size()), adj(n), value(value) {}

    /// Vertices reached from v without crossing edge v-p, p = -1 for v's whole tree
    vector<int> side(int v, int p) const{
        vector<int> seen(n, 0), res = {v};
        seen[v] = 1;
        if (p != -1) seen[p] = 1;
        for (int i = 0; i < (int)res.size(); i++){
            for (int w : adj[res[i]]){
                if (!seen[w]) seen[w] = 1, res.push_back(w);
            }
        }

        return res;
    }

    bool connected(int u, int v) const{
        vector<int> tree = side(u, -1);
        return find(tree.begin(), tree.end(), v) != tree.end();
    }

    long long subtree_sum(int v, int p) const{
        long long res = 0;
        for (int x : side(v, p)) res += value[x];
        return res;
    }
};

/// Random operations; cut picks a real edge most of the time and link often re-adds the edge just cut
void check_random(int n, int ops, long long max_value){
    vector<long long> values(n);
    for (auto& x : values) x = stress::rand_int(-max_value, max_value);
    bool from_values = stress::rand_int(0, 1);
    EulerTourTree<long long> ett(from_values ? values : vector<long long>(n));
    if (!from_values){
        for (int u = 0; u < n; u++) ett.set(u, values[u]);
    }
    NaiveForest ref(values);
    pair<int, int> last_cut = {0, 0};

    for (int op = 0; op < ops; op++){
        int kind = stress::rand_int(0, 9), u = stress::rand_int(0, n - 1), v = stress::rand_int(0, n - 1);
        long long x = stress::rand_int(-max_value, max_value);
        bool has_edge = !ref.adj[u].empty();
        if (kind >= 4 && kind <= 6 && has_edge){
            v = *next(ref.adj[u].begin(), stress::rand_int(0, ref.adj[u].size() - 1));
        }

        if (kind <= 1){
            if (kind == 1) tie(u, v) = last_cut;
            bool linked = !ref.connected(u, v);
            if (linked) ref.adj[u].insert(v), ref.adj[v].insert(u);
            assert(ett.link(u, v) == linked);
        }
        else if (kind <= 3){
            if (has_edge && kind == 3) v = *ref.adj[u].begin();
            bool exists = ref.adj[u].count(v);
            if (exists) ref.adj[u].erase(v), ref.adj[v].erase(u), last_cut = {u, v};
            assert(ett.cut(u, v) == exists);
        }
        else if (kind == 4 && has_edge){
            for (int w : ref.side(u, v)) ref.value[w] += x;
            ett.subtree_add(u, v, x);
        }
        else if (kind <= 6 && has_edge) assert(ett.subtree_query(u, v) == ref.subtree_sum(u, v));
        else if (kind == 7) assert(ett.connected(u, v) == ref.connected(u, v));
        else if (kind == 8){
            if (x & 1) ett.add(u, x), ref.value[u] += x;
            else ett.set(u, x), ref.value[u] = x;
        }
        else assert(ett.get(u) == ref.value[u]);
    }

    for (int u = 0; u < n; u++){
        assert(ett.get(u) == ref.value[u]);
        for (int w : ref.adj[u]) assert(ett.subtree_query(u, w) == ref.subtree_sum(u, w));
    }
}

/// A path 0-1-...-n-1 built in order: the deepest tour, checked against prefix sums after a range add
void check_long_path(int n){
    EulerTourTree<long long> ett(n);
    vector<long long> value(n), prefix(n + 1, 0);
    for (int i = 0; i < n; i++){
        value[i] = stress::rand_int(0, 1000000000);
        ett.set(i, value[i]);
    }
    for (int i = 1; i < n; i++) assert(ett.link(i - 1, i));

    int mid = n / 2;
    ett.subtree_add(mid, mid - 1, 7);
    for (int i = mid; i < n; i++) value[i] += 7;
    for (int i = 0; i < n; i++) prefix[i + 1] = prefix[i] + value[i];

    for (int i = 1; i < n; i++){
        assert(ett.subtree_query(i, i - 1) == prefix[n] - prefix[i]);
        assert(ett.subtree_query(i - 1, i) == prefix[i]);
    }

    assert(ett.cut(mid, mid - 1));
    assert(!ett.connected(0, n - 1));
    assert(ett.link(n - 1, 0));
    assert(ett.subtree_query(mid, mid + 1) == value[mid]);
    assert(ett.subtree_query(mid + 1, mid) == prefix[n] - value[mid]);
}

int main(){
    /// Without the assert a query on a non-edge would read arc -1 and return garbage
    auto query_non_edge = []{ EulerTourTree<long long>(3).subtree_query(0, 1); };
    stress::expect_assert("data_structures/euler_tour_tree.cpp", query_non_edge);

    for (long long it = 0; it < stress::scaled(2000); it++){
        int n = it < 400 ? it % 8 + 1 : stress::rand_int(1, 40);
        check_random(n, stress::rand_int(1, 300), it % 3 ? 1000000000 : 3);
    }
    for (int it = 0; it < 3; it++) check_random(1000, 20000, 1000000000);

    check_long_path(200000);

    return 0;
}
