#include "../common.h"

#define main library_main
#include "../../code_library/trees/static_top_tree.cpp"
#undef main

const long long MOD = 998244353;

/// Every call is counted, so the test can check the O(log n) calls per update and per value
long long calls = 0;

/// The sum over every w of P(r, w), P the composite of the edge maps b x + c from r towards w applied to a[w]
struct Composite{
    struct Path{ long long A = 1, B = 0, sum = 0, cnt = 0; };
    struct Point{ long long sum = 0, cnt = 0; };
    vector<long long> *a, *b, *c;

    long long mul(int e) const{ return e == -1 ? 1 : (*b)[e]; }
    long long add(int e) const{ return e == -1 ? 0 : (*c)[e]; }

    Path vertex(int v, int e) const{
        calls++;
        return {mul(e), add(e), (mul(e) * (*a)[v] + add(e)) % MOD, 1};
    }
    Path add_vertex(const Point& x, int v, int e) const{
        calls++;
        long long inside = ((*a)[v] + x.sum) % MOD, count = 1 + x.cnt;
        return {mul(e), add(e), (mul(e) * inside + add(e) * count) % MOD, count};
    }
    Point add_edge(const Path& p) const{
        calls++;
        return {p.sum, p.cnt};
    }
    Point rake(const Point& x, const Point& y) const{
        calls++;
        return {(x.sum + y.sum) % MOD, x.cnt + y.cnt};
    }
    Path compress(const Path& t, const Path& d) const{
        calls++;
        return {t.A * d.A % MOD, (t.A * d.B + t.B) % MOD, (t.sum + t.A * d.sum + t.B * (d.cnt % MOD)) % MOD, t.cnt + d.cnt};
    }
};

/// The farthest vertex from r, in max-plus arithmetic, so compress does not commute
struct Farthest{
    struct Path{ long long len = 0, far = 0; };
    struct Point{ long long far = 0; };
    vector<long long>* weight;

    long long w(int e) const{ return e == -1 ? 0 : (*weight)[e]; }
    Path vertex(int, int e) const{ return {w(e), w(e)}; }
    Path add_vertex(const Point& x, int, int e) const{ return {w(e), w(e) + x.far}; }
    Point add_edge(const Path& p) const{ return {p.far}; }
    Point rake(const Point& x, const Point& y) const{ return {max(x.far, y.far)}; }
    Path compress(const Path& t, const Path& d) const{ return {t.len + d.len, max(t.far, t.len + d.far)}; }
};

/// Random labeled tree of a chosen shape: random parents, path, star, caterpillar, or a broom (path then random)
vector<pair<int, int>> make_tree(int n, int shape){
    vector<int> label(n);
    iota(label.begin(), label.end(), 0);
    shuffle(label.begin(), label.end(), stress::rng());

    vector<pair<int, int>> edges;
    for (int i = 1; i < n; i++){
        int p;
        if (shape == 0) p = stress::rand_int(0, i - 1);
        else if (shape == 1) p = i - 1;
        else if (shape == 2) p = 0;
        else if (shape == 3) p = i % 2 ? i - 1 : max(0, i - 2);
        else p = i < n / 2 ? i - 1 : stress::rand_int(max(0, n / 2 - 3), i - 1);
        if (stress::rand_int(0, 1)) edges.push_back({label[i], label[p]});
        else edges.push_back({label[p], label[i]});
    }
    shuffle(edges.begin(), edges.end(), stress::rng());
    return edges;
}

vector<vector<pair<int, int>>> adjacency(int n, const vector<pair<int, int>>& edges){
    vector<vector<pair<int, int>>> adj(n);
    for (int e = 0; e < (int)edges.size(); e++){
        adj[edges[e].first].push_back({edges[e].second, e});
        adj[edges[e].second].push_back({edges[e].first, e});
    }

    return adj;
}

/// Walks from r keeping the composite map from r to each vertex
long long naive_composite(const vector<vector<pair<int, int>>>& adj, int r, const vector<long long>& a,
                          const vector<long long>& b, const vector<long long>& c){
    int n = adj.size();
    vector<long long> A(n, 1), B(n, 0);
    vector<int> order = {r}, seen(n, 0);
    seen[r] = 1;
    long long res = 0;
    for (int i = 0; i < (int)order.size(); i++){
        int u = order[i];
        res = (res + A[u] * a[u] + B[u]) % MOD;
        for (auto [v, e] : adj[u]){
            if (seen[v]) continue;
            seen[v] = 1, order.push_back(v);
            A[v] = A[u] * b[e] % MOD, B[v] = (A[u] * c[e] + B[u]) % MOD;
        }
    }

    return res;
}

long long naive_farthest(const vector<vector<pair<int, int>>>& adj, int r, const vector<long long>& weight){
    int n = adj.size();
    vector<long long> dist(n, -1);
    vector<int> order = {r};
    dist[r] = 0;
    for (int i = 0; i < (int)order.size(); i++){
        for (auto [v, e] : adj[order[i]]){
            if (dist[v] == -1) dist[v] = dist[order[i]] + weight[e], order.push_back(v);
        }
    }

    return *max_element(dist.begin(), dist.end());
}

void run_composite(int n, int shape, int rounds){
    vector<pair<int, int>> edges = make_tree(n, shape);
    vector<vector<pair<int, int>>> adj = adjacency(n, edges);
    vector<long long> a(n), b(n - 1), c(n - 1);
    for (auto& x : a) x = stress::rand_int(0, MOD - 1);
    for (int e = 0; e < n - 1; e++) b[e] = stress::rand_int(0, MOD - 1), c[e] = stress::rand_int(0, MOD - 1);

    int root = stress::rand_int(0, n - 1);
    StaticTopTree<Composite> tree(n, edges, Composite{&a, &b, &c}, root);
    for (int round = 0; round < rounds; round++){
        if (n > 1 && stress::rand_int(0, 1)){
            int e = stress::rand_int(0, n - 2);
            b[e] = stress::rand_int(0, MOD - 1), c[e] = stress::rand_int(0, MOD - 1);
            tree.update_edge(e);
        }
        else{
            int v = stress::rand_int(0, n - 1);
            a[v] = stress::rand_int(0, MOD - 1);
            tree.update(v);
        }

        assert(tree.value().sum == naive_composite(adj, root, a, b, c));
        int r = stress::rand_int(0, n - 1);
        Composite::Path from_r = tree.value(r);
        assert(from_r.sum == naive_composite(adj, r, a, b, c));
        assert(from_r.cnt == n);
    }
}

void run_farthest(int n, int shape, int rounds){
    vector<pair<int, int>> edges = make_tree(n, shape);
    vector<vector<pair<int, int>>> adj = adjacency(n, edges);
    vector<long long> weight(n - 1);
    for (auto& w : weight) w = stress::rand_int(0, 1000);

    int root = stress::rand_int(0, n - 1);
    StaticTopTree<Farthest> tree(n, edges, Farthest{&weight}, root);
    StaticTopTree<Farthest> fixed(n, edges, Farthest{&weight}, root, StaticTopTree<Farthest>::FIXED_ROOT);
    for (int round = 0; round < rounds; round++){
        if (n > 1){
            int e = stress::rand_int(0, n - 2);
            weight[e] = stress::rand_int(0, 1000);
            tree.update_edge(e), fixed.update_edge(e);
        }

        assert(tree.value().far == naive_farthest(adj, root, weight));
        assert(fixed.value().far == tree.value().far);
        for (int k = 0; k < 3; k++){
            int r = stress::rand_int(0, n - 1);
            assert(tree.value(r).far == naive_farthest(adj, r, weight));
        }
    }
}

/// Every vertex of a small tree as the root, after every update, catches a wrong side at any chain position
void exhaustive_roots(int n, int shape){
    vector<pair<int, int>> edges = make_tree(n, shape);
    vector<vector<pair<int, int>>> adj = adjacency(n, edges);
    vector<long long> a(n), b(n - 1), c(n - 1);
    for (auto& x : a) x = stress::rand_int(0, 50);
    for (int e = 0; e < n - 1; e++) b[e] = stress::rand_int(0, 50), c[e] = stress::rand_int(0, 50);

    StaticTopTree<Composite> tree(n, edges, Composite{&a, &b, &c}, stress::rand_int(0, n - 1));
    for (int round = 0; round < 4; round++){
        for (int r = 0; r < n; r++) assert(tree.value(r).sum == naive_composite(adj, r, a, b, c));
        int v = stress::rand_int(0, n - 1);
        a[v] = stress::rand_int(0, 50);
        tree.update(v);
    }
}

/// O(log n) calls per update and per value, on a deep path and a star where a linear walk would show
void call_counts(int shape){
    int n = 1 << 16;
    vector<pair<int, int>> edges = make_tree(n, shape);
    vector<long long> a(n, 1), b(n - 1, 1), c(n - 1, 1);
    StaticTopTree<Composite> tree(n, edges, Composite{&a, &b, &c});

    long long worst_update = 0, worst_value = 0;
    for (int k = 0; k < 2000; k++){
        int v = stress::rand_int(0, n - 1);
        calls = 0, tree.update(v), worst_update = max(worst_update, calls);
        calls = 0, tree.value(v), worst_value = max(worst_value, calls);
    }
    assert(worst_update <= 8 * 16 && worst_value <= 8 * 16);
}

int main(){
    run_composite(1, 0, 5);
    run_farthest(1, 0, 5);
    for (int shape = 0; shape < 5; shape++){
        run_composite(2, shape, 10);
        for (int n = 3; n <= 9; n++) exhaustive_roots(n, shape);
    }

    for (int it = 0; it < stress::scaled(300); it++){
        int n = stress::rand_int(1, 40), shape = stress::rand_int(0, 4);
        run_composite(n, shape, 30);
        run_farthest(n, shape, 30);
    }
    for (int shape = 0; shape < 5; shape++){
        run_composite(3000, shape, 60);
        run_farthest(3000, shape, 60);
    }

    call_counts(1);
    call_counts(2);
    call_counts(0);

    /// Not a tree: one edge too few, and a cycle that leaves a vertex unreached
    stress::expect_assert("trees/static_top_tree.cpp", []{
        vector<long long> w = {1};
        StaticTopTree<Farthest> broken(3, {{0, 1}}, Farthest{&w});
    });
    stress::expect_assert("trees/static_top_tree.cpp", []{
        vector<long long> w = {1, 1, 1};
        StaticTopTree<Farthest> broken(4, {{0, 1}, {1, 2}, {2, 0}}, Farthest{&w});
    });

    /// FIXED_ROOT keeps no reversed Paths, so value(r) must refuse
    stress::expect_assert("trees/static_top_tree.cpp", []{
        vector<long long> w = {1};
        StaticTopTree<Farthest> fixed(2, {{0, 1}}, Farthest{&w}, 0, StaticTopTree<Farthest>::FIXED_ROOT);
        fixed.value(1);
    });

    return 0;
}
