#include "../common.h"

#define main library_main
#include "../../code_library/numerical/dinkelbach.cpp"
#undef main

using Frac = pair<long long, long long>;

bool less_frac(Frac x, Frac y){
    return (__int128)x.first * y.second < (__int128)y.first * x.second;
}

bool same_frac(Frac x, Frac y){
    return (__int128)x.first * y.second == (__int128)y.first * x.second;
}

/// Random inputs take at most lg(n M) calls at STRESS_SCALE=20; more means convergence is no longer superlinear
void check_calls(long long calls, long long n, long long m){
    assert(calls <= 2 + __lg(max(2LL, n * m)));
}

vector<long long> random_values(int n, long long lo, long long hi, int mode){
    long long c = stress::rand_int(lo, hi);
    vector<long long> v(n);
    for (auto& x : v) x = mode == 0 ? c : mode == 1 ? stress::rand_int(lo, min(hi, lo + 2)) : stress::rand_int(lo, hi);
    return v;
}

/// Oracle for "keep k of n items": the k largest q a_i - p b_i, ties broken arbitrarily
template<typename T>
pair<T, T> top_k(const vector<long long>& a, const vector<long long>& b, int k, T p, T q, bool largest){
    int n = a.size();
    vector<int> order(n);
    iota(order.begin(), order.end(), 0);
    shuffle(order.begin(), order.end(), stress::rng());
    using K = conditional_t<is_integral<T>::value, __int128, T>;
    auto key = [&](int i){ return largest ? (K)q * a[i] - (K)p * b[i] : (K)p * b[i] - (K)q * a[i]; };
    nth_element(order.begin(), order.begin() + (k - 1), order.end(), [&](int i, int j){ return key(i) > key(j); });

    pair<T, T> res(0, 0);
    for (int t = 0; t < k; t++) res.first += a[order[t]], res.second += b[order[t]];

    return res;
}

void check_items_small(int n, int k, long long m){
    vector<long long> a = random_values(n, stress::rand_int(0, 1) ? -m : 0, m, stress::rand_int(0, 3));
    vector<long long> b = random_values(n, 1, m, stress::rand_int(0, 3));

    set<Frac> sums;
    Frac lo(LLONG_MAX, 1), hi(LLONG_MIN, 1);
    for (int mask = 0; mask < 1 << n; mask++){
        if (__builtin_popcount(mask) != k) continue;
        Frac s(0, 0);
        for (int i = 0; i < n; i++) if (mask >> i & 1) s.first += a[i], s.second += b[i];
        sums.insert(s);
        if (lo.first == LLONG_MAX || less_frac(s, lo)) lo = s;
        if (hi.first == LLONG_MIN || less_frac(hi, s)) hi = s;
    }
    vector<Frac> all(sums.begin(), sums.end());
    Frac feasible = all[stress::rand_int(0, all.size() - 1)];

    long long calls = 0;
    auto best_max = [&](long long p, long long q){ calls++; return top_k<long long>(a, b, k, p, q, true); };
    auto best_min = [&](long long p, long long q){ calls++; return top_k<long long>(a, b, k, p, q, false); };

    /// A feasible start returns a real x's sums; a loose bound may return itself only if it equals the optimum
    bool loose = stress::rand_int(0, 1);
    Frac got = max_ratio(best_max, loose ? -m * n - stress::rand_int(0, 5) : feasible.first, loose ? 1LL : feasible.second);
    assert(same_frac(got, hi));
    assert(sums.count(got));
    check_calls(calls, n, m);

    calls = 0;
    got = min_ratio(best_min, loose ? m * n + stress::rand_int(0, 5) : feasible.first, loose ? 1LL : feasible.second);
    assert(same_frac(got, lo));
    assert(sums.count(got));
    check_calls(calls, n, m);

    /// An oracle returning a random x of positive value instead of the best, like a negative cycle finder
    auto any_improving = [&](long long p, long long q){
        vector<Frac> up;
        for (Frac s : all) if ((__int128)q * s.first - (__int128)p * s.second > 0) up.push_back(s);
        return up.empty() ? all[0] : up[stress::rand_int(0, up.size() - 1)];
    };
    got = max_ratio(any_improving, feasible.first, feasible.second);
    assert(same_frac(got, hi));
    assert(sums.count(got));

    /// Starting exactly at the optimum: one call, and the start comes back
    calls = 0;
    got = max_ratio(best_max, hi.first, hi.second);
    assert(got == hi);
    assert(calls == 1);

    auto best_max_d = [&](double p, double q){ return top_k<double>(a, b, k, p, q, true); };
    auto best_min_d = [&](double p, double q){ return top_k<double>(a, b, k, p, q, false); };
    pair<double, double> gd = max_ratio(best_max_d, double(feasible.first), double(feasible.second));
    double want = double(hi.first) / hi.second;
    assert(abs(gd.first / gd.second - want) <= 1e-12 * max(1.0, abs(want)));
    gd = min_ratio(best_min_d, double(feasible.first), double(feasible.second));
    want = double(lo.first) / lo.second;
    assert(abs(gd.first / gd.second - want) <= 1e-12 * max(1.0, abs(want)));
}

/// Oracle for the spanning tree problem: Prim on the complete graph with weight q c - p l, minimal or maximal
template<typename T>
pair<T, T> prim(const vector<vector<long long>>& c, const vector<vector<long long>>& l, T p, T q, bool smallest){
    int n = c.size();
    auto w = [&](int u, int v){ T x = q * T(c[u][v]) - p * T(l[u][v]); return smallest ? x : -x; };
    vector<bool> in(n, false);
    vector<int> from(n, 0);
    pair<T, T> res(0, 0);
    in[0] = true;
    for (int it = 1; it < n; it++){
        int v = -1;
        for (int u = 0; u < n; u++){
            if (!in[u] && (v == -1 || w(from[u], u) < w(from[v], v))) v = u;
        }
        in[v] = true;
        res.first += c[from[v]][v], res.second += l[from[v]][v];
        for (int u = 0; u < n; u++) if (!in[u] && w(v, u) < w(from[u], u)) from[u] = v;
    }

    return res;
}

int find_root(vector<int>& parent, int x){
    while (parent[x] != x) x = parent[x];
    return x;
}

/// Costs in [0, m] and lengths in [1, m] on the complete graph, as in POJ 2728 Desert King
pair<vector<vector<long long>>, vector<vector<long long>>> random_graph(int n, long long m){
    vector<vector<long long>> c(n, vector<long long>(n)), l = c;
    vector<long long> cv = random_values(n * n, 0, m, stress::rand_int(0, 3));
    vector<long long> lv = random_values(n * n, 1, m, stress::rand_int(0, 3));
    for (int u = 0; u < n; u++){
        for (int v = u + 1; v < n; v++) c[u][v] = c[v][u] = cv[u * n + v], l[u][v] = l[v][u] = lv[u * n + v];
    }

    return {c, l};
}

void check_trees_small(int n, long long m){
    auto [c, l] = random_graph(n, m);
    vector<pair<int, int>> edges;
    for (int u = 0; u < n; u++){
        for (int v = u + 1; v < n; v++) edges.push_back({u, v});
    }

    int e = edges.size();
    set<Frac> sums;
    Frac lo(-1, 0), hi(-1, 0);
    for (int mask = 0; mask < 1 << e; mask++){
        if (__builtin_popcount(mask) != n - 1) continue;
        vector<int> parent(n);
        iota(parent.begin(), parent.end(), 0);
        bool tree = true;
        Frac s(0, 0);
        for (int i = 0; i < e && tree; i++){
            if (!(mask >> i & 1)) continue;
            int ru = find_root(parent, edges[i].first), rv = find_root(parent, edges[i].second);
            if (ru == rv) tree = false;
            parent[ru] = rv;
            s.first += c[edges[i].first][edges[i].second], s.second += l[edges[i].first][edges[i].second];
        }
        if (!tree) continue;
        sums.insert(s);
        if (lo.second == 0 || less_frac(s, lo)) lo = s;
        if (hi.second == 0 || less_frac(hi, s)) hi = s;
    }

    long long calls = 0;
    auto best_min = [&](long long p, long long q){ calls++; return prim<long long>(c, l, p, q, true); };
    auto best_max = [&](long long p, long long q){ calls++; return prim<long long>(c, l, p, q, false); };
    Frac got = min_ratio(best_min, m * n, 1LL);
    assert(same_frac(got, lo));
    assert(sums.count(got));
    check_calls(calls, n, m);

    calls = 0;
    got = max_ratio(best_max, 0LL, 1LL);
    assert(same_frac(got, hi));
    assert(sums.count(got) || got == Frac(0, 1));  /// every cost 0: the start is already optimal and comes back
    check_calls(calls, n, m);
}

/// Too many subsets to list: the answer p / q must make the best value of q a(x) - p b(x) exactly 0
void check_items_large(int n, long long m){
    int k = stress::rand_int(1, n);
    vector<long long> a = random_values(n, -m, m, stress::rand_int(0, 5) ? 2 : stress::rand_int(0, 1));
    vector<long long> b = random_values(n, 1, m, stress::rand_int(0, 5) ? 2 : stress::rand_int(0, 1));
    long long calls = 0;
    auto best = [&](long long p, long long q){ calls++; return top_k<long long>(a, b, k, p, q, true); };
    Frac got = max_ratio(best, -m, 1LL);
    check_calls(calls, n, m);

    vector<__int128> value(n);
    for (int i = 0; i < n; i++) value[i] = (__int128)got.second * a[i] - (__int128)got.first * b[i];
    sort(value.rbegin(), value.rend());
    __int128 top = 0;
    for (int i = 0; i < k; i++) top += value[i];
    assert(top == 0);

    auto best_d = [&](double p, double q){ return top_k<double>(a, b, k, p, q, true); };
    pair<double, double> gd = max_ratio(best_d, double(-m), 1.0);
    double want = double(got.first) / got.second;
    assert(abs(gd.first / gd.second - want) <= 1e-12 * max(1.0, abs(want)));
}

/// Desert King sized graphs: exact against double, and the exact answer zeroes the best q c - p l
void check_trees_large(int n, long long m){
    auto [c, l] = random_graph(n, m);
    long long calls = 0;
    auto best = [&](long long p, long long q){ calls++; return prim<long long>(c, l, p, q, true); };
    Frac got = min_ratio(best, m, 1LL);
    check_calls(calls, n, m);

    Frac tree = prim<long long>(c, l, got.first, got.second, true);
    assert(got.second * tree.first - got.first * tree.second == 0);

    auto best_d = [&](double p, double q){ return prim<double>(c, l, p, q, true); };
    pair<double, double> gd = min_ratio(best_d, double(m), 1.0);
    double want = double(got.first) / got.second;
    assert(abs(gd.first / gd.second - want) <= 1e-12 * max(1.0, abs(want)));
}

int main(){
    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(1, 10);
        check_items_small(n, stress::rand_int(1, n), stress::rand_int(0, 2) ? 20 : 1000000);
    }

    for (long long it = 0; it < stress::scaled(400); it++){
        check_trees_small(stress::rand_int(2, 6), stress::rand_int(0, 1) ? 10 : 1000000);
    }

    for (long long it = 0; it < stress::scaled(300); it++) check_items_large(stress::rand_int(1, 2000), 1000000);

    /// Sums near 2^51, so the cross products in max_ratio's stop test overflow 64 bits
    for (long long it = 0; it < stress::scaled(100); it++) check_items_large(stress::rand_int(1, 2000), 1000000000000LL);

    for (long long it = 0; it < stress::scaled(100); it++) check_trees_large(stress::rand_int(2, 100), 1000000);

    stress::expect_assert("numerical/dinkelbach.cpp", []{
        max_ratio([](long long, long long){ return Frac(1, 0); }, 0LL, 1LL);
    });

    return 0;
}
