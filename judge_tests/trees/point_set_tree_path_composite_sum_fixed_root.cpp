// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/point_set_tree_path_composite_sum_fixed_root
// competitive-verifier: TLE 10
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/trees/static_top_tree.cpp"
#undef main

const long long MOD = 998244353;

int read_int(){
    int c = getchar_unlocked(), x = 0;
    while (c < '0' || c > '9') c = getchar_unlocked();
    for (; c >= '0' && c <= '9'; c = getchar_unlocked()) x = x * 10 + c - '0';
    return x;
}

/// sum is the sum over every w in the cluster of P(outside, w), the maps of the edges from the outside end applied
struct Composite{
    struct Path{ long long A = 1, B = 0, sum = 0, cnt = 0; };
    struct Point{ long long sum = 0, cnt = 0; };
    vector<long long> *a, *b, *c;

    long long mul(int e) const{ return e == -1 ? 1 : (*b)[e]; }
    long long add(int e) const{ return e == -1 ? 0 : (*c)[e]; }

    Path vertex(int v, int e) const{
        return {mul(e), add(e), (mul(e) * (*a)[v] + add(e)) % MOD, 1};
    }
    Path add_vertex(const Point& x, int v, int e) const{
        long long inside = ((*a)[v] + x.sum) % MOD, count = 1 + x.cnt;
        return {mul(e), add(e), (mul(e) * inside + add(e) * count) % MOD, count};
    }
    Point add_edge(const Path& p) const{
        return {p.sum, p.cnt};
    }
    Point rake(const Point& x, const Point& y) const{
        return {(x.sum + y.sum) % MOD, x.cnt + y.cnt};
    }
    Path compress(const Path& t, const Path& d) const{
        return {t.A * d.A % MOD, (t.A * d.B + t.B) % MOD, (t.sum + t.A * d.sum + t.B * d.cnt) % MOD, t.cnt + d.cnt};
    }
};

int main(){
    int n = read_int(), q = read_int();
    vector<long long> a(n), b(n - 1), c(n - 1);
    for (auto& x : a) x = read_int();

    vector<pair<int, int>> edges(n - 1);
    for (int e = 0; e < n - 1; e++){
        edges[e].first = read_int(), edges[e].second = read_int();
        b[e] = read_int(), c[e] = read_int();
    }

    StaticTopTree<Composite> tree(n, edges, Composite{&a, &b, &c}, 0, StaticTopTree<Composite>::FIXED_ROOT);
    while (q--){
        if (read_int() == 0){
            int w = read_int();
            a[w] = read_int();
            tree.update(w);
        }
        else{
            int e = read_int();
            b[e] = read_int(), c[e] = read_int();
            tree.update_edge(e);
        }
        printf("%lld\n", tree.value().sum);
    }

    return 0;
}
