// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/vertex_set_path_composite
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/trees/hld.cpp"
#undef main

const long long MOD = 998244353;

struct Affine{
    long long a = 1, b = 0;
};

/// g after f
Affine compose(Affine f, Affine g){
    return {g.a * f.a % MOD, (g.a * f.b + g.b) % MOD};
}

/// Point-set tree keeping each range's composite in both walking directions, 1-indexed like hld.pos
struct CompositeTree{
    int size;
    vector<Affine> fwd, bwd;

    CompositeTree(int n) : size(1){
        while (size < n + 1) size *= 2;
        fwd.assign(2 * size, Affine()), bwd.assign(2 * size, Affine());
    }

    void set(int p, Affine f){
        p += size;
        fwd[p] = bwd[p] = f;
        for (p /= 2; p; p /= 2){
            fwd[p] = compose(fwd[2 * p], fwd[2 * p + 1]);
            bwd[p] = compose(bwd[2 * p + 1], bwd[2 * p]);
        }
    }

    /// Composite of positions l..r walked upward when forward, downward otherwise
    Affine query(int l, int r, bool forward){
        Affine left, right;
        for (l += size, r += size + 1; l < r; l /= 2, r /= 2){
            if (l & 1){
                left = forward ? compose(left, fwd[l]) : compose(bwd[l], left);
                l++;
            }
            if (r & 1){
                r--;
                right = forward ? compose(fwd[r], right) : compose(right, bwd[r]);
            }
        }
        return forward ? compose(left, right) : compose(right, left);
    }
};

int read_int(){
    int c = getchar_unlocked(), x = 0;
    while (c < '0' || c > '9') c = getchar_unlocked();
    for (; c >= '0' && c <= '9'; c = getchar_unlocked()) x = x * 10 + c - '0';
    return x;
}

int main(){
    int n = read_int(), q = read_int();
    vector<Affine> f(n);
    for (auto& [a, b] : f) a = read_int(), b = read_int();
    HLD hld(n);
    for (int i = 0; i + 1 < n; i++){
        int u = read_int(), v = read_int();
        hld.add_edge(u, v);
    }
    hld.build(0);

    CompositeTree seg(n);
    for (int v = 0; v < n; v++) seg.set(hld.pos[v], f[v]);

    while (q--){
        int t = read_int(), x = read_int(), y = read_int(), z = read_int();
        if (t == 0){
            seg.set(hld.pos[x], {y, z});
            continue;
        }
        Affine res;
        for (auto [a, b] : hld.ordered_path(x, y)) res = compose(res, a <= b ? seg.query(a, b, true) : seg.query(b, a, false));
        printf("%lld\n", (res.a * z + res.b) % MOD);
    }
    return 0;
}
