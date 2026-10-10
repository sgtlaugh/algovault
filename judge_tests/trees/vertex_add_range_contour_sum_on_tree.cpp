// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/vertex_add_range_contour_sum_on_tree
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/trees/centroid_decomposition.cpp"
#undef main
#define main fenwick_main
#include "../../code_library/data_structures/fenwick_tree.cpp"
#undef main

struct Ancestor{
    int centroid, dist, piece;
};

int read_int(){
    int c = getchar_unlocked(), x = 0;
    while (c != '-' && (c < '0' || c > '9')) c = getchar_unlocked();
    bool negative = c == '-';
    if (negative) c = getchar_unlocked();
    for (; c >= '0' && c <= '9'; c = getchar_unlocked()) x = x * 10 + c - '0';
    return negative ? -x : x;
}

int main(){
    int n = read_int(), q = read_int();
    vector<int> a(n);
    for (int& x : a) x = read_int();
    CentroidDecomposition cd(n);
    for (int i = 0; i + 1 < n; i++){
        int u = read_int(), v = read_int();
        cd.add_edge(u, v);
    }

    /// One Fenwick per centroid over distances in its component, one per piece to cancel pairs inside that piece
    vector<vector<Ancestor>> anc(n);
    vector<FenwickPointUpdate<long long>> by_centroid(n), by_piece;
    cd.build([&](int c, const CentroidDecomposition::Branches& branches){
        int max_dist = 0;
        anc[c].push_back({c, 0, -1});
        for (auto& branch : branches){
            int piece = by_piece.size(), piece_max = 0;
            for (auto [v, d] : branch){
                anc[v].push_back({c, (int)d, piece});
                piece_max = max(piece_max, (int)d);
            }
            by_piece.emplace_back(piece_max + 1);
            max_dist = max(max_dist, piece_max);
        }
        by_centroid[c] = FenwickPointUpdate<long long>(max_dist + 1);
    });

    auto add = [&](int p, long long x){
        for (auto [c, d, piece] : anc[p]){
            by_centroid[c].update(d + 1, x);
            if (piece != -1) by_piece[piece].update(d + 1, x);
        }
    };
    auto range_sum = [&](FenwickPointUpdate<long long>& fen, int lo, int hi){
        lo = max(lo, 0), hi = min(hi, fen.n - 1);
        return lo > hi ? 0LL : fen.query(lo + 1, hi + 1);
    };

    for (int v = 0; v < n; v++) add(v, a[v]);
    while (q--){
        int t = read_int(), p = read_int(), x = read_int();
        if (t == 0){
            add(p, x);
            continue;
        }
        int r = read_int();
        long long res = 0;
        for (auto [c, d, piece] : anc[p]){
            res += range_sum(by_centroid[c], x - d, r - 1 - d);
            if (piece != -1) res -= range_sum(by_piece[piece], x - d, r - 1 - d);
        }
        printf("%lld\n", res);
    }
    return 0;
}
