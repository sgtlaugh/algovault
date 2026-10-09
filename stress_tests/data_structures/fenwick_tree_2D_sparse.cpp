#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/fenwick_tree_2D_sparse.cpp"
#undef main

/// Point updates and rectangle queries against a list of every update, one instance because each preallocates ~270 MB
int main(){
    const int N = 1000000000, M = 1000000000;
    static FenwickSparse2D<long long> fen(N, M);
    vector<array<long long, 3>> points;  /// i, j, value

    auto coordinate = [&](int limit){
        int mode = stress::rand_int(0, 3);
        if (mode == 0) return (int)stress::rand_int(1, limit);
        if (mode == 1) return (int)stress::rand_int(1, 50);                       /// dense corner, many shared tree nodes
        if (mode == 2) return limit - (int)stress::rand_int(0, 3);                /// the grid edge, where update loops end
        return (int)min<long long>(limit, 1000003LL * stress::rand_int(1, 900));   /// the stride the old key scheme collided on
    };

    auto brute = [&](int i, int j, int k, int l){
        long long sum = 0;
        for (auto& p : points) if (i <= p[0] && p[0] <= k && j <= p[1] && p[1] <= l) sum += p[2];
        return sum;
    };

    /// At most 3000 points keep the brute force linear and the hash table at a low load, larger scales add queries
    for (long long it = 0; it < stress::scaled(3000); it++){
        int i = coordinate(N), j = coordinate(M);
        if (points.size() < 3000){
            long long v = stress::rand_int(-1000, 1000);
            fen.update(i, j, v);
            points.push_back({i, j, v});
        }

        int a = coordinate(N), b = coordinate(M), c = coordinate(N), d = coordinate(M);
        if (a > c) swap(a, c);
        if (b > d) swap(b, d);
        assert(fen.query(a, b, c, d) == brute(a, b, c, d));
        assert(fen.query(i, j, i, j) == brute(i, j, i, j));
        assert(fen.query(c, d) == brute(1, 1, c, d));
    }

    /// Small tables filled to about 60% load, so probe chains often run past the last slot and wrap to the first
    for (long long round = 0; round < stress::scaled(100); round++){
        const int n = stress::rand_int(100, 1000), m = stress::rand_int(100, 1000), lg = 10;  /// enough cells to reach the target load
        FenwickSparse2D<long long> small(n, m, lg);
        vector<array<long long, 3>> pts;
        set<pair<int, int>> nodes;

        while (nodes.size() < (1 << lg) * 6 / 10){
            int i = stress::rand_int(1, n), j = stress::rand_int(1, m);
            long long v = stress::rand_int(-1000, 1000);
            small.update(i, j, v);
            pts.push_back({i, j, v});
            for (int x = i; x <= n; x += x & -x) for (int y = j; y <= m; y += y & -y) nodes.insert({x, y});
        }

        for (int q = 0; q < 200; q++){
            int a = stress::rand_int(1, n), b = stress::rand_int(1, m), c = stress::rand_int(1, n), d = stress::rand_int(1, m);
            if (a > c) swap(a, c);
            if (b > d) swap(b, d);
            long long sum = 0;
            for (auto& p : pts) if (a <= p[0] && p[0] <= c && b <= p[1] && p[1] <= d) sum += p[2];
            assert(small.query(a, b, c, d) == sum);
        }
    }

    return 0;
}
