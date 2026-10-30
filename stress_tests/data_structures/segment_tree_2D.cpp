#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/segment_tree_2D.cpp"
#undef main

template <typename T, typename Combine>
T brute(const vector<vector<T>>& a, int x1, int y1, int x2, int y2, T identity, Combine combine){
    T res = identity;
    for (int x = x1; x <= x2; x++) for (int y = y1; y <= y2; y++) res = combine(res, a[x][y]);
    return res;
}

/// Grid built tree and update built tree against a plain grid; lo and hi bound the values
template <typename T, typename Combine>
void check(int n, int m, int ops, T lo, T hi, T identity, Combine combine){
    vector<vector<T>> a(n, vector<T>(m));
    for (auto& row : a) for (auto& v : row) v = stress::rand_int(lo, hi);

    SegmentTree2D built(a, identity, combine);
    SegmentTree2D grown(n, m, identity, combine);
    for (int x = 0; x < n; x++) for (int y = 0; y < m; y++) grown.update(x, y, a[x][y]);

    /// Every rectangle on small grids catches a node that one size or one boundary leaves stale
    if (n * m <= 30){
        for (int x1 = 0; x1 < n; x1++) for (int x2 = x1; x2 < n; x2++){
            for (int y1 = 0; y1 < m; y1++) for (int y2 = y1; y2 < m; y2++){
                T want = brute(a, x1, y1, x2, y2, identity, combine);
                assert(built.query(x1, y1, x2, y2) == want);
                assert(grown.query(x1, y1, x2, y2) == want);
            }
        }
    }

    for (int op = 0; op < ops; op++){
        int x1 = stress::rand_int(0, n - 1), x2 = stress::rand_int(0, n - 1);
        int y1 = stress::rand_int(0, m - 1), y2 = stress::rand_int(0, m - 1);
        if (x1 > x2) swap(x1, x2);
        if (y1 > y2) swap(y1, y2);
        int kind = stress::rand_int(0, 9);

        if (kind < 4){
            T v = stress::rand_int(lo, hi);
            a[x1][y1] = v;
            built.update(x1, y1, v);
            grown.update(x1, y1, v);
        }
        else if (kind == 4){
            assert(built.query(x1, y1) == a[x1][y1]);
            assert(built.query(x1, y1, x1, y1) == a[x1][y1]);
        }
        else{
            if (kind == 5) x1 = 0, y1 = 0, x2 = n - 1, y2 = m - 1;
            T want = brute(a, x1, y1, x2, y2, identity, combine);
            assert(built.query(x1, y1, x2, y2) == want);
            assert(grown.query(x1, y1, x2, y2) == want);
        }
    }
}

template <typename Call>
void check_rejects(Call call){
    stress::expect_assert("data_structures/segment_tree_2D.cpp", call);
}

int main(){
    auto mn = [](int a, int b){ return min(a, b); };
    auto mx = [](int a, int b){ return max(a, b); };
    auto gcd_ll = [](long long a, long long b){ return __gcd(a, b); };

    /// An index past n or m lands inside the next row's storage, so without the asserts it corrupts silently
    for (int n : {1, 5}){
        int m = n + 2;
        check_rejects([&]{ SegmentTree2D<int>(n, m).update(n, 0, 1); });
        check_rejects([&]{ SegmentTree2D<int>(n, m).update(0, m, 1); });
        check_rejects([&]{ SegmentTree2D<int>(n, m).update(-1, 0, 1); });
        check_rejects([&]{ SegmentTree2D<int>(n, m).query(0, m); });
        check_rejects([&]{ SegmentTree2D<int>(n, m).query(0, 0, n, 0); });
        check_rejects([&]{ SegmentTree2D<int>(n, m).query(0, 1, 0, 0); });
    }

    for (long long it = 0; it < stress::scaled(2000); it++){
        int big = it % 10 == 0 ? 40 : 9;
        int n = stress::rand_int(1, big), m = stress::rand_int(1, big);
        if (it % 7 == 0) n = 1;
        if (it % 11 == 0) m = 1;

        check(n, m, 200, INT_MIN, INT_MAX, INT_MAX, mn);
        check(n, m, 200, INT_MIN, INT_MAX, INT_MIN, mx);
        check<long long>(n, m, 200, -1000000000000LL, 1000000000000LL, 0, plus<long long>());
        check<long long>(n, m, 200, 0, 60, 0, gcd_ll);
    }

    /// Odd sizes give nodes whose range wraps past the row end, which a query must never pick; a power of two has none
    check(64, 64, 3000, -50, 50, INT_MAX, mn);
    check(37, 53, 3000, -50, 50, INT_MIN, mx);

    return 0;
}
