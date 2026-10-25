#include "../common.h"

#define main library_main
#include "../../code_library/dp/convex_hull_trick.cpp"
#undef main

/// Every line ever added, dominated ones included, so the hull's pruning is checked against the plain minimum
template <bool MAX>
long long brute(const vector<pair<long long, long long>>& lines, long long x){
    long long res = MAX ? LLONG_MIN : LLONG_MAX;
    for (auto [k, b] : lines) res = MAX ? max(res, k * x + b) : min(res, k * x + b);
    return res;
}

/// Slopes in the order add_line needs, few distinct values so equal slopes and parallel lines are common
template <bool MAX>
vector<long long> sorted_slopes(int n, long long range){
    vector<long long> k(n);
    for (auto& x : k) x = stress::rand_int(0, 2) ? stress::rand_int(-range, range) : stress::rand_int(-3, 3);
    sort(k.begin(), k.end());
    if (!MAX) reverse(k.begin(), k.end());
    return k;
}

template <bool MAX>
void check_random(int n, int queries, long long k_range, long long b_range, long long x_range){
    vector<long long> slopes = sorted_slopes<MAX>(n, k_range);
    vector<pair<long long, long long>> lines;
    ConvexHullTrick<long long, MAX> hull;

    for (int i = 0; i < n; i++){
        long long b = stress::rand_int(-b_range, b_range);
        hull.add_line(slopes[i], b);
        lines.push_back({slopes[i], b});
        for (int q = 0; q < queries; q++){
            long long x = stress::rand_int(-x_range, x_range);
            assert(hull.query(x) == brute<MAX>(lines, x));
        }
    }
}

/// Adds interleaved with query_increasing on non-decreasing x, the pointer must survive lines popped behind it
template <bool MAX>
void check_increasing(int n, long long k_range, long long b_range, long long x_range){
    vector<long long> slopes = sorted_slopes<MAX>(n, k_range);
    vector<pair<long long, long long>> lines;
    ConvexHullTrick<long long, MAX> hull;
    long long x = -x_range;

    for (int i = 0; i < n; i++){
        long long b = stress::rand_int(-b_range, b_range);
        hull.add_line(slopes[i], b);
        lines.push_back({slopes[i], b});

        int queries = stress::rand_int(0, 3);
        for (int q = 0; q < queries; q++){
            x = min(x_range, x + stress::rand_int(0, max(1LL, 2 * x_range / n)));
            assert(hull.query_increasing(x) == brute<MAX>(lines, x));
            assert(hull.query(x) == brute<MAX>(lines, x));
        }
    }
}

int main(){
    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(1, it % 10 ? 12 : 200);
        check_random<false>(n, 3, 6, 20, 15);
        check_random<true>(n, 3, 6, 20, 15);
        check_increasing<false>(n, 6, 20, 15);
        check_increasing<true>(n, 6, 20, 15);
    }

    /// The documented extremes: |k|, |x| <= 1e9 and |b| <= 1e18, where the crossing tests need __int128
    const long long K = 1000000000, B = 1000000000000000000LL;
    for (long long it = 0; it < stress::scaled(300); it++){
        int n = stress::rand_int(1, 300);
        check_random<false>(n, 2, K, B, K);
        check_random<true>(n, 2, K, B, K);
        check_increasing<false>(n, K, B, K);
        check_increasing<true>(n, K, B, K);
    }

    /// A slope out of order would silently break the hull, so it must trip the library's assert
    stress::expect_assert("dp/convex_hull_trick.cpp", []{
        ConvexHullTrick<long long> hull;
        hull.add_line(1, 0);
        hull.add_line(2, 0);
    });
    stress::expect_assert("dp/convex_hull_trick.cpp", []{
        ConvexHullTrick<long long, true> hull;
        hull.add_line(2, 0);
        hull.add_line(1, 0);
    });

    return 0;
}
