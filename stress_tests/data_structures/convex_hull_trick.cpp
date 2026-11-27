#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/convex_hull_trick.cpp"
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

double rand_real(double lo, double hi){
    return uniform_real_distribution<double>(lo, hi)(stress::rng());
}

/// Floating-point lines, near-parallel or nearly through one point, must stay within rounding of the exact optimum
template <typename T, bool MAX>
void check_floating(int mode){
    int n = stress::rand_int(1, 200);
    double range = mode == 2 ? 1e9 : 1e3, x0 = rand_real(-range, range), y0 = rand_real(-range, range);
    vector<T> k(n), b(n);
    for (auto& x : k) x = mode == 1 ? 1 + 1e-12 * stress::rand_int(0, 50) : rand_real(-range, range);
    sort(k.begin(), k.end());
    if (!MAX) reverse(k.begin(), k.end());
    for (int i = 0; i < n; i++) b[i] = mode == 3 ? y0 - k[i] * x0 + rand_real(-1e-9, 1e-9) : rand_real(-range, range) * range;

    ConvexHullTrick<T, MAX> hull, walk;
    for (int i = 0; i < n; i++) hull.add_line(k[i], b[i]), walk.add_line(k[i], b[i]);

    vector<T> xs(30);
    for (auto& x : xs) x = mode == 3 ? x0 + rand_real(-1e-6, 1e-6) : rand_real(-range, range);
    sort(xs.begin(), xs.end());
    for (T x : xs){
        long double best = (long double)k[0] * x + b[0], scale = 1;
        for (int i = 0; i < n; i++){
            long double v = (long double)k[i] * x + b[i];
            best = MAX ? max(best, v) : min(best, v);
            scale = max(scale, fabsl((long double)k[i] * x) + fabsl((long double)b[i]));
        }
        assert(fabsl(hull.query(x) - best) <= 1e-14 * scale);
        assert(fabsl(walk.query_increasing(x) - best) <= 1e-14 * scale);
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

    /// |b| near 9e18 with small k, x: every k * x + b fits, but intercept differences overflow long long
    const long long K2 = 10000000, B2 = 9000000000000000000LL;
    for (long long it = 0; it < stress::scaled(300); it++){
        int n = stress::rand_int(1, 300);
        check_random<false>(n, 2, K2, B2, K2);
        check_random<true>(n, 2, K2, B2, K2);
        check_increasing<false>(n, K2, B2, K2);
        check_increasing<true>(n, K2, B2, K2);
    }

    /// int lines whose k * x + b fit in int but whose intercepts differ by more than INT_MAX
    ConvexHullTrick<int> small;
    small.add_line(2, -1200000000);
    small.add_line(1, -300000000);
    small.add_line(0, 1200000000);
    for (int x : {0, 300000000, 500000000, 900000000}){
        long long want = min({2LL * x - 1200000000, 1LL * x - 300000000, 1200000000LL});
        assert(small.query(x) == want);
    }

    /// int slopes near 1.5e9 differ by more than INT_MAX, so the slope differences must widen too
    ConvexHullTrick<int> steep;
    steep.add_line(1500000000, 0);
    steep.add_line(0, 1);
    steep.add_line(-1500000000, 0);
    for (int x : {-1, 0, 1}) assert(steep.query(x) == min({1500000000LL * x, 1LL, -1500000000LL * x}));

    for (long long it = 0; it < stress::scaled(1000); it++){
        int mode = it % 4;
        check_floating<double, false>(mode);
        check_floating<double, true>(mode);
        check_floating<long double, false>(mode);
    }

    /// A slope out of order would silently break the hull, so it must trip the library's assert
    stress::expect_assert("data_structures/convex_hull_trick.cpp", []{
        ConvexHullTrick<long long> hull;
        hull.add_line(1, 0);
        hull.add_line(2, 0);
    });
    stress::expect_assert("data_structures/convex_hull_trick.cpp", []{
        ConvexHullTrick<long long, true> hull;
        hull.add_line(2, 0);
        hull.add_line(1, 0);
    });

    return 0;
}
