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

    for (long long it = 0; it < stress::scaled(1000); it++){
        int mode = it % 4;
        check_floating<double, false>(mode);
        check_floating<double, true>(mode);
        check_floating<long double, false>(mode);
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
