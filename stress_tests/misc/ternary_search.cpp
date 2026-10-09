#include "../common.h"

#define main library_main
#include "../../code_library/misc/ternary_search.cpp"
#undef main

const long long LIM = 4000000000000000000LL;

long long scan_argmax(const vector<long long>& v){
    return max_element(v.begin(), v.end()) - v.begin();
}

/// Strictly increasing to a random peak, then non-increasing with frequent plateaus
vector<long long> random_peak(int n){
    vector<long long> v(n);
    int p = stress::rand_int(0, n - 1);
    v[p] = stress::rand_int(-1000000, 1000000);
    for (int i = p - 1; i >= 0; i--) v[i] = v[i + 1] - stress::rand_int(1, 5);
    for (int i = p + 1; i < n; i++) v[i] = v[i - 1] - (stress::rand_int(0, 2) ? 0 : stress::rand_int(1, 5));
    return v;
}

bool searchable(const vector<long long>& v){
    int i = 0, n = v.size();
    while (i + 1 < n && v[i] < v[i + 1]) i++;
    for (; i + 1 < n; i++) if (v[i] < v[i + 1]) return false;
    return true;
}

void check_integer(const vector<long long>& v, long long lo){
    int n = v.size();
    long long calls = 0;
    long long got = ternary_search_max(lo, lo + n - 1, [&](long long i){ calls++; return v[i - lo]; });
    assert(got - lo == scan_argmax(v));
    assert(calls <= 2 * (64 - __builtin_clzll(n)));

    /// Mirror recipe: reversed v has plateaus on the rising side, the answer is its largest argmax
    vector<long long> u(v.rbegin(), v.rend());
    long long hi = lo + n - 1;
    long long mirrored = lo + hi - ternary_search_max(lo, hi, [&](long long i){ return u[(lo + hi - i) - lo]; });
    long long top = *max_element(u.begin(), u.end()), largest = n - 1;
    while (u[largest] != top) largest--;
    assert(mirrored - lo == largest);
}

/// golden_section_min on a V with a flat bottom [c, c + w], slopes k1 and k2: x must land in the bottom
void check_golden_v(double lo, double hi, double c){
    double w = stress::rand_int(0, 3) ? 0 : (hi - c) * stress::rand_int(0, 1000) / 1e3;
    double k1 = stress::rand_int(1, 1000) / 10.0, k2 = stress::rand_int(1, 1000) / 10.0;
    auto f = [&](double x){ return x < c ? k1 * (c - x) : x > c + w ? k2 * (x - c - w) : 0.0; };

    double eps = stress::rand_int(0, 1) ? 1e-9 : 1e-12;
    long long calls = 0;
    double x = golden_section_min(lo, hi, [&](double t){ calls++; return f(t); }, eps);
    /// The final interval holds x and a minimizer, so its own stop threshold bounds the miss, not the initial bounds
    double tol = 2 * eps * max(1.0, 2 * abs(x));
    assert(lo <= x && x <= hi);
    assert(c - tol <= x && x <= c + w + tol);

    /// Every interval the loop shrinks contains x, so its stop threshold is >= eps * max(1, |x|) and width shrinks by r
    double r = (sqrt(5.0) - 1) / 2, floor_tol = eps * max(1.0, abs(x));
    double steps = hi - lo > floor_tol ? log((hi - lo) / floor_tol) / log(1 / r) : 0;
    /// Re-placing a drifted point costs 1 call; drift starts at rounding and grows by 1/r per step, measured gap >= 76
    assert(calls <= 2 + 1 + ceil(steps) + 2 + ceil(steps / 64));
}

/// golden_section_min on smooth unimodal functions against a linear scan of a fine grid
void check_golden_smooth(){
    double lo = stress::rand_int(-1000, 1000), hi = lo + stress::rand_int(0, 2000);
    double c = lo + (hi - lo) * stress::rand_int(0, 1000) / 1e3 + stress::rand_int(-100, 100);
    double a = stress::rand_int(1, 100) / 10.0, d = stress::rand_int(-1000, 1000);
    int kind = stress::rand_int(0, 2);
    auto f = [&](double x){
        double t = x - c;
        if (kind == 0) return a * t * t + d;
        if (kind == 1) return a * t * t * t * t + t * t + d;
        return cosh(t / (hi - lo + 1)) * a + d;
    };

    double x = golden_section_min(lo, hi, f);
    const int grid = 20000;
    double step = (hi - lo) / grid, best = f(lo), best_x = lo;
    for (int i = 1; i <= grid; i++){
        double g = lo + step * i;
        if (f(g) < best) best = f(g), best_x = g;
    }
    assert(lo <= x && x <= hi);
    assert(abs(x - best_x) <= step + 1e-5 * max(1.0, abs(lo) + abs(hi)));

    /// f' vanishes at an interior minimum, so f(x) matches the scan to rounding there (not at a sloped endpoint)
    double xs = min(max(c, lo), hi);
    if (lo + 1 < xs && xs < hi - 1) assert(f(x) <= best + 1e-9 * max(1.0, abs(best)));

    /// Rounding hides f's change within about sqrt(ulp(f) / a) of a smooth minimum, so x is checked loosely
    if (kind < 2) assert(abs(x - xs) <= 1e-5 * max(1.0, abs(lo) + abs(hi)));
}

int main(){
    /// Every array of length <= 7 over values 0..3 that satisfies the precondition, at a few offsets
    for (int n = 1; n <= 7; n++){
        int total = 1 << (2 * n);
        for (int mask = 0; mask < total; mask++){
            vector<long long> v(n);
            for (int i = 0; i < n; i++) v[i] = (mask >> (2 * i)) & 3;
            if (!searchable(v)) continue;
            check_integer(v, 0);
            check_integer(v, -LIM);
            check_integer(v, LIM - n + 1);
        }
    }

    /// Random peaks with plateaus, against a linear scan
    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = it < 300 ? it + 1 : stress::rand_int(1, 3000);
        long long lo = stress::rand_int(-LIM, LIM - n + 1);
        check_integer(random_peak(n), lo);
    }

    /// Full claimed range: peak anywhere, falling side flat in blocks of 1000
    for (long long it = 0; it < stress::scaled(2000); it++){
        long long p = it < 3 ? vector<long long>{-LIM, LIM, 0}[it] : stress::rand_int(-LIM, LIM);
        long long calls = 0;
        long long got = ternary_search_max(-LIM, LIM, [&](long long i){
            calls++;
            return i <= p ? i - p : -((i - p + 999) / 1000);
        });
        assert(got == p);
        assert(calls <= 2 * 63);
    }

    for (long long it = 0; it < stress::scaled(3000); it++){
        int scale = stress::rand_int(0, 2);
        double span = scale == 0 ? 1 : scale == 1 ? 1e4 : 1e12;
        double lo = (stress::rand_int(-1000000, 1000000) / 1e6) * span;
        double hi = lo + (stress::rand_int(0, 1000000) / 1e6) * span;
        check_golden_v(lo, hi, lo + (hi - lo) * stress::rand_int(0, 1000000) / 1e6);
    }

    /// Kink near 0 in a wide interval: interior points keep the rounding error from when the interval was wide
    for (long long it = 0; it < stress::scaled(3000); it++){
        double span = stress::rand_int(0, 1) ? 1e9 : 1e12;
        double c = stress::rand_int(-1000000, 1000000) / 1e9;
        double lo = c - (stress::rand_int(0, 1000000) / 1e6) * span;
        double hi = c + (stress::rand_int(0, 1000000) / 1e6) * span;
        check_golden_v(lo, hi, c);
    }

    for (long long it = 0; it < stress::scaled(500); it++) check_golden_smooth();

    return 0;
}
