#include "../common.h"

#define main library_main
#include "../../code_library/dp/slope_trick.cpp"
#undef main

/// f as an array over x in [-D, D]; every breakpoint stays inside |x| < C (see random_ops), so after each operation
/// the cells beyond |x| > C are rebuilt as the straight-line continuation, exact for the shift windows to read
const int C = 100, D = C + 3;

struct Brute{
    vector<long long> f = vector<long long>(2 * D + 1, 0);

    template <typename F>
    void add(F g){
        for (int x = -D; x <= D; x++) f[x + D] += g(x);
        extend();
    }

    void prefix_min(){
        for (int i = 1; i <= 2 * D; i++) f[i] = min(f[i], f[i - 1]);
        extend();
    }

    void shift(int lo, int hi){
        vector<long long> g = f;
        for (int x = -C; x <= C; x++){
            g[x + D] = LLONG_MAX;
            for (int y = x - hi; y <= x - lo; y++) g[x + D] = min(g[x + D], f[y + D]);
        }

        f = g;
        extend();
    }

    void suffix_min(){
        for (int i = 2 * D - 1; i >= 0; i--) f[i] = min(f[i], f[i + 1]);
        extend();
    }

private:
    void extend(){
        long long left_slope = f[-C + 1 + D] - f[-C + D], right_slope = f[C + D] - f[C - 1 + D];
        for (int x = C + 1; x <= D; x++){
            f[x + D] = f[x - 1 + D] + right_slope;
            f[-x + D] = f[-x + 1 + D] - left_slope;
        }
    }
};

/// breakpoints are added in [-10, 10] and at most 20 shifts of |lo|, |hi| <= 3 move them, so they stay inside [-70, 70]
void random_ops(){
    SlopeTrick<long long> f;
    Brute b;
    int ops = stress::rand_int(0, 40), shifts = 0;

    for (int op = 0; op < ops; op++){
        int type = stress::rand_int(0, 7);
        long long a = stress::rand_int(-10, 10);
        if (type >= 6 && shifts == 20) type = 0;

        if (type == 0) f.add_abs(a), b.add([&](long long x){ return abs(x - a); });
        if (type == 1) f.add_x_minus_a(a), b.add([&](long long x){ return max(0LL, x - a); });
        if (type == 2) f.add_a_minus_x(a), b.add([&](long long x){ return max(0LL, a - x); });
        if (type == 3) f.add_const(a), b.add([&](long long){ return a; });
        if (type == 4) f.prefix_min(), b.prefix_min();
        if (type == 5) f.suffix_min(), b.suffix_min();
        if (type == 6){
            int lo = stress::rand_int(-3, 3), hi = stress::rand_int(lo, 3);
            f.shift(lo, hi), b.shift(lo, hi), shifts++;
        }
        if (type == 7){
            int d = stress::rand_int(-3, 3);
            f.shift(d), b.shift(d, d), shifts++;
        }

        long long mn = LLONG_MAX;
        for (int x = -C; x <= C; x++){
            assert(f.eval(x) == b.f[x + D]);
            mn = min(mn, b.f[x + D]);
        }
        assert(f.get_min() == mn);

        long long l = -C, r = C;
        while (b.f[l + D] != mn) l++;
        while (b.f[r + D] != mn) r--;
        auto [al, ar] = f.argmin();
        assert(al == (l == -C ? -SlopeTrick<long long>::INF : l));
        assert(ar == (r == C ? SlopeTrick<long long>::INF : r));
    }
}

/// SlopeTrick<int> against SlopeTrick<long long> (checked by random_ops) with contest-sized values: breakpoints up to
/// +-1e9 sit outside INT_MAX / 3, so a hinge added while the opposite heap is empty must not treat the sentinel as a breakpoint;
/// shifted runs keep breakpoints within 5e8 and both accumulated shifts within 5e8, so every breakpoint stays within the 1e9 bound
void int_ops(){
    const long long E9 = 1000000000, LL_INF = SlopeTrick<long long>::INF, INT_INF = SlopeTrick<int>::INF;
    SlopeTrick<int> f;
    SlopeTrick<long long> g;
    bool shifting = stress::rand_int(0, 1);
    int ops = stress::rand_int(1, 8), shifts = 0;

    for (int op = 0; op < ops; op++){
        int type = stress::rand_int(0, shifting ? 6 : 5), pick = stress::rand_int(0, 3);
        long long bound = shifting ? E9 / 2 : E9, a = pick == 0 ? -bound : pick == 1 ? bound : stress::rand_int(-bound, bound);
        if (type == 6 && shifts == 2) type = 0;

        SlopeTrick<long long> next = g;
        if (type == 0) next.add_abs(a);
        if (type == 1) next.add_x_minus_a(a);
        if (type == 2) next.add_a_minus_x(a);
        if (type == 3) next.prefix_min();
        if (type == 4) next.suffix_min();
        if (type == 5) next.add_const(a);
        long long lo = stress::rand_int(-E9 / 4, E9 / 4), hi = stress::rand_int(lo, E9 / 4);
        if (type == 6) next.shift(lo, hi), shifts++;
        if (next.get_min() < INT_MIN || next.get_min() > INT_MAX) return;
        g = next;

        if (type == 0) f.add_abs(a);
        if (type == 1) f.add_x_minus_a(a);
        if (type == 2) f.add_a_minus_x(a);
        if (type == 3) f.prefix_min();
        if (type == 4) f.suffix_min();
        if (type == 5) f.add_const(a);
        if (type == 6) f.shift(lo, hi);

        assert(f.get_min() == g.get_min());
        auto [l, r] = g.argmin();
        auto [fl, fr] = f.argmin();
        assert(fl <= fr && fl == (l == -LL_INF ? -INT_INF : l) && fr == (r == LL_INF ? INT_INF : r));
        for (long long x : {-E9, -E9 / 2, 0LL, a, E9 / 2, E9, l, r}){
            if (abs(x) > E9) continue;
            long long want = g.eval(x);
            if (want >= INT_MIN && want <= INT_MAX) assert(f.eval(x) == want);
        }
    }
}

/// the optimal b only takes values from a, so a DP over the sorted values is exact
long long brute_non_decreasing(const vector<long long>& a){
    vector<long long> v = a;
    sort(v.begin(), v.end());
    v.erase(unique(v.begin(), v.end()), v.end());

    vector<long long> dp(v.size(), 0);
    for (long long x : a){
        long long best = LLONG_MAX;
        for (size_t j = 0; j < v.size(); j++){
            best = min(best, dp[j]);
            dp[j] = best + abs(x - v[j]);
        }
    }
    return a.empty() ? 0 : *min_element(dp.begin(), dp.end());
}

int main(){
    {
        const int E9 = 1000000000;
        SlopeTrick<int> f, g, h;
        f.add_a_minus_x(E9);
        assert(f.get_min() == 0 && f.eval(E9) == 0 && f.eval(0) == E9 && f.argmin() == make_pair(E9, SlopeTrick<int>::INF));
        g.add_x_minus_a(-E9);
        assert(g.get_min() == 0 && g.argmin() == make_pair(-SlopeTrick<int>::INF, -E9));
        h.add_abs(5);
        h.prefix_min();
        h.add_a_minus_x(900000000);
        assert(h.get_min() == 0 && h.argmin() == make_pair(900000000, SlopeTrick<int>::INF));
    }
    for (long long it = 0; it < stress::scaled(3000); it++) int_ops();

    for (long long it = 0; it < stress::scaled(1500); it++) random_ops();

    for (long long it = 0; it < stress::scaled(1500); it++){
        int n = stress::rand_int(0, it % 10 ? 12 : 300);
        long long spread = it % 3 ? 5 : 1000000000;
        vector<long long> a(n);
        for (auto& x : a) x = stress::rand_int(-spread, spread);
        assert(min_cost_non_decreasing(a) == brute_non_decreasing(a));
    }

    /// header range boundary: breakpoints at +-1e18, n <= 4 keeps the answer below 8e18
    const long long BIG = 1e18;
    for (long long it = 0; it < stress::scaled(500); it++){
        int n = stress::rand_int(1, 4);
        vector<long long> a(n);
        for (auto& x : a){
            int pick = stress::rand_int(0, 3);
            x = pick == 0 ? -BIG : pick == 1 ? BIG : stress::rand_int(-BIG, BIG);
        }
        assert(min_cost_non_decreasing(a) == brute_non_decreasing(a));
    }

    /// shift(0, 1e18) leaves a non-increasing f unchanged and shift(-1e18, 0) a non-decreasing one, so these compute
    /// the non-decreasing and non-increasing costs; more than 9 rounds overflow unless prefix_min / suffix_min reset the offset
    for (long long it = 0; it < stress::scaled(500); it++){
        int n = stress::rand_int(1, 30);
        vector<long long> a(n);
        for (auto& x : a) x = stress::rand_int(-1000000000, 1000000000);

        SlopeTrick<long long> inc, dec;
        for (long long x : a){
            inc.prefix_min();
            inc.shift(0, BIG);
            inc.add_abs(x);
            dec.suffix_min();
            dec.shift(-BIG, 0);
            dec.add_abs(x);
        }

        assert(inc.get_min() == brute_non_decreasing(a));
        reverse(a.begin(), a.end());
        assert(dec.get_min() == brute_non_decreasing(a));
    }

    return 0;
}
