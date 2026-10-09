#include "../common.h"

#define main library_main
#include "../../code_library/algebra/min_plus_convolution.cpp"
#undef main

const long long LIMIT = (1LL << 62) - 1;

/// Every comparison spends one unit of a shared budget, so a quadratic regression aborts early instead of timing out
struct Counted{
    long long value = 0;
    long long* budget = nullptr;

    Counted operator+(const Counted& other) const{
        return {value + other.value, budget};
    }

    bool operator<(const Counted& other) const{
        assert(--*budget >= 0);
        return value < other.value;
    }
};

vector<long long> brute(const vector<long long>& a, const vector<long long>& b, bool maximize){
    if (a.empty() || b.empty()) return {};

    vector<__int128> c(a.size() + b.size() - 1);
    vector<bool> seen(c.size());
    for (size_t i = 0; i < a.size(); i++){
        for (size_t j = 0; j < b.size(); j++){
            __int128 s = (__int128)a[i] + b[j];
            if (!seen[i + j] || (maximize ? s > c[i + j] : s < c[i + j])) c[i + j] = s, seen[i + j] = true;
        }
    }

    vector<long long> res;
    for (auto x : c){
        assert(x >= LLONG_MIN && x <= LLONG_MAX);
        res.push_back((long long)x);
    }
    return res;
}

bool in_range(const vector<long long>& v){
    for (auto x : v) if (x < -LIMIT || x > LIMIT) return false;
    return true;
}

bool is_convex(const vector<long long>& v){
    for (size_t i = 1; i + 1 < v.size(); i++){
        if ((__int128)v[i + 1] - v[i] < (__int128)v[i] - v[i - 1]) return false;
    }
    return true;
}

/// Sorted random differences from [-d, d], d halved until every value is within the 2^62 claim
vector<long long> random_convex(int n, long long d){
    while (true){
        vector<long long> diffs(max(n - 1, 0));
        for (auto& x : diffs) x = stress::rand_int(-d, d);
        sort(diffs.begin(), diffs.end());

        vector<__int128> v(n);
        bool fits = true;
        for (int i = 0; i < n; i++){
            v[i] = i ? v[i - 1] + diffs[i - 1] : (__int128)stress::rand_int(-d, d);
            fits = fits && v[i] >= -LIMIT && v[i] <= LIMIT;
        }
        if (fits) return vector<long long>(v.begin(), v.end());
        d /= 2;
    }
}

vector<long long> random_any(int n, long long r){
    vector<long long> v(n);
    for (auto& x : v) x = stress::rand_int(-r, r);
    return v;
}

vector<long long> negated(vector<long long> v){
    for (auto& x : v) x = -x;
    return v;
}

void check(const vector<long long>& a, const vector<long long>& b, const vector<long long>& arbitrary){
    auto expected = brute(a, b, false);
    assert(min_plus_convex_convex(a, b) == expected);
    assert(min_plus_convex_arbitrary(a, b) == expected);
    assert(min_plus_convex_arbitrary(b, a) == expected);

    assert(max_plus_concave_concave(negated(a), negated(b)) == brute(negated(a), negated(b), true));
    assert(min_plus_convex_arbitrary(a, arbitrary) == brute(a, arbitrary, false));
}

/// Every sequence of up to 3 values drawn from the claim boundary that is convex
void check_extremes(){
    const vector<long long> pool = {-LIMIT, -LIMIT + 1, -1, 0, 1, LIMIT - 1, LIMIT};
    vector<vector<long long>> convex;
    for (int n = 1; n <= 3; n++){
        int total = 1;
        for (int i = 0; i < n; i++) total *= pool.size();
        for (int mask = 0; mask < total; mask++){
            vector<long long> v;
            for (int i = 0, x = mask; i < n; i++, x /= pool.size()) v.push_back(pool[x % pool.size()]);
            if (is_convex(v)) convex.push_back(v);
        }
    }

    for (const auto& a : convex){
        for (const auto& b : convex) check(a, b, b);
    }
}

/// Each recursion level scans at most m + (nodes on the level) candidates, over log2(n + m) + 1 levels
void check_cost(int n, int m){
    long long budget = 2LL * (n + m) * (__lg(n + m) + 2);
    vector<Counted> a, b;
    for (auto x : random_convex(n, 1000000)) a.push_back({x, &budget});
    for (auto x : random_any(m, 1000000000)) b.push_back({x, &budget});

    assert((int)min_plus_convex_arbitrary(a, b).size() == n + m - 1);
}

int main(){
    check_extremes();

    const long long spreads[] = {0, 1, 3, 1000, 1000000000, LIMIT};
    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(0, it % 50 ? 12 : 400), m = stress::rand_int(0, it % 50 ? 12 : 400);
        long long d = spreads[stress::rand_int(0, 5)], r = spreads[stress::rand_int(0, 5)];

        auto a = random_convex(n, d), b = random_convex(m, d);
        assert(is_convex(a) && is_convex(b) && in_range(a) && in_range(b));
        check(a, b, random_any(m, r));
    }

    check_cost(200000, 200000);
    check_cost(1, 400000);
    check_cost(400000, 1);
    check_cost(2, 400000);

    return 0;
}
