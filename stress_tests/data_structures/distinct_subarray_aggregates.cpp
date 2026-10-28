#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/distinct_subarray_aggregates.cpp"
#undef main

/// Segments after the push of a[i] against f[l] = op(a[l..i]) for every l <= i
template <typename T, typename Segments>
void verify(const Segments& segs, const vector<T>& f, int i, size_t max_segments){
    assert(!segs.empty() && segs.size() <= max_segments);
    assert(segs[0].left == 0 && segs.back().right == i);

    set<T> seen;
    for (size_t k = 0; k < segs.size(); k++){
        auto& s = segs[k];
        assert(s.left <= s.right && (k == 0 || s.left == segs[k - 1].right + 1));
        assert(seen.insert(s.value).second);
        for (int l = s.left; l <= s.right; l++) assert(f[l] == s.value);
    }
}

/// Runs verify after every push and compares the O(n^2) map of all subarray values
template <typename T, typename Op>
void check(const vector<T>& a, Op op, size_t max_segments){
    int n = a.size();
    DistinctSubarrayAggregates<T, Op> agg(op);
    map<T, long long> expected_cnt, cnt;
    vector<T> f(n);

    for (int i = 0; i < n; i++){
        agg.push(a[i]);

        f[i] = a[i];
        for (int l = i - 1; l >= 0; l--) f[l] = op(a[l], f[l + 1]);
        for (int l = 0; l <= i; l++) expected_cnt[f[l]]++;

        verify(agg.segments(), f, i, max_segments);
        for (auto& s : agg.segments()) cnt[s.value] += s.right - s.left + 1;
    }

    assert(cnt == expected_cnt);
}

/// Every array of length at most len over alphabet, each one pushed once onto a copy of its parent prefix
template <typename T, typename Op, typename Bound>
void exhaustive(const vector<T>& alphabet, int len, Op op, Bound max_segments){
    vector<T> a;
    auto dfs = [&](auto& self, const DistinctSubarrayAggregates<T, Op>& agg) -> void{
        for (T x : alphabet){
            auto next = agg;
            next.push(x);
            a.push_back(x);

            int i = a.size() - 1;
            vector<T> f(i + 1);
            f[i] = x;
            for (int l = i - 1; l >= 0; l--) f[l] = op(a[l], f[l + 1]);
            verify(next.segments(), f, i, max_segments(a));

            if (i + 1 < len) self(self, next);
            a.pop_back();
        }
    };
    dfs(dfs, DistinctSubarrayAggregates<T, Op>(op));
}

size_t gcd_bound(const vector<long long>& a){
    long long A = *max_element(a.begin(), a.end());
    return A == 0 ? 1 : __lg(A) + 2;
}

long long smooth(long long limit){
    static const int primes[] = {2, 3, 5, 7};
    long long x = 1;
    while (stress::rand_int(0, 5)){
        int p = primes[stress::rand_int(0, 3)];
        if (x > limit / p) break;
        x *= p;
    }

    return x;
}

int main(){
    auto gcd_op = [](long long x, long long y){ return gcd(x, y); };

    /// Every array of length at most 6
    exhaustive(vector<long long>{0, 1, 2, 3, 4, 6}, 6, gcd_op, gcd_bound);
    vector<unsigned long long> three_bits = {0, 1, 2, 3, 4, 5, 6, 7};
    exhaustive(three_bits, 6, bit_or<unsigned long long>(), [](auto&){ return size_t(4); });
    exhaustive(three_bits, 6, bit_and<unsigned long long>(), [](auto&){ return size_t(4); });

    for (long long it = 0; it < stress::scaled(600); it++){
        int n = it < 400 ? it % 40 + 1 : stress::rand_int(1, 300);
        int kind = it % 4;

        vector<long long> g(n);
        long long A = vector<long long>{6, 1000, 1000000000000000000LL, 1LL << 62}[kind];
        for (auto& x : g){
            int r = stress::rand_int(0, 9);
            x = r == 0 ? 0 : r < 7 ? smooth(A) : stress::rand_int(0, A);
        }
        if (it % 7 == 0) for (int i = 0; i < n; i++) g[i] = 1LL << max(0, 62 - (n - 1 - i));
        check(g, gcd_op, gcd_bound(g));

        int bits = vector<int>{3, 10, 40, 64}[kind];
        vector<unsigned long long> u(n);
        for (auto& x : u){
            x = 0;
            for (int b = 0; b < bits; b++) if (stress::rand_int(0, 7) == 0) x |= 1ULL << b;
        }
        check(u, bit_or<unsigned long long>(), bits + 1);
        for (auto& x : u) x = bits == 64 ? ~x : x ^ ((1ULL << bits) - 1);
        check(u, bit_and<unsigned long long>(), bits + 1);

        vector<long long> s(n);
        for (auto& x : s) x = stress::rand_int(0, 3) ? -stress::rand_int(1, 1LL << 20) : stress::rand_int(LLONG_MIN, LLONG_MAX);
        check(s, bit_and<long long>(), 65);
        check(s, bit_or<long long>(), 65);
    }

    vector<long long> ladder;
    for (int b = 0; b <= 62; b++) ladder.push_back(1LL << b);
    ladder.push_back(0);
    check(ladder, gcd_op, 64);

    /// 65 segments, then a push that merges while the list is at its widest
    vector<unsigned long long> bits_up;
    for (int b = 0; b < 64; b++) bits_up.push_back(1ULL << b);
    bits_up.push_back(0), bits_up.push_back(1);
    check(bits_up, bit_or<unsigned long long>(), 65);

    /// O(n log A) op calls at n = 2e5 with gcd chains as long as 1e18 allows
    long long calls = 0;
    auto counted = [&](long long x, long long y){ calls++; return gcd(x, y); };
    DistinctSubarrayAggregates<long long, decltype(counted)> agg(counted);
    int n = 200000;
    for (int i = 0; i < n; i++){
        agg.push(stress::rand_int(0, 3) ? smooth(1000000000000000000LL) : stress::rand_int(0, 1000000000000000000LL));
        assert(agg.segments().size() <= 61);
    }
    assert(calls <= 61LL * n);

    return 0;
}
