#include "../common.h"

#define main library_main
#include "../../code_library/dp/meet_in_the_middle.cpp"
#undef main

vector<long long> brute_sums(const vector<long long>& a){
    int n = a.size();
    vector<long long> sums;
    for (int mask = 0; mask < (1 << n); mask++){
        long long s = 0;
        for (int i = 0; i < n; i++) if (mask >> i & 1) s += a[i];
        sums.push_back(s);
    }
    sort(sums.begin(), sums.end());

    return sums;
}

long long brute_knapsack(const vector<long long>& values, const vector<long long>& weights, long long capacity){
    int n = values.size();
    long long best = 0;
    for (int mask = 0; mask < (1 << n); mask++){
        long long w = 0, v = 0;
        for (int i = 0; i < n; i++) if (mask >> i & 1) w += weights[i], v += values[i];
        if (w <= capacity) best = max(best, v);
    }

    return best;
}

vector<long long> random_vector(int n, long long lo, long long hi){
    vector<long long> a(n);
    for (long long& x : a) x = stress::rand_int(lo, hi);

    return a;
}

void check_small(const vector<long long>& a, const vector<long long>& values){
    vector<long long> sums = brute_sums(a);
    assert(subset_sums(a) == sums);

    /// Ranges ending on, between and beyond the sums, including lo > hi and the extreme long longs
    vector<long long> ends = {LLONG_MIN, LLONG_MAX, sums.front() - 1, sums.back() + 1};
    for (int it = 0; it < 4; it++){
        long long s = sums[stress::rand_int(0, sums.size() - 1)];
        ends.push_back(s), ends.push_back(s + stress::rand_int(-1, 1));
    }
    for (long long lo : ends){
        for (long long hi : ends){
            long long expected = upper_bound(sums.begin(), sums.end(), hi) - lower_bound(sums.begin(), sums.end(), lo);
            assert(count_subset_sums(a, lo, hi) == max(0LL, expected));
        }
    }

    /// a doubles as the weights, so negative weights are covered too
    long long top = 0;
    for (long long w : a) top += abs(w);
    bool int_values = all_of(values.begin(), values.end(), [](long long v){ return INT_MIN <= v && v <= INT_MAX; });
    for (int it = 0; it < 4; it++){
        long long capacity = it == 0 ? 0 : stress::rand_int(0, top);
        long long best = brute_knapsack(values, a, capacity);
        assert(knapsack_large_weights(values, a, capacity) == best);
        if (int_values) assert(knapsack_large_weights(vector<int>(values.begin(), values.end()), a, capacity) == best);
    }
}

int main(){
    /// Every vector over {-1, 0, 1, 2} up to length 5: duplicates, zeros and negatives in each half
    for (int n = 0; n <= 5; n++){
        int total = 1;
        for (int i = 0; i < n; i++) total *= 4;
        for (int code = 0; code < total; code++){
            vector<long long> a(n);
            for (int i = 0, c = code; i < n; i++, c /= 4) a[i] = c % 4 - 1;
            check_small(a, random_vector(n, -3, 5));
        }
    }

    for (int it = 0; it < stress::scaled(400); it++){
        int n = stress::rand_int(0, 14);
        long long big = it % 3 == 0 ? 10 : it % 3 == 1 ? 1000000 : (long long)1e17;
        bool negative = it % 2;
        check_small(random_vector(n, negative ? -big : 0, big), random_vector(n, negative ? -big : 0, big));
    }

    /// n = 36 and 35 with small numbers against DP over sums: counts reach 2^36, beyond 32 bits
    for (int it = 0; it < stress::scaled(2); it++){
        int n = 36 - it % 2;
        vector<long long> a = random_vector(n, 0, 60), values = random_vector(n, -5, 1000000000);
        long long total = accumulate(a.begin(), a.end(), 0LL);

        vector<long long> ways(total + 1, 0), best(total + 1, LLONG_MIN);
        ways[0] = 1, best[0] = 0;
        for (int i = 0; i < n; i++){
            for (long long s = total; s >= a[i]; s--){
                ways[s] += ways[s - a[i]];
                if (best[s - a[i]] != LLONG_MIN) best[s] = max(best[s], best[s - a[i]] + values[i]);
            }
        }

        assert(count_subset_sums(a, 0, total) == 1LL << n);
        for (int q = 0; q < 5; q++){
            long long lo = stress::rand_int(0, total), hi = stress::rand_int(lo, total), expected = 0;
            for (long long s = lo; s <= hi; s++) expected += ways[s];
            assert(count_subset_sums(a, lo, hi) == expected);

            long long capacity = stress::rand_int(0, total), optimum = 0;
            for (long long s = 0; s <= capacity; s++) optimum = max(optimum, best[s]);
            assert(knapsack_large_weights(values, a, capacity) == optimum);
        }
    }

    stress::expect_assert("dp/meet_in_the_middle.cpp", []{ knapsack_large_weights(vector<int>{1}, vector<int>{1}, -1); });

    return 0;
}
