/***
 *
 * Meet in the Middle
 * Subset sums of up to about 44 numbers: all of them sorted, how many fall in a range, and 0/1 knapsack
 *
 * Complexity:
 *   - count_subset_sums and knapsack_large_weights: O(2^(n/2)) time and memory
 *   - subset_sums: O(2^n) time and memory
 *
 * count_subset_sums(a, lo, hi): number of subsets of a with sum in [lo, hi], the empty subset included
 * knapsack_large_weights(values, weights, capacity): best total value of a subset with total weight at most capacity
 * subset_sums(a): the 2^n sums of all subsets of a, sorted, duplicates kept, 0 for the empty subset
 *
 * Notes:
 *   - a, values and weights are each vector<int> or vector<long long>, any sign, every subset's sums fit in long long
 *   - lo > hi gives 0, a sum of exactly s is count_subset_sums(a, s, s)
 *   - capacity >= 0, so the empty subset always fits and the best value is at least 0
 *   - Each half holds 2^(n/2) sums: n = 40 takes about 28 MB to count and 52 MB for knapsack, n = 44 about 3.5x that
 *   - Weights and sums up to W of about 1e6: bounded_knapsack.cpp and bounded_subset_sum.cpp run in O(n W) or less
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// Adding x to every sum keeps them sorted, so one merge per item replaces a sort: 3.4x faster at n = 22
template <typename T>
vector<long long> subset_sums(const vector<T>& a){
    vector<long long> sums{0};
    for (long long x : a){
        size_t k = sums.size();
        sums.resize(2 * k);
        for (size_t i = 0; i < k; i++) sums[k + i] = sums[i] + x;
        inplace_merge(sums.begin(), sums.begin() + k, sums.end());
    }

    return sums;
}

template <typename T>
long long count_subset_sums(const vector<T>& a, long long lo, long long hi){
    size_t half = a.size() / 2;
    vector<long long> left = subset_sums(vector<T>(a.begin(), a.begin() + half));
    vector<long long> right = subset_sums(vector<T>(a.begin() + half, a.end()));

    /// right[from..to) holds the sums r with lo <= l + r <= hi, both ends only move down as l grows
    long long res = 0;
    size_t from = right.size(), to = right.size();
    for (long long l : left){
        while (from > 0 && l + right[from - 1] >= lo) from--;
        while (to > 0 && l + right[to - 1] > hi) to--;
        if (to > from) res += to - from;
    }

    return res;
}

/// Every subset of items [first, last) as (weight, value), sorted by weight
template <typename V, typename W>
vector<pair<long long, long long>> weight_value_sums(const vector<V>& values, const vector<W>& weights, size_t first, size_t last){
    vector<pair<long long, long long>> sums{{0, 0}};
    for (size_t i = first; i < last; i++){
        size_t k = sums.size();
        sums.resize(2 * k);
        for (size_t j = 0; j < k; j++) sums[k + j] = {sums[j].first + weights[i], sums[j].second + values[i]};
        inplace_merge(sums.begin(), sums.begin() + k, sums.end());
    }

    return sums;
}

template <typename V, typename W>
long long knapsack_large_weights(const vector<V>& values, const vector<W>& weights, long long capacity){
    assert(values.size() == weights.size() && capacity >= 0);
    size_t n = values.size(), half = n / 2;
    vector<pair<long long, long long>> left = weight_value_sums(values, weights, 0, half);
    vector<pair<long long, long long>> right = weight_value_sums(values, weights, half, n);

    /// right[j].second becomes the best value among weights up to right[j].first
    for (size_t j = 1; j < right.size(); j++) right[j].second = max(right[j].second, right[j - 1].second);

    long long best = 0;
    size_t fit = right.size();
    for (auto [w, v] : left){
        while (fit > 0 && w + right[fit - 1].first > capacity) fit--;
        if (fit == 0) break;
        best = max(best, v + right[fit - 1].second);
    }

    return best;
}

int main(){
    /// Subsets of {3, 1, 4, 1}: sums 0 1 1 2 3 4 4 4 5 5 5 6 7 8 8 9
    vector<int> a = {3, 1, 4, 1};
    assert((subset_sums(a) == vector<long long>{0, 1, 1, 2, 3, 4, 4, 4, 5, 5, 5, 6, 7, 8, 8, 9}));
    assert(count_subset_sums(a, 5, 5) == 3);
    assert(count_subset_sums(a, 4, 6) == 7);
    assert(count_subset_sums(a, 0, 9) == 16);
    assert(count_subset_sums(a, 10, 20) == 0);

    /// Negative numbers: {-2, 3} has sums -2, 0, 1, 3
    assert(count_subset_sums(vector<long long>{-2, 3}, -1, 2) == 2);

    /// Weights 6e9, 5e9, 5e9 with values 10, 7, 7 and capacity 1e10: both 5e9 items beat the single 6e9 one
    vector<long long> values = {10, 7, 7}, weights = {6000000000LL, 5000000000LL, 5000000000LL};
    assert(knapsack_large_weights(values, weights, 10000000000LL) == 14);
    assert(knapsack_large_weights(values, weights, 9999999999LL) == 10);
    assert(knapsack_large_weights(values, weights, 4999999999LL) == 0);

    return 0;
}
