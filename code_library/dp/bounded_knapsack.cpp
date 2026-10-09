/***
 *
 * Bounded Knapsack
 * Maximum total value of items with total weight at most capacity, item i usable up to counts[i] times
 *
 * Complexity: O(n * capacity), O(capacity) memory, independent of the counts
 *
 * bounded_knapsack(values, weights, counts, capacity) returns the best total value, 0 for no items
 *
 * Capacity, weights and counts must be >= 0, values may be negative (such items are never taken)
 * The answer and every positive values[i] * (capacity / weights[i]) must fit in long long,
 * for zero-weight items values[i] * counts[i] must fit
 *
 * Per item, dp positions x = s, s + w, s + 2w, ... of one residue s mod w form a sequence where
 * new[x_k] = max over k - c <= j <= k of old[x_j] + (k - j) * v, a sliding window maximum of old[x_j] - j * v
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

long long bounded_knapsack(const vector<long long>& values, const vector<int>& weights, const vector<long long>& counts, int capacity){
    vector<long long> dp(capacity + 1, 0);
    vector<pair<int, long long>> window(capacity + 1);

    for (size_t i = 0; i < values.size(); i++){
        long long v = values[i], c = counts[i];
        int w = weights[i];
        /// Skipping items that never help also keeps x += w and dp[x] - k * v below from overflowing
        if (v <= 0 || w > capacity) continue;
        if (w == 0){
            for (auto& best : dp) best += v * c;
            continue;
        }

        /// In place is safe: position x reads old dp[x] before writing it, later positions only read old values
        for (int s = 0; s < w && s <= capacity; s++){
            int head = 0, tail = 0;
            for (int k = 0, x = s; x <= capacity; k++, x += w){
                long long cur = dp[x] - k * v;
                while (tail > head && window[tail - 1].second <= cur) tail--;
                window[tail++] = {k, cur};
                if (window[head].first < k - c) head++;
                dp[x] = window[head].second + k * v;
            }
        }
    }

    return dp[capacity];
}

int main(){
    assert(bounded_knapsack({3, 4, 1}, {2, 3, 1}, {2, 1, 5}, 7) == 10);
    assert(bounded_knapsack({10, 7}, {5, 4}, {1, 3}, 12) == 21);
    assert(bounded_knapsack({10, 7}, {5, 4}, {1, 1}, 12) == 17);
    assert(bounded_knapsack({60, 100, 120}, {10, 20, 30}, {1, 1, 1}, 50) == 220);

    assert(bounded_knapsack({}, {}, {}, 10) == 0);
    assert(bounded_knapsack({5}, {3}, {4}, 0) == 0);
    assert(bounded_knapsack({5}, {11}, {4}, 10) == 0);
    assert(bounded_knapsack({5}, {2}, {0}, 10) == 0);
    assert(bounded_knapsack({-5, 2}, {1, 3}, {10, 10}, 10) == 6);
    assert(bounded_knapsack({5, 2}, {0, 3}, {3, 10}, 7) == 19);
    assert(bounded_knapsack({-5}, {0}, {3}, 7) == 0);
    assert(bounded_knapsack({5}, {INT_MAX}, {3}, 10) == 0);
    assert(bounded_knapsack({(long long)4e18, -(long long)9e17}, {10, 1}, {1, 10}, 10) == (long long)4e18);

    assert(bounded_knapsack({1}, {1}, {(long long)1e18}, 100000) == 100000);
    assert(bounded_knapsack({(long long)1e13}, {1}, {(long long)1e18}, 100000) == (long long)1e18);
    assert(bounded_knapsack({LLONG_MAX / 100000}, {1}, {(long long)1e18}, 100000) == 9223372036854700000LL);
    assert(bounded_knapsack({(long long)1e9, 3}, {0, 2}, {(long long)9e9, 1}, 5) == 9000000000000000003LL);
    assert(bounded_knapsack({7, 3}, {3, 1}, {2, 100}, 1000) == 7 * 2 + 3 * 100);

    return 0;
}
