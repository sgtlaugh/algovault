/***
 *
 * Bounded Knapsack
 * Maximum total value of items with total weight at most capacity, item i usable up to counts[i] times
 *
 * Complexity: O(n * capacity), O(capacity) memory, independent of the counts
 *
 * bounded_knapsack(values, weights, counts, capacity) returns the best total value, 0 for no items
 * values, weights and counts share one integer type: vector<int> or vector<long long> for all three
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

template <typename T = long long>
long long bounded_knapsack(const vector<T>& values, const vector<T>& weights, const vector<T>& counts, int capacity){
    vector<long long> dp(capacity + 1, 0);
    vector<pair<int, long long>> window(capacity + 1);

    for (size_t i = 0; i < values.size(); i++){
        long long v = values[i], c = counts[i];
        /// Skipping items that never help also keeps x += w and dp[x] - k * v below from overflowing
        if (v <= 0 || weights[i] > capacity) continue;
        int w = weights[i];
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
    int n, W;
    scanf("%d %d", &n, &W);
    vector<int> values(n), weights(n), counts(n);
    for (int i = 0; i < n; i++) scanf("%d %d %d", &values[i], &weights[i], &counts[i]);
    printf("%lld\n", bounded_knapsack(values, weights, counts, W));
    return 0;
}
