#include "../common.h"

#define main library_main
#include "../../code_library/dp/bounded_knapsack.cpp"
#undef main

/// Each item expanded into min(count, capacity / weight) copies, then plain 0/1 knapsack
long long expanded_knapsack(const vector<long long>& values, const vector<int>& weights, const vector<long long>& counts, int capacity){
    vector<long long> dp(capacity + 1, 0);
    for (size_t i = 0; i < values.size(); i++){
        long long copies = weights[i] == 0 ? counts[i] : min(counts[i], (long long)capacity / weights[i]);
        for (long long t = 0; t < copies; t++){
            for (int x = capacity; x >= weights[i]; x--) dp[x] = max(dp[x], dp[x - weights[i]] + values[i]);
        }
    }
    return dp[capacity];
}

/// The original source (ShahjalalShohag code-library), positive weights only
int shohag_bounded_knapsack(vector<int> ps, vector<int> ws, vector<int> ms, int W){
    int n = ps.size();
    vector<vector<int>> dp(n + 1, vector<int>(W + 1));
    for (int i = 0; i < n; ++i){
        for (int s = 0; s < ws[i]; ++s){
            int alpha = 0;
            queue<int> que;
            deque<int> peek;
            for (int w = s; w <= W; w += ws[i]){
                alpha += ps[i];
                int a = dp[i][w] - alpha;
                que.push(a);
                while (!peek.empty() && peek.back() < a) peek.pop_back();
                peek.push_back(a);
                while ((int)que.size() > ms[i] + 1){
                    if (que.front() == peek.front()) peek.pop_front();
                    que.pop();
                }
                dp[i + 1][w] = peek.front() + alpha;
            }
        }
    }
    int ans = 0;
    for (int w = 0; w <= W; ++w) ans = max(ans, dp[n][w]);
    return ans;
}

int main(){
    for (long long it = 0; it < stress::scaled(3000); it++){
        bool big = it >= 2500;
        int n = stress::rand_int(0, big ? 40 : 5), capacity = stress::rand_int(0, big ? 2000 : 20);
        int max_weight = stress::rand_int(1, big ? 300 : 12);
        bool positive = stress::rand_int(0, 1), huge = stress::rand_int(0, 3) == 0;
        /// At most capacity positive-weight copies plus 6 per zero-weight item are taken, so the answer fits
        long long max_value = huge ? LLONG_MAX / (capacity + 6 * n + 1) : 1000;
        vector<long long> values(n), counts(n);
        vector<int> weights(n);
        for (int i = 0; i < n; i++){
            values[i] = positive ? stress::rand_int(1, max_value) : stress::rand_int(huge ? -max_value : -20, max_value);
            /// Negative values carry no bound in the header since they are never taken
            if (!positive && huge && stress::rand_int(0, 3) == 0) values[i] = stress::rand_int(LLONG_MIN, -1);
            weights[i] = positive ? stress::rand_int(1, max_weight) : stress::rand_int(0, max_weight);
            if (!positive && stress::rand_int(0, 9) == 0) weights[i] = INT_MAX - stress::rand_int(0, 5);
            counts[i] = stress::rand_int(0, 5) == 0 ? stress::rand_int(0, (long long)1e18) : stress::rand_int(0, 6);
            if (weights[i] == 0) counts[i] = stress::rand_int(0, 6);
        }

        long long got = bounded_knapsack(values, weights, counts, capacity);
        assert(got == expanded_knapsack(values, weights, counts, capacity));

        if (positive && !huge){
            vector<int> ps(values.begin(), values.end()), ms(n);
            for (int i = 0; i < n; i++) ms[i] = min(counts[i], (long long)capacity);
            assert(got == shohag_bounded_knapsack(ps, weights, ms, capacity));
        }
    }

    return 0;
}
