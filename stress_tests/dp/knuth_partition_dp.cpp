#include "../common.h"

#define main divide_conquer_main
#include "../../code_library/dp/divide_conquer_dp.cpp"
#undef main

#define main library_main
#include "../../code_library/dp/knuth_partition_dp.cpp"
#undef main

const long long INF = numeric_limits<long long>::max();

/// Sum of Monge parts: non-negative pair weights, convex in length, squared non-negative sum, f(j) + g(i) of any sign
struct MongeCost{
    vector<vector<long long>> pairs;
    vector<long long> by_length, sums, f, g;

    MongeCost(int n, long long spread){
        pairs.assign(n + 1, vector<long long>(n + 1, 0));
        bool use_pairs = stress::rand_int(0, 2), sparse = stress::rand_int(0, 1);
        for (int x = 0; x < n; x++){
            for (int y = 0; y < n; y++){
                long long w = (use_pairs && x < y && (!sparse || stress::rand_int(0, 9) == 0)) ? stress::rand_int(0, spread) : 0;
                pairs[x + 1][y + 1] = pairs[x][y + 1] + pairs[x + 1][y] - pairs[x][y] + w;
            }
        }

        by_length.assign(n + 1, 0);
        long long step = stress::rand_int(-spread, spread);
        bool use_length = stress::rand_int(0, 2);
        for (int len = 1; len <= n; len++){
            step += use_length ? stress::rand_int(0, spread) : 0;
            by_length[len] = by_length[len - 1] + step;
        }

        sums.assign(n + 1, 0);
        bool use_squares = stress::rand_int(0, 1);
        for (int x = 0; x < n; x++) sums[x + 1] = sums[x] + (use_squares ? stress::rand_int(0, min(spread, 1000LL)) : 0);

        f.resize(n + 1), g.resize(n + 1);
        for (int i = 0; i <= n; i++) f[i] = stress::rand_int(-spread, spread), g[i] = stress::rand_int(-spread, spread);
    }

    long long operator()(int j, int i) const{
        long long inside = pairs[i][i] - pairs[j][i] - pairs[i][j] + pairs[j][j];
        long long sum = sums[i] - sums[j];
        return inside + by_length[i - j] + sum * sum + f[j] + g[i];
    }
};

/// layers[g][i] by the naive O(k n^2) recurrence
vector<vector<long long>> brute(int n, int k, const MongeCost& cost){
    vector<vector<long long>> layers(k + 1, vector<long long>(n + 1, INF));
    layers[0][0] = 0;

    for (int g = 1; g <= k; g++){
        for (int i = 1; i <= n; i++){
            for (int j = 0; j < i; j++){
                if (layers[g - 1][j] != INF) layers[g][i] = min(layers[g][i], layers[g - 1][j] + cost(j, i));
            }
        }
    }

    return layers;
}

/// Asserts the call domain and the O(n^2 + n k) bound, which a scan ignoring either opt bound exceeds
vector<long long> checked_run(int n, int k, const MongeCost& cost){
    long long calls = 0;
    auto counted = [&](int j, int i){
        assert(0 <= j && j < i && i <= n);
        calls++;
        return cost(j, i);
    };

    auto dp = knuth_partition_dp<long long>(n, k, counted);
    assert(calls <= 1LL * n * n + 1LL * n * k);

    return dp;
}

int main(){
    for (long long it = 0; it < stress::scaled(1200); it++){
        int n = it < 600 ? it % 9 : stress::rand_int(9, 40);
        long long spread = it % 3 ? 3 : 1000000;
        MongeCost cost(n, spread);
        auto layers = brute(n, n, cost);

        for (int k = 0; k <= n + 1; k++){
            auto expected = k <= n ? layers[k] : vector<long long>(n + 1, INF);
            auto dp = checked_run(n, k, cost);
            assert(dp == expected);
            assert(dp == divide_conquer_dp<long long>(n, k, cost));
        }
    }

    for (long long it = 0; it < stress::scaled(12); it++){
        int n = it % 2 ? stress::rand_int(300, 500) : stress::rand_int(100, 200);
        int k = it % 3 == 0 ? stress::rand_int(1, 6) : it % 3 == 1 ? stress::rand_int(1, n) : stress::rand_int(n - 5, n);
        MongeCost cost(n, it % 3 ? 3 : 1000000);
        auto layers = brute(n, k, cost);

        auto dp = checked_run(n, k, cost);
        assert(dp == layers[k]);
        assert(dp == divide_conquer_dp<long long>(n, k, cost));
    }

    return 0;
}
