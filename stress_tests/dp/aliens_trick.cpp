#include "../common.h"

#define main library_main
#include "../../code_library/dp/aliens_trick.cpp"
#undef main

const long long INF = LLONG_MAX / 4;

enum Ties{FEWEST, MOST, HASHED, RANDOM};

/// RANDOM may report another count for the same lambda on every call, HASHED always the same one
template <typename T>
long long pick_tied(const vector<long long>& tied, Ties ties, T lambda, unsigned long long salt){
    if (ties == FEWEST) return tied.front();
    if (ties == MOST) return tied.back();
    if (ties == RANDOM) return tied[stress::rand_int(0, tied.size() - 1)];

    unsigned long long bits = 0;
    memcpy(&bits, &lambda, sizeof(lambda));
    return tied[(bits ^ salt) * 0x9E3779B97F4A7C15ULL % tied.size()];
}

/// Convex f on [first, n] given directly; the solver scans every count, so ties come from repeated slopes
void check_synthetic(Ties ties){
    int first = stress::rand_int(0, 3), n = first + stress::rand_int(0, 10);
    int range = stress::rand_int(1, 2) == 1 ? 3 : 1000;

    vector<long long> slopes;
    for (int c = first + 1; c <= n; c++) slopes.push_back(stress::rand_int(-range, range));
    sort(slopes.begin(), slopes.end());

    vector<long long> f(n + 1, INF);
    f[first] = stress::rand_int(-range, range) + (stress::rand_int(0, 1) ? 5000000000LL : 0);
    for (int c = first + 1; c <= n; c++) f[c] = f[c - 1] + slopes[c - first - 1];

    unsigned long long salt = stress::rand_int(0, LLONG_MAX);
    vector<long long> tied;
    auto solver = [&](long long lambda){
        long long best = LLONG_MAX;
        tied.clear();
        for (int c = first; c <= n; c++){
            long long value = f[c] + lambda * c;
            if (value < best){
                best = value;
                tied.clear();
            }
            if (value == best) tied.push_back(c);
        }
        return pair<long long, long long>{best, pick_tied(tied, ties, lambda, salt)};
    };

    long long lo = 0, hi = 0;
    if (!slopes.empty()) lo = -slopes.back(), hi = -slopes.front();
    if (stress::rand_int(0, 1)) lo -= stress::rand_int(0, 1000000000000LL), hi += stress::rand_int(0, 1000000000000LL);

    for (int k = first; k <= n; k++) assert(aliens_trick(k, lo, hi, solver) == f[k]);

    /// int bounds with an int64 solver must still return f(k) beyond 2^31
    if (lo == int(lo) && hi == int(hi) && hi - lo <= INT_MAX){
        for (int k = first; k <= n; k++) assert(aliens_trick(k, int(lo), int(hi), solver) == f[k]);
    }

    /// Shifting lambda by s solves f(c) - s c; bounds from INT_MIN catch a search that evaluates lambda = lo - 1
    if (hi - lo <= INT_MAX){
        long long shift = INT_MIN - lo;
        auto shifted = [&](long long lambda){ return solver(lambda - shift); };
        for (int k = first; k <= n; k++) assert(aliens_trick(k, INT_MIN, int(INT_MIN + (hi - lo)), shifted) == f[k] - shift * k);
    }
}

/// Split into k segments minimizing the sum of squared segment sums, against the O(n^2 k) DP
void check_squared_split(int n, int value_max, Ties ties){
    vector<long long> a(n), pre(n + 1, 0);
    for (auto& x : a) x = stress::rand_int(0, value_max);
    for (int i = 0; i < n; i++) pre[i + 1] = pre[i] + a[i];

    vector<vector<long long>> exact(n + 1, vector<long long>(n + 1, INF));
    exact[0][0] = 0;
    for (int c = 1; c <= n; c++){
        for (int j = 1; j <= n; j++){
            for (int i = 0; i < j; i++){
                if (exact[c - 1][i] == INF) continue;
                long long s = pre[j] - pre[i];
                exact[c][j] = min(exact[c][j], exact[c - 1][i] + s * s);
            }
        }
    }

    /// Counting each segment as -1 makes the pair minimum report the most segments among the tied
    int step = ties == MOST ? -1 : 1;
    auto solver = [&](long long lambda){
        vector<pair<long long, int>> dp(n + 1, {LLONG_MAX, 0});
        dp[0] = {0, 0};
        for (int j = 1; j <= n; j++){
            for (int i = 0; i < j; i++){
                long long s = pre[j] - pre[i];
                dp[j] = min(dp[j], {dp[i].first + s * s + lambda, dp[i].second + step});
            }
        }
        return pair<long long, int>{dp[n].first, dp[n].second * step};
    };

    for (int k = 1; k <= n; k++) assert(aliens_trick(k, 0LL, pre[n] * pre[n], solver) == exact[k][n]);
}

/// Maximum total of exactly k non-overlapping non-empty subarrays, against the O(n^2 k) DP
void check_max_k_subarrays(int n, int value_max, Ties ties){
    vector<long long> a(n), pre(n + 1, 0);
    for (auto& x : a) x = stress::rand_int(-value_max, value_max);
    for (int i = 0; i < n; i++) pre[i + 1] = pre[i] + a[i];

    vector<vector<long long>> exact(n + 1, vector<long long>(n + 1, -INF));
    for (int j = 0; j <= n; j++) exact[0][j] = 0;
    for (int c = 1; c <= n; c++){
        for (int j = 1; j <= n; j++){
            exact[c][j] = exact[c][j - 1];
            for (int i = 0; i < j; i++){
                if (exact[c - 1][i] != -INF) exact[c][j] = max(exact[c][j], exact[c - 1][i] + pre[j] - pre[i]);
            }
        }
    }

    long long bound = 1;
    for (long long x : a) bound += abs(x);

    int step = ties == MOST ? -1 : 1;
    auto solver = [&](long long lambda){
        pair<long long, int> out = {0, 0}, in = {INF, 0};
        for (long long x : a){
            in = min(in, {out.first + lambda, out.second + step});
            in.first -= x;
            out = min(out, in);
        }
        return pair<long long, int>{out.first, out.second * step};
    };

    for (int k = 1; k <= n; k++) assert(-aliens_trick(k, -bound, bound, solver) == exact[k][n]);
}

double rand_real(double lo, double hi){
    return uniform_real_distribution<double>(lo, hi)(stress::rng());
}

/// Real convex f where each step adds scale items, so k reaches 1e12 and repeated slopes make ties
void check_synthetic_real(Ties ties){
    int first = stress::rand_int(0, 3), n = first + stress::rand_int(0, 10);
    long long scale = stress::rand_int(0, 1) ? 1 : stress::rand_int(1, 100000000000LL);
    double range = stress::rand_int(0, 1) ? 1 : 1e6;

    vector<double> slopes;
    for (int c = first + 1; c <= n; c++){
        slopes.push_back(stress::rand_int(0, 2) ? rand_real(-range, range) : (slopes.empty() ? 0 : slopes.back()));
    }
    sort(slopes.begin(), slopes.end());

    vector<double> f(n + 1, 0);
    f[first] = rand_real(-range, range);
    for (int c = first + 1; c <= n; c++) f[c] = f[c - 1] + slopes[c - first - 1];

    unsigned long long salt = stress::rand_int(0, LLONG_MAX);
    vector<long long> tied;
    auto solver = [&](double lambda){
        double best = 1e300;
        tied.clear();
        for (int c = first; c <= n; c++){
            double value = f[c] + lambda * c * scale;
            if (value < best){
                best = value;
                tied.clear();
            }
            if (value == best) tied.push_back(c * scale);
        }
        return pair<double, long long>{best, pick_tied(tied, ties, lambda, salt)};
    };

    double lo = -range / scale - 1, hi = range / scale + 1;
    for (int c = first; c <= n; c++){
        double err = fabs(aliens_trick(c * scale, lo, hi, solver) - f[c]);
        assert(err <= 1e-9 * (range * (n + 1) + 1) + 1e-14 * max(-lo, hi) * c * scale);
    }
}

/// Split k workers over floors minimizing the sum of a[i] / c[i] with every c[i] >= 1, against an exact DP
void check_split_work(int n, int k, int a_max){
    vector<double> a(n);
    for (auto& x : a) x = stress::rand_int(1, a_max);

    vector<double> exact(k + 1, 1e300);
    exact[0] = 0;
    for (int i = 0; i < n; i++){
        vector<double> next(k + 1, 1e300);
        for (int used = 0; used <= k; used++){
            for (int c = 1; used + c <= k; c++) next[used + c] = min(next[used + c], exact[used] + a[i] / c);
        }
        exact = next;
    }

    auto solver = [&](double lambda){
        pair<double, long long> total = {0, 0};
        for (double x : a){
            pair<double, long long> best = {1e300, 0};
            long long guess = max(1LL, (long long)sqrt(x / max(lambda, 1e-18)));
            for (long long c = max(1LL, guess - 2); c <= min((long long)k, guess + 2); c++) best = min(best, {x / c + lambda * c, c});
            if (best.second == 0) best = {x / k + lambda * k, k};
            total.first += best.first, total.second += best.second;
        }
        return total;
    };

    assert(fabs(aliens_trick(k, 0.0, (double)a_max, solver) - exact[k]) <= 1e-9 * (exact[k] + 1));
}

int main(){
    for (long long it = 0; it < stress::scaled(4000); it++) check_synthetic(Ties(it % 4));

    for (long long it = 0; it < stress::scaled(1500); it++){
        int n = stress::rand_int(1, 12), value_max = it % 3 == 0 ? 2 : (it % 3 == 1 ? 10 : 1000);
        check_squared_split(n, value_max, it % 2 ? MOST : FEWEST);
        check_max_k_subarrays(n, value_max, it % 2 ? MOST : FEWEST);
    }

    for (long long it = 0; it < stress::scaled(4000); it++) check_synthetic_real(Ties(it % 4));

    for (long long it = 0; it < stress::scaled(1500); it++){
        int n = stress::rand_int(1, 6);
        check_split_work(n, n + stress::rand_int(0, 30), it % 2 ? 10 : 1000000);
    }

    for (int it = 0; it < 3; it++){
        check_squared_split(70, 1000000, it == 0 ? FEWEST : MOST);
        check_max_k_subarrays(70, 1000000000, it == 0 ? FEWEST : MOST);
    }

    return 0;
}
