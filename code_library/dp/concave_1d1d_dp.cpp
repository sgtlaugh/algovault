/***
 *
 * Concave 1D1D DP (quadrangle inequality / Monge cost)
 * dp[x] = min over 0 <= i < x of dp[i] + w(i, x), for x = 1..n, with dp[0] given
 *
 * Complexity: O(n log n) evaluations of w, O(n) memory
 *
 * concave_1d1d_dp(n, dp0, w, &opt): returns dp[0..n], opt (optional) gets opt[x] = an argmin i for x >= 1
 * w(i, x) is called only with 0 <= i < x <= n, T is any type with +, < and copy (int64, double, ...)
 *
 * Requirement, the quadrangle inequality: w(a, c) + w(b, d) <= w(a, d) + w(b, c) for all a <= b <= c <= d
 * It holds for g(x - i) with g convex, for -A[x] * B[i] with A and B nondecreasing, for any term that
 * depends on only i or only x, and for sums of these. Checking a few values by brute force is the practical test
 *
 * Under it the optimal i is nondecreasing in x, and a newer candidate that beats an older one at some x
 * keeps beating it at every larger x. The deque holds (start, candidate) segments covering x+1..n:
 * a new candidate pops the segments it wins entirely, then binary searches where it overtakes the last one
 *
 * Example: split n items into segments, a segment of length L costing L^2 + C
 *   auto dp = concave_1d1d_dp(n, 0LL, [&](int i, int x){ return 1LL * (x - i) * (x - i) + C; });
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T, typename Cost>
vector<T> concave_1d1d_dp(int n, T dp0, Cost w, vector<int>* opt = nullptr){
    vector<T> dp(n + 1, dp0);
    if (opt) opt->assign(n + 1, -1);
    auto value = [&](int i, int x){ return dp[i] + w(i, x); };

    deque<pair<int, int>> segments = {{1, 0}};
    for (int x = 1; x <= n; x++){
        while (segments.size() > 1 && segments[1].first <= x) segments.pop_front();
        int k = segments.front().second;
        dp[x] = value(k, x);
        if (opt) (*opt)[x] = k;

        while (segments.back().first > x){
            auto [start, old] = segments.back();
            if (value(old, start) < value(x, start)) break;
            segments.pop_back();
        }

        int old = segments.back().second, lo = max(segments.back().first, x + 1), hi = n + 1;
        while (lo < hi){
            int mid = lo + (hi - lo) / 2;
            if (!(value(old, mid) < value(x, mid))) hi = mid;
            else lo = mid + 1;
        }
        if (lo <= n) segments.push_back({lo, x});
    }

    return dp;
}

int main(){
    auto squares = [](long long C){
        return [C](int i, int x){ return 1LL * (x - i) * (x - i) + C; };
    };

    vector<int> opt;
    assert((concave_1d1d_dp(6, 0LL, squares(5), &opt) == vector<long long>{0, 6, 9, 14, 18, 23, 27}));
    assert(opt[6] == 4 && opt[4] == 2 && opt[2] == 0);

    assert((concave_1d1d_dp(4, 0LL, squares(100)) == vector<long long>{0, 101, 104, 109, 116}));
    assert((concave_1d1d_dp(0, 7LL, squares(1)) == vector<long long>{7}));
    assert((concave_1d1d_dp(1, -3LL, squares(1)) == vector<long long>{-3, -1}));

    int n = 5, s = 1;
    vector<long long> T = {0, 1, 3, 4, 2, 1}, F = {0, 3, 2, 3, 3, 4};
    for (int i = 1; i <= n; i++) T[i] += T[i - 1], F[i] += F[i - 1];
    auto batch = [&](int i, int x){ return s * (F[n] - F[i]) + T[x] * (F[x] - F[i]); };
    assert(concave_1d1d_dp(n, 0LL, batch)[n] == 153);

    auto halves = concave_1d1d_dp(4, 0.0, [](int i, int x){ return (x - i) * (x - i) * 0.5 + 1.25; });
    assert(abs(halves[4] - 6.5) < 1e-9);

    return 0;
}
