/***
 *
 * Aliens Trick (WQS Binary Search, Lagrangian Relaxation)
 * Minimum cost of a solution with exactly k items, when that cost f(k) is convex in k
 *
 * Complexity: O(log(hi - lo)) calls to the solver, plus one
 *
 * aliens_trick(k, lo, hi, solver) returns f(k)
 *   solver(lambda) returns {value, count}: value = min over all solutions of cost + lambda * items,
 *   count = the MINIMUM number of items among the solutions reaching that value
 *   (a DP over pairs {cost, items} compared lexicographically gives exactly this)
 *
 * Requirements:
 *   - integer costs, f convex on its domain and k inside the domain
 *   - [lo, hi] contains every marginal -(f(c) - f(c - 1)), so M = max |f(c) - f(c - 1)| gives [-M, M]
 *   - the bounds' type is wide enough for hi - lo, the solver's value type for cost + max(|lo|, |hi|) * (max items)
 *   - f(k) comes back in the solver's value type, so int bounds with an int64 solver still return int64
 *
 * Ties: several counts can be optimal for one lambda, so a lambda with count == k may not exist
 * The driver takes the smallest lambda whose minimum optimal count is <= k; k is then optimal there too,
 * and f(k) = value - lambda * k (never value - lambda * count)
 *
 * Maximizing a concave f: minimize the negated gains and negate the result
 *
 * Example: split a non-negative array into k segments minimizing the sum of squared segment sums
 *   aliens_trick(k, 0LL, S * S, solver), solver(lambda) runs the O(n^2) split DP with + lambda per segment
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T, typename Solver>
auto aliens_trick(int k, T lo, T hi, Solver solver){
    while (lo < hi){
        T mid = lo + (hi - lo) / 2;
        if (solver(mid).second <= k) hi = mid;
        else lo = mid + 1;
    }

    auto value = solver(lo).first;
    return value - decltype(value)(lo) * k;
}

struct Line{
    long long slope, intercept;
    int count;

    pair<long long, int> at(long long x) const{
        return {slope * x + intercept, count};
    }
};

bool useless(const Line& a, const Line& b, const Line& c){
    return (__int128)(c.intercept - a.intercept) * (a.slope - b.slope) < (__int128)(b.intercept - a.intercept) * (a.slope - c.slope);
}

long long take_photos(int n, int m, int k, vector<int> r, vector<int> c){
    vector<pair<int, int>> segs(n);
    for (int i = 0; i < n; i++) segs[i] = {min(r[i], c[i]), -max(r[i], c[i])};
    sort(segs.begin(), segs.end());

    vector<long long> lo = {0}, hi = {-1};
    for (auto& s : segs){
        if (-s.second > hi.back()) lo.push_back(s.first), hi.push_back(-s.second);
    }
    int cnt = lo.size() - 1;

    auto sq = [](long long v){ return v * v; };
    auto solver = [&](long long lambda){
        vector<Line> hull;
        size_t head = 0;
        pair<long long, int> dp = {0, 0};
        for (int i = 1; i <= cnt; i++){
            long long a = lo[i] - 1, overlap = i > 1 ? sq(max(0LL, hi[i - 1] - a)) : 0;
            Line line = {-2 * a, dp.first - overlap + a * a, dp.second};
            while (hull.size() - head >= 2 && useless(hull[hull.size() - 2], hull.back(), line)) hull.pop_back();
            hull.push_back(line);

            while (hull.size() - head >= 2 && hull[head + 1].at(hi[i]) <= hull[head].at(hi[i])) head++;
            dp = hull[head].at(hi[i]);
            dp.first += sq(hi[i]) + lambda, dp.second++;
        }

        return dp;
    };

    return aliens_trick(min(k, cnt), 0LL, (long long)m * m, solver);
}

#ifndef EVAL
int main(){
    int n, m, k;

    scanf("%d %d %d", &n, &m, &k);
    vector<int> r(n), c(n);
    for (int i = 0; i < n; i++) scanf("%d %d", &r[i], &c[i]);
    printf("%lld\n", take_photos(n, m, k, r, c));
    return 0;
}
#endif
