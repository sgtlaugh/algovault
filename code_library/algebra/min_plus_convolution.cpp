/***
 *
 * Min-Plus Convolution
 * c[k] = min over i + j = k of a[i] + b[j], fast when one or both of a and b are convex
 *
 * Complexity:
 *   - min_plus_convex_convex, max_plus_concave_concave: O(n + m)
 *   - min_plus_convex_arbitrary: O((n + m) log(n + m))
 *
 * max_plus_concave_concave(a, b): c of size n + m - 1 with c[k] = max over i + j = k of a[i] + b[j], both concave
 * min_plus_convex_arbitrary(a, b): c of size n + m - 1, a convex and b arbitrary
 * min_plus_convex_convex(a, b): c of size n + m - 1 for a of size n and b of size m, both convex
 *
 * Convex means a[i + 1] - a[i] >= a[i] - a[i - 1] for every inner i, so sizes 1 and 2 are always convex
 * Convexity is not checked, and a non-convex input gives wrong results
 * An empty input gives an empty result
 * Max-plus with a concave and b arbitrary: negate both, call min_plus_convex_arbitrary and negate the result
 * With T = long long every |a[i]| and |b[j]| must be below 2^62, so sums and adjacent differences fit
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T>
vector<T> min_plus_convex_convex(const vector<T>& a, const vector<T>& b){
    int n = a.size(), m = b.size();
    if (!n || !m) return {};

    vector<T> c(n + m - 1);
    c[0] = a[0] + b[0];
    /// Merges the difference sequences of a and b in increasing order, as in merge sort
    for (int i = 0, j = 0; i + j < n + m - 2; ){
        if (j == m - 1 || (i < n - 1 && a[i + 1] - a[i] < b[j + 1] - b[j])) i++;
        else j++;
        c[i + j] = a[i] + b[j];
    }

    return c;
}

template <typename T>
vector<T> max_plus_concave_concave(vector<T> a, vector<T> b){
    for (auto& x : a) x = -x;
    for (auto& x : b) x = -x;

    auto c = min_plus_convex_convex(a, b);
    for (auto& x : c) x = -x;
    return c;
}

template <typename T>
vector<T> min_plus_convex_arbitrary(const vector<T>& a, const vector<T>& b){
    int n = a.size(), m = b.size();
    if (!n || !m) return {};

    vector<T> c(n + m - 1);
    /// a convex keeps the leftmost optimal j for c[k] from decreasing as k grows, so each half searches one side of best
    auto solve = [&](auto&& self, int lo, int hi, int jlo, int jhi) -> void {
        if (lo > hi) return;

        int k = (lo + hi) / 2, best = max(jlo, k - n + 1);
        for (int j = best + 1; j <= min(jhi, k); j++){
            if (b[j] + a[k - j] < b[best] + a[k - best]) best = j;
        }
        c[k] = b[best] + a[k - best];

        self(self, lo, k - 1, jlo, best);
        self(self, k + 1, hi, best, jhi);
    };

    solve(solve, 0, n + m - 2, 0, m - 1);
    return c;
}

int main(){
    /// a = {0, 1, 4} is convex, c[2] = min(0 + 7, 1 + 0, 4 + 5) = 1 against the arbitrary b
    vector<long long> a = {0, 1, 4};
    assert((min_plus_convex_convex(a, vector<long long>{2, 2}) == vector<long long>{2, 2, 3, 6}));
    assert((min_plus_convex_arbitrary(a, vector<long long>{5, 0, 7, 2}) == vector<long long>{5, 0, 1, 2, 3, 6}));

    /// Both concave, c[2] = max(3 + 2, 4 + 0) = 5
    assert((max_plus_concave_concave<long long>({0, 3, 4}, {0, 2}) == vector<long long>{0, 3, 5, 6}));

    return 0;
}
