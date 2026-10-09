/***
 *
 * Min-Plus Convolution
 * c[k] = min over i + j = k of a[i] + b[j], for a of size n and b of size m, c has size n + m - 1
 *
 * Complexity:
 *   min_plus_convex_convex(a, b)      O(n + m)                 a and b convex
 *   max_plus_concave_concave(a, b)    O(n + m)                 a and b concave, max instead of min
 *   min_plus_convex_arbitrary(a, b)   O((n + m) log(n + m))    a convex, b arbitrary
 *
 * Convex means a[i + 1] - a[i] >= a[i] - a[i - 1] for every inner i, sizes 1 and 2 are always convex
 * Convexity is not checked, a non-convex input gives wrong results
 * An empty input gives an empty result
 * Max-plus with a concave and b arbitrary: negate both, call min_plus_convex_arbitrary, negate the result
 *
 * With T = long long every |a[i]|, |b[j]| must be below 2^62, so sums and adjacent differences fit
 *
 * convex_convex merges the two difference sequences in sorted order, as in merge sort
 * convex_arbitrary uses monotone minima: the leftmost optimal j for c[k] never decreases as k grows
 *
 * Example:
 *   min_plus_convex_convex<long long>({0, 1, 4}, {2, 2})               // {2, 2, 3, 6}
 *   min_plus_convex_arbitrary<long long>({0, 1, 4}, {5, 0, 7, 2})      // {5, 0, 1, 2, 3, 6}
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
    const long long L = (1LL << 62) - 1;

    assert((min_plus_convex_convex<long long>({0, 1, 4}, {2, 2}) == vector<long long>{2, 2, 3, 6}));
    assert((min_plus_convex_convex<long long>({5}, {3, 1, 0}) == vector<long long>{8, 6, 5}));
    assert((min_plus_convex_convex<long long>({4, 1, 0, 1, 4}, {0, -1, 0}) == vector<long long>{4, 1, 0, -1, 0, 1, 4}));
    assert((min_plus_convex_convex<long long>({L, -L, L}, {-L, L}) == vector<long long>{0, -2 * L, 0, 2 * L}));
    assert((min_plus_convex_convex<long long>({}, {1}).empty()));
    assert((min_plus_convex_convex<long long>({1}, {}).empty()));

    assert((max_plus_concave_concave<long long>({0, 3, 4}, {0, 2}) == vector<long long>{0, 3, 5, 6}));
    assert((max_plus_concave_concave<long long>({-L, L, -L}, {L, -L}) == vector<long long>{0, 2 * L, 0, -2 * L}));

    assert((min_plus_convex_arbitrary<long long>({0, 1, 4}, {5, 0, 7, 2}) == vector<long long>{5, 0, 1, 2, 3, 6}));
    assert((min_plus_convex_arbitrary<long long>({4, 1, 0, 1, 4}, {0, -1, 0}) == vector<long long>{4, 1, 0, -1, 0, 1, 4}));
    assert((min_plus_convex_arbitrary<long long>({7}, {3, -2, 9}) == vector<long long>{10, 5, 16}));
    assert((min_plus_convex_arbitrary<long long>({L, -L, L}, {L, -L}) == vector<long long>{2 * L, 0, -2 * L, 0}));
    assert((min_plus_convex_arbitrary<long long>({}, {1}).empty()));

    return 0;
}
