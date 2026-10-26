/***
 *
 * Convex Hull Trick
 * Minimum or maximum of lines y = k * x + b at a point, for lines added in slope order
 *
 * Complexity: O(1) amortized per add_line and query_increasing, O(log n) per query
 *
 * ConvexHullTrick<T> hull: no lines yet, answers the minimum
 * ConvexHullTrick<T, true> hull: the same hull answering the maximum
 * hull.add_line(k, b): add y = k * x + b, k at most every earlier slope for the minimum, at least for the maximum
 * hull.query(x): the minimum over the lines at x, any x in any order, at least one line
 * hull.query_increasing(x): the same, for x never smaller than in the previous query_increasing call
 *
 * Notes:
 *   - Every k * x + b must fit in T, like |k|, |x| <= 1e9 and |b| <= 1e18 for long long, and needs __int128
 *   - T may be double or long double, every answer then within rounding of the exact optimum
 *   - Slopes in any order or lines only on part of the range: li_chao_tree.cpp
 *   - dp[i] = min over j < i of dp[j] + a[j] * x[i], a non-increasing: query(x[i]), then add_line(a[i], dp[i])
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T, bool MAX = false>
struct ConvexHullTrick{
    /// The lower hull with strictly decreasing slopes; the maximum hull is stored negated
    vector<T> k, b;
    int ptr = 0;

    void add_line(T slope, T intercept){
        if (MAX) slope = -slope, intercept = -intercept;
        assert(k.empty() || slope <= k.back());
        if (!k.empty() && slope == k.back()){
            if (intercept >= b.back()) return;
            k.pop_back(), b.pop_back();
        }

        while (k.size() >= 2 && useless(k.size() - 2, k.size() - 1, slope, intercept)) k.pop_back(), b.pop_back();
        k.push_back(slope), b.push_back(intercept);
        ptr = min(ptr, (int)k.size() - 1);
    }

    T query(T x) const{
        assert(!k.empty());
        int lo = 0, hi = k.size() - 1;
        while (lo < hi){
            int mid = (lo + hi) / 2;
            if (value(mid, x) >= value(mid + 1, x)) lo = mid + 1;
            else hi = mid;
        }

        return MAX ? -value(lo, x) : value(lo, x);
    }

    T query_increasing(T x){
        assert(!k.empty());
        while (ptr + 1 < (int)k.size() && value(ptr + 1, x) <= value(ptr, x)) ptr++;
        return MAX ? -value(ptr, x) : value(ptr, x);
    }

    T value(int i, T x) const{
        return k[i] * x + b[i];
    }

    /// Line j between lines i and new is never the minimum once new crosses i no later than j does
    bool useless(int i, int j, T slope, T intercept) const{
        /// A cast to __int128 would drop the fractions of floating-point lines
        if constexpr (is_floating_point_v<T>){
            return (long double)(intercept - b[i]) * (k[i] - k[j]) <= (long double)(b[j] - b[i]) * (k[i] - slope);
        }
        return (__int128)(intercept - b[i]) * (k[i] - k[j]) <= (__int128)(b[j] - b[i]) * (k[i] - slope);
    }
};

int main(){
    /// Lines y = 3x, y = x + 4, y = -x + 9 and y = -2x + 20, the last never the minimum below x = 11
    ConvexHullTrick<long long> lower;
    lower.add_line(3, 0);
    lower.add_line(1, 4);
    lower.add_line(-1, 9);
    lower.add_line(-2, 20);
    assert(lower.query(0) == 0);       /// 3 * 0
    assert(lower.query(3) == 6);       /// min(9, 7, 6, 14)
    assert(lower.query(5) == 4);       /// -5 + 9
    assert(lower.query(12) == -4);     /// -24 + 20 beats -12 + 9
    assert(lower.query_increasing(0) == 0);
    assert(lower.query_increasing(3) == 6);
    assert(lower.query_increasing(12) == -4);

    /// The maximum takes slopes in non-decreasing order
    ConvexHullTrick<long long, true> upper;
    upper.add_line(-1, 5);
    upper.add_line(0, 2);
    upper.add_line(2, -3);
    assert(upper.query(-4) == 9);      /// 4 + 5
    assert(upper.query(2) == 3);       /// max(3, 2, 1)
    assert(upper.query(10) == 17);     /// 20 - 3

    return 0;
}
