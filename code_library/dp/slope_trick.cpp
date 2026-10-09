/***
 *
 * Slope Trick
 * Maintains a convex piecewise linear function f(x) with integer slopes as two heaps of breakpoints
 *
 * Complexity: O(log n) per add, O(1) per shift, prefix_min, suffix_min and min query, O(n) per eval, n = breakpoints added so far
 *
 * SlopeTrick<T> f: starts as f(x) = 0, T is a signed integer type
 * add_const(c): f(x) += c
 * add_x_minus_a(a): f(x) += max(0, x - a)
 * add_a_minus_x(a): f(x) += max(0, a - x)
 * add_abs(a): f(x) += |x - a|
 * prefix_min(): f(x) = min of f(y) over y <= x, drops every breakpoint right of the minimum
 * suffix_min(): f(x) = min of f(y) over y >= x
 * shift(lo, hi): f(x) = min of f(y) over x - hi <= y <= x - lo, needs lo <= hi
 * shift(d): f(x) = f(x - d), moves the graph right by d
 * get_min(): the minimum of f
 * argmin(): [l, r], the range where f is minimal, -INF or INF when it is unbounded on that side
 * eval(x): f(x), without changing f
 *
 * For long long keep every breakpoint and each side's accumulated shift within [-1e18, 1e18],
 * and the minimum and any evaluated value within the range of T
 * prefix_min() resets the accumulated shift of the right side, suffix_min() the one of the left side
 *
 * Min cost to make a non-decreasing: for each a[i], f.prefix_min() then f.add_abs(a[i]), answer f.get_min()
 * Strictly increasing: run it on a[i] - i
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T>
struct SlopeTrick{
    static constexpr T INF = numeric_limits<T>::max() / 3;

    T min_f = 0, add_left = 0, add_right = 0;
    vector<T> left, right;  /// max-heap and min-heap, every left breakpoint <= every right one, stored without the lazy add

    void add_a_minus_x(T a){
        min_f += max(T(0), a - top_right());
        push_right(a);
        push_left(pop_right());
    }

    void add_abs(T a){
        add_x_minus_a(a);
        add_a_minus_x(a);
    }

    void add_const(T c){
        min_f += c;
    }

    void add_x_minus_a(T a){
        min_f += max(T(0), top_left() - a);
        push_left(a);
        push_right(pop_left());
    }

    pair<T, T> argmin() const{
        return {top_left(), top_right()};
    }

    T eval(T x) const{
        T res = min_f;
        for (T l : left) res += max(T(0), l + add_left - x);
        for (T r : right) res += max(T(0), x - r - add_right);
        return res;
    }

    T get_min() const{
        return min_f;
    }

    void prefix_min(){
        right.clear();
        add_right = 0;
    }

    void shift(T lo, T hi){
        assert(lo <= hi);
        add_left += lo, add_right += hi;
    }

    void shift(T d){
        shift(d, d);
    }

    void suffix_min(){
        left.clear();
        add_left = 0;
    }

private:
    T pop_left(){
        T x = top_left();
        pop_heap(left.begin(), left.end());
        left.pop_back();
        return x;
    }

    T pop_right(){
        T x = top_right();
        pop_heap(right.begin(), right.end(), greater<T>());
        right.pop_back();
        return x;
    }

    void push_left(T x){
        left.push_back(x - add_left);
        push_heap(left.begin(), left.end());
    }

    void push_right(T x){
        right.push_back(x - add_right);
        push_heap(right.begin(), right.end(), greater<T>());
    }

    T top_left() const{
        return left.empty() ? -INF : left[0] + add_left;
    }

    T top_right() const{
        return right.empty() ? INF : right[0] + add_right;
    }
};

/// min of sum |a[i] - b[i]| over non-decreasing b
long long min_cost_non_decreasing(const vector<long long>& a){
    SlopeTrick<long long> f;
    for (long long x : a){
        f.prefix_min();
        f.add_abs(x);
    }
    return f.get_min();
}

int main(){
    const long long INF = SlopeTrick<long long>::INF, BIG = 1e18;

    assert(min_cost_non_decreasing({}) == 0);
    assert(min_cost_non_decreasing({-7}) == 0);
    assert(min_cost_non_decreasing({1, 2, 2, 7}) == 0);
    assert(min_cost_non_decreasing({3, 1, 2}) == 2);
    assert(min_cost_non_decreasing({5, 4, 3, 2, 1}) == 6);
    assert(min_cost_non_decreasing({-BIG, BIG, -BIG, BIG}) == 2 * BIG);

    vector<long long> sonya = {2, 1, 5, 11, 5, 9, 11}, falling = {5, 4, 3, 2, 1};
    for (int i = 0; i < 7; i++) sonya[i] -= i;
    for (int i = 0; i < 5; i++) falling[i] -= i;
    assert(min_cost_non_decreasing(sonya) == 9);
    assert(min_cost_non_decreasing(falling) == 12);

    SlopeTrick<long long> f;
    assert(f.get_min() == 0 && f.argmin() == make_pair(-INF, INF) && f.eval(42) == 0);

    f.add_abs(2);
    f.add_abs(5);
    assert(f.get_min() == 3 && f.argmin() == make_pair(2LL, 5LL));
    assert(f.eval(0) == 7 && f.eval(6) == 5);

    f.add_x_minus_a(4);
    assert(f.get_min() == 3 && f.argmin() == make_pair(2LL, 4LL));
    assert(f.eval(5) == 4 && f.eval(10) == 19);

    f.prefix_min();
    assert(f.get_min() == 3 && f.argmin() == make_pair(2LL, INF));
    assert(f.eval(0) == 7 && f.eval(10) == 3);

    f.shift(1, 3);
    assert(f.get_min() == 3 && f.argmin() == make_pair(3LL, INF));
    assert(f.eval(0) == 9 && f.eval(3) == 3 && f.eval(100) == 3);

    f.add_const(10);
    assert(f.get_min() == 13 && f.eval(0) == 19);

    SlopeTrick<long long> g;
    g.add_a_minus_x(4);
    assert(g.get_min() == 0 && g.argmin() == make_pair(4LL, INF));

    g.add_x_minus_a(6);
    assert(g.argmin() == make_pair(4LL, 6LL) && g.eval(8) == 2 && g.eval(1) == 3);

    g.suffix_min();
    assert(g.argmin() == make_pair(-INF, 6LL) && g.eval(1) == 0 && g.eval(8) == 2);

    g.shift(-2);
    assert(g.argmin() == make_pair(-INF, 4LL) && g.eval(7) == 3);

    SlopeTrick<long long> h;
    h.add_abs(-BIG);
    h.add_abs(BIG);
    assert(h.get_min() == 2 * BIG && h.argmin() == make_pair(-BIG, BIG) && h.eval(0) == 2 * BIG);

    SlopeTrick<long long> s;
    s.shift(-BIG, BIG);
    s.add_abs(BIG);
    s.add_abs(-BIG);
    assert(s.get_min() == 2 * BIG && s.argmin() == make_pair(-BIG, BIG) && s.eval(0) == 2 * BIG);

    SlopeTrick<long long> w;
    w.add_abs(0);
    w.shift(-BIG, BIG);
    assert(w.get_min() == 0 && w.argmin() == make_pair(-BIG, BIG) && w.eval(0) == 0 && w.eval(-BIG) == 0);

    /// a non-increasing f is unchanged by shift(0, hi), so this is min_cost_non_decreasing of 12, 11, ..., 1:
    /// b = 6 everywhere costs 5 + 4 + ... + 0 + 1 + ... + 6 = 36, and 12 shifts of 1e18 would overflow without the reset
    SlopeTrick<long long> r;
    for (long long x = 12; x >= 1; x--){
        r.prefix_min();
        r.shift(0, BIG);
        r.add_abs(x);
    }
    assert(r.get_min() == 36 && r.argmin() == make_pair(6LL, 7LL));

    return 0;
}
