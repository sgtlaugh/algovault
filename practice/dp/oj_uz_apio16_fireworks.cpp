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
 * argmin(): [l, r], the range where f is minimal, -INF or INF when it is unbounded on that side, INF = max of T / 2
 * eval(x): f(x), without changing f
 *
 * Keep every breakpoint, each side's accumulated shift and every eval argument within [-1e18, 1e18] for long long
 * or [-1e9, 1e9] for int, and the minimum and any evaluated value within the range of T
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
    /// above every allowed breakpoint, so an empty heap's sentinel adds 0 to min_f in add_*, and INF - (-INF) still fits in T
    static constexpr T INF = numeric_limits<T>::max() / 2;

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

void merge_into(SlopeTrick<long long>& to, SlopeTrick<long long>& from){
    if (to.left.size() + to.right.size() < from.left.size() + from.right.size()) swap(to, from);
    to.add_const(from.get_min());
    for (long long l : from.left) to.add_a_minus_x(l + from.add_left);
    for (long long r : from.right) to.add_x_minus_a(r + from.add_right);
    from = SlopeTrick<long long>();
}

int main(){
    int n, m;

    scanf("%d %d", &n, &m);
    vector<int> parent(n + m + 1);
    vector<long long> len(n + m + 1);
    for (int i = 2; i <= n + m; i++) scanf("%d %lld", &parent[i], &len[i]);

    vector<SlopeTrick<long long>> f(n + m + 1);
    for (int i = n + m; i >= 2; i--){
        if (i > n) f[i].add_abs(0);

        long long c = len[i], l = f[i].pop_left(), r = f[i].top_right();
        f[i].push_left(l + c);
        f[i].prefix_min();
        f[i].push_right(r + c);
        merge_into(f[parent[i]], f[i]);
    }

    printf("%lld\n", f[1].get_min());
    return 0;
}
