/***
 *
 * Dinkelbach Fractional Programming
 * Largest or smallest ratio a(x) / b(x) over choices x with b(x) > 0, given an oracle for a(x) - lambda b(x)
 *
 * Complexity: O(log(n M)) calls of best when a(x) and b(x) are sums of n integers of absolute value at most M
 *
 * max_ratio(best, a, b): the (a(x), b(x)) of an x maximizing a(x) / b(x), starting from any a / b <= the maximum
 * min_ratio(best, a, b): the (a(x), b(x)) of an x minimizing a(x) / b(x), starting from any a / b >= the minimum
 *
 * Notes:
 *   - best(p, q) returns the pair (a(x), b(x)) of an x maximizing q a(x) - p b(x), for min_ratio minimizing it
 *   - T is long long for exact answers or double, and best returns exactly pair<T, T>
 *   - long long: the template compares in __int128, best must keep q a(x) - p b(x) from overflowing itself
 *   - best may return any x with q a(x) - p b(x) > 0 (< 0 for min_ratio) when one exists: still exact, but slower
 *   - An oracle that gives only the value max(q a(x) - p b(x)): binary search p / q on its sign instead
 *   - With a feasible start the answer is always some x's (a(x), b(x)), otherwise it can be the start itself
 *   - Dinkelbach, On Nonlinear Fractional Programming, Management Science 13(7), 1967
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template<typename T, typename F>
pair<T, T> max_ratio(F best, T a, T b){
    /// A literal start like (0, 1) would deduce T = int and silently truncate the oracle's sums
    static_assert(is_same<decltype(best(a, b)), pair<T, T>>::value, "best must return pair<T, T>");
    using W = conditional_t<is_integral<T>::value, __int128, T>;
    assert(b > 0);

    while (true){
        auto [na, nb] = best(a, b);
        assert(nb > 0);
        /// No eps: stopping once the ratio fails to strictly increase also ends the loop when doubles round
        if ((W)na * b <= (W)a * nb) return {a, b};
        a = na, b = nb;
    }
}

template<typename T, typename F>
pair<T, T> min_ratio(F best, T a, T b){
    auto negated = [&](T p, T q){
        auto res = best(-p, q);
        res.first = -res.first;
        return res;
    };
    auto [na, nb] = max_ratio(negated, -a, b);

    return {-na, nb};
}

int main(){
    /// Keep 2 of 3 tests maximizing total score / total weight (POJ 2976): scores {5, 0, 2}, weights {5, 1, 6}
    vector<long long> score = {5, 0, 2}, weight = {5, 1, 6};
    auto best_two = [&](long long p, long long q){
        vector<int> order = {0, 1, 2};
        sort(order.begin(), order.end(), [&](int i, int j){
            return q * score[i] - p * weight[i] > q * score[j] - p * weight[j];
        });
        return make_pair(score[order[0]] + score[order[1]], weight[order[0]] + weight[order[1]]);
    };
    auto [a, b] = max_ratio(best_two, 0LL, 1LL);
    assert(a == 5);
    assert(b == 6);

    /// Spanning tree of a triangle minimizing cost / length: drop the edge where q cost - p length is largest
    vector<double> cost = {3, 1, 2}, length = {1, 4, 2};
    auto best_tree = [&](double p, double q){
        int drop = 0;
        for (int e = 1; e < 3; e++){
            if (q * cost[e] - p * length[e] > q * cost[drop] - p * length[drop]) drop = e;
        }
        pair<double, double> res(0, 0);
        for (int e = 0; e < 3; e++){
            if (e != drop) res.first += cost[e], res.second += length[e];
        }

        return res;
    };
    auto [c, l] = min_ratio(best_tree, 6.0, 1.0);
    assert(abs(c / l - 0.5) < 1e-12);  /// edges 1 and 2: 3 / 6

    return 0;
}
