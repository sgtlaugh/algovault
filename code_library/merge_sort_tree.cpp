/***
 *
 * Merge Sort Tree
 * Static array, counts values below a threshold in a range and finds the k-th smallest in a range
 *
 * Complexity: O(n log n) to build and memory, O(log^2 n) per query
 *
 * MergeSortTree<T> tree(values), 0-based positions, ranges [l, r] inclusive
 * tree.count_less(l, r, x), tree.count_less_equal(l, r, x): how many values in [l, r] are < x or <= x
 * tree.kth(l, r, k): the k-th smallest value in [l, r], 0-based, needs 0 <= k <= r - l
 *
 * kth descends a second tree built over positions sorted by value, avoiding a binary search on the answer
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T>
struct MergeSortTree{
    int n;
    vector<T> sorted_values;
    vector<vector<T>> by_position;  /// node covering positions [a, b] keeps their values sorted
    vector<vector<int>> by_rank;    /// node covering ranks [a, b] keeps the positions of those values sorted

    MergeSortTree(const vector<T>& values) : n(values.size()), by_position(4 * max(1, (int)values.size())), by_rank(4 * max(1, (int)values.size())){
        vector<int> order(n);
        iota(order.begin(), order.end(), 0);
        stable_sort(order.begin(), order.end(), [&](int i, int j){ return values[i] < values[j]; });
        for (int i : order) sorted_values.push_back(values[i]);
        if (n) build(1, 0, n - 1, values, order);
    }

    void build(int node, int a, int b, const vector<T>& values, const vector<int>& order){
        if (a == b){
            by_position[node] = {values[a]};
            by_rank[node] = {order[a]};
            return;
        }
        int m = (a + b) / 2;
        build(2 * node, a, m, values, order);
        build(2 * node + 1, m + 1, b, values, order);
        merge(by_position[2 * node].begin(), by_position[2 * node].end(), by_position[2 * node + 1].begin(), by_position[2 * node + 1].end(), back_inserter(by_position[node]));
        merge(by_rank[2 * node].begin(), by_rank[2 * node].end(), by_rank[2 * node + 1].begin(), by_rank[2 * node + 1].end(), back_inserter(by_rank[node]));
    }

    template <typename Bound>
    int count(int node, int a, int b, int l, int r, const T& x, Bound bound) const{
        if (r < a || b < l) return 0;
        if (l <= a && b <= r) return bound(by_position[node].begin(), by_position[node].end(), x) - by_position[node].begin();
        int m = (a + b) / 2;
        return count(2 * node, a, m, l, r, x, bound) + count(2 * node + 1, m + 1, b, l, r, x, bound);
    }

    int count_less(int l, int r, const T& x) const{
        assert(0 <= l && l <= r && r < n);
        return count(1, 0, n - 1, l, r, x, [](auto first, auto last, const T& v){ return lower_bound(first, last, v); });
    }

    int count_less_equal(int l, int r, const T& x) const{
        assert(0 <= l && l <= r && r < n);
        return count(1, 0, n - 1, l, r, x, [](auto first, auto last, const T& v){ return upper_bound(first, last, v); });
    }

    T kth(int l, int r, int k) const{
        assert(0 <= l && l <= r && r < n && 0 <= k && k <= r - l);
        int node = 1, a = 0, b = n - 1;
        while (a != b){
            const vector<int>& left = by_rank[2 * node];
            int inside = upper_bound(left.begin(), left.end(), r) - lower_bound(left.begin(), left.end(), l);
            int m = (a + b) / 2;
            if (k < inside) node = 2 * node, b = m;
            else k -= inside, node = 2 * node + 1, a = m + 1;
        }
        return sorted_values[a];
    }
};

int main(){
    MergeSortTree<int> tree({5, 1, 4, 1, 3, 9, 2});
    assert(tree.count_less(0, 6, 4) == 4);
    assert(tree.count_less_equal(0, 6, 4) == 5);
    assert(tree.count_less(1, 3, 1) == 0 && tree.count_less_equal(1, 3, 1) == 2);
    assert(tree.count_less_equal(5, 5, 8) == 0 && tree.count_less_equal(5, 5, 9) == 1);

    assert(tree.kth(0, 6, 0) == 1 && tree.kth(0, 6, 1) == 1 && tree.kth(0, 6, 2) == 2 && tree.kth(0, 6, 6) == 9);
    assert(tree.kth(2, 4, 0) == 1 && tree.kth(2, 4, 1) == 3 && tree.kth(2, 4, 2) == 4);
    assert(tree.kth(5, 5, 0) == 9);

    MergeSortTree<long long> negative({-5, LLONG_MIN, LLONG_MAX, 0});
    assert(negative.kth(0, 3, 0) == LLONG_MIN && negative.kth(0, 3, 3) == LLONG_MAX);
    assert(negative.count_less(0, 3, 0) == 2);

    MergeSortTree<int> empty({});
    assert(empty.n == 0);

    return 0;
}
