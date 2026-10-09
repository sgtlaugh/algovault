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
 * The tree is built over ranks (positions sorted by value), each node keeping its positions sorted
 * Counts binary search x in the sorted values, then descend over that rank prefix
 * kth descends the same tree, avoiding a binary search on the answer
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T>
struct MergeSortTree{
    int n;
    vector<T> sorted_values;
    vector<vector<int>> by_rank;  /// node covering ranks [a, b] keeps the positions of those values sorted

    MergeSortTree(const vector<T>& values) : n(values.size()), by_rank(4 * max(1, (int)values.size())){
        vector<int> order(n);
        iota(order.begin(), order.end(), 0);
        sort(order.begin(), order.end(), [&](int i, int j){ return values[i] < values[j]; });
        for (int i : order) sorted_values.push_back(values[i]);
        if (n) build(1, 0, n - 1, order);
    }

    void build(int node, int a, int b, const vector<int>& order){
        if (a == b){
            by_rank[node] = {order[a]};
            return;
        }

        int m = (a + b) / 2;
        build(2 * node, a, m, order);
        build(2 * node + 1, m + 1, b, order);
        merge(by_rank[2 * node].begin(), by_rank[2 * node].end(), by_rank[2 * node + 1].begin(), by_rank[2 * node + 1].end(), back_inserter(by_rank[node]));
    }

    int count_in(int node, int l, int r) const{
        const vector<int>& v = by_rank[node];
        return upper_bound(v.begin(), v.end(), r) - lower_bound(v.begin(), v.end(), l);
    }

    /// how many of the ranks [0, p) sit at positions in [l, r]
    int count_prefix(int p, int l, int r) const{
        int res = 0, node = 1, a = 0, b = n - 1;
        while (a < p){
            if (b < p) return res + count_in(node, l, r);
            int m = (a + b) / 2;
            if (m < p) res += count_in(2 * node, l, r), node = 2 * node + 1, a = m + 1;
            else node = 2 * node, b = m;
        }
        return res;
    }

    int count_less(int l, int r, const T& x) const{
        assert(0 <= l && l <= r && r < n);
        return count_prefix(lower_bound(sorted_values.begin(), sorted_values.end(), x) - sorted_values.begin(), l, r);
    }

    int count_less_equal(int l, int r, const T& x) const{
        assert(0 <= l && l <= r && r < n);
        return count_prefix(upper_bound(sorted_values.begin(), sorted_values.end(), x) - sorted_values.begin(), l, r);
    }

    T kth(int l, int r, int k) const{
        assert(0 <= l && l <= r && r < n && 0 <= k && k <= r - l);
        int node = 1, a = 0, b = n - 1;
        while (a != b){
            int inside = count_in(2 * node, l, r), m = (a + b) / 2;
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
