/***
 *
 * Merge Sort Tree and Wavelet Matrix
 * Counts values below a threshold or equal to a value and finds the k-th smallest in a range of a static array
 *
 * Complexity:
 *   - MergeSortTree: O(n log n) to build and memory, O(log^2 n) per query
 *   - WaveletMatrix: O(n log n) to build, O(n log sigma) memory, O(log sigma) per query, sigma = number of distinct values
 *
 * MergeSortTree<T> tree(values): positions 0..n - 1
 * WaveletMatrix<T> tree(values): positions 0..n - 1, the same queries plus count_equal
 * tree.count_less(l, r, x): how many of values[l..r] are < x, 0 <= l <= r < n
 * tree.count_less_equal(l, r, x): how many of values[l..r] are <= x, 0 <= l <= r < n
 * tree.kth(l, r, k): the k-th smallest of values[l..r], 0-based, 0 <= k <= r - l
 * tree.count_equal(l, r, x): WaveletMatrix only, how many of values[l..r] equal x
 *
 * WaveletMatrix takes 4 bytes per element per level, about 80 MB for n = 1e6 distinct values
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

    /// The tree is over ranks, not positions: a count binary searches x once in sorted_values, then descends over
    /// that rank prefix, and kth descends the same tree instead of binary searching on the answer
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
};

template <typename T>
struct WaveletMatrix{
    int n, levels = 0;
    vector<T> distinct;
    vector<vector<int>> zeros_before;  /// zeros_before[d][i]: elements among the first i of level d whose bit is 0
    vector<int> zeros;

    /// Values become codes [0, sigma) split on one bit per level, high bit first; each level stably moves the 0-bit
    /// elements before the 1-bit ones, so a range maps to one range per level
    WaveletMatrix(const vector<T>& values) : n(values.size()), distinct(values){
        sort(distinct.begin(), distinct.end());
        distinct.erase(unique(distinct.begin(), distinct.end()), distinct.end());
        while ((1 << levels) < (int)distinct.size()) levels++;

        vector<int> codes(n);
        for (int i = 0; i < n; i++) codes[i] = code_of(values[i]);

        zeros_before.assign(levels, vector<int>(n + 1, 0));
        zeros.assign(levels, 0);
        for (int d = 0; d < levels; d++){
            int bit = levels - 1 - d;
            for (int i = 0; i < n; i++) zeros_before[d][i + 1] = zeros_before[d][i] + !(codes[i] >> bit & 1);
            zeros[d] = zeros_before[d][n];
            stable_partition(codes.begin(), codes.end(), [&](int c){ return !(c >> bit & 1); });
        }
    }

    int count_less(int l, int r, const T& x) const{
        assert(0 <= l && l <= r && r < n);
        return count_code_less(l, r + 1, code_of(x));
    }

    int count_less_equal(int l, int r, const T& x) const{
        assert(0 <= l && l <= r && r < n);
        return count_code_less(l, r + 1, upper_bound(distinct.begin(), distinct.end(), x) - distinct.begin());
    }

    int count_equal(int l, int r, const T& x) const{
        assert(0 <= l && l <= r && r < n);
        int c = code_of(x);
        if (c == (int)distinct.size() || distinct[c] != x) return 0;
        return count_code_less(l, r + 1, c + 1) - count_code_less(l, r + 1, c);
    }

    T kth(int l, int r, int k) const{
        assert(0 <= l && l <= r && r < n && 0 <= k && k <= r - l);
        int code = 0, lo = l, hi = r + 1;
        for (int d = 0; d < levels; d++){
            int zl = zeros_before[d][lo], zr = zeros_before[d][hi];
            if (k < zr - zl){
                lo = zl, hi = zr;
                continue;
            }
            k -= zr - zl, code |= 1 << (levels - 1 - d);
            lo = zeros[d] + lo - zl, hi = zeros[d] + hi - zr;
        }
        return distinct[code];
    }

    int code_of(const T& x) const{
        return lower_bound(distinct.begin(), distinct.end(), x) - distinct.begin();
    }

    /// how many codes below p sit at positions [l, r)
    int count_code_less(int l, int r, int p) const{
        if (p >= (int)distinct.size()) return r - l;

        int res = 0;
        for (int d = 0; d < levels; d++){
            int zl = zeros_before[d][l], zr = zeros_before[d][r];
            if (p >> (levels - 1 - d) & 1) res += zr - zl, l = zeros[d] + l - zl, r = zeros[d] + r - zr;
            else l = zl, r = zr;
        }
        return res;
    }
};

int main(){
    int n, m;
    scanf("%d %d", &n, &m);
    vector<int> a(n);
    for (int i = 0; i < n; i++) scanf("%d", &a[i]);

    MergeSortTree<int> tree(a);
    while (m--){
        int i, j, k;
        scanf("%d %d %d", &i, &j, &k);
        printf("%d\n", tree.kth(i - 1, j - 1, k - 1));
    }
    return 0;
}
