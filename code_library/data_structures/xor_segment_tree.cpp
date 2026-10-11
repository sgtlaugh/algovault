/***
 *
 * XOR Segment Tree
 * Fold of a[p ^ x] over l <= p <= r for any x under any associative merge, with point add for sums
 *
 * Complexity: O(n log n) time and memory to build, O(log n) per query and per add
 *
 * The size n must be a power of two (pad by hand, what a padded a[p ^ x] means depends on the problem)
 * Queries need 0 <= l <= r < n and 0 <= x < n
 * Every node over [b, b + len) stores t[node][y] = merge of a[p ^ y] over its range in p order, for all 0 <= y < len
 * For sums every t[node][y] equals the block sum, since xor by y < len only permutes the block
 * The layout only pays off for order-dependent merges such as max subarray sum (CF 1716E)
 * The query maps p -> p ^ x level by level: a set bit of x swaps the halves the query range lands in
 *
 * XorSegmentTree<T, Merge = plus<T>> st(a, merge), Merge needs a const operator() and need not be commutative
 * add only exists for Merge = plus<T>, any other merge builds a static tree
 *
 * Usage:
 *   XorSegmentTree<long long> st(a);       // a.size() is a power of two
 *   st.add(i, v);                          // a[i] += v
 *   st.query(l, r, x);                     // sum of a[p ^ x] for l <= p <= r
 *
 *   struct Node{ long long sum, pre, suf, best; };
 *   struct Join{ Node operator()(const Node& a, const Node& b) const{ ... } };
 *   XorSegmentTree<Node, Join> best(nodes); // CF 1716E: best.query(0, n - 1, x).best
 *
 * Related problem: https://codeforces.com/contest/1654/problem/F
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T, typename Merge = plus<T>>
struct XorSegmentTree{
    int n;
    Merge merge;
    vector<vector<T>> t;
    vector<T> added;

    XorSegmentTree(const vector<T>& a, Merge merge = Merge()) : n((int)a.size()), merge(merge), t(2 * a.size()), added(2 * a.size()){
        assert(n > 0 && (n & (n - 1)) == 0);
        build(1, 0, n, a);
    }

    void add(int i, T v){
        static_assert(is_same_v<Merge, plus<T>>, "add relies on every node value being the plain block sum");
        assert(0 <= i && i < n);
        add(1, 0, n, i, v);
    }

    T query(int l, int r, int x) const{
        assert(0 <= l && l <= r && r < n && 0 <= x && x < n);
        return query(1, 0, n, l, r, x);
    }

private:
    void add(int node, int b, int len, int i, T v){
        added[node] += v;
        if (len == 1) return;

        int half = len >> 1;
        if (i < b + half) add(node << 1, b, half, i, v);
        else add(node << 1 | 1, b + half, half, i, v);
    }

    void build(int node, int b, int len, const vector<T>& a){
        if (len == 1){
            t[node] = {a[b]};
            return;
        }

        int half = len >> 1, lc = node << 1, rc = lc | 1;
        build(lc, b, half, a);
        build(rc, b + half, half, a);

        t[node].resize(len);
        for (int y = 0; y < half; y++){
            t[node][y] = merge(t[lc][y], t[rc][y]);
            t[node][y + half] = merge(t[rc][y], t[lc][y]);
        }
    }

    /// x < len, b <= l <= r < b + len; folds a[p ^ x] over p in [l, r], never visiting an empty range
    T query(int node, int b, int len, int l, int r, int x) const{
        if (l == b && r == b + len - 1){
            if constexpr (is_same_v<Merge, plus<T>>) return added[node] + t[node][x];
            else return t[node][x];
        }

        /// a set bit half of x sends positions [b, mid) to the right child and [mid, b + len) to the left one
        int half = len >> 1, mid = b + half, sw = x & half ? half : 0;
        int lo = (node << 1) | (sw > 0), hi = lo ^ 1;
        x &= half - 1;
        if (r < mid) return query(lo, b + sw, half, l + sw, r + sw, x);
        if (l >= mid) return query(hi, b + half - sw, half, l - sw, r - sw, x);
        return merge(query(lo, b + sw, half, l + sw, mid - 1 + sw, x), query(hi, b + half - sw, half, mid - sw, r - sw, x));
    }
};

int main(){
    XorSegmentTree<long long> st({5, 1, 4, 1, 5, 9, 2, 6});
    assert(st.query(3, 3, 0) == 1);             /// x = 0 is a plain range sum
    assert(st.query(2, 5, 5) == 14);            /// p in 2..5 reads a[7], a[6], a[1], a[0]: 6 + 2 + 1 + 5
    assert(st.query(0, 7, 3) == 33);            /// the whole array, any x

    st.add(6, 10);                              /// a[6] = 12
    assert(st.query(0, 0, 6) == 12);
    assert(st.query(4, 7, 4) == 11);            /// a[0..3]

    XorSegmentTree concat(vector<string>{"a", "b", "c", "d"}, [](const string& x, const string& y){ return x + y; });
    assert(concat.query(0, 3, 1) == "badc");
    assert(concat.query(1, 2, 2) == "da");
    assert(concat.query(1, 3, 3) == "cba");
    assert(concat.query(2, 2, 0) == "c");
    return 0;
}
