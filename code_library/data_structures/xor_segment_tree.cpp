/***
 *
 * XOR Segment Tree
 * Fold of a[p ^ x] over l <= p <= r for any x, with point add when the fold is a sum
 *
 * Complexity: O(n log n) time and memory to build, O(log n) per query and per add
 *
 * XorSegmentTree<T> st(a): sums of a, a.size() a power of two
 * XorSegmentTree<T, Combine> st(a, combine): folds with combine, which must be associative and have a const operator() but need not commute
 * st.query(l, r, x): combine of a[p ^ x] over l <= p <= r in p order, 0 <= l <= r < n, 0 <= x < n
 * st.add(i, v): a[i] += v, only with the default Combine = plus<T>, every other combine builds a static tree
 *
 * Pad a to a power of two by hand, the padding values depend on what a padded a[p ^ x] should mean
 * For sums every stored fold of a node is the same block sum, the n log n layout pays off for order-dependent folds like max subarray sum (CF 1716E)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T, typename Combine = plus<T>>
struct XorSegmentTree{
    int n;
    Combine combine;
    vector<vector<T>> t;  /// t[node][y] = combine of a[p ^ y] in p order over the node's block, 0 <= y < block length, adds excluded
    vector<T> added;      /// added[node] = sum of the add() values at positions in the node's block

    XorSegmentTree(const vector<T>& a, Combine combine = Combine()) : n((int)a.size()), combine(combine), t(2 * a.size()), added(2 * a.size()){
        assert(n > 0 && (n & (n - 1)) == 0);
        build(1, 0, n, a);
    }

    void add(int i, T v){
        static_assert(is_same_v<Combine, plus<T>>, "add relies on every node value being the plain block sum");
        assert(0 <= i && i < n);
        add(1, 0, n, i, v);
    }

    T query(int l, int r, int x) const{
        assert(0 <= l && l <= r && r < n && 0 <= x && x < n);
        return query(1, 0, n, l, r, x);
    }

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
            t[node][y] = combine(t[lc][y], t[rc][y]);
            t[node][y + half] = combine(t[rc][y], t[lc][y]);
        }
    }

    /// x < len, b <= l <= r < b + len; folds a[p ^ x] over p in [l, r], never visiting an empty range
    T query(int node, int b, int len, int l, int r, int x) const{
        if (l == b && r == b + len - 1){
            if constexpr (is_same_v<Combine, plus<T>>) return added[node] + t[node][x];
            else return t[node][x];
        }

        /// a set bit half of x sends positions [b, mid) to the right child and [mid, b + len) to the left one
        int half = len >> 1, mid = b + half, shift = x & half ? half : 0;
        int left_src = (node << 1) | (shift > 0), right_src = left_src ^ 1;  /// the children read by [b, mid) and by [mid, b + len)
        x &= half - 1;
        if (r < mid) return query(left_src, b + shift, half, l + shift, r + shift, x);
        if (l >= mid) return query(right_src, b + half - shift, half, l - shift, r - shift, x);
        return combine(query(left_src, b + shift, half, l + shift, mid - 1 + shift, x),
                     query(right_src, b + half - shift, half, mid - shift, r - shift, x));
    }
};

int main(){
    /// Sums of a[p ^ x], then a point add
    XorSegmentTree<long long> st({5, 1, 4, 1, 5, 9, 2, 6});
    assert(st.query(2, 5, 5) == 14);            /// p in 2..5 reads a[7], a[6], a[1], a[0]: 6 + 2 + 1 + 5
    assert(st.query(0, 7, 3) == 33);            /// the whole array, any x
    st.add(6, 10);                              /// a[6] = 12
    assert(st.query(0, 0, 6) == 12);
    assert(st.query(4, 7, 4) == 11);            /// a[0..3]

    /// Concatenation does not commute, so the fold must follow p order
    XorSegmentTree concat(vector<string>{"a", "b", "c", "d"}, [](const string& x, const string& y){ return x + y; });
    assert(concat.query(0, 3, 1) == "badc");
    assert(concat.query(1, 2, 2) == "da");
    assert(concat.query(1, 3, 3) == "cba");

    return 0;
}
