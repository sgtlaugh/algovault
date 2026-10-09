/***
 *
 * XOR Segment Tree
 * Range sum of a[p ^ x] over l <= p <= r for any x, with point add
 *
 * Complexity: O(n log n) time and memory to build, O(log n) per query and per add
 *
 * The size n must be a power of two (pad with zeros if needed), 0 <= l, r < n and 0 <= x < n
 * Every node over [b, b + len) stores t[node][y] = merge of a[p ^ y] over its range in p order, for all 0 <= y < len
 * For sums every t[node][y] equals the block sum, since xor by y < len only permutes the block
 * The layout only pays off for order-dependent merges such as max subarray sum (CF 1716E)
 * The query maps p -> p ^ x level by level: a set bit of x swaps the halves the query range lands in
 *
 * Usage:
 *   XorSegmentTree<long long> st(a);       // a.size() is a power of two
 *   st.add(i, v);                          // a[i] += v
 *   st.query(l, r, x);                     // sum of a[p ^ x] for l <= p <= r
 *
 * Related problem: https://codeforces.com/contest/1654/problem/F
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T>
struct XorSegmentTree{
    int n;
    vector<vector<T>> t;
    vector<T> added;

    XorSegmentTree(const vector<T>& a) : n((int)a.size()), t(2 * a.size()), added(2 * a.size()){
        assert(n > 0 && (n & (n - 1)) == 0);
        build(1, 0, n, a);
    }

    void add(int i, T v){
        assert(0 <= i && i < n);
        add(1, 0, n, i, v);
    }

    T query(int l, int r, int x){
        assert(0 <= l && l < n && 0 <= r && r < n && 0 <= x && x < n);
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
            t[node][y] = t[lc][y] + t[rc][y];
            t[node][y + half] = t[rc][y] + t[lc][y];
        }
    }

    /// x < len; sums a[p ^ x] for p in [l, r] intersected with [b, b + len)
    T query(int node, int b, int len, int l, int r, int x){
        int e = b + len - 1;
        if (l > r || r < b || l > e) return 0;
        if (l <= b && e <= r) return added[node] + t[node][x];

        int half = len >> 1, mid = b + half, lc = node << 1, rc = lc | 1;
        if (!(x & half)) return query(lc, b, half, l, r, x) + query(rc, mid, half, l, r, x);

        /// p ^ half moves the left part of the range into the right child and vice versa
        x ^= half;
        T from_left = query(rc, mid, half, max(l, b) + half, min(r, mid - 1) + half, x);
        T from_right = query(lc, b, half, max(l, mid) - half, min(r, e) - half, x);
        return from_left + from_right;
    }
};

int main(){
    XorSegmentTree<long long> small({1, 2, 3, 4});
    assert(small.query(0, 1, 1) == 3);
    assert(small.query(1, 2, 3) == 5);
    assert(small.query(0, 0, 2) == 3);
    for (int x = 0; x < 4; x++) assert(small.query(0, 3, x) == 10);

    XorSegmentTree<long long> st({5, 1, 4, 1, 5, 9, 2, 6});
    assert(st.query(2, 5, 5) == 14);
    assert(st.query(1, 6, 7) == 22);
    assert(st.query(3, 3, 0) == 1);

    st.add(6, 10);
    assert(st.query(1, 6, 7) == 32);
    assert(st.query(0, 0, 6) == 12);
    assert(st.query(4, 7, 4) == 11);

    XorSegmentTree<long long> single({7});
    assert(single.query(0, 0, 0) == 7);
    single.add(0, -3);
    assert(single.query(0, 0, 0) == 4);

    XorSegmentTree<long long> big({1000000000000000000LL, 1000000000000000000LL});
    assert(big.query(0, 1, 1) == 2000000000000000000LL);
    assert(big.query(1, 1, 1) == 1000000000000000000LL);

    return 0;
}
