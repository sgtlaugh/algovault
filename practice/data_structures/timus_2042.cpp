/***
 *
 * Segment Tree with Lazy Propagation
 * Range update, range query and binary search for any associative combine and lazy update
 *
 * Complexity: O(n) to build, O(log n) per update, query, max_right and min_left
 *
 * SegmentTree<T, Ops = AddSum<T>, L = T> st(n): n copies of T(), range add and range sum by default
 * SegmentTree<T, Ops, L> st(ar): built from ar
 * st.update(l, r, lz): apply lz to a[l..r], 0 <= l <= r < n
 * st.query(l, r): combine over a[l..r], 0 <= l <= r < n
 * st.max_right(l, pred): largest r with pred(query(l, r)) true, l - 1 if none, 0 <= l <= n
 * st.min_left(r, pred): smallest l with pred(query(l, r)) true, r + 1 if none, -1 <= r < n
 *
 * pred must be monotone: once false it stays false as the range grows
 * pred is never called on an empty range, so combine needs no identity anywhere
 *
 * Ops holds three static functions, L is the lazy type when it differs from T:
 *   combine(a, b): fold of two adjacent ranges, a on the left
 *   apply(val, lz, len): value of a range of len elements after lz
 *   compose(old, new): one lazy meaning old then new, a lazy that changes nothing must return old
 *
 * Arithmetic progression add (x, x + d, ... on [l, r]) with range sum: T = {sum, cnt, idx_sum}, L = {c0, c1}
 *   leaves are {ar[i], 1, i}, combine and compose add componentwise, apply adds c0 * cnt + c1 * idx_sum
 *   the update is update(l, r, {x - d * l, d})
 *   sums reach d * n^2 / 2, so large inputs need a modulus or __int128
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template<typename T>
struct AddSum{
    static T combine(T a, T b){ return a + b; }
    static T apply(T val, T lz, int len){ return val + lz * len; }
    static T compose(T old_lz, T new_lz){ return old_lz + new_lz; }
};

/// Inherits Ops so seg.combine(a, b) stays callable from outside, hld.cpp's query_path relies on it
template<typename T, typename Ops = AddSum<T>, typename L = T>
struct SegmentTree : Ops{
    int n;
    vector<T> tree;
    vector<L> lazy;
    vector<char> has_lazy;

    SegmentTree(int n) : SegmentTree(vector<T>(n)) {}

    SegmentTree(const vector<T>& ar) : n(ar.size()), tree(4 * n), lazy(4 * n), has_lazy(4 * n){
        if (n) build(ar, 1, 0, n - 1);
    }

    void update(int l, int r, L val){
        assert(0 <= l && l <= r && r < n);
        update(1, 0, n - 1, l, r, val);
    }

    T query(int l, int r){
        assert(0 <= l && l <= r && r < n);
        return query(1, 0, n - 1, l, r);
    }

    template<typename Pred>
    int max_right(int l, Pred pred){
        assert(0 <= l && l <= n);
        T acc{};
        bool empty = true;
        int pos = l < n ? first_fail(1, 0, n - 1, l, pred, acc, empty) : -1;
        return pos == -1 ? n - 1 : pos - 1;
    }

    template<typename Pred>
    int min_left(int r, Pred pred){
        assert(-1 <= r && r < n);
        T acc{};
        bool empty = true;
        int pos = r >= 0 ? last_fail(1, 0, n - 1, r, pred, acc, empty) : -1;
        return pos == -1 ? 0 : pos + 1;
    }

    void build(const vector<T>& ar, int idx, int a, int b){
        if (a == b){
            tree[idx] = ar[a];
            return;
        }

        int p = idx << 1, q = p | 1, c = (a + b) >> 1;
        build(ar, p, a, c);
        build(ar, q, c + 1, b);
        tree[idx] = Ops::combine(tree[p], tree[q]);
    }

    void push_lazy(int idx, const L& val){
        lazy[idx] = has_lazy[idx] ? Ops::compose(lazy[idx], val) : val;
        has_lazy[idx] = 1;
    }

    void propagate(int idx, int a, int b){
        if (!has_lazy[idx]) return;
        tree[idx] = Ops::apply(tree[idx], lazy[idx], b - a + 1);
        if (a != b){
            push_lazy(idx << 1, lazy[idx]);
            push_lazy(idx << 1 | 1, lazy[idx]);
        }
        has_lazy[idx] = 0;
    }

    void update(int idx, int a, int b, int l, int r, const L& val){
        if (a == l && b == r){
            push_lazy(idx, val);
            propagate(idx, a, b);
            return;
        }
        propagate(idx, a, b);

        /// the untouched child gets propagated too, so the combine below reads its true value
        int p = idx << 1, q = p | 1, c = (a + b) >> 1;
        if (r <= c){
            propagate(q, c + 1, b);
            update(p, a, c, l, r, val);
        }
        else if (l > c){
            propagate(p, a, c);
            update(q, c + 1, b, l, r, val);
        }
        else{
            update(p, a, c, l, c, val);
            update(q, c + 1, b, c + 1, r, val);
        }
        tree[idx] = Ops::combine(tree[p], tree[q]);
    }

    T query(int idx, int a, int b, int l, int r){
        propagate(idx, a, b);
        if (a == l && b == r) return tree[idx];

        int p = idx << 1, q = p | 1, c = (a + b) >> 1;
        if (r <= c) return query(p, a, c, l, r);
        if (l > c) return query(q, c + 1, b, l, r);
        return Ops::combine(query(p, a, c, l, c), query(q, c + 1, b, c + 1, r));
    }

    /// first position in [max(a, l), b] where pred(acc merged with values l..pos) fails, or -1, acc absorbs every passing node
    /// acc holds nothing until empty turns false, so the descent never needs an identity
    template<typename Pred>
    int first_fail(int idx, int a, int b, int l, Pred& pred, T& acc, bool& empty){
        if (b < l) return -1;
        propagate(idx, a, b);
        if (l <= a){
            T nxt = empty ? tree[idx] : Ops::combine(acc, tree[idx]);
            if (pred(nxt)){
                acc = nxt, empty = false;
                return -1;
            }
            if (a == b) return a;
        }

        int c = (a + b) >> 1;
        int res = first_fail(idx << 1, a, c, l, pred, acc, empty);
        return res != -1 ? res : first_fail(idx << 1 | 1, c + 1, b, l, pred, acc, empty);
    }

    /// mirror of first_fail scanning leftwards from r, combine order kept as values pos..r for non-commutative combines
    template<typename Pred>
    int last_fail(int idx, int a, int b, int r, Pred& pred, T& acc, bool& empty){
        if (a > r) return -1;
        propagate(idx, a, b);
        if (b <= r){
            T nxt = empty ? tree[idx] : Ops::combine(tree[idx], acc);
            if (pred(nxt)){
                acc = nxt, empty = false;
                return -1;
            }
            if (a == b) return a;
        }

        int c = (a + b) >> 1;
        int res = last_fail(idx << 1 | 1, c + 1, b, r, pred, acc, empty);
        return res != -1 ? res : last_fail(idx << 1, a, c, r, pred, acc, empty);
    }
};

const int MAXK = 50;

int K;

/// a palindrome of length <= K crossing a node boundary only sees K - 1 chars on each side, so nodes keep just those borders
struct Node{
    int cnt, len;
    char pre[MAXK], suf[MAXK];
};

int border(int len){
    return min(len, K - 1);
}

/// palindromes of length <= K in a + b that use chars from both sides, Manacher over the joined borders
int count_crossing(const Node& a, const Node& b){
    static char w[2 * MAXK];
    static int d1[2 * MAXK], d2[2 * MAXK];
    int sa = border(a.len), sb = border(b.len), m = sa + sb;
    memcpy(w, a.suf + (K - 1 - sa), sa);
    memcpy(w + sa, b.pre, sb);

    for (int i = 0, l = 0, r = -1; i < m; i++){
        int t = i > r ? 1 : min(d1[l + r - i], r - i + 1);
        while (i - t >= 0 && i + t < m && w[i - t] == w[i + t]) t++;
        d1[i] = t;
        if (i + t - 1 > r) l = i - t + 1, r = i + t - 1;
    }
    for (int i = 0, l = 0, r = -1; i < m; i++){
        int t = i > r ? 0 : min(d2[l + r - i + 1], r - i + 1);
        while (i - t - 1 >= 0 && i + t < m && w[i - t - 1] == w[i + t]) t++;
        d2[i] = t;
        if (i + t - 1 > r) l = i - t, r = i + t - 1;
    }

    int res = 0;
    for (int c = 0; c < m; c++){
        int lo = max(0, max(c - sa + 1, sa - c)), hi = min(d1[c] - 1, (K - 1) / 2);
        if (hi >= lo) res += hi - lo + 1;
        lo = max(1, max(c - sa + 1, sa - c + 1)), hi = min(d2[c], K / 2);
        if (hi >= lo) res += hi - lo + 1;
    }
    return res;
}

/// suf is right-aligned in its buffer so the last K - 1 chars always end at suf[K - 2]
struct PalindromeOps{
    static Node combine(const Node& a, const Node& b){
        Node res;
        res.len = a.len + b.len;
        res.cnt = a.cnt + b.cnt + count_crossing(a, b);

        int sa = border(a.len), pb = border(b.len), bl = border(res.len);
        int take = min(sa, bl);
        memcpy(res.pre, a.pre, take);
        memcpy(res.pre + take, b.pre, bl - take);

        take = min(pb, bl);
        memcpy(res.suf + (K - 1 - take), b.suf + (K - 1 - take), take);
        memcpy(res.suf + (K - 1 - bl), a.suf + (K - 1 - (bl - take)), bl - take);
        return res;
    }
    static Node apply(Node val, char c, int len){
        val.cnt = 0;
        for (int t = 1; t <= min(K, len); t++) val.cnt += len - t + 1;
        memset(val.pre, c, sizeof(val.pre));
        memset(val.suf, c, sizeof(val.suf));
        return val;
    }
    static char compose(char, char c){ return c; }
};

int main(){
    static char s[100005];
    scanf("%s %d", s, &K);
    int n = strlen(s);

    vector<Node> ar(n);
    for (int i = 0; i < n; i++){
        ar[i].len = 1;
        ar[i] = PalindromeOps::apply(ar[i], s[i], 1);
    }
    SegmentTree<Node, PalindromeOps, char> st(ar);
    /// the tree alone is 43 MB of the 64 MB limit
    vector<Node>().swap(ar);

    int m;
    scanf("%d", &m);
    while (m--){
        int t, l, r;
        scanf("%d %d %d", &t, &l, &r);
        if (t == 1){
            char c[2];
            scanf("%s", c);
            st.update(l - 1, r - 1, c[0]);
        }
        else printf("%d\n", st.query(l - 1, r - 1).cnt);
    }
    return 0;
}
