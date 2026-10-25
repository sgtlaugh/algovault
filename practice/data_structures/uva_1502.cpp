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

private:
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

struct ChmaxMax{
    static int combine(int a, int b){ return max(a, b); }
    static int apply(int val, int lz, int){ return max(val, lz); }
    static int compose(int old_lz, int new_lz){ return max(old_lz, new_lz); }
};

int main(){
    int t;

    scanf("%d", &t);
    static char buf[300010];
    for (int tc = 1; tc <= t; tc++){
        int n;
        scanf("%d", &n);
        vector<string> words(n);
        vector<int> w(n);
        for (int i = 0; i < n; i++){
            scanf("%s %d", buf, &w[i]);
            words[i] = buf;
        }

        vector<array<int, 26>> trie(1);
        trie[0].fill(0);
        for (auto& s : words){
            int cur = 0;
            for (char ch : s){
                int c = ch - 'a';
                if (!trie[cur][c]){
                    trie[cur][c] = trie.size();
                    trie.emplace_back();
                    trie.back().fill(0);
                }
                cur = trie[cur][c];
            }
        }

        int m = trie.size();
        vector<int> fail(m, 0), order = {0};
        vector<vector<int>> children(m);
        for (int f = 0; f < (int)order.size(); f++){
            int u = order[f];
            for (int c = 0; c < 26; c++){
                int& v = trie[u][c];
                if (!v) v = u ? trie[fail[u]][c] : 0;
                else{
                    fail[v] = u ? trie[fail[u]][c] : 0;
                    children[fail[v]].push_back(v);
                    order.push_back(v);
                }
            }
        }

        vector<int> tin(m), tout(m), stack = {0};
        for (int timer = 0; !stack.empty(); ){
            int u = stack.back();
            if (u >= 0){
                tin[u] = timer++;
                stack.back() = ~u;
                for (int v : children[u]) stack.push_back(v);
            }
            else{
                tout[~u] = timer - 1;
                stack.pop_back();
            }
        }

        SegmentTree<int, ChmaxMax> best(m);
        int res = 0;
        for (int i = 0; i < n; i++){
            if (w[i] <= 0) continue;
            int cur = 0, prev = 0;
            for (char ch : words[i]){
                cur = trie[cur][ch - 'a'];
                prev = max(prev, best.query(tin[cur], tin[cur]));
            }
            int dp = prev + w[i];
            res = max(res, dp);
            best.update(tin[cur], tout[cur], dp);
        }

        printf("Case #%d: %d\n", tc, res);
    }
    return 0;
}
