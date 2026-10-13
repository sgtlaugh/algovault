/***
 *
 * Run-Length Implicit Treap
 * Implicit treap whose nodes are runs of equal (or arithmetic progression) values, for sequences up to 1e18 long
 *
 * Complexity: expected O(log m) per operation, m = number of runs, at most the initial runs plus 3 per operation
 *
 * RunLengthImplicitTreap<Ops = SumRuns>: long long positions, 0-based, ranges [l, r] inclusive
 *   RunLengthImplicitTreap(len, v): len copies of v, RunLengthImplicitTreap(runs): from {value, len} pairs in O(n)
 *   insert(pos, v, len): len copies of v before position pos, pos = size() appends, erase(l, r), size()
 *   assign(l, r, v): [l, r] becomes a single run of v
 *   apply(l, r, lz), reverse(l, r), query(l, r): aggregate of [l, r], get(pos): the value at pos as a run of length 1
 *   max_right(l, pred), min_left(r, pred): as in segment_tree.cpp, pred monotone and also tried on parts of a run
 *   for_each_run(l, r, f): f(v, len) for the runs covering [l, r] in order, cut at l and r
 *   The total length must stay below 4.6e18
 *
 * Ops: V run value, S aggregate (default constructible), L lazy, derive from RunOps for the defaults
 *   summarize(v, len): aggregate of a run of len values described by v
 *   merge(a, b): aggregate of two adjacent ranges, a on the left
 *   cut(v, k): {first k values, the rest} of a run, default {v, v} for equal values
 *   apply_run(v, lz, len), apply_agg(s, lz, len), compose(old, new): a lazy as in segment_tree.cpp, it must
 *       commute with reverse, and compose must return old for a lazy that changes nothing, default no lazy
 *   reverse_run(v, len), reverse_agg(s): the run and the aggregate read backwards, default unchanged
 * SumRuns: equal long long values with range add and range sum
 * Arithmetic progression runs: V = {a, d}, cut(v, k) = {v, {a + k * d, d}}, reverse_run(v, len) = {a + (len - 1) * d, -d}
 *
 * For fixed positions known in advance, a segment tree over compressed coordinates is simpler and faster,
 * this treap is for positions that shift (insert, erase, reverse) or arrive online
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// Defaults for Ops without a lazy whose runs and aggregates read the same backwards
struct RunOps{
    using L = char;

    template<typename V>
    static pair<V, V> cut(const V& v, long long){ return {v, v}; }

    template<typename V>
    static V apply_run(const V& v, L, long long){ return v; }

    template<typename S>
    static S apply_agg(const S& s, L, long long){ return s; }

    static L compose(L, L lz){ return lz; }

    template<typename V>
    static V reverse_run(const V& v, long long){ return v; }

    template<typename S>
    static S reverse_agg(const S& s){ return s; }
};

struct SumRuns : RunOps{
    using V = long long;
    using S = long long;
    using L = long long;

    static S summarize(V v, long long len){ return v * len; }
    static S merge(S a, S b){ return a + b; }
    static V apply_run(V v, L lz, long long){ return v + lz; }
    static S apply_agg(S s, L lz, long long len){ return s + lz * len; }
    static L compose(L old_lz, L new_lz){ return old_lz + new_lz; }
};

template<typename Ops = SumRuns>
struct RunLengthImplicitTreap{
    using V = typename Ops::V;
    using S = typename Ops::S;
    using L = typename Ops::L;

    struct Node{
        V v;
        S agg;
        L lz;
        long long len, size;
        unsigned int priority;
        int l, r;
        bool has_lz, rev;
    };

    vector<Node> nodes;
    int root = -1;
    mt19937 rng{(unsigned int)chrono::steady_clock::now().time_since_epoch().count()};

    RunLengthImplicitTreap(long long len, const V& v){
        assert(0 <= len && len <= LLONG_MAX / 2);
        if (len > 0) root = new_node(v, len);
    }

    RunLengthImplicitTreap(const vector<pair<V, long long>>& runs){
        /// Cartesian tree on random priorities with a stack, O(n)
        vector<int> stack;
        for (auto& [v, len] : runs){
            assert(len >= 1);
            int t = new_node(v, len), last = -1;
            while (!stack.empty() && nodes[stack.back()].priority < nodes[t].priority){
                last = stack.back();
                stack.pop_back();
                pull(last);
            }
            nodes[t].l = last;
            if (!stack.empty()) nodes[stack.back()].r = t;
            stack.push_back(t);
        }

        if (!stack.empty()) root = stack[0];
        while (!stack.empty()){
            pull(stack.back());
            stack.pop_back();
        }
        assert(size() <= LLONG_MAX / 2);
    }

    long long size() const{
        return size(root);
    }

    void insert(long long pos, const V& v, long long len){
        assert(0 <= pos && pos <= size() && len >= 1 && len <= LLONG_MAX / 2 - size());
        int a, b;
        split(root, pos, a, b);
        root = merge(merge(a, new_node(v, len)), b);
    }

    void erase(long long l, long long r){
        on_range(l, r, [&](int& t){ t = -1; });
    }

    void assign(long long l, long long r, const V& v){
        on_range(l, r, [&](int& t){ t = new_node(v, r - l + 1); });
    }

    void apply(long long l, long long r, const L& lz){
        on_range(l, r, [&](int& t){ apply_lazy(t, lz); });
    }

    void reverse(long long l, long long r){
        on_range(l, r, [&](int& t){ apply_reverse(t); });
    }

    S query(long long l, long long r){
        S res{};
        on_range(l, r, [&](int& t){ res = nodes[t].agg; });
        return res;
    }

    V get(long long pos){
        V res{};
        on_range(pos, pos, [&](int& t){ res = nodes[t].v; });
        return res;
    }

    template<typename F>
    void for_each_run(long long l, long long r, F f){
        on_range(l, r, [&](int& t){
            vector<int> stack;
            for (int x = t; x != -1 || !stack.empty(); ){
                for (; x != -1; x = nodes[x].l) push(x), stack.push_back(x);
                x = stack.back();
                stack.pop_back();
                f(nodes[x].v, nodes[x].len);
                x = nodes[x].r;
            }
        });
    }

    /// Largest r in [l - 1, size() - 1] with pred(query(l, r)) true, pred is never called on an empty range
    template<typename Pred>
    long long max_right(long long l, Pred pred){
        assert(0 <= l && l <= size());
        int a, b;
        split(root, l, a, b);
        long long taken = descend(b, pred, false);
        root = merge(a, b);
        return l + taken - 1;
    }

    /// Smallest l in [0, r + 1] with pred(query(l, r)) true, pred is never called on an empty range
    template<typename Pred>
    long long min_left(long long r, Pred pred){
        assert(-1 <= r && r < size());
        int a, b;
        split(root, r + 1, a, b);
        long long taken = descend(a, pred, true);
        root = merge(a, b);
        return r - taken + 1;
    }

    long long size(int t) const{
        return t == -1 ? 0 : nodes[t].size;
    }

    int new_node(const V& v, long long len){
        nodes.push_back({v, Ops::summarize(v, len), L(), len, len, (unsigned int)rng(), -1, -1, false, false});
        return nodes.size() - 1;
    }

    void apply_lazy(int t, const L& lz){
        if (t == -1) return;
        Node& u = nodes[t];
        u.v = Ops::apply_run(u.v, lz, u.len);
        u.agg = Ops::apply_agg(u.agg, lz, u.size);
        u.lz = u.has_lz ? Ops::compose(u.lz, lz) : lz;
        u.has_lz = true;
    }

    void apply_reverse(int t){
        if (t == -1) return;
        Node& u = nodes[t];
        swap(u.l, u.r);
        u.v = Ops::reverse_run(u.v, u.len);
        u.agg = Ops::reverse_agg(u.agg);
        u.rev ^= 1;
    }

    void push(int t){
        Node& u = nodes[t];
        if (u.has_lz){
            apply_lazy(u.l, u.lz);
            apply_lazy(u.r, u.lz);
            u.has_lz = false;
        }
        if (u.rev){
            apply_reverse(u.l);
            apply_reverse(u.r);
            u.rev = false;
        }
    }

    void pull(int t){
        Node& u = nodes[t];
        u.size = u.len + size(u.l) + size(u.r);
        S s = Ops::summarize(u.v, u.len);
        if (u.l != -1) s = Ops::merge(nodes[u.l].agg, s);
        if (u.r != -1) s = Ops::merge(s, nodes[u.r].agg);
        u.agg = s;
    }

    /// First k elements into a, the rest into b, cutting the run that straddles position k
    /// Children go through locals: new_node may reallocate nodes, so no reference into it is held across a call
    void split(int t, long long k, int& a, int& b){
        if (t == -1){
            a = b = -1;
            return;
        }

        push(t);
        long long left = size(nodes[t].l), len = nodes[t].len;
        if (k <= left){
            int x, y;
            split(nodes[t].l, k, x, y);
            nodes[t].l = y;
            a = x, b = t;
        }
        else if (k >= left + len){
            int x, y;
            split(nodes[t].r, k - left - len, x, y);
            nodes[t].r = x;
            a = t, b = y;
        }
        else{
            /// t keeps the first part of its run, the rest becomes a new node merged in front of t's right subtree
            /// The new node draws a fresh priority: inheriting t's would give every piece of a run the same one,
            /// and equal priorities chain up (200000 cuts of one run overflowed the stack)
            auto [first, rest] = Ops::cut(nodes[t].v, k - left);
            int u = new_node(rest, len - (k - left)), right = nodes[t].r;
            nodes[t].r = -1;
            nodes[t].v = first, nodes[t].len = k - left;
            a = t, b = merge(u, right);
        }
        pull(t);
    }

    int merge(int a, int b){
        if (a == -1 || b == -1) return a == -1 ? b : a;
        if (nodes[a].priority > nodes[b].priority){
            push(a);
            int r = merge(nodes[a].r, b);
            nodes[a].r = r;
            pull(a);
            return a;
        }
        push(b);
        int l = merge(a, nodes[b].l);
        nodes[b].l = l;
        pull(b);
        return b;
    }

    /// Applies f to the treap holding exactly [l, r], f may replace it, then joins everything back
    template<typename F>
    void on_range(long long l, long long r, F f){
        assert(0 <= l && l <= r && r < size());
        int a, b, c;
        split(root, l, a, b);
        split(b, r - l + 1, b, c);
        f(b);
        root = merge(merge(a, b), c);
    }

    /// Number of elements of t, from its left end (or its right end when backwards), whose fold passes pred
    template<typename Pred>
    long long descend(int t, Pred& pred, bool backwards){
        long long taken = 0;
        bool empty = true;
        S acc{};
        auto extend = [&](const S& s){ return empty ? s : backwards ? Ops::merge(s, acc) : Ops::merge(acc, s); };
        while (t != -1){
            push(t);
            int near = backwards ? nodes[t].r : nodes[t].l, far = backwards ? nodes[t].l : nodes[t].r;
            if (near != -1){
                S s = extend(nodes[near].agg);
                if (!pred(s)){
                    t = near;
                    continue;
                }
                acc = s, empty = false, taken += nodes[near].size;
            }

            const V v = nodes[t].v;
            long long len = nodes[t].len;
            S whole = extend(Ops::summarize(v, len));
            if (pred(whole)){
                acc = whole, empty = false, taken += len;
                t = far;
                continue;
            }

            /// The answer ends inside this run: binary search the longest part of it that still passes
            long long lo = 0, hi = len - 1;
            while (lo < hi){
                long long mid = lo + (hi - lo + 1) / 2;
                auto parts = Ops::cut(v, backwards ? len - mid : mid);
                if (pred(extend(Ops::summarize(backwards ? parts.second : parts.first, mid)))) lo = mid;
                else hi = mid - 1;
            }
            return taken + lo;
        }
        return taken;
    }
};

int main(){
    RunLengthImplicitTreap<> seq(10, 0);                  /// 0 0 0 0 0 0 0 0 0 0
    seq.assign(2, 5, 3);                                  /// 0 0 3 3 3 3 0 0 0 0
    seq.insert(4, 7, 2);                                  /// 0 0 3 3 7 7 3 3 0 0 0 0
    seq.apply(5, 9, 1);                                   /// 0 0 3 3 7 8 4 4 1 1 0 0
    assert(seq.query(0, 11) == 31);
    seq.reverse(2, 6);                                    /// 0 0 4 8 7 3 3 4 1 1 0 0
    assert(seq.get(3) == 8);
    assert(seq.query(2, 4) == 19);

    RunLengthImplicitTreap<> big(1000000000000000000LL, 1);
    big.assign(10, 19, 5);                                /// ten 1s, ten 5s, then 1s
    big.erase(0, 4);
    assert(big.size() == 999999999999999995LL);
    assert(big.query(0, 24) == 5 + 50 + 10);
    assert(big.max_right(0, [](long long s){ return s <= 20; }) == 7);   /// 1 1 1 1 1 5 5 5 sums to 20
    return 0;
}
