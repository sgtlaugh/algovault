/***
 *
 * Dynamic Segment Tree
 * Lazy segment tree over a huge index range that creates nodes only where updates reach, for online queries
 *
 * Complexity: O(log C) per update and query for C = hi - lo + 1, about 4 log2(C) new nodes per update, none per query
 *
 * DynamicSegmentTree<T, Ops = DynamicAddSum<T>, L = T> st(lo, hi, fill = T()): positions lo..hi, each holding fill
 * st.update(l, r, lz): apply lz to a[l..r], lo <= l <= r <= hi
 * st.query(l, r): combine over a[l..r], lo <= l <= r <= hi
 *
 * Notes:
 *   - Ops is SegmentTree's contract (combine, apply, compose), except that apply takes len as long long
 *   - combine(fill, fill) must equal fill: 0 for sums, an identity for composites, +INF for min
 *   - A start other than such a fill, like all ones under sums, is one update(lo, hi, ...) before the first query
 *   - Positions are long long in [-4e18, 4e18], and a sum over such lengths overflows fast: use a modulus or unsigned
 *   - Memory is st.nodes.size() * sizeof(Node), and st.nodes.reserve(k) avoids the peak of vector doubling
 *   - No max_right or min_left: compress the coordinates into a SegmentTree when binary search is needed offline
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template<typename T>
struct DynamicAddSum{
    static T combine(T a, T b){ return a + b; }
    static T apply(T val, T lz, long long len){ return val + lz * len; }
    static T compose(T old_lz, T new_lz){ return old_lz + new_lz; }
};

/// Converts only to long long: an Ops::apply written for SegmentTree's int len fails to compile instead of truncating
struct RangeLength{
    long long len;

    operator long long() const{ return len; }

    template<typename U>
    operator U() const = delete;
};

template<typename T, typename Ops = DynamicAddSum<T>, typename L = T>
struct DynamicSegmentTree : Ops{
    /// Child 0 is absent: every position under it holds fill, apart from lazies still pending at its parent
    struct Node{
        T val;
        L lz;
        int left = 0, right = 0;
        bool has_lazy = false;
    };

    long long lo, hi;
    T fill;
    vector<Node> nodes;

    DynamicSegmentTree(long long lo, long long hi, T fill = T()) : lo(lo), hi(hi), fill(fill), nodes(2){
        assert(lo <= hi);
        nodes[1].val = fill;
    }

    void update(long long l, long long r, L lz){
        assert(lo <= l && l <= r && r <= hi);
        update(1, lo, hi, l, r, lz);
    }

    T query(long long l, long long r){
        assert(lo <= l && l <= r && r <= hi);
        return query(1, lo, hi, l, r);
    }

    void update(int idx, long long a, long long b, long long l, long long r, const L& lz){
        if (l <= a && b <= r){
            push_lazy(idx, a, b, lz);
            return;
        }
        propagate(idx, a, b);

        long long c = a + (b - a) / 2;
        if (l <= c) update(child(idx, false), a, c, l, r, lz);
        if (r > c) update(child(idx, true), c + 1, b, l, r, lz);
        nodes[idx].val = Ops::combine(value(nodes[idx].left), value(nodes[idx].right));
    }

    /// Never pushes, so queries create no nodes: a node's own lazy is applied to the part of the range below it
    T query(int idx, long long a, long long b, long long l, long long r){
        if (!idx) return fill;
        if (l <= a && b <= r) return nodes[idx].val;

        long long c = a + (b - a) / 2;
        T res;
        if (r <= c) res = query(nodes[idx].left, a, c, l, r);
        else if (l > c) res = query(nodes[idx].right, c + 1, b, l, r);
        else res = Ops::combine(query(nodes[idx].left, a, c, l, r), query(nodes[idx].right, c + 1, b, l, r));

        if (nodes[idx].has_lazy) res = Ops::apply(res, nodes[idx].lz, RangeLength{min(b, r) - max(a, l) + 1});
        return res;
    }

    int child(int idx, bool right){
        int id = right ? nodes[idx].right : nodes[idx].left;
        if (id) return id;

        nodes.push_back({fill, L()});
        id = nodes.size() - 1;
        (right ? nodes[idx].right : nodes[idx].left) = id;
        return id;
    }

    T value(int idx) const{
        return idx ? nodes[idx].val : fill;
    }

    /// nodes[idx].val already includes nodes[idx].lz, which is still owed to the children
    void push_lazy(int idx, long long a, long long b, const L& lz){
        nodes[idx].val = Ops::apply(nodes[idx].val, lz, RangeLength{b - a + 1});
        if (a == b) return;
        nodes[idx].lz = nodes[idx].has_lazy ? Ops::compose(nodes[idx].lz, lz) : lz;
        nodes[idx].has_lazy = true;
    }

    void propagate(int idx, long long a, long long b){
        if (!nodes[idx].has_lazy) return;

        L lz = nodes[idx].lz;
        long long c = a + (b - a) / 2;
        push_lazy(child(idx, false), a, c, lz);
        push_lazy(child(idx, true), c + 1, b, lz);
        nodes[idx].has_lazy = false;
    }
};

int main(){
    /// Range add and range sum over positions 0..10^18, starting at 0
    DynamicSegmentTree<long long> st(0, 1000000000000000000LL);
    st.update(10, 19, 3);                                               /// ten 3s
    st.update(15, 1000000000000000000LL, 1);
    assert(st.query(0, 9) == 0);
    assert(st.query(10, 14) == 15);
    assert(st.query(15, 16) == 8);                                      /// 4 + 4
    assert(st.query(1000000000000000000LL, 1000000000000000000LL) == 1);

    /// Range assign and range min through an Ops struct, every position starting at 100
    struct MinSet{
        static int combine(int a, int b){ return min(a, b); }
        static int apply(int, int lz, long long){ return lz; }
        static int compose(int, int lz){ return lz; }
    };

    DynamicSegmentTree<int, MinSet> mn(-1000000000, 1000000000, 100);
    mn.update(-5, 5, 7);
    mn.update(0, 2, 9);                                                 /// -5..-1 and 3..5 hold 7, 0..2 hold 9
    assert(mn.query(0, 2) == 9);
    assert(mn.query(-1000000000, 1000000000) == 7);
    assert(mn.query(6, 1000000000) == 100);

    return 0;
}
