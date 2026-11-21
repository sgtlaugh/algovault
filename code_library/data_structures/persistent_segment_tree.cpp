/***
 *
 * Persistent Segment Tree
 * Updates that each create a new version, point or range with lazy propagation, with range queries on any version
 *
 * Complexity:
 *   - build: O(n)
 *   - add, set: O(log n) time and O(log n) new nodes
 *   - get, query, kth: O(log n)
 *   - RangeKth: O(n log n) build time and memory, O(log n) per kth
 *   - PersistentLazySegmentTree: O(log n) per call, update and copy add at most 4 ceil(log2 n) + 4 nodes, query none
 *
 * PersistentSegmentTree<T> tree(values): version 0 holds values, ranges fold with +
 * PersistentSegmentTree<T, Combine> tree(values, combine): ranges fold with combine, associative, const operator()
 * tree.add(v, i, x): id of a new version equal to version v with x added at position i
 * tree.set(v, i, x): id of a new version equal to version v with position i replaced by x
 * tree.get(v, i): the value at position i in version v
 * tree.kth(lo, hi, k): smallest p whose prefix sum [0, p] of version hi minus version lo exceeds k, Combine = plus only
 * tree.query(v, l, r): combine over positions l..r in version v, 0 <= l <= r < n
 * tree.versions(): number of versions, the ids are 0..versions() - 1
 * RangeKth<T> rk(values): k-th smallest queries over values
 * rk.kth(l, r, k): the k-th smallest of values[l..r], 0-based, 0 <= k <= r - l
 * PersistentLazySegmentTree<T, Ops = PersistentAddSum<T>, L = T> seg(values): version 0 holds values, n >= 1
 * seg.copy(v, src, l, r): id of a new version equal to version v with positions l..r taken from version src
 * seg.update(v, l, r, lz): id of a new version equal to version v with lz applied to positions l..r
 * seg.query(v, l, r): combine over positions l..r in version v, 0 <= l <= r < n
 * seg.versions(): number of versions, the ids are 0..versions() - 1
 *
 * Notes:
 *   - kth needs 0 <= k < the total difference and no negative position difference, as with counts on prefix versions
 *   - A persistent array uses only set and get, with a combine that cannot overflow T, like FirstCombine in main()
 *   - PersistentLazySegmentTree takes the Ops of segment_tree.cpp, range add and range sum by default:
 *       - combine(a, b): fold of two adjacent ranges, a on the left
 *       - apply(val, lz, len): value of a range of len elements after lz
 *       - compose(old, new): one lazy meaning old then new
 *   - L must be default constructible, and a version may update or copy from any older one, so history can branch
 *   - A node takes 9 + sizeof(T) + sizeof(L) bytes, about 240 MB after 1e5 range affine updates at n = 1e5
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template<typename T>
struct PersistentAddSum{
    static T combine(T a, T b){ return a + b; }
    static T apply(T val, T lz, int len){ return val + lz * len; }
    static T compose(T old_lz, T new_lz){ return old_lz + new_lz; }
};

template<typename T, typename Ops = PersistentAddSum<T>, typename L = T>
struct PersistentLazySegmentTree{
    int n;
    vector<int> roots, lc, rc;
    vector<T> val;
    vector<L> lazy;
    vector<char> has_lazy;

    PersistentLazySegmentTree(const vector<T>& values) : n(values.size()){
        assert(n >= 1);
        roots.push_back(build(values, 0, n - 1));
    }

    int copy(int version, int source, int l, int r){
        assert(0 <= version && version < versions() && 0 <= source && source < versions());
        assert(0 <= l && l <= r && r < n);
        roots.push_back(copy(roots[version], nullptr, roots[source], nullptr, 0, n - 1, l, r));
        return roots.size() - 1;
    }

    int update(int version, int l, int r, const L& lz){
        assert(0 <= version && version < versions() && 0 <= l && l <= r && r < n);
        roots.push_back(update(roots[version], nullptr, 0, n - 1, l, r, lz));
        return roots.size() - 1;
    }

    T query(int version, int l, int r) const{
        assert(0 <= version && version < versions() && 0 <= l && l <= r && r < n);
        return query(roots[version], nullptr, 0, n - 1, l, r);
    }

    int versions() const{
        return roots.size();
    }

    int build(const vector<T>& values, int a, int b){
        if (a == b) return new_node(-1, -1, values[a]);

        int m = (a + b) / 2;
        int l = build(values, a, m), r = build(values, m + 1, b);
        return new_node(l, r, Ops::combine(val[l], val[r]));
    }

    /// Nodes are never written after creation, a lazy owed from above travels down as pending instead of copying
    int copy(int node, const L* pending, int src, const L* src_pending, int a, int b, int l, int r){
        if (b < l || r < a) return with_lazy(node, a, b, pending);
        if (l <= a && b <= r) return with_lazy(src, a, b, src_pending);

        L below, src_below;
        const L* down = pending_below(node, pending, below);
        const L* src_down = pending_below(src, src_pending, src_below);
        int m = (a + b) / 2;
        int left = copy(lc[node], down, lc[src], src_down, a, m, l, r);
        int right = copy(rc[node], down, rc[src], src_down, m + 1, b, l, r);
        return new_node(left, right, Ops::combine(val[left], val[right]));
    }

    int update(int node, const L* pending, int a, int b, int l, int r, const L& lz){
        if (b < l || r < a) return with_lazy(node, a, b, pending);
        if (l <= a && b <= r){
            int fresh = with_lazy(node, a, b, pending);
            if (fresh == node) fresh = clone(node);
            apply_lazy(fresh, a, b, lz);
            return fresh;
        }

        L below;
        const L* down = pending_below(node, pending, below);
        int m = (a + b) / 2;
        int left = update(lc[node], down, a, m, l, r, lz);
        int right = update(rc[node], down, m + 1, b, l, r, lz);
        return new_node(left, right, Ops::combine(val[left], val[right]));
    }

    T query(int node, const L* pending, int a, int b, int l, int r) const{
        if (l <= a && b <= r) return pending ? Ops::apply(val[node], *pending, b - a + 1) : val[node];

        L below;
        const L* down = pending_below(node, pending, below);
        int m = (a + b) / 2;
        if (r <= m) return query(lc[node], down, a, m, l, r);
        if (l > m) return query(rc[node], down, m + 1, b, l, r);
        return Ops::combine(query(lc[node], down, a, m, l, r), query(rc[node], down, m + 1, b, l, r));
    }

    /// val[node] already includes lazy[node] but its children still owe it, then whatever is pending from above
    const L* pending_below(int node, const L* pending, L& below) const{
        if (!has_lazy[node]) return pending;
        below = pending ? Ops::compose(lazy[node], *pending) : lazy[node];
        return &below;
    }

    int with_lazy(int node, int a, int b, const L* pending){
        if (!pending) return node;
        int fresh = clone(node);
        apply_lazy(fresh, a, b, *pending);
        return fresh;
    }

    /// Only on a node made by the current call, so no older version sees the write
    void apply_lazy(int node, int a, int b, const L& lz){
        val[node] = Ops::apply(val[node], lz, b - a + 1);
        if (a == b) return;
        lazy[node] = has_lazy[node] ? Ops::compose(lazy[node], lz) : lz;
        has_lazy[node] = 1;
    }

    int clone(int node){
        int fresh = new_node(lc[node], rc[node], val[node]);
        lazy[fresh] = lazy[node], has_lazy[fresh] = has_lazy[node];
        return fresh;
    }

    int new_node(int l, int r, const T& x){
        lc.push_back(l), rc.push_back(r), val.push_back(x), lazy.emplace_back(), has_lazy.push_back(0);
        return val.size() - 1;
    }
};

template<typename T, typename Combine = plus<T>>
struct PersistentSegmentTree{
    int n;
    Combine combine;
    vector<int> roots, lc, rc;
    vector<T> val;

    PersistentSegmentTree(const vector<T>& values, Combine combine = Combine()) : n(values.size()), combine(combine){
        lc.reserve(2 * n), rc.reserve(2 * n), val.reserve(2 * n);
        roots.push_back(n ? build(values, 0, n - 1) : -1);
    }

    int add(int version, int i, T x){
        assert(0 <= version && version < versions() && 0 <= i && i < n);
        roots.push_back(update(roots[version], 0, n - 1, i, x, true));
        return roots.size() - 1;
    }

    int set(int version, int i, T x){
        assert(0 <= version && version < versions() && 0 <= i && i < n);
        roots.push_back(update(roots[version], 0, n - 1, i, x, false));
        return roots.size() - 1;
    }

    T get(int version, int i) const{
        assert(0 <= version && version < versions() && 0 <= i && i < n);
        int node = roots[version], a = 0, b = n - 1;
        while (a != b){
            int m = (a + b) / 2;
            if (i <= m) node = lc[node], b = m;
            else node = rc[node], a = m + 1;
        }

        return val[node];
    }

    int kth(int lo_version, int hi_version, T k) const{
        static_assert(is_same_v<Combine, plus<T>>, "kth descends by subtracting subtree sums, only valid for Combine = plus");
        assert(0 <= lo_version && lo_version < versions() && 0 <= hi_version && hi_version < versions());
        int lo = roots[lo_version], hi = roots[hi_version], a = 0, b = n - 1;
        assert(n > 0 && 0 <= k && k < val[hi] - val[lo]);

        while (a != b){
            T inside = val[lc[hi]] - val[lc[lo]];
            int m = (a + b) / 2;
            if (k < inside) lo = lc[lo], hi = lc[hi], b = m;
            else k -= inside, lo = rc[lo], hi = rc[hi], a = m + 1;
        }

        return a;
    }

    T query(int version, int l, int r) const{
        assert(0 <= version && version < versions() && 0 <= l && l <= r && r < n);
        return query(roots[version], 0, n - 1, l, r);
    }

    T query(int node, int a, int b, int l, int r) const{
        if (l <= a && b <= r) return val[node];

        int m = (a + b) / 2;
        if (r <= m) return query(lc[node], a, m, l, r);
        if (l > m) return query(rc[node], m + 1, b, l, r);
        return combine(query(lc[node], a, m, l, r), query(rc[node], m + 1, b, l, r));
    }

    int versions() const{
        return roots.size();
    }

    int build(const vector<T>& values, int a, int b){
        if (a == b) return new_node(-1, -1, values[a]);

        int m = (a + b) / 2;
        int l = build(values, a, m), r = build(values, m + 1, b);
        return new_node(l, r, combine(val[l], val[r]));
    }

    /// every node of the old path is copied, never written, so older versions stay intact
    int update(int node, int a, int b, int i, T x, bool add){
        if (a == b) return new_node(-1, -1, add ? val[node] + x : x);

        int m = (a + b) / 2, l = lc[node], r = rc[node];
        if (i <= m) l = update(l, a, m, i, x, add);
        else r = update(r, m + 1, b, i, x, add);
        return new_node(l, r, combine(val[l], val[r]));
    }

    int new_node(int l, int r, T x){
        lc.push_back(l), rc.push_back(r), val.push_back(x);
        return val.size() - 1;
    }
};

template<typename T>
struct RangeKth{
    vector<T> sorted_values;
    PersistentSegmentTree<int> counts;

    /// Version i counts the ranks of values[0..i - 1], so [l, r] is version r + 1 minus version l
    RangeKth(const vector<T>& values) : sorted_values(distinct(values)), counts(vector<int>(sorted_values.size(), 0)){
        for (int i = 0; i < (int)values.size(); i++){
            int rank = lower_bound(sorted_values.begin(), sorted_values.end(), values[i]) - sorted_values.begin();
            counts.add(i, rank, 1);
        }
    }

    T kth(int l, int r, int k) const{
        assert(0 <= l && l <= r && r + 1 < counts.versions() && 0 <= k && k <= r - l);
        return sorted_values[counts.kth(l, r + 1, k)];
    }

    static vector<T> distinct(vector<T> values){
        sort(values.begin(), values.end());
        values.erase(unique(values.begin(), values.end()), values.end());
        return values;
    }
};

int main(){
    struct MinCombine{
        long long operator()(long long a, long long b) const{
            return min(a, b);
        }
    };

    struct FirstCombine{
        long long operator()(long long a, long long) const{
            return a;
        }
    };

    /// Range sums over versions
    PersistentSegmentTree<long long> tree({3, 1, 4, 1, 5});
    int v1 = tree.add(0, 2, 10);        /// 3 1 14 1 5
    int v2 = tree.set(v1, 4, 0);        /// 3 1 14 1 0
    assert(tree.versions() == 3);
    assert(tree.query(0, 0, 4) == 14);  /// version 0 is untouched
    assert(tree.query(v1, 0, 4) == 24);
    assert(tree.query(v2, 2, 4) == 15);
    assert(tree.get(v2, 2) == 14);
    assert(tree.kth(0, v1, 3) == 2);    /// version v1 minus version 0 is 10 at position 2

    /// Range minimum
    PersistentSegmentTree<long long, MinCombine> mins({5, 2, 8, 6});
    int m1 = mins.set(0, 1, 9);
    assert(mins.query(0, 0, 3) == 2);
    assert(mins.query(m1, 0, 3) == 5);

    /// Persistent array: FirstCombine never adds, so values at the long long extremes are safe
    PersistentSegmentTree<long long, FirstCombine> arr({LLONG_MAX, 7, LLONG_MIN});
    int a1 = arr.set(0, 1, LLONG_MAX);
    assert(arr.get(0, 1) == 7);
    assert(arr.get(a1, 1) == LLONG_MAX);

    /// k-th smallest in a range
    RangeKth<int> rk({5, 1, 4, 1, 3, 9, 2});
    assert(rk.kth(0, 6, 2) == 2);  /// sorted 1 1 2 3 4 5 9
    assert(rk.kth(2, 4, 1) == 3);  /// sorted 1 3 4

    /// Range add over versions, then a range copied from one version into another
    PersistentLazySegmentTree<long long> seg({1, 2, 3, 4, 5});
    int s1 = seg.update(0, 1, 3, 10);      /// 1 12 13 14 5
    int s2 = seg.copy(0, s1, 2, 4);        /// 1 2 13 14 5
    assert(seg.query(s1, 0, 4) == 45);
    assert(seg.query(0, 0, 4) == 15);      /// version 0 is untouched
    assert(seg.query(s2, 0, 4) == 35);
    assert(seg.query(s2, 1, 2) == 15);

    /// Range set and range min through an Ops struct
    struct MinSet{
        static int combine(int a, int b){ return min(a, b); }
        static int apply(int, int lz, int){ return lz; }
        static int compose(int, int lz){ return lz; }
    };

    PersistentLazySegmentTree<int, MinSet> mn({5, 3, 8});
    int t1 = mn.update(0, 0, 1, 9);        /// 9 9 8
    assert(mn.query(t1, 0, 2) == 8);
    assert(mn.query(0, 0, 2) == 3);

    return 0;
}
