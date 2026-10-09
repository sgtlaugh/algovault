/***
 *
 * Persistent Segment Tree
 * Point updates that each create a new version, range queries and k-th search on any version
 *
 * Complexity: O(n) to build, O(log n) time and O(log n) new nodes per update, O(log n) per query
 * Memory: O(n + u log n) after u updates, nodes live in index pools (vectors), no raw new
 *
 * PersistentSegmentTree<T, Merge = plus<T>> tree(values), 0-based positions, ranges [l, r] inclusive
 * Version 0 is the initial array, add and set return the id of the version they create
 * tree.add(v, i, x): new version equal to version v with x added at position i
 * tree.set(v, i, x): new version equal to version v with position i replaced by x
 * tree.get(v, i): the value at position i in version v
 * tree.query(v, l, r): Merge over [l, r] in version v, the sum by default
 * tree.kth(lo, hi, k): smallest position p with sum over [0, p] of (version hi - version lo) > k
 *   needs Merge = plus and nonnegative differences, as with counts on prefix versions
 *
 * Merge is any associative functor with a const operator(), min or max work too
 * Persistent array: use set and get only, with a Merge that cannot overflow T for arbitrary values
 *   struct FirstMerge{ long long operator()(long long a, long long) const{ return a; } };
 *   PersistentSegmentTree<long long, FirstMerge> arr(values);
 *
 * RangeKth<T> rk(values): rk.kth(l, r, k) is the k-th smallest value in [l, r], 0-based k <= r - l
 *   version i counts the compressed ranks of values[0..i-1], so [l, r] is version r + 1 minus version l
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template<typename T, typename Merge = plus<T>>
struct PersistentSegmentTree{
    int n;
    Merge merge;
    vector<int> roots, lc, rc;
    vector<T> val;

    PersistentSegmentTree(const vector<T>& values, Merge merge = Merge()) : n(values.size()), merge(merge){
        lc.reserve(2 * n), rc.reserve(2 * n), val.reserve(2 * n);
        roots.push_back(n ? build(values, 0, n - 1) : -1);
    }

    int versions() const{
        return roots.size();
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

    T query(int version, int l, int r) const{
        assert(0 <= version && version < versions() && 0 <= l && l <= r && r < n);
        return query(roots[version], 0, n - 1, l, r);
    }

    int kth(int lo_version, int hi_version, T k) const{
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

private:
    int new_node(int l, int r, T x){
        lc.push_back(l), rc.push_back(r), val.push_back(x);
        return val.size() - 1;
    }

    int build(const vector<T>& values, int a, int b){
        if (a == b) return new_node(-1, -1, values[a]);

        int m = (a + b) / 2;
        int l = build(values, a, m), r = build(values, m + 1, b);
        return new_node(l, r, merge(val[l], val[r]));
    }

    /// every node of the old path is copied, never written, so older versions stay intact
    int update(int node, int a, int b, int i, T x, bool add){
        if (a == b) return new_node(-1, -1, add ? val[node] + x : x);

        int m = (a + b) / 2, l = lc[node], r = rc[node];
        if (i <= m) l = update(l, a, m, i, x, add);
        else r = update(r, m + 1, b, i, x, add);
        return new_node(l, r, merge(val[l], val[r]));
    }

    T query(int node, int a, int b, int l, int r) const{
        if (l <= a && b <= r) return val[node];

        int m = (a + b) / 2;
        if (r <= m) return query(lc[node], a, m, l, r);
        if (l > m) return query(rc[node], m + 1, b, l, r);
        return merge(query(lc[node], a, m, l, r), query(rc[node], m + 1, b, l, r));
    }
};

template<typename T>
struct RangeKth{
    vector<T> sorted_values;
    PersistentSegmentTree<int> counts;

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

private:
    static vector<T> distinct(vector<T> values){
        sort(values.begin(), values.end());
        values.erase(unique(values.begin(), values.end()), values.end());
        return values;
    }
};

struct MinMerge{
    long long operator()(long long a, long long b) const{
        return min(a, b);
    }
};

struct FirstMerge{
    long long operator()(long long a, long long) const{
        return a;
    }
};

int main(){
    PersistentSegmentTree<long long> tree({3, 1, 4, 1, 5});
    int v1 = tree.add(0, 2, 10);
    int v2 = tree.set(0, 0, -7);
    int v3 = tree.set(v1, 4, 0);
    assert(v1 == 1 && v2 == 2 && v3 == 3 && tree.versions() == 4);
    assert(tree.query(0, 0, 4) == 14 && tree.query(0, 1, 3) == 6 && tree.query(0, 4, 4) == 5);
    assert(tree.query(v1, 0, 4) == 24 && tree.query(v1, 2, 2) == 14 && tree.query(v1, 0, 1) == 4);
    assert(tree.query(v2, 0, 4) == 4 && tree.query(v2, 0, 0) == -7 && tree.get(v2, 2) == 4);
    assert(tree.query(v3, 0, 4) == 19 && tree.get(v3, 4) == 0 && tree.get(v3, 2) == 14);
    assert(tree.get(0, 0) == 3 && tree.get(0, 2) == 4 && tree.get(0, 4) == 5);

    PersistentSegmentTree<int> counts({0, 0, 0, 0});
    counts.add(0, 2, 1), counts.add(1, 0, 1), counts.add(2, 2, 1), counts.add(3, 3, 1);
    assert(counts.kth(0, 4, 0) == 0 && counts.kth(0, 4, 1) == 2 && counts.kth(0, 4, 2) == 2 && counts.kth(0, 4, 3) == 3);
    assert(counts.kth(1, 3, 0) == 0 && counts.kth(1, 3, 1) == 2 && counts.kth(2, 4, 1) == 3);

    PersistentSegmentTree<long long, MinMerge> mins({5, 2, 8, 6});
    int m1 = mins.set(0, 1, 9);
    assert(mins.query(0, 0, 3) == 2 && mins.query(m1, 0, 3) == 5 && mins.query(m1, 1, 3) == 6 && mins.query(0, 2, 3) == 6);

    PersistentSegmentTree<long long, FirstMerge> arr({LLONG_MAX, LLONG_MAX, 0});
    int a1 = arr.set(0, 2, LLONG_MIN), a2 = arr.set(a1, 0, -1);
    assert(arr.get(0, 2) == 0 && arr.get(a1, 2) == LLONG_MIN && arr.get(a1, 0) == LLONG_MAX && arr.get(a2, 0) == -1 && arr.get(a2, 1) == LLONG_MAX);
    assert(arr.query(0, 1, 2) == LLONG_MAX && arr.query(a2, 0, 2) == -1 && arr.query(a1, 2, 2) == LLONG_MIN);

    PersistentSegmentTree<int> single({42});
    int s1 = single.add(0, 0, -50);
    assert(single.get(0, 0) == 42 && single.get(s1, 0) == -8 && single.query(s1, 0, 0) == -8);

    RangeKth<int> rk({5, 1, 4, 1, 3, 9, 2});
    assert(rk.kth(0, 6, 0) == 1 && rk.kth(0, 6, 1) == 1 && rk.kth(0, 6, 2) == 2 && rk.kth(0, 6, 6) == 9);
    assert(rk.kth(2, 4, 0) == 1 && rk.kth(2, 4, 1) == 3 && rk.kth(2, 4, 2) == 4);
    assert(rk.kth(5, 5, 0) == 9 && rk.kth(0, 0, 0) == 5);

    RangeKth<long long> extremes({-5, LLONG_MIN, LLONG_MAX, 0});
    assert(extremes.kth(0, 3, 0) == LLONG_MIN && extremes.kth(0, 3, 3) == LLONG_MAX && extremes.kth(2, 3, 0) == 0);

    RangeKth<int> empty({});
    PersistentSegmentTree<int> none({});
    assert(empty.counts.versions() == 1 && none.versions() == 1);

    return 0;
}
