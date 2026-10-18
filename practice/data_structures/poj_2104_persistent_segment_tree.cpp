/***
 *
 * Persistent Segment Tree
 * Point updates that each create a new version, with range queries and k-th search on any version
 *
 * Complexity:
 *   - build: O(n)
 *   - add, set: O(log n) time and O(log n) new nodes
 *   - get, query, kth: O(log n)
 *   - RangeKth: O(n log n) build time and memory, O(log n) per kth
 *
 * PersistentSegmentTree<T> tree(values): version 0 holds values, ranges fold with +
 * PersistentSegmentTree<T, Combine> tree(values, combine): ranges fold with combine, any associative functor with a const operator()
 * tree.add(v, i, x): id of a new version equal to version v with x added at position i
 * tree.set(v, i, x): id of a new version equal to version v with position i replaced by x
 * tree.get(v, i): the value at position i in version v
 * tree.query(v, l, r): combine over positions l..r in version v, 0 <= l <= r < n
 * tree.versions(): number of versions, the ids are 0..versions() - 1
 * tree.kth(lo, hi, k): smallest position p whose prefix sum [0, p] of version hi minus version lo exceeds k, Combine = plus only
 * RangeKth<T> rk(values): k-th smallest queries over values
 * rk.kth(l, r, k): the k-th smallest of values[l..r], 0-based, 0 <= k <= r - l
 *
 * kth needs 0 <= k < the total difference and every position's difference nonnegative, as with counts on prefix versions
 * A persistent array uses only set and get, with a combine that cannot overflow T, like FirstCombine in main()
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

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

    int new_node(int l, int r, T x){
        lc.push_back(l), rc.push_back(r), val.push_back(x);
        return val.size() - 1;
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

    T query(int node, int a, int b, int l, int r) const{
        if (l <= a && b <= r) return val[node];

        int m = (a + b) / 2;
        if (r <= m) return query(lc[node], a, m, l, r);
        if (l > m) return query(rc[node], m + 1, b, l, r);
        return combine(query(lc[node], a, m, l, r), query(rc[node], m + 1, b, l, r));
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

int main(){
    int n, m;
    scanf("%d %d", &n, &m);
    vector<int> a(n);
    for (int i = 0; i < n; i++) scanf("%d", &a[i]);

    RangeKth<int> rk(a);
    while (m--){
        int i, j, k;
        scanf("%d %d %d", &i, &j, &k);
        printf("%d\n", rk.kth(i - 1, j - 1, k - 1));
    }
    return 0;
}
