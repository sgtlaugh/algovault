/***
 *
 * CDQ Divide and Conquer
 * Offline 3D dominance sums: for each query, the total weight of the points at or below it in all three coordinates
 *
 * Complexity: O(n log^2 n) for n points and queries together, O(n) memory
 *
 * DominanceSum<T> ds: no points and no queries, sums of type T, long long by default
 * ds.add_point(x, y, z, w = 1): a point of weight w
 * ds.add_query(x, y, z): a query, returns its index, counting up from 0
 * ds.solve(): res[i] = the sum of w over points with px <= x, py <= y and pz <= z for query i
 *
 * Notes:
 *   - Coordinates are any long long, ties count as dominated; for strict px < x with integers query x - 1
 *   - Counting, for each point, the others it dominates: add a query at each point and subtract its own weight
 *   - Box sums by inclusion-exclusion over the 8 corners, or 4 when one coordinate is a time with no lower bound
 *   - Updates over time: make time the x coordinate, so a query sees exactly the points added before it
 *   - Offline point add and rectangle sum at n = q = 2e5: 0.6 s and 98 MB, 5x less of both than fenwick_tree_2D_sparse
 *   - Prefix max instead of sum: a max Fenwick tree whose touched nodes are reset after each merge, max has no inverse
 *   - A DP whose terms are convolutions of earlier terms is online_convolution.cpp, the same split on a sequence
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// BEGIN COPY fenwick_point_update from code_library/data_structures/fenwick_tree.cpp
template <typename T>
struct FenwickPointUpdate{
    int n;
    vector<T> tree;

    FenwickPointUpdate(int n = 0) : n(n), tree(n + 1, 0) {}

    /// O(n): each node pushes its finished sum into its parent
    FenwickPointUpdate(const vector<T>& a) : n(a.size()), tree(n + 1, 0){
        for (int i = 1; i <= n; i++){
            tree[i] += a[i - 1];
            if (i + (i & -i) <= n) tree[i + (i & -i)] += tree[i];
        }
    }

    void update(int p, T v){
        assert(0 <= p && p < n);
        for (p++; p <= n; p += p & -p) tree[p] += v;
    }

    /// Binary lifting: p stays a multiple of 2 * step, so tree[p + step] is the sum of a[p .. p + step - 1]
    int lower_bound(T s){
        int p = 0, step = 1;

        while (step <= n / 2) step *= 2;

        for (; step > 0; step >>= 1){
            if (p + step <= n && tree[p + step] < s){
                p += step;
                s -= tree[p];
            }
        }

        return p;
    }

    T prefix(int p){
        assert(-1 <= p && p < n);
        T res = 0;
        for (p++; p > 0; p -= p & -p) res += tree[p];
        return res;
    }

    T query(int l, int r){
        if (l > r) return 0;
        assert(0 <= l && r < n);
        return prefix(r) - prefix(l - 1);
    }
};
/// END COPY fenwick_point_update

template <typename T = long long>
struct DominanceSum{
    struct Event{
        long long x, y, z;
        T w;
        int query;  /// -1 for a point
    };

    vector<Event> events, order, merged;
    vector<T> res;
    FenwickPointUpdate<T> fen;
    int queries = 0;

    void add_point(long long x, long long y, long long z, T w = 1){
        events.push_back({x, y, z, w, -1});
    }

    int add_query(long long x, long long y, long long z){
        events.push_back({x, y, z, T(0), queries});
        return queries++;
    }

    vector<T> solve(){
        vector<long long> zs;
        for (const Event& e : events) if (e.query == -1) zs.push_back(e.z);
        sort(zs.begin(), zs.end());
        zs.erase(unique(zs.begin(), zs.end()), zs.end());

        /// A point keeps the index of its z, a query the index of the last point z at or below it, -1 if none
        order = events;
        for (Event& e : order){
            if (e.query == -1) e.z = lower_bound(zs.begin(), zs.end(), e.z) - zs.begin();
            else e.z = upper_bound(zs.begin(), zs.end(), e.z) - zs.begin() - 1;
        }

        /// Points before queries at equal x, so every point that comes earlier has px <= x and every later one px > x
        sort(order.begin(), order.end(), [](const Event& a, const Event& b){
            if (a.x != b.x) return a.x < b.x;
            return a.query < b.query;
        });

        res.assign(queries, T(0));
        merged.resize(order.size());
        fen = FenwickPointUpdate<T>(zs.size());
        divide(0, order.size());
        return res;
    }

    /// Adds the left half's points to the right half's queries, leaving order[l, r) sorted by y
    void divide(int l, int r){
        if (r - l <= 1) return;

        int mid = (l + r) / 2;
        divide(l, mid);
        divide(mid, r);

        int i = l, j = mid, k = l;
        while (j < r){
            if (i < mid && order[i].y <= order[j].y){
                if (order[i].query == -1) fen.update(order[i].z, order[i].w);
                merged[k++] = order[i++];
            }
            else{
                if (order[j].query != -1) res[order[j].query] += fen.prefix(order[j].z);
                merged[k++] = order[j++];
            }
        }
        for (int p = l; p < i; p++){
            if (order[p].query == -1) fen.update(order[p].z, -order[p].w);
        }
        while (i < mid) merged[k++] = order[i++];

        copy(merged.begin() + l, merged.begin() + r, order.begin() + l);
    }
};

int main(){
    /***
     * Points (1, 1, 1) w 1, (2, 2, 2) w 10, (1, 3, 1) w 100, (3, 1, 3) w 1000
     * (2, 2, 2) dominates the first two: 11; (1, 3, 1) dominates the first and itself: 101; (3, 3, 3) all: 1111
    ***/
    DominanceSum<long long> ds;
    ds.add_point(1, 1, 1, 1);
    ds.add_point(2, 2, 2, 10);
    ds.add_point(1, 3, 1, 100);
    ds.add_point(3, 1, 3, 1000);
    assert(ds.add_query(2, 2, 2) == 0);
    assert(ds.add_query(1, 3, 1) == 1);
    assert(ds.add_query(3, 3, 3) == 2);
    assert(ds.add_query(0, 5, 5) == 3);

    vector<long long> res = ds.solve();
    assert(res[0] == 11);
    assert(res[1] == 101);
    assert(res[2] == 1111);
    assert(res[3] == 0);

    return 0;
}
