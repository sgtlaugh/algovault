/***
 *
 * Range LIS
 * Length of the longest increasing subsequence of any subarray, all queries answered in one batch
 *
 * Complexity:
 *   - RangeLIS: O(n log^2 n) to build, O((n + q) log n) for a batch of q queries, O(n) memory
 *   - UnitMongeProduct: O(n log n) per product of two permutations of size n
 *
 * RangeLIS lis(a): a of any type with <, strictly increasing subsequences
 * RangeLIS lis(a, true): non-decreasing subsequences instead
 * lis.query(queries): res[i] = the LIS length of a[l..r] for queries[i] = {l, r}, 0 when l > r
 * UnitMongeProduct product: the permutation product RangeLIS is built on, callable on its own
 * product(a, b): the permutation c with S(c)(i, k) = min over j of S(a)(i, j) + S(b)(j, k)
 *
 * Notes:
 *   - a, b and c are permutations of 0..n-1, and S(p)(i, j) = #{i' >= i : p[i'] < j} for i, j in [0, n]
 *   - n = 10^5 builds in about 0.5 s at -O2
 *   - Source: the seaweed braid of A. Tiskin, Semi-local string comparison (arXiv 0707.3619)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct UnitMongeProduct{
    vector<int> pool;
    size_t top = 0;

    vector<int> operator()(const vector<int>& a, const vector<int>& b){
        int n = a.size();
        assert((int)b.size() == n);
        if (pool.size() < 12 * (size_t)n + 256) pool.resize(12 * (size_t)n + 256);

        vector<int> c(n);
        multiply(a.data(), b.data(), c.data(), n);
        return c;
    }

    /// Splitting the middle index j in two halves, each row keeps its half's point except where the boundary steps
    void multiply(const int* a, const int* b, int* c, int n){
        if (n <= 1){
            if (n) c[0] = 0;
            return;
        }

        int mid = n / 2;
        size_t saved = top;
        int *rows = take(n), *sub_a = take(n), *cols = take(n), *rank = take(n), *sub_b = take(n), *sub_c = take(n);
        int lo = 0, hi = mid;
        for (int i = 0; i < n; i++){
            if (a[i] < mid) rows[lo] = i, sub_a[lo++] = a[i];
            else rows[hi] = i, sub_a[hi++] = a[i] - mid;
        }
        fill(rank, rank + n, 1);
        for (int j = 0; j < mid; j++) rank[b[j]] = 0;
        lo = 0, hi = mid;
        for (int k = 0; k < n; k++){
            if (rank[k] == 0) rank[k] = lo, cols[lo++] = k;
            else rank[k] = hi - mid, cols[hi++] = k;
        }
        for (int j = 0; j < n; j++) sub_b[j] = rank[b[j]];
        multiply(sub_a, sub_b, sub_c, mid);
        multiply(sub_a + mid, sub_b + mid, sub_c + mid, n - mid);

        int *row_of = rank, *in_hi = sub_a, *bnd = take(n + 1);
        for (int t = 0; t < n; t++){
            c[rows[t]] = cols[t < mid ? sub_c[t] : mid + sub_c[t]];
            in_hi[rows[t]] = t >= mid;
        }
        for (int i = 0; i < n; i++) row_of[c[i]] = i;

        /// d(i, k) = #hi points [row < i, col < k] - #lo points [row >= i, col >= k]; bnd[i] = first k with d >= 0
        bnd[n] = 0;
        int d = 0;
        for (int i = n - 1; i >= 0; i--){
            int k = bnd[i + 1];
            d -= in_hi[i] ? c[i] < k : c[i] >= k;
            while (d < 0){
                int r = row_of[k++];
                d += in_hi[r] ? r < i : r >= i;
            }
            bnd[i] = k;
        }
        for (int i = 0; i < n; i++){
            if (bnd[i] > bnd[i + 1]) c[i] = bnd[i] - 1;
        }

        top = saved;
    }

    /// Scratch from one buffer in stack order: a vector per call made the range LIS build 5x slower
    int* take(int len){
        assert(top + len <= pool.size());
        int* res = pool.data() + top;
        top += len;
        return res;
    }
};

struct RangeLIS{
    int n;
    vector<int> start_at;  /// start_at[e]: the top column whose seaweed leaves the bottom at column e, -1 for none
    UnitMongeProduct product;

    template <typename T>
    RangeLIS(const vector<T>& a, bool allow_equal = false) : n(a.size()), start_at(n, -1){
        /// Equal values never rise strictly, and always rise weakly: ranks tie-break by index accordingly
        vector<int> order(n), p(n);
        iota(order.begin(), order.end(), 0);
        sort(order.begin(), order.end(), [&](int i, int j){
            if (a[i] < a[j] || a[j] < a[i]) return a[i] < a[j];
            return allow_equal ? i < j : i > j;
        });
        for (int r = 0; r < n; r++) p[order[r]] = r;

        vector<int> out = braid(p);
        for (int s = 0; s < n; s++){
            if (out[n + s] < n) start_at[out[n + s]] = s;
        }
        product.pool = {};
    }

    /// LIS(a[l..r]) = (r - l + 1) - #{seaweeds that start at top column >= l and leave the bottom at column <= r}
    vector<int> query(const vector<pair<int, int>>& queries) const{
        int q = queries.size();
        vector<int> order(q), res(q, 0), tree(n + 1, 0);
        iota(order.begin(), order.end(), 0);
        sort(order.begin(), order.end(), [&](int i, int j){ return queries[i].second < queries[j].second; });

        int e = 0, added = 0;
        for (int i : order){
            auto [l, r] = queries[i];
            assert(0 <= l && l <= n && -1 <= r && r < n);
            for (; e <= r; e++){
                if (start_at[e] == -1) continue;
                for (int x = start_at[e] + 1; x <= n; x += x & -x) tree[x]++;
                added++;
            }

            int below = 0;
            for (int x = l; x > 0; x -= x & -x) below += tree[x];
            if (l <= r) res[i] = (r - l + 1) - (added - below);
        }

        return res;
    }

    /// p's braid on 2k strands, rows by value: in at left then top, out at bottom then right, each from bottom-left
    vector<int> braid(const vector<int>& p){
        int k = p.size();
        if (k <= 1) return k ? vector<int>{0, 1} : vector<int>{};

        int h = k / 2;
        vector<int> lo_cols, hi_cols, p_lo, p_hi;
        for (int j = 0; j < k; j++){
            if (p[j] < h) lo_cols.push_back(j), p_lo.push_back(p[j]);
            else hi_cols.push_back(j), p_hi.push_back(p[j] - h);
        }
        vector<int> lo = expand(braid(p_lo), lo_cols, k), hi = expand(braid(p_hi), hi_cols, k);

        /// The low rows lie above the high rows; each block passes the other block's left or right strands straight on
        vector<int> above(2 * k), below(2 * k);
        for (int i = 0; i < k - h; i++) above[i] = i;
        for (int i = 0; i < k + h; i++) above[k - h + i] = k - h + lo[i];
        for (int i = 0; i < 2 * k - h; i++) below[i] = hi[i];
        for (int i = 2 * k - h; i < 2 * k; i++) below[i] = i;
        return product(above, below);
    }

    /// A child's braid over the k parent columns: a column the child does not own passes straight down
    static vector<int> expand(const vector<int>& b, const vector<int>& cols, int k){
        int rows = cols.size();
        auto parent_out = [&](int o){ return o < rows ? cols[o] : k + o - rows; };

        vector<int> res(rows + k);
        for (int r = 0; r < rows; r++) res[r] = parent_out(b[r]);
        for (int c = 0; c < k; c++) res[rows + c] = c;
        for (int t = 0; t < rows; t++) res[rows + cols[t]] = parent_out(b[rows + t]);
        return res;
    }
};

int main(){
    /***
     * a = 3 1 4 1 5 9 2 6
     * a[0..7]: 1 4 5 9, or 1 4 5 6, has length 4; a[1..3] = 1 4 1 gives 1 4; a[5..7] = 9 2 6 gives 2 6
     * Non-decreasing on a[1..4] = 1 4 1 5: 1 4 5 or 1 1 5, still 3; on a[1..3] the two 1s give 1 1 or 1 4, so 2
    ***/
    vector<int> a = {3, 1, 4, 1, 5, 9, 2, 6};
    RangeLIS lis(a);
    assert((lis.query({{0, 7}, {1, 3}, {5, 7}, {4, 4}, {3, 2}}) == vector<int>{4, 2, 2, 1, 0}));

    RangeLIS weak(a, true);
    assert((weak.query({{1, 4}, {1, 3}}) == vector<int>{3, 2}));

    /***
     * a = b = 1 0: S(a) = S(b) = [[0, 1, 2], [0, 1, 1], [0, 0, 0]], and the min over j gives that matrix again
     * So a crossing times itself stays 1 0, where composing the permutations would give 0 1
    ***/
    UnitMongeProduct product;
    assert((product({1, 0}, {1, 0}) == vector<int>{1, 0}));
    assert((product({1, 0}, {0, 1}) == vector<int>{1, 0}));

    return 0;
}
