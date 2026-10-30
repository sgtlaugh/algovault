/***
 *
 * 2D Segment Tree
 * Point assignment and rectangle min, max, sum or gcd on a grid, for any commutative combine with an identity
 *
 * Complexity:
 *   - build: O(n m)
 *   - update, rectangle query: O(log n log m)
 *   - cell query: O(1)
 *
 * SegmentTree2D<T, Combine = plus<T>> st(n, m, identity = T(), combine = Combine()): n x m grid of identity
 * SegmentTree2D<T, Combine> st(grid, identity = T(), combine = Combine()): built from a rectangular grid
 * st.update(x, y, v): a[x][y] = v, 0 <= x < n, 0 <= y < m
 * st.query(x, y): a[x][y], 0 <= x < n, 0 <= y < m
 * st.query(x1, y1, x2, y2): combine over rows x1..x2 and columns y1..y2, 0 <= x1 <= x2 < n, 0 <= y1 <= y2 < m
 *
 * Notes:
 *   - combine must be associative and commutative, as a rectangle's cells have no order, and combine(identity, a) = a
 *   - Point add on sums is st.update(x, y, st.query(x, y) + v)
 *   - Memory is 4 n m values, use fenwick_tree_2D.cpp for plain sums in a quarter of it
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T, typename Combine = plus<T>>
struct SegmentTree2D{
    int n, m;
    T identity;
    Combine combine;
    /// at(r, c) combines the rows of row node r with the columns of column node c, cell (x, y) is at(n + x, m + y)
    vector<T> tree;

    SegmentTree2D(int n, int m, T identity = T(), Combine combine = Combine())
        : n(n), m(m), identity(identity), combine(combine), tree(4LL * n * m, identity) {}

    SegmentTree2D(const vector<vector<T>>& grid, T identity = T(), Combine combine = Combine())
        : SegmentTree2D(grid.size(), grid.empty() ? 0 : grid[0].size(), identity, combine){
        for (int i = 0; i < n; i++){
            for (int j = 0; j < m; j++) at(n + i, m + j) = grid[i][j];
            for (int c = m - 1; c > 0; c--) at(n + i, c) = combine(at(n + i, 2 * c), at(n + i, 2 * c + 1));
        }

        for (int r = n - 1; r > 0; r--){
            for (int c = 1; c < 2 * m; c++) at(r, c) = combine(at(2 * r, c), at(2 * r + 1, c));
        }
    }

    void update(int x, int y, T v){
        assert(0 <= x && x < n && 0 <= y && y < m);
        int r = n + x;

        at(r, m + y) = v;
        for (int c = (m + y) >> 1; c > 0; c >>= 1) at(r, c) = combine(at(r, 2 * c), at(r, 2 * c + 1));

        for (r >>= 1; r > 0; r >>= 1){
            for (int c = m + y; c > 0; c >>= 1) at(r, c) = combine(at(2 * r, c), at(2 * r + 1, c));
        }
    }

    T query(int x, int y){
        assert(0 <= x && x < n && 0 <= y && y < m);
        return at(n + x, m + y);
    }

    T query(int x1, int y1, int x2, int y2){
        assert(0 <= x1 && x1 <= x2 && x2 < n && 0 <= y1 && y1 <= y2 && y2 < m);
        T res = identity;
        int lo = n + x1, hi = n + x2 + 1;

        while (lo < hi){
            if (lo & 1) res = combine(res, query_row(lo++, y1, y2));
            if (hi & 1) res = combine(res, query_row(--hi, y1, y2));
            lo >>= 1, hi >>= 1;
        }

        return res;
    }

    T& at(int r, int c){
        return tree[(size_t)r * 2 * m + c];
    }

    T query_row(int r, int y1, int y2){
        T res = identity;
        int lo = m + y1, hi = m + y2 + 1;

        while (lo < hi){
            if (lo & 1) res = combine(res, at(r, lo++));
            if (hi & 1) res = combine(res, at(r, --hi));
            lo >>= 1, hi >>= 1;
        }

        return res;
    }
};

int main(){
    /// Rectangle sums on a 3 x 4 grid, then a point assignment
    vector<vector<long long>> grid = {{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}};
    SegmentTree2D<long long> sum(grid);
    assert(sum.query(0, 0, 2, 3) == 78);
    assert(sum.query(1, 1, 2, 2) == 34);                /// 6 + 7 + 10 + 11

    sum.update(1, 2, 100);
    assert(sum.query(1, 1, 2, 2) == 127);               /// the 7 became 100
    assert(sum.query(1, 2) == 100);

    /// Rectangle minimum with a lambda over a grid that starts at the identity, the largest int
    SegmentTree2D mn(3, 4, INT_MAX, [](int a, int b){ return min(a, b); });
    mn.update(0, 3, 7);
    mn.update(2, 0, 4);
    assert(mn.query(0, 0, 2, 3) == 4);
    assert(mn.query(0, 1, 1, 3) == 7);
    assert(mn.query(1, 1, 1, 1) == INT_MAX);

    return 0;
}
