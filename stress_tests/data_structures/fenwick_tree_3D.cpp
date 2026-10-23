#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/fenwick_tree_3D.cpp"
#undef main

typedef vector<vector<vector<long long>>> Grid;

long long box_sum(const Grid& g, int x1, int y1, int z1, int x2, int y2, int z2){
    long long s = 0;
    for (int x = x1; x <= x2; x++) for (int y = y1; y <= y2; y++) for (int z = z1; z <= z2; z++) s += g[x][y][z];
    return s;
}

template <typename Call>
void check_rejects(Call call){
    stress::expect_assert("data_structures/fenwick_tree_3D.cpp", call);
}

/// All three variants against a plain 3D grid, boxes touching the borders included
int main(){
    /// Without the asserts a coordinate of -1 loops forever at i += i & -i = 0, and one past the edge silently drops the update
    for (int n : {1, 3}){
        int m = n + 1, r = n + 2;
        check_rejects([&]{ FenwickPointUpdate3D<int>(n, m, r).update(n, 0, 0, 1); });
        check_rejects([&]{ FenwickPointUpdate3D<int>(n, m, r).update(0, 0, r, 1); });
        check_rejects([&]{ FenwickPointUpdate3D<int>(n, m, r).prefix(0, -2, 0); });
        check_rejects([&]{ FenwickPointUpdate3D<int>(n, m, r).query(0, 0, 0, 0, m, 0); });
        check_rejects([&]{ FenwickRangeUpdate3D<int>(n, m, r).update(0, 0, 0, 0, 0, r, 1); });
        check_rejects([&]{ FenwickRangeUpdate3D<int>(n, m, r).update(-1, 0, 0, 0, 0, 0, 1); });
        check_rejects([&]{ FenwickRangeUpdate3D<int>(n, m, r).query(0, m, 0); });
        check_rejects([&]{ FenwickFull3D<long long>(n, m, r).update(0, -1, 0, 0, 0, 0, 1); });
        check_rejects([&]{ FenwickFull3D<long long>(n, m, r).update(0, 0, 0, n, 0, 0, 1); });
        check_rejects([&]{ FenwickFull3D<long long>(n, m, r).prefix(0, 0, r); });
        check_rejects([&]{ FenwickFull3D<long long>(n, m, r).query(0, 0, 0, 0, 0, r); });
    }

    for (long long it = 0; it < stress::scaled(1000); it++){
        int n = stress::rand_int(1, it % 10 ? 5 : 12), m = stress::rand_int(1, it % 10 ? 5 : 12), r = stress::rand_int(1, it % 10 ? 5 : 12);
        FenwickPointUpdate3D<long long> point(n, m, r);
        FenwickRangeUpdate3D<long long> range(n, m, r);
        FenwickFull3D<long long> full(n, m, r);
        Grid a(n, vector<vector<long long>>(m, vector<long long>(r))), b = a;

        for (int op = 0; op < 60; op++){
            int c[2][3], lim[3] = {n, m, r};
            for (int d = 0; d < 3; d++){
                c[0][d] = stress::rand_int(0, lim[d] - 1), c[1][d] = stress::rand_int(0, lim[d] - 1);
                if (op % 5) tie(c[0][d], c[1][d]) = make_pair(min(c[0][d], c[1][d]), max(c[0][d], c[1][d]));  /// otherwise often empty, which must be a no-op
            }
            auto [x1, y1, z1] = make_tuple(c[0][0], c[0][1], c[0][2]);
            auto [x2, y2, z2] = make_tuple(c[1][0], c[1][1], c[1][2]);
            long long v = stress::rand_int(-1000000, 1000000);

            if (stress::rand_int(0, 1)){
                point.update(x1, y1, z1, v);
                a[x1][y1][z1] += v;
                range.update(x1, y1, z1, x2, y2, z2, v);
                full.update(x1, y1, z1, x2, y2, z2, v);
                for (int x = x1; x <= x2; x++) for (int y = y1; y <= y2; y++) for (int z = z1; z <= z2; z++) b[x][y][z] += v;
            }
            else{
                assert(point.query(x1, y1, z1, x2, y2, z2) == box_sum(a, x1, y1, z1, x2, y2, z2));
                assert(range.query(x1, y1, z1) == b[x1][y1][z1]);
                assert(full.query(x1, y1, z1, x2, y2, z2) == box_sum(b, x1, y1, z1, x2, y2, z2));
                assert(point.prefix(x1 - 1, y2, z2) == box_sum(a, 0, 0, 0, x1 - 1, y2, z2));
                assert(full.prefix(x2, y1 - 1, z1) == box_sum(b, 0, 0, 0, x2, y1 - 1, z1));
            }
        }
    }

    return 0;
}
