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

/// All three variants against a plain 3D grid, boxes touching the borders included
int main(){
    for (long long it = 0; it < stress::scaled(1000); it++){
        int n = stress::rand_int(1, it % 10 ? 5 : 12), m = stress::rand_int(1, it % 10 ? 5 : 12), r = stress::rand_int(1, it % 10 ? 5 : 12);
        FenwickPointUpdate3D<long long> point(n, m, r);
        FenwickRangeUpdate3D<long long> range(n, m, r);
        FenwickFull3D<long long> full(n, m, r);
        Grid a(n + 1, vector<vector<long long>>(m + 1, vector<long long>(r + 1))), b = a;

        for (int op = 0; op < 60; op++){
            int c[2][3], lim[3] = {n, m, r};
            for (int d = 0; d < 3; d++){
                c[0][d] = stress::rand_int(1, lim[d]), c[1][d] = stress::rand_int(1, lim[d]);
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
            }
        }
    }

    return 0;
}
