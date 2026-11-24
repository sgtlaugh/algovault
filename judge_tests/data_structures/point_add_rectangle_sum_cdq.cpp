// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/point_add_rectangle_sum
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/cdq_divide_conquer.cpp"
#undef main

int read_int(){
    int c = getchar();
    while (c < '0' || c > '9') c = getchar();

    int x = 0;
    for (; c >= '0' && c <= '9'; c = getchar()) x = x * 10 + (c - '0');
    return x;
}

int main(){
    int n = read_int(), q = read_int();

    /// Time is the first coordinate, so a sum sees exactly the points added before it
    DominanceSum<long long> ds;
    for (int i = 0; i < n; i++){
        int x = read_int(), y = read_int(), w = read_int();
        ds.add_point(0, x, y, w);
    }

    /// x in [l, r) and y in [d, u) by inclusion-exclusion over the 4 corners
    vector<int> first;
    for (int t = 1; t <= q; t++){
        if (read_int() == 0){
            int x = read_int(), y = read_int(), w = read_int();
            ds.add_point(t, x, y, w);
            continue;
        }

        int l = read_int(), d = read_int(), r = read_int(), u = read_int();
        first.push_back(ds.add_query(t, r - 1, u - 1));
        ds.add_query(t, l - 1, u - 1);
        ds.add_query(t, r - 1, d - 1);
        ds.add_query(t, l - 1, d - 1);
    }

    vector<long long> res = ds.solve();
    string out;
    for (int i : first) out += to_string(res[i] - res[i + 1] - res[i + 2] + res[i + 3]) + '\n';
    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
