// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/1068
// competitive-verifier: TLE 6
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/segment_tree_2D.cpp"
#undef main

int read_int(){
    int c = getchar();
    while (c != EOF && (c < '0' || c > '9')) c = getchar();

    int x = 0;
    for (; c >= '0' && c <= '9'; c = getchar()) x = x * 10 + c - '0';
    return x;
}

struct Min{
    int operator()(int a, int b) const{
        return min(a, b);
    }
};

int main(){
    string out;
    while (true){
        int r = read_int(), c = read_int(), q = read_int();
        if (r == 0 && c == 0 && q == 0) break;

        vector<vector<int>> grid(r, vector<int>(c));
        for (auto& row : grid){
            for (int& x : row) x = read_int();
        }

        SegmentTree2D<int, Min> st(grid, INT_MAX);
        while (q--){
            int r1 = read_int(), c1 = read_int(), r2 = read_int(), c2 = read_int();
            out += to_string(st.query(r1, c1, r2, c2)) + '\n';
        }
    }

    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
