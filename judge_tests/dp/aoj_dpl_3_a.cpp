// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/DPL_3_A
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/dp/maximum_square.cpp"
#undef main

static char buf[1 << 23];
int buf_len, buf_pos;

int read_int(){
    while (buf_pos < buf_len && (buf[buf_pos] < '0' || buf[buf_pos] > '9')) buf_pos++;
    int x = 0;
    while (buf_pos < buf_len && buf[buf_pos] >= '0' && buf[buf_pos] <= '9') x = x * 10 + (buf[buf_pos++] - '0');
    return x;
}

int main(){
    buf_len = fread(buf, 1, sizeof(buf), stdin);

    int h = read_int(), w = read_int();
    vector<vector<int>> g(h, vector<int>(w));
    for (auto& row : g){
        for (auto& cell : row) cell = read_int() == 0;
    }

    int best = 0;
    for (const auto& row : square_sizes(g)){
        for (int side : row) best = max(best, side);
    }

    printf("%d\n", best * best);
    return 0;
}
