// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/two_edge_connected_components
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/online_bridges.cpp"
#undef main

static char in_buf[1 << 16];
static size_t in_len = 0, in_pos = 0;

int read_char(){
    if (in_pos == in_len){
        in_len = fread(in_buf, 1, sizeof(in_buf), stdin), in_pos = 0;
        if (in_len == 0) return -1;
    }
    return in_buf[in_pos++];
}

int read_int(){
    int c = read_char();
    while (c != -1 && (c < '0' || c > '9')) c = read_char();

    int x = 0;
    while (c >= '0' && c <= '9') x = x * 10 + (c - '0'), c = read_char();
    return x;
}

int main(){
    int n = read_int(), m = read_int();
    OnlineBridges g(n);
    for (int i = 0; i < m; i++){
        int a = read_int(), b = read_int();
        g.add_edge(a, b);
    }

    vector<vector<int>> groups(n);
    for (int v = 0; v < n; v++) groups[g.component(v)].push_back(v);

    int count = 0;
    string body;
    for (auto& group : groups){
        if (group.empty()) continue;
        count++;
        body += to_string(group.size());
        for (int v : group) body += ' ' + to_string(v);
        body += '\n';
    }

    string out = to_string(count) + '\n' + body;
    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
