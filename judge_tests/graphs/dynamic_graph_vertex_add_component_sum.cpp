// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/dynamic_graph_vertex_add_component_sum
// competitive-verifier: TLE 10
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/online_dynamic_connectivity.cpp"
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
    int n = read_int(), q = read_int();
    vector<long long> a(n);
    for (auto& x: a) x = read_int();

    OnlineDynamicConnectivity g(a);
    string out;
    for (int i = 0; i < q; i++){
        int type = read_int(), u = read_int();
        if (type == 0) g.add_edge(u, read_int());
        else if (type == 1) g.remove_edge(u, read_int());
        else if (type == 2) g.add(u, read_int());
        else out += to_string(g.sum(u)) + '\n';
    }

    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
