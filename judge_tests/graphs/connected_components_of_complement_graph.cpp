// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/connected_components_of_complement_graph
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/complement_graph_bfs.cpp"
#undef main

static char buf[1 << 25];
int buf_len = 0, buf_pos = 0;

inline int read_int(){
    while (buf[buf_pos] < '0') buf_pos++;
    int x = 0;
    while (buf[buf_pos] >= '0') x = x * 10 + buf[buf_pos++] - '0';
    return x;
}

int main(){
    buf_len = fread(buf, 1, sizeof(buf) - 1, stdin);
    buf[buf_len] = 0;
    int n = read_int(), m = read_int();

    ComplementGraph g(n);
    for (int i = 0; i < m; i++){
        int u = read_int(), v = read_int();
        g.add_edge(u, v);
    }

    auto comp = g.components();
    int k = *max_element(comp.begin(), comp.end()) + 1;
    vector<vector<int>> groups(k);
    for (int v = 0; v < n; v++) groups[comp[v]].push_back(v);

    string out = to_string(k) + "\n";
    for (auto& group : groups){
        out += to_string(group.size());
        for (int v : group) out += " " + to_string(v);
        out += "\n";
    }
    fputs(out.c_str(), stdout);
    return 0;
}
