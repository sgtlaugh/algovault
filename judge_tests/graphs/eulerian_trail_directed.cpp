// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/eulerian_trail_directed
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/euler_path.cpp"
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

void write_list(string& out, const vector<int>& values){
    for (int i = 0; i < (int)values.size(); i++){
        if (i) out += ' ';
        out += to_string(values[i]);
    }
    out += '\n';
}

int main(){
    int t = read_int();
    string out;
    while (t--){
        int n = read_int(), m = read_int();
        EulerDirected g(n);
        for (int i = 0; i < m; i++){
            int u = read_int(), v = read_int();
            g.add_edge(u, v);
        }

        auto path = g.solve(-1);
        if (path.empty()){
            out += "No\n";
            continue;
        }

        out += "Yes\n";
        write_list(out, path);
        write_list(out, g.walk_edges);
    }
    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
