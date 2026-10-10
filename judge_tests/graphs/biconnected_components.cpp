// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/biconnected_components
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/articulation_points.cpp"
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

void write_int(string& out, long long x, char end){
    char digits[24];
    int len = 0;
    if (x < 0) out += '-', x = -x;
    do digits[len++] = '0' + x % 10, x /= 10; while (x);
    while (len) out += digits[--len];
    out += end;
}

int main(){
    int n = read_int(), m = read_int();
    Graph g(n);
    for (int i = 0; i < m; i++){
        int a = read_int(), b = read_int();
        g.add_edge(a, b);
    }

    auto blocks = g.get_blocks();

    /// A vertex outside every block (isolated, or touching only self loops) is a component of its own
    vector<vector<int>> comps;
    vector<int> last(n, -1);
    for (int b = 0; b < (int)blocks.size(); b++){
        comps.emplace_back();
        for (int id : blocks[b]){
            for (int u : {g.edges[id].first, g.edges[id].second}){
                if (last[u] == b) continue;
                last[u] = b;
                comps.back().push_back(u);
            }
        }
    }
    for (int u = 0; u < n; u++){
        if (last[u] == -1) comps.push_back({u});
    }

    string out;
    write_int(out, comps.size(), '\n');
    for (auto& comp : comps){
        write_int(out, comp.size(), comp.empty() ? '\n' : ' ');
        for (int i = 0; i < (int)comp.size(); i++) write_int(out, comp[i], i + 1 == (int)comp.size() ? '\n' : ' ');
    }
    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
