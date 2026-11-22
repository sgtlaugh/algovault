// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/tree_diameter
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/trees/leaf_peeling.cpp"
#undef main

int read_int(){
    int c = getchar();
    while (c < '0' || c > '9') c = getchar();

    int x = 0;
    for (; c >= '0' && c <= '9'; c = getchar()) x = x * 10 + c - '0';
    return x;
}

int main(){
    int n = read_int();
    vector<pair<int, int>> edges(n - 1);
    vector<int> weight(n - 1);
    for (int i = 0; i < n - 1; i++){
        edges[i].first = read_int(), edges[i].second = read_int();
        weight[i] = read_int();
    }

    PeeledTree t = peel_leaves(n, edges);
    vector<long long> up(n, 0);
    for (int i = 0; i < n - 1; i++){
        auto [u, v] = edges[i];
        up[t.parent[u] == v ? u : v] = weight[i];
    }

    vector<long long> best(n, 0), second(n, 0);
    vector<int> best_child(n, -1), second_child(n, -1);
    for (int v : t.order){
        int p = t.parent[v];
        if (p == -1) continue;

        long long len = best[v] + up[v];
        if (best_child[p] == -1 || len > best[p]){
            second[p] = best[p], second_child[p] = best_child[p];
            best[p] = len, best_child[p] = v;
        }
        else if (second_child[p] == -1 || len > second[p]) second[p] = len, second_child[p] = v;
    }

    int top = 0;
    for (int v = 0; v < n; v++){
        if (best[v] + second[v] > best[top] + second[top]) top = v;
    }

    vector<int> path;
    for (int v = best_child[top]; v != -1; v = best_child[v]) path.push_back(v);
    reverse(path.begin(), path.end());
    path.push_back(top);
    for (int v = second_child[top]; v != -1; v = best_child[v]) path.push_back(v);

    string out = to_string(best[top] + second[top]) + ' ' + to_string(path.size()) + '\n';
    for (int i = 0; i < (int)path.size(); i++) out += to_string(path[i]) + (i + 1 < (int)path.size() ? ' ' : '\n');
    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
