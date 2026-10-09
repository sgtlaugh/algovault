#include "../common.h"

#define main library_main
#include "../../code_library/graphs/hopcroft_karp.cpp"
#undef main

/// The original Library implementation, unchanged, as a second reference
#define clr(ar) memset(ar, 0, sizeof(ar))
namespace hc{
    const int MAX = 100010;
    bool visited[MAX];
    vector <int> adj[MAX];
    int n, L[MAX], R[MAX], Q[MAX], len[MAX], dis[MAX], parent[MAX];
    inline void init(int nodes){
        n = nodes, clr(len);
        for (int i = 0; i < MAX; i++) adj[i].clear();
    }
    inline void add_edge(int u, int v){
        len[u]++;
        adj[u].push_back(v);
    }
    bool dfs(int i){
        for (int j = 0; j < len[i]; j++){
            int x = adj[i][j];
            if (L[x] == -1 || (parent[L[x]] == i)){
                if (L[x] == -1 || dfs(L[x])){
                    L[x] = i, R[i] = x;
                    return true;
                }
            }
        }
        return false;
    }
    bool bfs(){
        clr(visited);
        int i, j, x, d, f = 0, l = 0;
        for (i = 0; i < n; i++){
            if (R[i] == -1){
                visited[i] = true;
                Q[l++] = i, dis[i] = 0;
            }
        }
        while (f < l){
            i = Q[f++];
            for (j = 0; j < len[i]; j++){
                x = adj[i][j], d = L[x];
                if (d == -1) return true;
                else if (!visited[d]){
                    Q[l++] = d;
                    parent[d] = i, visited[d] = true, dis[d] = dis[i] + 1;
                }
            }
        }
        return false;
    }
    int hopcroft_karp(){
        int res = 0;
        memset(L, -1, sizeof(L));
        memset(R, -1, sizeof(R));
        while (bfs()){
            for (int i = 0; i < n; i++){
                if (R[i] == -1 && dfs(i)) res++;
            }
        }
        return res;
    }
}

int kuhn(int n_left, int n_right, const vector<vector<int>>& adj){
    vector<int> match(n_right, -1);
    vector<int> seen;
    function<bool(int)> try_kuhn = [&](int u){
        for (int v : adj[u]){
            if (seen[v]) continue;
            seen[v] = 1;
            if (match[v] == -1 || try_kuhn(match[v])){
                match[v] = u;
                return true;
            }
        }
        return false;
    };
    int res = 0;
    for (int u = 0; u < n_left; u++){
        seen.assign(n_right, 0);
        res += try_kuhn(u);
    }
    return res;
}

/// Exact maximum matching by DP over subsets of used right vertices
int brute(int n_left, int n_right, const vector<vector<int>>& adj){
    vector<int> best(1 << n_right, -1);
    best[0] = 0;
    for (int u = 0; u < n_left; u++){
        vector<int> next = best;
        for (int mask = 0; mask < (1 << n_right); mask++){
            if (best[mask] < 0) continue;
            for (int v : adj[u]){
                if (!(mask >> v & 1)) next[mask | 1 << v] = max(next[mask | 1 << v], best[mask] + 1);
            }
        }
        best = next;
    }
    return *max_element(best.begin(), best.end());
}

int run(int n_left, int n_right, const vector<vector<int>>& adj){
    HopcroftKarp hk(n_left, n_right);
    for (int u = 0; u < n_left; u++){
        for (int v : adj[u]) hk.add_edge(u, v);
    }
    int res = hk.max_matching();

    int pairs = 0;
    for (int u = 0; u < n_left; u++){
        int v = hk.match_left[u];
        if (v == -1) continue;
        pairs++;
        assert(hk.match_right[v] == u);
        assert(find(adj[u].begin(), adj[u].end(), v) != adj[u].end());
    }
    for (int v = 0; v < n_right; v++){
        if (hk.match_right[v] != -1) assert(hk.match_left[hk.match_right[v]] == v);
    }
    assert(pairs == res);
    return res;
}

vector<vector<int>> random_graph(int n_left, int n_right, int m){
    vector<vector<int>> adj(n_left);
    for (int i = 0; i < m && n_left && n_right; i++) adj[stress::rand_int(0, n_left - 1)].push_back(stress::rand_int(0, n_right - 1));
    return adj;
}

int main(){
    for (long long it = 0; it < stress::scaled(4000); it++){
        int n_left = stress::rand_int(0, 7), n_right = stress::rand_int(0, 7);
        auto adj = random_graph(n_left, n_right, stress::rand_int(0, 20));
        assert(run(n_left, n_right, adj) == brute(n_left, n_right, adj));
    }

    for (long long it = 0; it < stress::scaled(400); it++){
        int n_left = stress::rand_int(1, 300), n_right = stress::rand_int(1, 300);
        auto adj = random_graph(n_left, n_right, stress::rand_int(0, 1500));
        assert(run(n_left, n_right, adj) == kuhn(n_left, n_right, adj));
    }

    for (long long it = 0; it < stress::scaled(6); it++){
        int n_left = stress::rand_int(20000, 50000), n_right = stress::rand_int(20000, 50000);
        auto adj = random_graph(n_left, n_right, stress::rand_int(n_left, 3 * n_left));
        hc::init(max(n_left, n_right));
        for (int u = 0; u < n_left; u++){
            for (int v : adj[u]) hc::add_edge(u, v);
        }
        assert(run(n_left, n_right, adj) == hc::hopcroft_karp());
    }

    /// Staircase: left i -> right i and i + 1, with the greedy-friendly order reversed, forces long augmenting paths
    int n = 100000;
    vector<vector<int>> stairs(n);
    for (int i = 0; i < n; i++){
        if (i + 1 < n) stairs[i].push_back(i + 1);
        stairs[i].push_back(i);
    }
    assert(run(n, n, stairs) == n);

    return 0;
}
