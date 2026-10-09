#include "../common.h"
#include <sys/wait.h>
#include <unistd.h>

#define main library_main
#include "../../code_library/trees/rerooting.cpp"
#undef main

const unsigned long long MOD = 1000000007ULL;

using Info = pair<long long, long long>;

/// Random tree of a given shape with random labels and weights
vector<array<long long, 3>> make_tree(int n, int shape, long long max_w){
    vector<int> label(n);
    iota(label.begin(), label.end(), 0);
    shuffle(label.begin(), label.end(), stress::rng());

    vector<array<long long, 3>> edges;
    for (int i = 1; i < n; i++){
        int p;
        if (shape == 0) p = stress::rand_int(0, i - 1);
        else if (shape == 1) p = i - 1;
        else if (shape == 2) p = 0;
        else if (shape == 3) p = i % 2 ? i - 1 : max(0, i - 2);
        else p = i < n / 2 ? i - 1 : stress::rand_int(max(0, n / 2 - 3), i - 1);
        edges.push_back({label[i], label[p], stress::rand_int(0, max_w)});
    }
    shuffle(edges.begin(), edges.end(), stress::rng());
    return edges;
}

/// The plain tree DP rooted at root: children before parents in reverse BFS order
template<typename T, typename Merge, typename Apply>
T naive_dp(const vector<vector<pair<int, long long>>>& adj, int root, Merge merge, T identity, Apply apply_edge){
    int n = adj.size();
    vector<int> parent(n, -1), order = {root};
    vector<bool> seen(n, false);
    seen[root] = true;
    for (int i = 0; i < (int)order.size(); i++){
        for (auto [v, w] : adj[order[i]]){
            if (!seen[v]) seen[v] = true, parent[v] = order[i], order.push_back(v);
        }
    }

    vector<T> dp(n, identity);
    for (int i = n - 1; i >= 0; i--){
        int u = order[i];
        for (auto [v, w] : adj[u]){
            if (v != parent[u]) dp[u] = merge(dp[u], apply_edge(dp[v], v, u, w));
        }
    }
    return dp[root];
}

template<typename T, typename E = long long, typename Merge, typename Apply>
void compare(int n, const vector<array<long long, 3>>& edges, Merge merge, T identity, Apply apply_edge){
    Rerooting<T, E> tree(n);
    vector<vector<pair<int, long long>>> adj(n);
    for (auto [u, v, w] : edges){
        adj[u].push_back({(int)v, w}), adj[v].push_back({(int)u, w});
        tree.add_edge(u, v, (E)w);
    }
    vector<T> res = tree.solve(merge, identity, apply_edge);

    assert((int)res.size() == n);
    for (int r = 0; r < n; r++) assert(res[r] == naive_dp(adj, r, merge, identity, apply_edge));
}

/// Labelled tree from a Pruefer sequence of length n - 2, by the O(n^2) smallest-leaf definition
vector<array<long long, 3>> pruefer_tree(int n, const vector<int>& code, long long max_w){
    vector<int> degree(n, 1);
    for (int x : code) degree[x]++;

    vector<array<long long, 3>> edges;
    for (int x : code){
        int leaf = find(degree.begin(), degree.end(), 1) - degree.begin();
        edges.push_back({leaf, x, stress::rand_int(0, max_w)});
        degree[leaf]--, degree[x]--;
    }

    int a = find(degree.begin(), degree.end(), 1) - degree.begin();
    int b = find(degree.begin() + a + 1, degree.end(), 1) - degree.begin();
    if (n >= 2) edges.push_back({a, b, stress::rand_int(0, max_w)});
    return edges;
}

/// Four DPs against the DP run from every root: distance sums, eccentricities (long long and int weights), a node-valued hash
void check_edges(int n, const vector<array<long long, 3>>& edges){
    vector<unsigned long long> value(n);
    for (auto& x : value) x = stress::rand_int(0, MOD - 1);

    auto add = [](const Info& a, const Info& b){ return Info{a.first + b.first, a.second + b.second}; };
    auto lift = [](const Info& x, int, int, long long w){ return Info{x.first + 1, x.second + (x.first + 1) * w}; };
    compare(n, edges, add, Info{0, 0}, lift);

    auto longest = [](long long a, long long b){ return max(a, b); };
    auto extend = [](long long x, int, int, long long w){ return x + w; };
    compare(n, edges, longest, 0LL, extend);

    auto extend_int = [](long long x, int, int, int w){ return x + w; };
    compare<long long, int>(n, edges, longest, 0LL, extend_int);

    auto product = [](long long a, long long b){ return (long long)((unsigned long long)a * b % MOD); };
    auto hash_edge = [&](long long x, int from, int to, long long w){
        return (long long)(((unsigned long long)x * 31 + value[from] * 7 + value[to] + (unsigned long long)w) % MOD);
    };
    compare(n, edges, product, 1LL, hash_edge);
}

void check(int n, int shape, long long max_w){
    check_edges(n, make_tree(n, shape, max_w));
}

/// Every labelled tree on n nodes, one per Pruefer sequence
void check_all_trees(int n, long long max_w){
    vector<int> code(max(0, n - 2), 0);
    while (true){
        check_edges(n, pruefer_tree(n, code, max_w));

        int i = 0;
        while (i < (int)code.size() && code[i] == n - 1) code[i++] = 0;
        if (i == (int)code.size()) break;
        code[i]++;
    }
}

/// Eccentricity against the farthest of the two diameter endpoints, for trees too large to root n times
void check_large(int n, int shape){
    auto edges = make_tree(n, shape, 1000000000LL);
    vector<vector<pair<int, long long>>> adj(n);
    Rerooting<long long> tree(n);
    for (auto [u, v, w] : edges){
        adj[u].push_back({(int)v, w}), adj[v].push_back({(int)u, w});
        tree.add_edge(u, v, w);
    }

    auto distances = [&](int src){
        vector<long long> dist(n, -1);
        vector<int> order = {src};
        dist[src] = 0;
        for (int i = 0; i < (int)order.size(); i++){
            for (auto [v, w] : adj[order[i]]){
                if (dist[v] == -1) dist[v] = dist[order[i]] + w, order.push_back(v);
            }
        }
        return dist;
    };
    auto da = distances(0);
    int a = max_element(da.begin(), da.end()) - da.begin();
    da = distances(a);
    int b = max_element(da.begin(), da.end()) - da.begin();
    auto db = distances(b);

    auto longest = [](long long x, long long y){ return max(x, y); };
    auto extend = [](long long x, int, int, long long w){ return x + w; };
    vector<long long> res = tree.solve(longest, 0LL, extend);
    for (int v = 0; v < n; v++) assert(res[v] == max(da[v], db[v]));
}

/// Unit path 0 - 1 - ... - n-1: node i has distance sum i(i+1)/2 + (n-1-i)(n-i)/2
void check_path(int n){
    Rerooting<Info> tree(n);
    for (int i = 1; i < n; i++) tree.add_edge(i - 1, i, 1);

    auto add = [](const Info& a, const Info& b){ return Info{a.first + b.first, a.second + b.second}; };
    auto lift = [](const Info& x, int, int, long long w){ return Info{x.first + 1, x.second + (x.first + 1) * w}; };
    vector<Info> res = tree.solve(add, Info{0, 0}, lift);
    for (long long i = 0; i < n; i++) assert(res[i].first == n - 1 && res[i].second == i * (i + 1) / 2 + (n - 1 - i) * (n - i) / 2);
}

/// solve must abort on a non-tree, the alarm turns an endless traversal into a failure instead of a hang
void check_rejects(int n, const vector<pair<int, int>>& edges){
    pid_t pid = fork();
    assert(pid >= 0);
    if (pid == 0){
        alarm(5);
        assert(freopen("/dev/null", "w", stderr));
        Rerooting<long long> tree(n);
        for (auto [u, v] : edges) tree.add_edge(u, v);
        tree.solve([](long long a, long long b){ return a + b; }, 0LL, [](long long x, int, int, long long){ return x + 1; });
        _exit(0);
    }

    int status;
    assert(waitpid(pid, &status, 0) == pid);
    assert(WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT);
}

int main(){
    check_rejects(4, {{0, 1}, {1, 2}, {2, 0}});          /// cycle plus an isolated node, n - 1 edges
    check_rejects(5, {{0, 1}, {1, 2}, {2, 3}, {3, 1}});  /// cycle below the root plus an isolated node
    check_rejects(3, {{0, 1}, {1, 2}, {2, 0}});          /// cycle, n edges
    check_rejects(3, {{0, 1}, {0, 1}, {1, 2}});          /// multi-edge
    check_rejects(3, {{0, 1}});                          /// disconnected
    check_rejects(2, {{0, 0}, {0, 1}});                  /// self-loop at the root, only the edge count sees it
    check_rejects(0, {});

    for (int n = 1; n <= 6; n++) check_all_trees(n, n % 2 ? 1000000000LL : 3);  /// 1296 trees at n = 6, n = 7 adds 16807 for little new structure

    for (long long it = 0; it < stress::scaled(2000); it++){
        int n = it < 200 ? it % 40 + 1 : stress::rand_int(1, 80);
        check(n, it % 5, it % 3 ? 1000000000LL : 3);
    }

    for (long long it = 0; it < stress::scaled(5); it++) check(stress::rand_int(300, 600), it % 5, 1000000000LL);

    for (int shape = 0; shape < 5; shape++) check_large(200000, shape);  /// the path shape is deep enough to overflow an 8 MB stack with recursion

    for (int n : {1, 2, 3, 200000}) check_path(n);

    return 0;
}
