#include "../common.h"
#include <sys/wait.h>
#include <unistd.h>

#define main library_main
#include "../../code_library/trees/prufer_code.cpp"
#undef main

/// Textbook O(n^2) encoding: repeatedly scan for the smallest leaf, record its neighbor, delete it
vector<int> naive_encode(int n, const vector<pair<int, int>>& edges){
    vector<set<int>> adj(n);
    for (auto [u, v] : edges) adj[u].insert(v), adj[v].insert(u);

    vector<int> code;
    for (int step = 0; step < n - 2; step++){
        int leaf = 0;
        while (adj[leaf].size() != 1) leaf++;
        int next = *adj[leaf].begin();
        code.push_back(next);
        adj[next].erase(leaf), adj[leaf].clear();
    }
    return code;
}

/// Asserts the edges form a tree and every edge is {child, parent} when rooted at n - 1
void check_rooted_tree(int n, const vector<pair<int, int>>& edges){
    assert((int)edges.size() == n - 1);
    vector<vector<int>> adj(n);
    for (auto [u, v] : edges){
        assert(0 <= u && u < n && 0 <= v && v < n);
        adj[u].push_back(v), adj[v].push_back(u);
    }

    vector<int> parent(n, -2), queue = {n - 1};
    parent[n - 1] = -1;
    for (int i = 0; i < (int)queue.size(); i++){
        for (int v : adj[queue[i]]){
            if (parent[v] != -2) continue;
            parent[v] = queue[i];
            queue.push_back(v);
        }
    }
    assert((int)queue.size() == n);

    for (auto [child, par] : edges) assert(parent[child] == par);
}

vector<pair<int, int>> normalized(vector<pair<int, int>> edges){
    for (auto& [u, v] : edges){
        if (u > v) swap(u, v);
    }
    sort(edges.begin(), edges.end());
    return edges;
}

vector<int> encode(int n, const vector<pair<int, int>>& edges){
    PruferCode tree(n);
    for (auto [u, v] : edges) tree.add_edge(u, v);
    return tree.encode();
}

/// Random labeled tree with a chosen shape: random attachment, path, star, caterpillar, broom
vector<pair<int, int>> random_tree(int n){
    vector<int> label(n);
    iota(label.begin(), label.end(), 0);
    shuffle(label.begin(), label.end(), stress::rng());

    int shape = stress::rand_int(0, 4), spine = max(1, n / 2);
    vector<pair<int, int>> edges;
    for (int i = 1; i < n; i++){
        int p;
        if (shape == 0) p = stress::rand_int(0, i - 1);
        else if (shape == 1) p = i - 1;
        else if (shape == 2) p = 0;
        else if (shape == 3) p = i < spine ? i - 1 : stress::rand_int(0, spine - 1);
        else p = i < spine ? i - 1 : spine - 1;
        if (stress::rand_int(0, 1)) edges.push_back({label[i], label[p]});
        else edges.push_back({label[p], label[i]});
    }
    shuffle(edges.begin(), edges.end(), stress::rng());
    return edges;
}

/// Every code of length n - 2 decodes to a valid tree whose textbook code is that code, all trees distinct
void check_exhaustive(int n){
    vector<int> code(n - 2, 0);
    set<vector<pair<int, int>>> trees;
    while (true){
        auto edges = prufer_decode(code);
        check_rooted_tree(n, edges);
        assert(naive_encode(n, edges) == code);
        assert(encode(n, edges) == code);
        trees.insert(normalized(edges));

        int i = 0;
        while (i < n - 2 && code[i] == n - 1) code[i++] = 0;
        if (i == n - 2) break;
        code[i]++;
    }

    long long cayley = 1;
    for (int i = 0; i < n - 2; i++) cayley *= n;
    assert((long long)trees.size() == cayley);
}

void check_random(int n, bool with_naive){
    auto edges = random_tree(n);
    auto code = encode(n, edges);
    assert((int)code.size() == n - 2);
    if (with_naive) assert(code == naive_encode(n, edges));

    vector<int> degree(n, -1);
    for (auto [u, v] : edges) degree[u]++, degree[v]++;
    for (int v : code) degree[v]--;
    for (int v = 0; v < n; v++) assert(degree[v] == 0);

    auto decoded = prufer_decode(code);
    check_rooted_tree(n, decoded);
    assert(normalized(decoded) == normalized(edges));
}

const int ABORT_EXIT_CODE = 77;

/// encode must abort on a non-tree or n < 2, decode on a value outside [0, n)
/// The child turns SIGABRT into a plain exit, a core dump per fork makes the enumeration crawl
template <typename Run>
void check_aborts(Run run){
    pid_t pid = fork();
    assert(pid >= 0);
    if (pid == 0){
        alarm(5);
        signal(SIGABRT, [](int){ _exit(ABORT_EXIT_CODE); });
        assert(freopen("/dev/null", "w", stderr));
        run();
        _exit(0);
    }

    int status;
    assert(waitpid(pid, &status, 0) == pid);
    assert(WIFEXITED(status) && WEXITSTATUS(status) == ABORT_EXIT_CODE);
}

/// Every multiset of n - 1 edges over [0, n), self-loops and repeats included: trees encode, the rest abort
void check_all_edge_multisets(int n){
    vector<pair<int, int>> pairs;
    for (int u = 0; u < n; u++){
        for (int v = u; v < n; v++) pairs.push_back({u, v});
    }

    int m = pairs.size();
    vector<int> pick(n - 1, 0);
    while (true){
        vector<pair<int, int>> edges;
        for (int i : pick) edges.push_back(pairs[i]);

        vector<int> root(n);
        iota(root.begin(), root.end(), 0);
        auto find = [&](int v){
            while (root[v] != v) v = root[v];
            return v;
        };
        bool is_tree = true;
        for (auto [u, v] : edges){
            int ru = find(u), rv = find(v);
            if (ru == rv) is_tree = false;
            root[ru] = rv;
        }

        if (is_tree) assert(encode(n, edges) == naive_encode(n, edges));
        else check_aborts([&]{ encode(n, edges); });

        int i = n - 2;
        while (i >= 0 && pick[i] == m - 1) i--;
        if (i < 0) break;
        pick[i]++;
        for (int j = i + 1; j < n - 1; j++) pick[j] = pick[i];
    }
}

void check_rejects(){
    for (int n = 2; n <= 4; n++) check_all_edge_multisets(n);

    vector<pair<int, vector<pair<int, int>>>> bad = {
        {3, {{0, 3}, {1, 2}}},
        {3, {{-1, 0}, {1, 2}}},
        {4, {{0, 1}, {1, 2}, {2, 0}}},
        {4, {{3, 1}, {1, 2}, {2, 3}}},
        {5, {{0, 1}, {1, 2}, {2, 3}, {3, 1}}},
        {6, {{0, 1}, {2, 3}, {3, 4}, {4, 5}, {5, 2}}},
        {3, {{0, 1}, {1, 2}, {2, 0}}},
        {3, {{0, 1}, {0, 1}}},
        {3, {{0, 0}, {1, 2}}},
        {3, {{0, 1}}},
        {1, {}},
    };
    for (auto& [n, edges] : bad) check_aborts([&]{ encode(n, edges); });

    check_aborts([]{ prufer_decode({3, 4}); });
    check_aborts([]{ prufer_decode({-1}); });
}

int main(){
    for (int n = 2; n <= 7; n++) check_exhaustive(n);

    for (long long it = 0; it < stress::scaled(2000); it++){
        int n = it < 500 ? it % 60 + 2 : stress::rand_int(2, 200);
        check_random(n, true);
    }

    for (long long it = 0; it < stress::scaled(100); it++) check_random(stress::rand_int(2, 5000), false);

    for (long long it = 0; it < stress::scaled(4); it++) check_random(stress::rand_int(100000, 200000), false);

    check_rejects();
    return 0;
}
