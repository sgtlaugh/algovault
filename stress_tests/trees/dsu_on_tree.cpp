#include "../common.h"

#define main library_main
#include "../../code_library/trees/dsu_on_tree.cpp"
#undef main

vector<pair<int, int>> make_tree(int n, int shape){
    vector<int> label(n);
    iota(label.begin(), label.end(), 0);
    shuffle(label.begin(), label.end(), stress::rng());
    vector<pair<int, int>> edges;
    for (int i = 1; i < n; i++){
        int p = shape == 0 ? stress::rand_int(0, i - 1) : shape == 1 ? i - 1 : shape == 2 ? 0 : (i - 1) / 2;
        edges.push_back({label[i], label[p]});
    }
    return edges;
}

/// Distinct colors and color sum per subtree, and the exact window contents at each answer, against a DFS
void check(int n, int shape, int colors){
    auto edges = make_tree(n, shape);
    int root = stress::rand_int(0, n - 1);
    DsuOnTree tree(n);
    for (auto [u, v] : edges) tree.add_edge(u, v);

    vector<int> color(n);
    for (auto& c : color) c = stress::rand_int(0, colors - 1);
    vector<int> freq(colors, 0), in_window(n, 0), distinct(n, -1);
    vector<long long> sum(n, -1);
    int count = 0, window = 0;
    long long total = 0;
    tree.run(root,
        [&](int v){ assert(!in_window[v]); in_window[v] = 1, window++, total += color[v], count += freq[color[v]]++ == 0; },
        [&](int v){ assert(in_window[v]); in_window[v] = 0, window--, total -= color[v], count -= --freq[color[v]] == 0; },
        [&](int v){ distinct[v] = count, sum[v] = total; assert(window == tree.size[v]); });
    assert(window == 0 && count == 0 && total == 0);

    vector<vector<int>> children(n);
    for (int v = 0; v < n; v++){
        if (v != root) children[tree.parent[v]].push_back(v);
    }
    for (int v = 0; v < n; v++){
        set<int> seen;
        long long s = 0;
        vector<int> stack = {v};
        int nodes = 0;
        while (!stack.empty()){
            int u = stack.back();
            stack.pop_back();
            seen.insert(color[u]), s += color[u], nodes++;
            for (int c : children[u]) stack.push_back(c);
        }
        assert(distinct[v] == (int)seen.size() && sum[v] == s && nodes == tree.size[v]);
        assert(tree.tout[v] - tree.tin[v] + 1 == nodes);
        for (int i = tree.tin[v]; i <= tree.tout[v]; i++){
            int u = tree.order[i];
            while (u != v && u != -1) u = tree.parent[u];
            assert(u == v);
        }
    }
}

int main(){
    for (long long it = 0; it < stress::scaled(3000); it++) check(stress::rand_int(1, 40), it % 4, it % 2 ? 3 : 30);
    check(2000, 0, 50);
    check(1500, 1, 50);

    /// Running twice, from different roots, must give the same answers as fresh objects
    for (long long it = 0; it < stress::scaled(300); it++){
        int n = stress::rand_int(1, 30);
        auto edges = make_tree(n, it % 4);
        DsuOnTree reused(n);
        for (auto [u, v] : edges) reused.add_edge(u, v);
        for (int round = 0; round < 3; round++){
            int root = stress::rand_int(0, n - 1);
            DsuOnTree fresh(n);
            for (auto [u, v] : edges) fresh.add_edge(u, v);
            vector<int> a(n), b(n);
            int wa = 0, wb = 0;
            reused.run(root, [&](int){ wa++; }, [&](int){ wa--; }, [&](int v){ a[v] = wa; });
            fresh.run(root, [&](int){ wb++; }, [&](int){ wb--; }, [&](int v){ b[v] = wb; });
            assert(a == b && wa == 0);
        }
    }

    /// A comb: spine edges first, so a wrong heavy child choice costs O(n^2) adds instead of O(n log n)
    for (int order = 0; order < 2; order++){
        int spine = 2000, n = 2 * spine;
        DsuOnTree comb(n);
        vector<pair<int, int>> edges;
        for (int i = 1; i < spine; i++) edges.push_back({i - 1, i});
        for (int i = 0; i < spine; i++) edges.push_back({i, spine + i});
        if (order) reverse(edges.begin(), edges.end());
        for (auto [u, v] : edges) comb.add_edge(u, v);
        long long adds = 0;
        comb.run(0, [&](int){ adds++; }, [&](int){}, [&](int){});
        assert(adds <= (long long)n * (__lg(n) + 1));
    }

    /// A path of 200000 nodes: deep enough to overflow an 8 MB stack with recursion
    int n = 200000;
    DsuOnTree path(n);
    for (int i = 1; i < n; i++) path.add_edge(i - 1, i);
    long long adds = 0, answers = 0;
    path.run(0, [&](int){ adds++; }, [&](int){}, [&](int){ answers++; });
    assert(answers == n && adds == n);
    return 0;
}
