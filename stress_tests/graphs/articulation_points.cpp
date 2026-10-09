#include "../common.h"

#define main library_main
#include "../../code_library/graphs/articulation_points.cpp"
#undef main

/// Component label of every vertex once `removed` is deleted, the removed vertex gets -1
vector<int> labels(int n, const vector<pair<int, int>>& edges, int removed){
    vector<int> parent(n);
    iota(parent.begin(), parent.end(), 0);
    auto find = [&](int x){
        while (parent[x] != x) x = parent[x] = parent[parent[x]];
        return x;
    };

    for (auto [u, v] : edges){
        if (u != removed && v != removed) parent[find(u)] = find(v);
    }

    vector<int> res(n);
    for (int v = 0; v < n; v++) res[v] = v == removed ? -1 : find(v);
    return res;
}

int components(int n, const vector<pair<int, int>>& edges, int removed){
    auto label = labels(n, edges, removed);
    int res = 0;
    for (int v = 0; v < n; v++) res += label[v] == v;
    return res;
}

/// Same class number for edges in the same block, numbered by first appearance, -1 for self loops
vector<int> canonical(const vector<int>& block_of){
    map<int, int> renumber;
    vector<int> res;
    for (int b : block_of) res.push_back(b < 0 ? -1 : renumber.emplace(b, renumber.size()).first->second);
    return res;
}

/// Two edges share a block iff no vertex x leaves them in different components of G - x
vector<int> separator_blocks(int n, const vector<pair<int, int>>& edges){
    int m = edges.size();
    vector<int> cls(m, 0);
    for (int e = 0; e < m; e++){
        if (edges[e].first == edges[e].second) cls[e] = -1;
    }

    for (int x = 0; x < n; x++){
        auto label = labels(n, edges, x);
        map<pair<int, int>, int> refined;
        for (int e = 0; e < m; e++){
            if (cls[e] < 0) continue;
            auto [a, b] = edges[e];
            int side = label[a == x ? b : a];
            cls[e] = refined.emplace(make_pair(cls[e], side), refined.size()).first->second;
        }
    }
    return canonical(cls);
}

/// Two edges share a block iff some simple cycle passes through both, by enumerating every simple cycle
vector<int> cycle_blocks(int n, const vector<pair<int, int>>& edges){
    int m = edges.size();
    vector<vector<pair<int, int>>> adj(n);
    for (int e = 0; e < m; e++){
        adj[edges[e].first].push_back({edges[e].second, e});
        adj[edges[e].second].push_back({edges[e].first, e});
    }

    vector<vector<bool>> share(m, vector<bool>(m, false));
    vector<bool> on_path(n, false);
    vector<int> path;
    int s = 0;
    function<void(int)> extend = [&](int u){
        for (auto [w, e] : adj[u]){
            if (w == s && !path.empty() && e != path[0]){
                path.push_back(e);
                for (int a : path) for (int b : path) share[a][b] = true;
                path.pop_back();
            }
            if (w > s && !on_path[w]){
                on_path[w] = true;
                path.push_back(e);
                extend(w);
                path.pop_back();
                on_path[w] = false;
            }
        }
    };
    for (s = 0; s < n; s++){
        on_path[s] = true;
        extend(s);
        on_path[s] = false;
    }

    vector<int> cls(m, -1);
    int count = 0;
    for (int e = 0; e < m; e++){
        if (edges[e].first == edges[e].second || cls[e] >= 0) continue;
        for (int f = e; f < m; f++){
            if (f == e || share[e][f]) cls[f] = count;
        }
        count++;
    }
    return cls;
}

/// Every non-loop edge sits in exactly one block, self loops in none
vector<int> block_of(int m, const vector<pair<int, int>>& edges, const vector<vector<int>>& blocks){
    vector<int> res(m, -1);
    for (int b = 0; b < (int)blocks.size(); b++){
        assert(!blocks[b].empty());
        for (int e : blocks[b]){
            assert(0 <= e && e < m && res[e] == -1 && edges[e].first != edges[e].second);
            res[e] = b;
        }
    }
    for (int e = 0; e < m; e++) assert((res[e] == -1) == (edges[e].first == edges[e].second));
    return canonical(res);
}

void check_tree(int n, const vector<pair<int, int>>& edges, const vector<vector<int>>& blocks,
                const vector<vector<int>>& tree, const vector<int>& cuts, bool check_paths){
    int nodes = n + blocks.size();
    assert((int)tree.size() == nodes);

    long long tree_edges = 0;
    for (int b = 0; b < (int)blocks.size(); b++){
        set<int> touched;
        for (int e : blocks[b]) touched.insert(edges[e].first), touched.insert(edges[e].second);
        auto linked = tree[n + b];
        sort(linked.begin(), linked.end());
        assert(linked == vector<int>(touched.begin(), touched.end()));
        for (int v : touched) assert(count(tree[v].begin(), tree[v].end(), n + b) == 1);
        tree_edges += touched.size();
    }

    vector<bool> is_cut(n, false);
    for (int v : cuts) is_cut[v] = true;
    long long vertex_side = 0;
    for (int v = 0; v < n; v++){
        assert(((int)tree[v].size() >= 2) == is_cut[v]);
        vertex_side += tree[v].size();
    }
    assert(vertex_side == tree_edges);

    /// A forest: edges = nodes - trees, and tree connectivity matches graph connectivity
    vector<int> comp(nodes, -1), parent(nodes, -1);
    int trees = 0;
    for (int s = 0; s < nodes; s++){
        if (comp[s] != -1) continue;
        vector<int> stack = {s};
        comp[s] = trees;
        while (!stack.empty()){
            int u = stack.back();
            stack.pop_back();
            for (int v : tree[u]){
                if (comp[v] == -1) comp[v] = trees, stack.push_back(v);
            }
        }
        trees++;
    }
    assert(tree_edges == nodes - trees);

    auto label = labels(n, edges, -1);
    for (int u = 0; u < n; u++){
        for (int w = 0; w < n; w++) assert((label[u] == label[w]) == (comp[u] == comp[w]));
    }
    if (!check_paths) return;

    /// x separates u from w in the graph iff x is a cut vertex strictly inside the tree path from u to w
    vector<vector<int>> removed_labels;
    for (int x = 0; x < n; x++) removed_labels.push_back(labels(n, edges, x));
    for (int u = 0; u < n; u++){
        fill(parent.begin(), parent.end(), -1);
        vector<int> queue = {u};
        parent[u] = u;
        for (size_t i = 0; i < queue.size(); i++){
            for (int v : tree[queue[i]]){
                if (parent[v] == -1) parent[v] = queue[i], queue.push_back(v);
            }
        }

        for (int w = 0; w < n; w++){
            if (w == u || label[u] != label[w]) continue;
            vector<bool> inside(n, false);
            for (int x = parent[w]; x != u; x = parent[x]){
                if (x < n) inside[x] = true;
            }
            for (int x = 0; x < n; x++){
                if (x == u || x == w) continue;
                assert(inside[x] == (removed_labels[x][u] != removed_labels[x][w]));
            }
        }
    }
}

vector<pair<int, int>> random_graph(int n, bool simple){
    set<pair<int, int>> edges;
    int density = stress::rand_int(0, n > 12 ? 3 : 100), shape = stress::rand_int(0, 3);
    auto add = [&](int u, int v){ if (u != v) edges.insert({min(u, v), max(u, v)}); };

    if (shape == 3 && n >= 2){
        /// Short cycles glued at shared nodes, the bowties where back edge handling decides a cut
        for (int c = stress::rand_int(1, n); c; c--){
            vector<int> cycle(stress::rand_int(2, min(n, 5)));
            for (auto& x : cycle) x = stress::rand_int(0, n - 1);
            for (size_t i = 0; i < cycle.size(); i++) add(cycle[i], cycle[(i + 1) % cycle.size()]);
        }
    }
    for (int u = 0; u < n && shape != 3; u++){
        for (int v = u + 1; v < n; v++){
            bool take = stress::rand_int(0, 99) < density;
            if (shape == 1) take = take && v == u + 1;                     /// paths with gaps, every inner node is a cut
            if (shape == 2) take = v == u + 1 || (take && density < 15);  /// a path with a few chords
            if (take) add(u, v);
        }
    }

    vector<pair<int, int>> res(edges.begin(), edges.end());
    if (!simple){
        /// Doubled edges turn bridges into 2-edge blocks, self loops must change nothing
        int doubles = stress::rand_int(0, res.size()), loops = stress::rand_int(0, 2);
        for (int i = 0; i < doubles; i++) res.push_back(res[stress::rand_int(0, res.size() - 1)]);
        for (int i = 0; i < loops; i++){
            int v = stress::rand_int(0, n - 1);
            res.push_back({v, v});
        }
    }
    for (auto& [u, v] : res){
        if (stress::rand_int(0, 1)) swap(u, v);
    }
    shuffle(res.begin(), res.end(), stress::rng());

    /// Relabel so structure does not follow node order
    vector<int> label(n);
    iota(label.begin(), label.end(), 0);
    shuffle(label.begin(), label.end(), stress::rng());
    for (auto& [u, v] : res) u = label[u], v = label[v];
    return res;
}

vector<pair<int, int>> tiny_graph(int n){
    vector<pair<int, int>> res(stress::rand_int(0, 10));
    for (auto& [u, v] : res){
        u = stress::rand_int(0, n - 1);
        v = stress::rand_int(0, n - 1);
        if (u == v && stress::rand_int(0, 3)) v = (u + 1) % n;
    }
    return res;
}

void check(int n, const vector<pair<int, int>>& edges, bool by_cycles){
    int m = edges.size();
    Graph g(n);

    /// Query halfway through the edges too, a second call must not keep results from the first
    int half = stress::rand_int(0, m);
    for (int i = 0; i < m; i++){
        if (i == half) g.get_cuts(), g.get_blocks();
        assert(g.add_edge(edges[i].first, edges[i].second) == i);
    }

    vector<int> expected;
    int base = components(n, edges, -1);
    for (int v = 0; v < n; v++){
        if (components(n, edges, v) > base) expected.push_back(v);
    }
    assert(g.get_cuts() == expected);

    auto blocks = g.get_blocks();
    auto got = block_of(m, edges, blocks);
    assert(got == separator_blocks(n, edges));
    if (by_cycles) assert(got == cycle_blocks(n, edges));

    check_tree(n, edges, blocks, g.get_block_cut_tree(), expected, n <= 40);
}

/// Cuts against vertex removal, blocks against vertex separators and against simple cycle enumeration,
/// the block-cut tree against both
int main(){
    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(1, it % 20 ? 12 : 250);
        check(n, random_graph(n, it % 2), false);
    }

    for (long long it = 0; it < stress::scaled(4000); it++){
        int n = stress::rand_int(1, 7);
        check(n, it % 3 ? tiny_graph(n) : random_graph(n, it % 2), true);
    }

    /// Past the old fixed MAX = 100010 arrays
    {
        int n = 150000;
        Graph g(n);
        for (int v = 1; v < n; v++) g.add_edge(0, v);
        assert(g.get_cuts() == vector<int>({0}));
    }

    /// DFS depth n, a recursive search overflows the 8 MB stack here
    {
        int n = 200000;
        Graph path(n), cycle(n);
        for (int v = 0; v + 1 < n; v++) path.add_edge(v, v + 1), cycle.add_edge(v, v + 1);
        cycle.add_edge(n - 1, 0);

        auto cuts = path.get_cuts();
        assert((int)cuts.size() == n - 2 && cuts.front() == 1 && cuts.back() == n - 2);
        assert((int)path.get_blocks().size() == n - 1);
        auto tree = path.get_block_cut_tree();
        assert((int)tree.size() == 2 * n - 1 && tree[0].size() == 1u && tree[n / 2].size() == 2u);

        assert(cycle.get_cuts().empty());
        auto blocks = cycle.get_blocks();
        assert(blocks.size() == 1u && (int)blocks[0].size() == n);
    }
    return 0;
}
