/***
 *
 * Dynamic Tree Diameter
 * Maintains the diameter of a weighted tree while edge weights change, 0-based nodes
 *
 * Complexity: O(n) to build, O(log n) per weight update, O(1) per diameter query, O(n) memory
 *
 * DynamicDiameter tree(n);           n >= 1
 * int id = tree.add_edge(u, v, w);   edge ids are 0, 1, ..., n - 2 in insertion order
 * tree.build();                      after all n - 1 edges
 * tree.update(id, w);                sets the weight of edge id to w
 * tree.diameter();                   longest path length, 0 for a single node
 *
 * Weights must be non-negative and their sum must stay at most 2e18
 *
 * With f[v] the distance from the root, the diameter is max f[x] - 2 f[y] + f[z] over positions x <= y <= z
 * of the Euler tour: the minimum f between two tour positions is at the lca, which needs non-negative weights
 * A weight change adds the same delta to f over the child's subtree, one contiguous tour range,
 * so a lazy segment tree keeps the five partial maxima of that expression per range
 * Traversals are iterative so deep trees are safe
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct DynamicDiameter{
    /// x: max f[x], y: max -2 f[y], xy: max f[x] - 2 f[y], yz: max -2 f[y] + f[z], xyz: max f[x] - 2 f[y] + f[z]
    struct Node{
        long long x, y, xy, yz, xyz;
    };

    int n, m = 0;
    vector<vector<pair<int, int>>> adj;
    vector<long long> weight, lazy;
    vector<int> child, first, last;
    vector<Node> seg;

    DynamicDiameter(int n) : n(n), adj(n){}

    int add_edge(int u, int v, long long w){
        assert(w >= 0);
        int id = weight.size();
        adj[u].push_back({v, id});
        adj[v].push_back({u, id});
        weight.push_back(w);
        return id;
    }

    void build(){
        assert((int)weight.size() == n - 1);
        child.assign(n - 1, -1), first.assign(n, -1), last.assign(n, -1);
        vector<int> parent_edge(n, -1), next_edge(n, 0), stack = {0};
        vector<long long> dist(n, 0), tour = {0};
        first[0] = 0;

        /// Tour lists a node on entry and again after each child, so the subtree of v is exactly [first[v], last[v]]
        while (!stack.empty()){
            int u = stack.back();
            if (next_edge[u] < (int)adj[u].size()){
                auto [v, id] = adj[u][next_edge[u]++];
                if (id == parent_edge[u]) continue;
                assert(first[v] == -1);
                parent_edge[v] = id, child[id] = v, dist[v] = dist[u] + weight[id];
                first[v] = tour.size(), tour.push_back(dist[v]);
                stack.push_back(v);
                continue;
            }
            last[u] = tour.size() - 1;
            stack.pop_back();
            if (!stack.empty()) tour.push_back(dist[stack.back()]);
        }
        assert((int)tour.size() == 2 * n - 1);

        m = tour.size();
        seg.assign(4 * m, Node{}), lazy.assign(4 * m, 0);
        build_seg(1, 0, m - 1, tour);
    }

    long long diameter() const{
        return seg[1].xyz;
    }

    void update(int id, long long w){
        assert(w >= 0);
        int v = child[id];
        add(1, 0, m - 1, first[v], last[v], w - weight[id]);
        weight[id] = w;
    }

    void add(int node, int b, int e, int l, int r, long long delta){
        if (r < b || e < l) return;
        if (l <= b && e <= r){
            apply(node, delta);
            return;
        }

        if (lazy[node]){
            apply(2 * node, lazy[node]), apply(2 * node + 1, lazy[node]);
            lazy[node] = 0;
        }

        int mid = (b + e) / 2;
        add(2 * node, b, mid, l, r, delta), add(2 * node + 1, mid + 1, e, l, r, delta);
        seg[node] = merge(seg[2 * node], seg[2 * node + 1]);
    }

    /// xyz is unchanged: f[x] - 2 f[y] + f[z] cancels a delta shared by the whole range
    void apply(int node, long long delta){
        Node& t = seg[node];
        t.x += delta, t.y -= 2 * delta, t.xy -= delta, t.yz -= delta;
        lazy[node] += delta;
    }

    void build_seg(int node, int b, int e, const vector<long long>& f){
        if (b == e){
            seg[node] = {f[b], -2 * f[b], -f[b], -f[b], 0};
            return;
        }

        int mid = (b + e) / 2;
        build_seg(2 * node, b, mid, f), build_seg(2 * node + 1, mid + 1, e, f);
        seg[node] = merge(seg[2 * node], seg[2 * node + 1]);
    }

    static Node merge(const Node& l, const Node& r){
        Node res;
        res.x = max(l.x, r.x);
        res.y = max(l.y, r.y);
        res.xy = max({l.xy, r.xy, l.x + r.y});
        res.yz = max({l.yz, r.yz, l.y + r.x});
        res.xyz = max({l.xyz, r.xyz, l.xy + r.x, l.x + r.yz});
        return res;
    }
};

int main(){
    /***
     *          0
     *     3  /   \  1
     *       1     2
     *    2 / \ 7   \ 4
     *     3   4     5
     *  5 /
     *   6
    ***/
    DynamicDiameter tree(7);
    vector<array<int, 3>> edges = {{0, 1, 3}, {0, 2, 1}, {1, 3, 2}, {1, 4, 7}, {2, 5, 4}, {3, 6, 5}};
    for (int i = 0; i < 6; i++) assert(tree.add_edge(edges[i][0], edges[i][1], edges[i][2]) == i);
    tree.build();
    assert(tree.diameter() == 15);

    tree.update(3, 20);
    assert(tree.diameter() == 28);
    tree.update(1, 0);
    assert(tree.diameter() == 27);
    tree.update(3, 0);
    assert(tree.diameter() == 14);
    for (int i = 0; i < 6; i++) tree.update(i, 0);
    assert(tree.diameter() == 0);
    tree.update(4, 9);
    assert(tree.diameter() == 9);

    DynamicDiameter single(1);
    single.build();
    assert(single.diameter() == 0);

    DynamicDiameter pair_tree(2);
    pair_tree.add_edge(1, 0, 5);
    pair_tree.build();
    assert(pair_tree.diameter() == 5);
    pair_tree.update(0, 0);
    assert(pair_tree.diameter() == 0);

    const long long big = 1000000000000000000LL;
    DynamicDiameter star(4);
    star.add_edge(0, 1, big), star.add_edge(0, 2, big), star.add_edge(0, 3, 0);
    star.build();
    assert(star.diameter() == 2 * big);
    star.update(0, 0), star.update(2, big);
    assert(star.diameter() == 2 * big);
    star.update(1, 0);
    assert(star.diameter() == big);

    DynamicDiameter path(5);
    for (int i = 0; i < 4; i++) path.add_edge(i, i + 1, big / 2);
    path.build();
    assert(path.diameter() == 2 * big);
    path.update(0, 0);
    assert(path.diameter() == 3 * big / 2);
    path.update(0, big / 2), path.update(3, 0);
    assert(path.diameter() == 3 * big / 2);
    path.update(3, big / 2);
    assert(path.diameter() == 2 * big);

    return 0;
}
