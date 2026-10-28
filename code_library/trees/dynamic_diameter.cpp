/***
 *
 * Dynamic Tree Diameter
 * Maintains the diameter of a weighted tree while edge weights change, 0-based nodes
 *
 * Complexity: O(n) to build, O(log n) per weight update, O(1) per diameter query, O(n) memory, about 315 MB at n = 1e6
 *
 * DynamicDiameter tree(n): an empty tree on n >= 1 nodes
 * tree.add_edge(u, v, w): adds the edge u - v of weight w, returns its id, ids are 0, 1, ..., n - 2 in insertion order
 * tree.build(): call after all n - 1 edges
 * tree.update(id, w): sets the weight of edge id to w
 * tree.diameter(): the longest path length, 0 for a single node
 *
 * Notes:
 *   - Weights must be non-negative and their sum must stay at most 2e18
 *   - Traversals are iterative so deep trees are safe
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct DynamicDiameter{
    /// g(i): tour depth at i minus the depth before the range; xyz = max g(i) - 2 g(j) + g(k) over i <= j <= k
    struct Node{
        long long sum, x, y, xy, yz, xyz;
    };

    int n, leaves = 0;
    vector<vector<pair<int, int>>> adj;
    vector<long long> weight;
    vector<int> enter, leave;
    vector<Node> seg;

    DynamicDiameter(int n) : n(n), adj(n) {}

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
        enter.assign(n - 1, -1), leave.assign(n - 1, -1);
        vector<int> parent_edge(n, -1), next_edge(n, 0), stack = {0};
        vector<long long> dist(n, 0), tour = {0};

        /// Tour lists a node on entry and again after each child, so edge id covers [enter[id], leave[id] - 1]
        while (!stack.empty()){
            int u = stack.back();
            if (next_edge[u] < (int)adj[u].size()){
                auto [v, id] = adj[u][next_edge[u]++];
                if (id == parent_edge[u]) continue;
                assert(enter[id] == -1);
                parent_edge[v] = id, dist[v] = dist[u] + weight[id];
                enter[id] = tour.size(), tour.push_back(dist[v]);
                stack.push_back(v);
                continue;
            }
            stack.pop_back();
            if (stack.empty()) break;
            leave[parent_edge[u]] = tour.size(), tour.push_back(dist[stack.back()]);
        }
        assert((int)tour.size() == 2 * n - 1);

        /// A power of two keeps every leaf on one level, so seg[1] combines the whole tour left to right
        leaves = 1;
        while (leaves < (int)tour.size()) leaves *= 2;
        seg.assign(2 * leaves, Node{});
        for (int i = 1; i < (int)tour.size(); i++) seg[leaves + i] = leaf(tour[i] - tour[i - 1]);
        for (int i = leaves - 1; i >= 1; i--) seg[i] = combine(seg[2 * i], seg[2 * i + 1]);
    }

    /// Both leaves change before any parent is recombined, so no node ever sees a half-applied update
    void update(int id, long long w){
        assert(w >= 0);
        long long delta = w - weight[id];
        weight[id] = w;
        int l = leaves + enter[id], r = leaves + leave[id];
        seg[l] = leaf(seg[l].sum + delta), seg[r] = leaf(seg[r].sum - delta);
        while (l > 1){
            l /= 2, r /= 2;
            seg[l] = combine(seg[2 * l], seg[2 * l + 1]);
            if (r != l) seg[r] = combine(seg[2 * r], seg[2 * r + 1]);
        }
    }

    long long diameter() const{
        return seg[1].xyz;
    }

    static Node leaf(long long d){
        return {d, d, -2 * d, -d, -d, 0};
    }

    /// r's values are relative to the end of l, so each term shifts by l.sum once per g it holds
    static Node combine(const Node& l, const Node& r){
        long long s = l.sum;
        Node res;
        res.sum = s + r.sum;
        res.x = max(l.x, r.x + s);
        res.y = max(l.y, r.y - 2 * s);
        res.xy = max({l.xy, r.xy - s, l.x + (r.y - 2 * s)});
        res.yz = max({l.yz, r.yz - s, l.y + (r.x + s)});
        res.xyz = max({l.xyz, r.xyz, l.xy + (r.x + s), l.x + (r.yz - s)});
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
    for (auto [u, v, w] : edges) tree.add_edge(u, v, w);
    tree.build();
    assert(tree.diameter() == 15);  /// 6 3 1 0 2 5

    tree.update(3, 20);             /// edge ids follow insertion order, 3 is 1 - 4
    assert(tree.diameter() == 28);  /// 4 1 0 2 5
    tree.update(1, 0);              /// 0 - 2
    assert(tree.diameter() == 27);  /// 4 1 0 2 5 and 6 3 1 4 tie
    return 0;
}
