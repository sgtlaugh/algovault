/***
 *
 * Directed minimum spanning tree (minimum arborescence), Tarjan / Gabow version of Chu–Liu/Edmonds
 * Constructs a rooted tree of minimum total weight from the root node, with the parent of every node
 *
 * Complexity: O(m log m) time, O(n + m) memory
 *
 * DirectedMST g(n); g.add_edge(u, v, w); auto res = g.solve(root);
 *     0-based nodes, long long weights (negative allowed), n * max|w| <= 1e18
 *     res.exists is false when some node is unreachable from root (weight 0, parent empty)
 *     otherwise res.weight is the minimum total and res.parent[v] the tail of v's tree edge, -1 for root
 *     solve is const, so one graph answers several roots
 *
 * Each node keeps a skew heap of its incoming edges with a lazy add, cycles are contracted with a DSU
 * by melding their heaps, so no round ever rescans the whole edge list
 * Every contracted cycle gets a fresh id, and the contraction forest gives the parents back in O(n):
 * the edge chosen by a contracted node enters one leaf, every node on the way down to that leaf drops
 * its own choice, every other node keeps it
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct Edge{
    int u, v;
    long long w;

    Edge(){}
    Edge(int u, int v, long long w) : u(u), v(v), w(w) {}
};

struct Arborescence{
    bool exists;
    long long weight;
    vector<int> parent;
};

struct DirectedMST{
    int n;
    vector<Edge> edges;

    DirectedMST(int n) : n(n) {}

    void add_edge(int u, int v, long long w){
        edges.push_back(Edge(u, v, w));
    }

    Arborescence solve(int root) const{
        struct Node{
            long long key, lazy;
            int l, r;
        };

        int m = edges.size();
        vector<Node> t(m);
        vector<int> heap(2 * n, -1), dsu(2 * n), up(2 * n, -1), in(2 * n, -1), seen(2 * n, -1), path(n), spine;
        iota(dsu.begin(), dsu.end(), 0);

        auto find = [&](int x){
            while (dsu[x] != x) x = dsu[x] = dsu[dsu[x]];
            return x;
        };

        auto push = [&](int a){
            if (t[a].l != -1) t[t[a].l].lazy += t[a].lazy;
            if (t[a].r != -1) t[t[a].r].lazy += t[a].lazy;
            t[a].key += t[a].lazy, t[a].lazy = 0;
        };

        /// Iterative because a single skew heap merge can walk an O(m) right spine
        auto merge = [&](int a, int b){
            spine.clear();
            while (a != -1 && b != -1){
                push(a), push(b);
                if (t[a].key > t[b].key) swap(a, b);
                spine.push_back(a);
                a = t[a].r;
            }

            int res = a != -1 ? a : b;
            while (!spine.empty()){
                int x = spine.back();
                spine.pop_back();
                t[x].r = t[x].l, t[x].l = res, res = x;
            }
            return res;
        };

        for (int i = 0; i < m; i++){
            t[i] = {edges[i].w, 0, -1, -1};
            heap[edges[i].v] = merge(heap[edges[i].v], i);
        }

        long long weight = 0;
        int nodes = n;
        seen[root] = root;
        for (int s = 0; s < n; s++){
            int u = s, len = 0;
            while (seen[u] == -1){
                if (heap[u] == -1) return {false, 0, {}};

                int e = heap[u];
                push(e);
                weight += t[e].key;
                heap[u] = merge(t[e].l, t[e].r);
                if (heap[u] != -1) t[heap[u]].lazy -= t[e].key;

                in[u] = e, path[len++] = u, seen[u] = s;
                u = find(edges[e].u);
                if (seen[u] != s) continue;

                /// A self loop of a contracted node is no cycle, popping it again just takes its next edge
                if (u == path[len - 1]){
                    len--, seen[u] = -1;
                    continue;
                }

                int cycle = nodes++, x;
                do{
                    x = path[--len];
                    heap[cycle] = merge(heap[cycle], heap[x]);
                    dsu[x] = up[x] = cycle;
                } while (x != u);
                u = cycle;
            }
        }

        /// Parents are built before children, so a node is known to be dropped before it is reached
        vector<int> parent(n, -1);
        vector<bool> dropped(nodes, false);
        for (int x = nodes - 1; x >= 0; x--){
            if (x == root || dropped[x]) continue;
            int e = in[x], v = edges[e].v;
            parent[v] = edges[e].u;
            for (int y = v; y != x; y = up[y]) dropped[y] = true;
        }

        return {true, weight, parent};
    }
};

int main(){
    DirectedMST g(5);
    g.add_edge(0, 1, 5);
    g.add_edge(1, 2, 10);
    g.add_edge(1, 3, 5);
    g.add_edge(3, 2, 1);
    g.add_edge(1, 4, 2);
    g.add_edge(2, 3, 1);
    g.add_edge(0, 2, 1);
    g.add_edge(4, 0, 15);

    Arborescence res = g.solve(0);
    assert(res.exists && res.weight == 9 && (res.parent == vector<int>{-1, 0, 0, 2, 1}));
    res = g.solve(1);
    assert(res.exists && res.weight == 19 && (res.parent == vector<int>{4, -1, 0, 2, 1}));
    res = g.solve(4);
    assert(res.exists && res.weight == 22 && (res.parent == vector<int>{4, 0, 0, 2, -1}));
    assert(!g.solve(2).exists && !g.solve(3).exists);

    DirectedMST single(1);
    res = single.solve(0);
    assert(res.exists && res.weight == 0 && (res.parent == vector<int>{-1}));

    DirectedMST loops(2);
    loops.add_edge(1, 1, -5);
    loops.add_edge(0, 1, 3);
    loops.add_edge(0, 0, -7);
    res = loops.solve(0);
    assert(res.exists && res.weight == 3 && (res.parent == vector<int>{-1, 0}));
    assert(!loops.solve(1).exists);

    DirectedMST negative(3);
    negative.add_edge(0, 1, -3);
    negative.add_edge(1, 2, 2);
    negative.add_edge(0, 2, 5);
    negative.add_edge(2, 1, -4);
    res = negative.solve(0);
    assert(res.exists && res.weight == -1 && (res.parent == vector<int>{-1, 0, 1}));

    return 0;
}
