/***
 *
 * Directed MST
 * Constructs a rooted tree of minimum total weight from the root node, with the parent of every node
 *
 * Complexity: O(m log m)
 *
 * DirectedMST g(n): nodes are 0 to n - 1
 * g.add_edge(u, v, w): directed edge u -> v with long long weight w, negative allowed
 * g.solve(root): {exists, weight, parent}, parent[v] the tail of v's tree edge and -1 at root
 *
 * Notes:
 *   - Requires n max|w| <= 1e18
 *   - exists is false when some node is unreachable from root, then weight is 0 and parent empty
 *   - solve is const, so one graph answers several roots
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct Arborescence{
    bool exists;
    long long weight;
    vector<int> parent;
};

struct DirectedMST{
    struct Edge{
        int u, v;
        long long w;
    };

    int n;
    vector<Edge> edges;

    DirectedMST(int n) : n(n) {}

    void add_edge(int u, int v, long long w){
        edges.push_back({u, v, w});
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
            if (!t[a].lazy) return;  /// most nodes carry no pending add, skipping saves two scattered child writes
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
    /// From 0 the tree takes 0 -> 1 -> 4 and 0 -> 2 -> 3, while from 2 only 3 is reachable
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
    assert(res.exists);
    assert(res.weight == 9);
    assert((res.parent == vector<int>{-1, 0, 0, 2, 1}));

    res = g.solve(1);
    assert(res.weight == 19);
    assert((res.parent == vector<int>{4, -1, 0, 2, 1}));

    assert(!g.solve(2).exists);
    return 0;
}
