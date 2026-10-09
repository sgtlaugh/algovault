/***
 *
 * Directed minimum spanning tree (minimum arborescence), Tarjan / Gabow version of Chu–Liu/Edmonds
 * Constructs a rooted tree of minimum total weight from the root node
 * Returns -1 if no solution from root
 * 0 based index for nodes, non-negative weights since -1 is the unreachable sentinel
 *
 * Each node keeps a skew heap of its incoming edges with a lazy add, cycles are contracted with a DSU
 * by melding their heaps, so no round ever rescans the whole edge list
 *
 * Complexity: O(m log n)
 *
***/

#include <stdio.h>
#include <bits/stdtr1c++.h>

using namespace std;

struct Edge{
    int u, v, w;

    Edge(){}
    Edge(int u, int v, int w) : u(u), v(v), w(w) {}
};

long long directed_mst(int n, int root, const vector <Edge>& edges){
    for (auto edge: edges) assert(edge.w >= 0);

    struct Node{
        long long key, lazy;
        int l, r;
    };

    int m = edges.size();
    vector <Node> t(m);
    vector <int> heap(n, -1), dsu(n), seen(n, -1), path(n), spine;
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

    long long res = 0;
    seen[root] = root;
    for (int s = 0; s < n; s++){
        int u = s, len = 0;
        while (seen[u] == -1){
            if (heap[u] == -1) return -1;

            int e = heap[u];
            push(e);
            res += t[e].key;
            heap[u] = merge(t[e].l, t[e].r);
            if (heap[u] != -1) t[heap[u]].lazy -= t[e].key;

            path[len++] = u, seen[u] = s;
            u = find(edges[e].u);
            if (seen[u] == s){
                int cycle = -1, x;
                do{
                    x = path[--len];
                    cycle = merge(cycle, heap[x]);
                    if (x != u) dsu[x] = u;
                } while (x != u);
                heap[u] = cycle, seen[u] = -1;
            }
        }
    }

    return res;
}

int main(){
    vector <Edge> edges = {
        Edge(0, 1, 5),
        Edge(1, 2, 10),
        Edge(1, 3, 5),
        Edge(3, 2, 1),
        Edge(1, 4, 2),
        Edge(2, 3, 1),
        Edge(0, 2, 1),
        Edge(4, 0, 15),
    };

    assert(directed_mst(5, 0, edges) == 9);
    assert(directed_mst(5, 1, edges) == 19);
    assert(directed_mst(5, 2, edges) == -1);
    assert(directed_mst(5, 3, edges) == -1);
    assert(directed_mst(5, 4, edges) == 22);

    return 0;
}
