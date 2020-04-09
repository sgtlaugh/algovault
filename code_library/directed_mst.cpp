#include <stdio.h>
#include <bits/stdtr1c++.h>

using namespace std;

struct Edge{
    int u, v, w;

    Edge(){}
    Edge(int u, int v, int w) : u(u), v(v), w(w) {}
};

/***
 *
 * directed minimum spanning tree in O(n * m)
 * constructs a rooted tree of minimum total weight from the root node
 * returns -1 if no solution from root
 * 0 based index for nodes
 *
***/

long long directed_mst(int n, int root, vector <Edge> edges){
    const int INF = INT_MAX;
    for (auto edge: edges) assert(edge.w >= 0 && edge.w < INF);

    int i, k, x, y;
    long long res = 0;
    vector <int> cost(n), parent(n), label(n), comp(n);

    for (; ;){
        for (i = 0; i < n; i++) cost[i] = INF;
        for (auto edge: edges){
            if (edge.u != edge.v && cost[edge.v] > edge.w){
                cost[edge.v] = edge.w;
                parent[edge.v] = edge.u;
            }
        }

        cost[root] = 0;
        for (i = 0; i < n && cost[i] != INF; i++){};
        if (i != n) return -1;

        for (i = 0, k = 0; i < n; i++) res += cost[i];
        for (i = 0; i < n; i++) label[i] = comp[i] = -1;

        for (i = 0; i < n; i++){
            for (x = i; x != root && comp[x] == -1; x = parent[x]) comp[x] = i;
            if (x != root && comp[x] == i){
                for (k++; label[x] == -1; x = parent[x]) label[x] = k - 1;
            }
        }

        if (k == 0) break;
        for (i = 0; i < n; i++){
            if (label[i] == -1) label[i] = k++;
        }

        for (auto &edge: edges){
            x = label[edge.u], y = label[edge.v];
            if (x != y) edge.w -= cost[edge.v];
            edge.u = x, edge.v = y;
        }
        root = label[root], n = k;
    }
    return res;
}

int main(){
    vector <Edge> edges = {
        Edge(0, 1, 5), Edge(1, 2, 10), Edge(1, 3, 5), Edge(3, 2, 1),
        Edge(1, 4, 2), Edge(2, 3, 1), Edge(0, 2, 1), Edge(4, 0, 15)
    };

    printf("%lld\n", directed_mst(5, 0, edges));  /// 9
    printf("%lld\n", directed_mst(5, 1, edges));  /// 19
    printf("%lld\n", directed_mst(5, 2, edges));  /// -1
    printf("%lld\n", directed_mst(5, 3, edges));  /// -1
    printf("%lld\n", directed_mst(5, 4, edges));  /// 22

    return 0;
}
