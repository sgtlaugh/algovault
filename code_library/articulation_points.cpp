#include <stdio.h>
#include <bits/stdtr1c++.h>

#define MAX 100010


using namespace std;

/***
 *
 * finds articulation points in simple graphs in O(N + M)
 * an articulation point or cut vertex is defined as a vertex which, when removed along with associated edges, makes the graph disconnected
 * or more precisely, increases the number of connected components in the graph
 * nodes are numbered from 0 to N-1
 *
***/

struct Graph{
    vector <int> adj[MAX];
    bool visited[MAX], is_cut[MAX];
    int n, disc_t, discover[MAX], low[MAX];

    Graph() {}
    Graph(int nodes){
        n = nodes;
        for (int i = 0; i < MAX; i++) adj[i].clear();
    }

    void dfs(int u, int p){
        visited[u] = true;
        discover[u] = low[u] = ++disc_t;

        int children = 0;
        for (auto v: adj[u]){
            if (!visited[v]){
                children++;
                dfs(v, u);
                low[u] = min(low[u], low[v]);

                if ((low[v] >= discover[u]) && !(p == -1 && children < 2)) {
                    is_cut[u] = true;
                }
            }
            else if (v != p) low[u] = min(low[u], discover[v]);
        }
    }

    /// adds undirected edge from u to v
    void add_edge(int u, int v){
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    vector <int> get_cuts(){
        memset(is_cut, 0, sizeof(is_cut));
        memset(visited, 0, sizeof(visited));

        vector <int> cuts;
        for (int i = 0; i < n; i++){
            if (!visited[i]){
                disc_t = 0;
                dfs(i, -1);
            }
            if (is_cut[i]) cuts.push_back(i);
        }

        return cuts;
    }
};


int main(){
    auto graph = Graph(10);

    graph.add_edge(0, 1);
    graph.add_edge(1, 2);
    graph.add_edge(2, 0);
    graph.add_edge(5, 6);
    graph.add_edge(6, 7);
    graph.add_edge(7, 6);
    graph.add_edge(2, 3);
    graph.add_edge(3, 4);
    graph.add_edge(2, 4);
    graph.add_edge(4, 5);
    graph.add_edge(8, 9);

    auto cuts = graph.get_cuts();
    printf("%d articulation points found\n\n", (int)cuts.size());
    for (auto node: cuts){
        printf("Node %d is an articulation point\n", node);
    }

    /***
     4 articulation points found

     Node 2 is an articulation point
     Node 4 is an articulation point
     Node 5 is an articulation point
     Node 6 is an articulation point

    ***/

    return 0;
}
