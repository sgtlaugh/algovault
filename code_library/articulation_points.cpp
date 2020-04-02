#include <bits/stdtr1c++.h>

#define MAX 100010

using namespace std;

/***
 *
 * Finds articulation points in O(N + M)
 * An articulation point or cut vertex is defined as a vertex which, when removed along with associated edges, makes the graph disconnected
 * Or more precisely, increases the number of connected components in the graph
 * Nodes are numbered from 0 to N-1
 *
***/

struct ArticulationPoint{
    vector <int> adj[MAX];
    bool visited[MAX], is_cut[MAX];
    int n, disc_t, discover[MAX], low[MAX];

    ArticulationPoint() {}
    ArticulationPoint(int nodes){
        n = nodes;
        for (int i = 0; i < MAX; i++) adj[i].clear();
    }

    void dfs(int i, int p){
        visited[i] = true;
        discover[i] = low[i] = ++disc_t;
        int j, x, children = 0, len = adj[i].size();

        for (j = 0; j < len; j++){
            x = adj[i][j];
            if (!visited[x]){
                children++;
                dfs(x, i);
                low[i] = min(low[i], low[x]);

                if ((low[x] >= discover[i])){
                    if (!(p == -1 && children < 2)){
                        is_cut[i] = true;
                    }
                }
            }
            else if (x != p) low[i] = min(low[i], discover[x]);
        }
    }

    void add_edge(int a, int b){
        adj[a].push_back(b);
        adj[b].push_back(a);
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
    auto ap = ArticulationPoint(10);
    ap.add_edge(0, 1);
    ap.add_edge(1, 2);
    ap.add_edge(2, 0);
    ap.add_edge(5, 6);
    ap.add_edge(6, 7);
    ap.add_edge(7, 6);
    ap.add_edge(2, 3);
    ap.add_edge(3, 4);
    ap.add_edge(2, 4);
    ap.add_edge(4, 5);
    ap.add_edge(8, 9);

    auto cuts = ap.get_cuts();
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
