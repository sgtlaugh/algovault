/***
 *
 * Online Dynamic Connectivity (Holm, de Lichtenberg, Thorup)
 * Connectivity, component sizes and component value sums of an undirected graph under edge insertions and deletions
 *
 * Complexity:
 *   - add_edge, remove_edge: amortized O(log^2 n)
 *   - every other call: amortized O(log n)
 *   - memory: O((n + m) log n) for m add_edge calls, about 430 MB at n = m = 3e5 for a path cut at random points
 *
 * OnlineDynamicConnectivity g(n): an empty graph on vertices 0 to n - 1, every value 0
 * OnlineDynamicConnectivity g(values): an empty graph on vertices 0 to n - 1 with the given values
 * g.add(v, x): adds x to the value of v
 * g.add_edge(u, v): adds a copy of the edge (u, v), parallel edges and self loops allowed
 * g.remove_edge(u, v): deletes one live copy of (u, v) in either orientation, the copy must exist
 * g.connected(u, v): u and v are in the same connected component
 * g.size(v): number of vertices in the component of v
 * g.sum(v): sum of the values in the component of v
 * g.components: number of connected components
 *
 * Notes:
 *   - Fully online; when all operations are known in advance, offline_dynamic_connectivity.cpp is simpler and faster
 *   - Everything is iterative, deep trees are safe
 *   - Algorithm from J. Holm, K. de Lichtenberg, M. Thorup, J. ACM 48(4), 2001, with Euler tour trees on splay trees
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct OnlineDynamicConnectivity{
    enum Mark : unsigned char{TREE = 1, NONTREE = 2};

    /// Level i keeps a spanning forest of the edges of level >= i, each tree as its cyclic Euler tour in one splay tree
    struct Node{
        int ch[2], par;
        int id;  /// a vertex slot if >= 0, else one of the two arcs of edge ~id
        int size;
        long long sum;
        unsigned char mark, any;  /// TREE on one arc of each tree edge of this level, NONTREE on a vertex with a list
    };

    struct Edge{
        int u, v, level;  /// level -1 once deleted
        vector<int> arcs;  /// arcs[i] and arcs[i] + 1 are its arcs in level i's forest
    };

    int n, components;
    vector<long long> value;
    vector<Node> t;  /// nodes 0..n - 1 are the vertices of level 0, the only nodes holding values
    vector<int> vertex_of;  /// per vertex slot, its vertex; slot k is the k-th vertex node made, slot v < n is node v
    vector<vector<int>> nontree;  /// per vertex slot, edges once non-tree at its level; stale when their level moved on
    vector<vector<int>> vertex_at;  /// per level, the node of each vertex, -1 until it gets a tree edge there
    vector<Edge> edges;
    vector<int> spare_arcs;
    map<pair<int, int>, vector<int>> live;

    OnlineDynamicConnectivity(int n) : OnlineDynamicConnectivity(vector<long long>(n)) {}

    OnlineDynamicConnectivity(const vector<long long>& values)
        : n(values.size()), components(n), value(values), vertex_at(1, vector<int>(n, -1)){
        for (int v = 0; v < n; v++) node(0, v);
    }

    void add(int v, long long x){
        splay(v);
        value[v] += x;
        pull(v);
    }

    void add_edge(int u, int v){
        int e = edges.size();
        edges.push_back({u, v, 0, {}});
        live[minmax(u, v)].push_back(e);
        if (u == v) return;

        if (same_tree(u, v)) add_nontree(e);
        else{
            link(e);
            components--;
        }
    }

    void remove_edge(int u, int v){
        auto it = live.find(minmax(u, v));
        assert(it != live.end());
        int e = it->second.back();
        it->second.pop_back();
        if (it->second.empty()) live.erase(it);

        int level = edges[e].level;
        edges[e].level = -1;
        if (edges[e].arcs.empty()) return;

        for (int a: edges[e].arcs) cut(a);
        edges[e].arcs.clear();
        for (int i = level; i >= 0; i--){
            if (replace(i, u, v)) return;
        }
        components++;
    }

    bool connected(int u, int v){
        return same_tree(u, v);
    }

    int size(int v){
        splay(v);
        return t[v].size;
    }

    long long sum(int v){
        splay(v);
        return t[v].sum;
    }

    int node(int i, int v){
        if (i == (int)vertex_at.size()) vertex_at.emplace_back(n, -1);
        if (vertex_at[i][v] == -1){
            vertex_at[i][v] = t.size();
            t.push_back({{-1, -1}, -1, (int)vertex_of.size(), 1, i == 0 ? value[v] : 0, 0, 0});
            vertex_of.push_back(v);
            nontree.emplace_back();
        }
        return vertex_at[i][v];
    }

    bool same_tree(int a, int b){
        splay(a);
        splay(b);
        return a == b || t[a].par != -1;
    }

    void add_nontree(int e){
        int i = edges[e].level;
        for (int w: {edges[e].u, edges[e].v}){
            int x = node(i, w);
            nontree[t[x].id].push_back(e);
            set_mark(x, NONTREE);
        }
    }

    /// Splices e into the forests of levels arcs.size()..level, its endpoints must be in different trees there
    void link(int e){
        while ((int)edges[e].arcs.size() <= edges[e].level){
            int i = edges[e].arcs.size(), a = new_arcs(e, i == edges[e].level ? TREE : 0);
            edges[e].arcs.push_back(a);

            /// v's tour A v B becomes A v (v->u) (u's tour from u) (u->v) B
            int ru = make_first(node(i, edges[e].u)), v = node(i, edges[e].v);
            splay(v);
            int after = t[v].ch[1];
            if (after != -1) t[after].par = -1;
            t[v].ch[1] = -1;
            join(join(join(v, a), ru), join(a + 1, after));
        }
    }

    void cut(int a){
        /// The two arcs split the cyclic tour into the part between them and two outer parts that close into one tour
        splay(a);
        splay(a + 1, a);
        int d = t[a].ch[1] == a + 1;
        int between = t[a + 1].ch[!d], outer_a = t[a].ch[!d], outer_b = t[a + 1].ch[d];
        for (int c: {between, outer_a, outer_b}){
            if (c != -1) t[c].par = -1;
        }
        join(outer_a, outer_b);
        spare_arcs.push_back(a);
    }

    /// After cutting a tree edge u-v of level >= i: searches level i for a replacement, raising what it scans to i + 1
    bool replace(int i, int u, int v){
        int a = node(i, u), b = node(i, v);
        splay(a);
        splay(b);
        if (t[a].size > t[b].size) swap(a, b);

        /// Raising the smaller tree's edges keeps every tree of level i + 1 at most half a tree of level i
        for (int x; (x = marked(a, TREE)) != -1; a = x){
            set_mark(x, 0);
            int e = ~t[x].id;
            edges[e].level++;
            link(e);
        }

        for (int x; (x = marked(a, NONTREE)) != -1; a = x){
            int slot = t[x].id, found = -1;
            while (found == -1 && !nontree[slot].empty()){
                int e = nontree[slot].back();
                nontree[slot].pop_back();
                /// Tree edges of level i rise in the TREE pass first, so a level check alone spots stale entries
                if (edges[e].level != i) continue;

                int y = node(i, edges[e].u ^ edges[e].v ^ vertex_of[slot]);
                if (!same_tree(x, y)) found = e;
                else{
                    edges[e].level++;
                    add_nontree(e);
                }
            }
            set_mark(x, nontree[slot].empty() ? 0 : NONTREE);

            if (found != -1){
                link(found);
                return true;
            }
        }

        return false;
    }

    /// A node of x's tree with the mark bit, splayed to the root, or -1
    int marked(int x, int bit){
        splay(x);
        if (!(t[x].any & bit)) return -1;

        while (!(t[x].mark & bit)){
            int l = t[x].ch[0];
            x = l != -1 && (t[l].any & bit) ? l : t[x].ch[1];
        }
        splay(x);
        return x;
    }

    void set_mark(int x, unsigned char mark){
        splay(x);
        t[x].mark = mark;
        pull(x);
    }

    int new_arcs(int e, unsigned char mark){
        int a = t.size();
        if (!spare_arcs.empty()){
            a = spare_arcs.back();
            spare_arcs.pop_back();
        }
        else{
            t.resize(a + 2);
        }

        t[a] = {{-1, -1}, -1, ~e, 0, 0, mark, mark};
        t[a + 1] = {{-1, -1}, -1, ~e, 0, 0, 0, 0};
        return a;
    }

    /// The tour is cyclic, so rotating it to start at any node keeps it a valid tour; returns the new root
    int make_first(int x){
        splay(x);
        int left = t[x].ch[0];
        if (left == -1) return x;

        t[left].par = -1;
        t[x].ch[0] = -1;
        pull(x);
        return join(x, left);
    }

    /// Joins the tours rooted at a and b, a first; returns the new root
    int join(int a, int b){
        if (a == -1 || b == -1) return a == -1 ? b : a;
        int x = a;
        while (t[x].ch[1] != -1) x = t[x].ch[1];

        splay(x);
        t[x].ch[1] = b;
        t[b].par = x;
        pull(x);
        return x;
    }

    void splay(int x, int goal = -1){
        while (t[x].par != goal){
            int p = t[x].par, g = t[p].par;
            if (g != goal) rotate((t[g].ch[1] == p) == (t[p].ch[1] == x) ? p : x);
            rotate(x);
        }
        pull(x);
    }

    /// Pulls only p, x stays stale until splay pulls it once at the end
    void rotate(int x){
        int p = t[x].par, g = t[p].par, d = t[p].ch[1] == x, child = t[x].ch[!d];

        if (g != -1) t[g].ch[t[g].ch[1] == p] = x;
        t[x].par = g;
        t[p].ch[d] = child;
        if (child != -1) t[child].par = p;
        t[x].ch[!d] = p;
        t[p].par = x;

        pull(p);
    }

    void pull(int x){
        Node& c = t[x];
        c.size = c.id >= 0;
        c.sum = x < n ? value[x] : 0;
        c.any = c.mark;
        for (int y: c.ch){
            if (y == -1) continue;
            c.size += t[y].size;
            c.sum += t[y].sum;
            c.any |= t[y].any;
        }
    }
};

int main(){
    /// A cycle 0-1-2-3 with a tail 3-4: deleting cycle edges finds replacements until the cycle is gone
    OnlineDynamicConnectivity g(vector<long long>{5, 3, 8, 1, 7, 2});
    g.add_edge(0, 1), g.add_edge(1, 2), g.add_edge(2, 3), g.add_edge(3, 0), g.add_edge(3, 4);
    assert(g.components == 2);
    assert(g.size(4) == 5);
    assert(g.sum(0) == 24);

    g.remove_edge(1, 0);                            /// either orientation, 0-3-2-1 still holds
    assert(g.connected(0, 1));
    g.remove_edge(2, 3);
    assert(!g.connected(0, 1));
    assert(g.connected(0, 4));
    assert(g.sum(1) == 11);                         /// {1, 2}: 3 + 8
    assert(g.sum(0) == 13);                         /// {0, 3, 4}: 5 + 1 + 7
    assert(g.components == 3);

    g.add(2, 10);
    assert(g.sum(1) == 21);
    assert(g.sum(5) == 2);

    /// Parallel copies: removing one leaves the other live, self loops change nothing
    g.add_edge(4, 5), g.add_edge(5, 4), g.add_edge(5, 5);
    g.remove_edge(4, 5);
    assert(g.connected(0, 5));
    assert(g.size(5) == 4);
    g.remove_edge(5, 5);
    g.remove_edge(5, 4);
    assert(!g.connected(0, 5));
    assert(g.components == 3);

    return 0;
}
