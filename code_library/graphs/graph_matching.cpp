/***
 *
 * General Graph Matching
 * Maximum matching in a general graph: Blossom (Edmonds) gives the pairs, TutteMatching (randomized) only the size
 *
 * Complexity:
 *   - Blossom: O(n (n^2 + m)), so O(n^3) without repeated edges, O(n + m) memory
 *   - TutteMatching: O(n^3) even on sparse graphs, O(n^2) memory
 *
 * Blossom b(n): vertices 0 to n - 1, no edges
 * b.add_edge(u, v): undirected edge, u != v, so skip self loops when reading the input
 * b.maximum_matching(): the size of a maximum matching, callable again after adding edges
 * b.mate[v]: after maximum_matching(), the vertex matched to v or -1
 *
 * TutteMatching g(n): vertices 0 to n - 1, no edges
 * g.add_edge(u, v): undirected edge, u != v
 * g.maximum_matching(): the maximum matching size with probability >= 1 - n / mod, a wrong answer is too small
 *
 * Notes:
 *   - Use Blossom: exact and faster in practice, TutteMatching is a size-only cross-check
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct TutteMatching{
    static constexpr unsigned int mod = 1073750017; /// prime, constexpr so % compiles to a multiply

    int n;
    vector<vector<bool>> adj;
    vector<vector<unsigned int>> tutte_matrix;
    mt19937 rng = mt19937(chrono::steady_clock::now().time_since_epoch().count());

    TutteMatching(int n): n(n){
        adj = vector<vector<bool>>(n, vector<bool>(n, 0));
        tutte_matrix = vector<vector<unsigned int>>(n, vector<unsigned int>(n, 0));
    }

    void add_edge(int u, int v){
        assert(u != v && u >= 0 && u < n && v >= 0 && v < n);
        adj[u][v] = adj[v][u] = 1;
    }

    int maximum_matching(){
        build_matrix();
        return get_matrix_rank() / 2;
    }

    void build_matrix(){
        for (auto &v : tutte_matrix) fill(v.begin(), v.end(), 0);
        for (int i = 0; i < n; i++){
            for (int j = i + 1; j < n; j++){
                if (adj[i][j]){
                    unsigned int v = rng() % (mod - 1) + 1;
                    tutte_matrix[i][j] = v, tutte_matrix[j][i] = mod - v;
                }
            }
        }
    }

    int get_matrix_rank(){
        int r = 0;
        vector<int> nonzero;

        for (int j = 0; j < n; j++){
            int k = r;
            while (k < n && !tutte_matrix[k][j]) k++;
            if (k == n) continue;

            swap(tutte_matrix[k], tutte_matrix[r]);
            auto& pivot = tutte_matrix[r];
            unsigned long long inv = expo(pivot[j], mod - 2);

            nonzero.clear();
            for (int v = j + 1; v < n; v++){
                if (pivot[v]){
                    pivot[v] = inv * pivot[v] % mod;
                    nonzero.push_back(v);
                }
            }

            for (int u = r + 1; u < n; u++){
                auto& row = tutte_matrix[u];
                if (!row[j]) continue;
                unsigned long long f = mod - row[j];
                for (int v : nonzero) row[v] = (row[v] + f * pivot[v]) % mod;
            }
            r++;
        }

        return r;
    }

    unsigned int expo(unsigned long long x, unsigned int n){
        unsigned long long res = 1;

        while (n){
            if (n & 1) res = res * x % mod;
            x = x * x % mod;
            n >>= 1;
        }

        return res;
    }
};

struct Blossom{
    int n;
    vector<vector<int>> adj;
    vector<int> mate, parent, base, seen_stamp;
    vector<bool> in_tree, in_blossom;
    int stamp = 0;

    Blossom(int n) : n(n), adj(n), mate(n, -1), parent(n), base(n), seen_stamp(n, 0), in_tree(n), in_blossom(n) {}

    void add_edge(int u, int v){
        assert(u != v && u >= 0 && u < n && v >= 0 && v < n);
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    int maximum_matching(){
        fill(mate.begin(), mate.end(), -1);
        for (int u = 0; u < n; u++){
            for (int v : adj[u]){
                if (mate[u] == -1 && mate[v] == -1) mate[u] = v, mate[v] = u;
            }
        }

        for (int root = 0; root < n; root++){
            if (mate[root] != -1) continue;
            for (int v = find_augmenting_path(root); v != -1;){
                int next = mate[parent[v]];
                mate[v] = parent[v], mate[parent[v]] = v;
                v = next;
            }
        }

        int size = 0;
        for (int v = 0; v < n; v++) size += mate[v] > v;
        return size;
    }

    /// Returns the free vertex ending an augmenting path from root, parent[] leads back along it, or -1 if none exists
    int find_augmenting_path(int root){
        fill(parent.begin(), parent.end(), -1);
        fill(in_tree.begin(), in_tree.end(), false);
        iota(base.begin(), base.end(), 0);
        /// A search contracts at most n / 2 blossoms, so resetting here keeps stamp far from overflow
        stamp = 0;
        fill(seen_stamp.begin(), seen_stamp.end(), 0);

        vector<int> queue = {root};
        in_tree[root] = true;
        for (size_t i = 0; i < queue.size(); i++){
            int v = queue[i];
            for (int to : adj[v]){
                if (base[v] == base[to] || mate[v] == to) continue;

                /// to is an even (outer) vertex of the tree, the edge closes an odd cycle
                if (to == root || (mate[to] != -1 && parent[mate[to]] != -1)){
                    int lca = lowest_common_base(v, to);
                    fill(in_blossom.begin(), in_blossom.end(), false);
                    mark_path(v, lca, to);
                    mark_path(to, lca, v);
                    for (int u = 0; u < n; u++){
                        if (!in_blossom[base[u]]) continue;
                        base[u] = lca;
                        if (!in_tree[u]) in_tree[u] = true, queue.push_back(u);
                    }
                }
                else if (parent[to] == -1){
                    parent[to] = v;
                    if (mate[to] == -1) return to;
                    in_tree[mate[to]] = true;
                    queue.push_back(mate[to]);
                }
            }
        }

        return -1;
    }

    int lowest_common_base(int u, int v){
        stamp++;
        while (true){
            u = base[u];
            seen_stamp[u] = stamp;
            if (mate[u] == -1) break;
            u = parent[mate[u]];
        }

        while (true){
            v = base[v];
            if (seen_stamp[v] == stamp) return v;
            v = parent[mate[v]];
        }
    }

    /// Walks from v up to base lca, flagging the bases it crosses and pointing odd vertices back across the cycle
    void mark_path(int v, int lca, int child){
        while (base[v] != lca){
            in_blossom[base[v]] = in_blossom[base[mate[v]]] = true;
            parent[v] = child;
            child = mate[v];
            v = parent[mate[v]];
        }
    }
};

int main(){
    /// A triangle 0-1-2 with the path 2-3-4-5: only 0-1, 2-3, 4-5 match every vertex
    Blossom b(6);
    b.add_edge(0, 1);
    b.add_edge(1, 2);
    b.add_edge(2, 0);
    b.add_edge(2, 3);
    b.add_edge(3, 4);
    b.add_edge(4, 5);
    assert(b.maximum_matching() == 3);
    assert((b.mate == vector<int>{1, 0, 3, 2, 5, 4}));

    TutteMatching g(6);
    g.add_edge(0, 1);
    g.add_edge(1, 2);
    g.add_edge(2, 0);
    g.add_edge(2, 3);
    g.add_edge(3, 4);
    g.add_edge(4, 5);
    assert(g.maximum_matching() == 3);

    return 0;
}
