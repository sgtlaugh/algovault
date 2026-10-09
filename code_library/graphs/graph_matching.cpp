/***
 *
 * Maximum matching in general graphs
 * Two structs with the same API: Graph (randomized Tutte matrix rank, size only) and Blossom (Edmonds, the pairs too)
 *
 * Complexity: O(n ^ 3) for Graph even on sparse graphs, O(n ^ 2) memory
 *             O(n * (n ^ 2 + m)) for Blossom, so O(n ^ 3) without repeated edges, O(n + m) memory
 *
 * Graph g(n), g.add_edge(u, v), g.maximum_matching(): the rank of the Tutte matrix is twice the size of the
 * maximum matching with probability >= 1 - n / mod (Schwartz-Zippel), a wrong answer is always too small
 *
 * Blossom b(n), b.add_edge(u, v), b.maximum_matching(): size of a maximum matching, always exact
 * b.mate[v]: vertex matched to v or -1, valid after maximum_matching(), which can be called again after adding edges
 * Each search grows an alternating tree from one free vertex and contracts odd cycles (blossoms) by relabelling base[],
 * at most n / 2 contractions of O(n) each per search and n searches, which is the O(n ^ 3) bound
 * The search is a BFS, so long augmenting paths cannot overflow the stack
 *
**/

#include <stdio.h>
#include <bits/stdtr1c++.h>

using namespace std;

struct Graph{
    static constexpr unsigned int mod = 1073750017; /// prime, constexpr so % compiles to a multiply

    int n;
    vector<vector<bool>> adj;
    vector<vector<unsigned int>> tutte_matrix;
    mt19937 rng = mt19937(chrono::steady_clock::now().time_since_epoch().count());

    Graph() {}
    Graph(int n): n(n){
        adj = vector<vector<bool>>(n, vector<bool>(n, 0));
        tutte_matrix = vector<vector<unsigned int>>(n, vector<unsigned int>(n, 0));
    }

    void add_edge(int u, int v){
        assert(u != v && u >= 0 && u < n && v >= 0 && v < n);
        adj[u][v] = adj[v][u] = 1;
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

    unsigned int expo(unsigned long long x, unsigned int n){
        unsigned long long res = 1;

        while (n){
            if (n & 1) res = res * x % mod;
            x = x * x % mod;
            n >>= 1;
        }

        return res;
    }

    int get_matrix_rank(){
        int j, k, u, r = 0;
        vector<int> nonzero;

        for (j = 0; j < n; j++){
            for (k = r; k < n && !tutte_matrix[k][j]; k++) {}
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

            for (u = r + 1; u < n; u++){
                auto& row = tutte_matrix[u];
                if (!row[j]) continue;
                unsigned long long f = mod - row[j];
                for (int v : nonzero) row[v] = (row[v] + f * pivot[v]) % mod;
            }
            r++;
        }

        return r;
    }

    int maximum_matching(){
        build_matrix();
        return get_matrix_rank() / 2;
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

    /// Returns the free vertex ending an augmenting path from root, parent[] then leads back along it, or -1 if none exists
    int find_augmenting_path(int root){
        fill(parent.begin(), parent.end(), -1);
        fill(in_tree.begin(), in_tree.end(), false);
        iota(base.begin(), base.end(), 0);
        /// A search contracts at most n / 2 blossoms, so resetting here keeps stamp far from overflow on long-lived instances
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

    /// Walks from v up to the blossom base lca, flagging the bases it crosses and pointing odd vertices back across the cycle
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
    auto is_matching = [](const Blossom& b, int size){
        int pairs = 0;
        for (int v = 0; v < b.n; v++){
            int u = b.mate[v];
            if (u == -1) continue;
            if (u < 0 || u >= b.n || b.mate[u] != v || find(b.adj[v].begin(), b.adj[v].end(), u) == b.adj[v].end()) return false;
            pairs += u > v;
        }
        return pairs == size;
    };

    auto g = Graph(12);
    auto b = Blossom(12);
    auto add_edge = [&](int u, int v){
        g.add_edge(u, v);
        b.add_edge(u, v);
    };
    auto matches = [&](int want){
        return g.maximum_matching() == want && b.maximum_matching() == want && is_matching(b, want);
    };
    assert(matches(0));

    add_edge(0, 1);
    assert(matches(1));
    add_edge(1, 2);
    assert(matches(1));
    add_edge(2, 0);
    assert(matches(1));

    add_edge(2, 3);
    assert(matches(2));
    add_edge(3, 4);
    assert(matches(2));

    add_edge(4, 5);
    assert(matches(3));
    add_edge(5, 6);
    assert(matches(3));

    add_edge(7, 8);
    assert(matches(4));
    add_edge(7, 9);
    assert(matches(4));
    add_edge(7, 10);
    assert(matches(4));

    add_edge(4, 10);
    assert(matches(5));

    add_edge(9, 11);
    assert(matches(6));
    add_edge(8, 11);
    assert(matches(6));

    auto cycle = Blossom(5);
    for (int i = 0; i < 5; i++) cycle.add_edge(i, (i + 1) % 5);
    assert(cycle.maximum_matching() == 2 && is_matching(cycle, 2));

    auto petersen = Blossom(10);
    for (int i = 0; i < 5; i++){
        petersen.add_edge(i, (i + 1) % 5);
        petersen.add_edge(i, i + 5);
        petersen.add_edge(i + 5, (i + 2) % 5 + 5);
    }
    assert(petersen.maximum_matching() == 5 && is_matching(petersen, 5));

    /// Greedy pairs 0-1, 2-3, 4-5, the only augmenting path 6-0-1-5-4-3-2-7 leaves the blossom 1-2-3-4-5 from 2,
    /// which a search without contraction labels odd and never expands
    auto flower = Blossom(8);
    for (auto [u, v] : vector<pair<int, int>>{{0, 1}, {2, 3}, {4, 5}, {1, 2}, {3, 4}, {5, 1}, {6, 0}, {2, 7}}) flower.add_edge(u, v);
    assert(flower.maximum_matching() == 4 && is_matching(flower, 4));
    assert((flower.mate == vector<int>{6, 5, 7, 4, 3, 1, 0, 2}));

    auto lonely = Blossom(1);
    assert(lonely.maximum_matching() == 0 && lonely.mate[0] == -1);

    srand(666);
    clock_t start = clock();
    int n = 2000, m = 5000;
    g = Graph(n);
    b = Blossom(n);

    for (int i = 0; i < m; i++){
        int u = 0, v = 0;
        while (u == v) u = rand() % n, v = rand() % n;
        add_edge(u, v);
    }

    assert(g.maximum_matching() == 992);
    fprintf(stderr, "Time taken = %0.5f\n", (clock() - start) / (double)CLOCKS_PER_SEC);  /// 1.35343 s

    start = clock();
    assert(b.maximum_matching() == 992 && is_matching(b, 992));
    fprintf(stderr, "Blossom time taken = %0.5f\n", (clock() - start) / (double)CLOCKS_PER_SEC);

    /// A matching never exceeds n / 2 and is_matching certifies the pairs, so 750 is the answer without trusting Blossom
    n = 1500;
    auto dense = Blossom(n);
    for (int u = 0; u < n; u++){
        for (int v = u + 1; v < n; v++){
            if (rand() % 2) dense.add_edge(u, v);
        }
    }

    start = clock();
    assert(dense.maximum_matching() == n / 2 && is_matching(dense, n / 2));
    fprintf(stderr, "Dense Blossom time taken = %0.5f\n", (clock() - start) / (double)CLOCKS_PER_SEC);

    return 0;
}
