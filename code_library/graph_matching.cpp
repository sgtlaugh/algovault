/***
 *
 * Maximum matching in general graphs
 * Uses a randomized algorithm to compute the rank of the Tutte matrix
 * The rank of the Tutte matrix is equal to twice the size of the maximum matching with high probability
 *
 * Complexity: O(n ^ 3), even on sparse graphs, since it eliminates the dense n x n Tutte matrix
 *
**/

#include <stdio.h>
#include <bits/stdtr1c++.h>

using namespace std;

struct Graph{
    static constexpr unsigned int mod = 1073750017; /// prime, constexpr so % compiles to a multiply

    int n;
    vector <vector<bool>> adj;
    vector <vector<unsigned int>> tutte_matrix;
    mt19937 rng = mt19937(chrono::steady_clock::now().time_since_epoch().count());

    Graph() {}
    Graph(int n): n(n) {
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
        vector <int> nonzero;

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


int main(){
    auto g = Graph(12);
    assert(g.maximum_matching() == 0);

    g.add_edge(0, 1);
    assert(g.maximum_matching() == 1);
    g.add_edge(1, 2);
    assert(g.maximum_matching() == 1);
    g.add_edge(2, 0);
    assert(g.maximum_matching() == 1);


    g.add_edge(2, 3);
    assert(g.maximum_matching() == 2);
    g.add_edge(3, 4);
    assert(g.maximum_matching() == 2);

    g.add_edge(4, 5);
    assert(g.maximum_matching() == 3);
    g.add_edge(5, 6);
    assert(g.maximum_matching() == 3);

    g.add_edge(7, 8);
    assert(g.maximum_matching() == 4);
    g.add_edge(7, 9);
    assert(g.maximum_matching() == 4);
    g.add_edge(7, 10);
    assert(g.maximum_matching() == 4);


    g.add_edge(4, 10);
    assert(g.maximum_matching() == 5);

    g.add_edge(9, 11);
    assert(g.maximum_matching() == 6);
    g.add_edge(8, 11);
    assert(g.maximum_matching() == 6);

    srand(666);
    clock_t start = clock();
    int n = 2000, m = 5000;
    g = Graph(n);

    for (int i = 0; i < m; i++){
        int u = 0, v = 0;
        while (u == v) u = rand() % n, v = rand() % n;
        g.add_edge(u, v);
    }

    assert(g.maximum_matching() == 992);
    fprintf(stderr, "Time taken = %0.5f\n", (clock() - start) / (double)CLOCKS_PER_SEC);  /// 1.35343 s

    return 0;
}
