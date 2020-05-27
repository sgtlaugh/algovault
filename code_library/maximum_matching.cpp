/***
 *
 * Maximum matching in general graphs
 * Uses a randomized algorithm to compute the rank of the Tutte matrix
 * The rank of the Tutte matrix is equal to twice the size of the maximum matching with high probability
 *
 * Complexity: O(n ^ 3) worst case, O(n * m) with high constant factor on average
 *
**/

#include <stdio.h>
#include <bits/stdtr1c++.h>

using namespace std;

struct Graph{
    int n, mod; /// mod must be a prime

    vector <vector<bool>> adj;
    vector <vector<int>> tutte_matrix;
    mt19937 rng = mt19937(chrono::steady_clock::now().time_since_epoch().count());

    Graph() {}
    Graph(int n, int mod=1073750017): n(n), mod(mod) {
        adj = vector<vector<bool>>(n, vector<bool>(n, 0));
        tutte_matrix = vector<vector<int>>(n, vector<int>(n, 0));
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

    int expo(long long x, int n){
        long long res = 1;

        while (n){
            if (n & 1) res = res * x % mod;
            x = x * x % mod;
            n >>= 1;
        }

        return res % mod;
    }

    int get_matrix_rank(){
        int i, j, k, u, v, x, r = 0;

        for (j = 0; j < n; j++){
            for (k = r; k < n && !tutte_matrix[k][j]; k++) {}
            if (k == n) continue;

            long long inv = expo(tutte_matrix[k][j], mod - 2);
            for (i = 0; i < n; i++){
                x = tutte_matrix[k][i];
                tutte_matrix[k][i] = tutte_matrix[r][i];
                tutte_matrix[r][i] = inv * x % mod;
            }

            for (u = r + 1; u < n; u++){
                for (v = j + 1; v < n && tutte_matrix[u][j]; v++){
                    if (tutte_matrix[r][v]){
                        tutte_matrix[u][v] = (tutte_matrix[u][v] - (long long)tutte_matrix[r][v] * tutte_matrix[u][j]) % mod;
                        if (tutte_matrix[u][v] < 0) tutte_matrix[u][v] += mod;
                    }
                }
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
