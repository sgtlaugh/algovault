/***
 *
 * Integer Determinant
 * Determinant of an integer matrix modulo any m, and the exact determinant when it is not too large
 *
 * Complexity: O(n^3 + n^2 log m)
 *
 * determinant_mod(a, m): det(a) mod m in [0, m) for any modulus 1 <= m < 2^62, prime or composite
 *     row reduction by repeated division, like the Euclidean algorithm, so no modular inverses are needed
 * determinant(a): the exact signed determinant, correct whenever |det(a)| <= 2258704744122758558 (~2.26e18)
 *     computed modulo the prime 4517409488245517117 and mapped back to the nearest signed value
 *
 * Entries may be any long long, negative included, a 0 x 0 matrix has determinant 1
 * Requires __int128 (64-bit GCC or Clang)
 *
 * Matrix tree theorem on top of determinant_mod, 0 based nodes, O(n^2) memory, count() is O(n^3 + n^2 log m)
 * SpanningTrees (Kirchhoff): number of spanning trees of an undirected multigraph mod m
 *     0 when disconnected, 1 for a single node, parallel edges count separately, self loops are ignored
 * Arborescences (Tutte): number of spanning trees of a directed multigraph rooted at root mod m,
 *     with every edge pointing away from root, so every other node has exactly one incoming edge
 *     for trees pointing towards root, add every edge reversed
 * Counts are exact when below the modulus, pass m = (1LL << 62) - 1 for counts up to ~4.6e18
 *
 *     SpanningTrees g(n);  g.add_edge(u, v);  g.count(MOD);
 *     Arborescences d(n);  d.add_edge(u, v);  d.count(root, MOD);
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

long long determinant_mod(vector<vector<long long>> a, long long m){
    int n = a.size();
    assert(1 <= m && m < (1LL << 62));
    for (auto& row : a){
        assert((int)row.size() == n);
        for (auto& x : row) x %= m;
    }

    __int128 res = 1;
    for (int i = 0; i < n; i++){
        for (int j = i + 1; j < n; j++){
            while (a[j][i] != 0){
                long long t = a[i][i] / a[j][i];
                if (t){
                    for (int k = i; k < n; k++) a[i][k] = (long long)((a[i][k] - (__int128)a[j][k] * t) % m);
                }
                swap(a[i], a[j]);
                res = -res;
            }
        }
        res = res * a[i][i] % m;
        if (res == 0) return 0;
    }

    res %= m;
    return (long long)(res < 0 ? res + m : res);
}

long long determinant(const vector<vector<long long>>& a){
    const long long P = 4517409488245517117LL;
    long long r = determinant_mod(a, P);
    return r > P / 2 ? r - P : r;
}

/// self loops cancel out in the Laplacian, so add_edge needs no u == v special case
long long laplacian_minor_det(const vector<vector<long long>>& lap, int skip, long long m){
    int n = lap.size();
    vector<vector<long long>> minor;
    minor.reserve(n - 1);

    for (int i = 0; i < n; i++){
        if (i == skip) continue;
        vector<long long> row;
        row.reserve(n - 1);
        for (int j = 0; j < n; j++){
            if (j != skip) row.push_back(lap[i][j]);
        }
        minor.push_back(row);
    }

    return determinant_mod(minor, m);
}

struct SpanningTrees{
    int n;
    vector<vector<long long>> lap;

    SpanningTrees(int n) : n(n), lap(n, vector<long long>(n, 0)){
        assert(n >= 1);
    }

    void add_edge(int u, int v){
        assert(0 <= u && u < n && 0 <= v && v < n);
        lap[u][u]++, lap[v][v]++;
        lap[u][v]--, lap[v][u]--;
    }

    long long count(long long m) const{
        return laplacian_minor_det(lap, 0, m);
    }
};

struct Arborescences{
    int n;
    vector<vector<long long>> lap;

    Arborescences(int n) : n(n), lap(n, vector<long long>(n, 0)){
        assert(n >= 1);
    }

    void add_edge(int u, int v){
        assert(0 <= u && u < n && 0 <= v && v < n);
        lap[v][v]++;
        lap[u][v]--;
    }

    long long count(int root, long long m) const{
        assert(0 <= root && root < n);
        return laplacian_minor_det(lap, root, m);
    }
};

const long long MOD = 998244353;

long long power(long long b, long long e){
    long long r = 1;
    for (b %= MOD; e; e >>= 1, b = b * b % MOD){
        if (e & 1) r = r * b % MOD;
    }

    return r;
}

int main(){
    int t;

    scanf("%d", &t);
    while (t--){
        int n, m;
        long long k;
        scanf("%d %d %lld", &n, &m, &k);
        long long v = n * k % MOD;

        vector<vector<long long>> a(n, vector<long long>(n, 0));
        for (int i = 0; i < n; i++) a[i][i] = v;
        while (m--){
            int x, y;
            scanf("%d %d", &x, &y);
            x--, y--;
            a[x][x]--, a[y][y]--;
            a[x][y]++, a[y][x]++;
        }

        long long inv = power(v, MOD - 2);
        printf("%lld\n", power(determinant_mod(a, MOD), k) * inv % MOD * inv % MOD);
    }

    return 0;
}
