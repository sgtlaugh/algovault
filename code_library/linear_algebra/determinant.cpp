/***
 *
 * Integer Determinant
 * Determinant of an integer matrix modulo any m, and the exact determinant when it is not too large
 *
 * Complexity: O(n^3 + n^2 log m) per determinant or count
 *
 * determinant(a): the exact signed determinant, correct whenever |det(a)| <= 2258704744122758558 (~2.26e18)
 * determinant_mod(a, m): det(a) mod m in [0, m) for any modulus 1 <= m < 2^62, prime or composite
 *
 * SpanningTrees g(n): undirected multigraph on nodes 0 to n - 1, counted by Kirchhoff's matrix tree theorem
 * g.add_edge(u, v): adds the edge u - v, parallel edges count separately and self loops are ignored
 * g.count(m): the number of spanning trees mod m, 0 when disconnected and 1 for a single node
 *
 * Arborescences d(n): directed multigraph on nodes 0 to n - 1, counted by Tutte's matrix tree theorem
 * d.add_edge(u, v): adds the edge u -> v, parallel edges count separately and self loops are ignored
 * d.count(root, m): the number of spanning trees mod m with every edge pointing away from root
 *
 * Notes:
 *   - Entries may be any long long, a 0 x 0 matrix has determinant 1
 *   - Counts are exact when below m, pass m = (1LL << 62) - 1 for counts up to ~4.6e18
 *   - For trees pointing towards root, add every edge reversed
 *   - Requires __int128 (64-bit GCC or Clang)
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

int main(){
    assert(determinant({{1, 2}, {3, 4}}) == -2);
    assert(determinant({{6, 1, 1}, {4, -2, 5}, {2, 8, 7}}) == -306);
    assert(determinant_mod({{1, 2}, {3, 4}}, 7) == 5);  /// -2 mod 7
    assert(determinant_mod({{2, 1}, {1, 3}}, 6) == 5);  /// composite modulus, 2 and 3 have no inverse mod 6

    const long long MOD = 998244353;
    SpanningTrees k4(4);
    for (int u = 0; u < 4; u++){
        for (int v = u + 1; v < 4; v++) k4.add_edge(u, v);
    }
    assert(k4.count(MOD) == 16);  /// Cayley, 4^(4 - 2)

    SpanningTrees forest(4);
    forest.add_edge(0, 1), forest.add_edge(2, 3);
    assert(forest.count(MOD) == 0);  /// disconnected

    Arborescences path(3);
    path.add_edge(0, 1), path.add_edge(1, 2);
    assert(path.count(0, MOD) == 1);
    assert(path.count(2, MOD) == 0);  /// edges must point away from the root
    return 0;
}
