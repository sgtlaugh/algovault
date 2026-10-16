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

typedef vector<vector<long long>> matrix;

long long X;

matrix identity(int n){
    matrix r(n, vector<long long>(n, 0));
    for (int i = 0; i < n; i++) r[i][i] = 1 % X;
    return r;
}

matrix add(const matrix& a, const matrix& b){
    int n = a.size();
    matrix r = a;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++) r[i][j] = (r[i][j] + b[i][j]) % X;
    }
    return r;
}

/// entries are below X < 2^23, so a row of 100 products fits in 64 bits before the single modulo
matrix mul(const matrix& a, const matrix& b){
    int n = a.size();
    matrix r(n, vector<long long>(n, 0));
    for (int i = 0; i < n; i++){
        for (int k = 0; k < n; k++){
            if (!a[i][k]) continue;
            for (int j = 0; j < n; j++) r[i][j] += a[i][k] * b[k][j];
        }
        for (int j = 0; j < n; j++) r[i][j] %= X;
    }
    return r;
}

/// sets sum = G^0 + ... + G^(k-1) and pw = G^k, all mod X
void geometric(const matrix& g, long long k, matrix& sum, matrix& pw){
    int n = g.size();
    if (k == 0){
        sum = matrix(n, vector<long long>(n, 0));
        pw = identity(n);
    }
    else if (k & 1){
        geometric(g, k - 1, sum, pw);
        sum = add(identity(n), mul(g, sum));
        pw = mul(g, pw);
    }
    else{
        geometric(g, k / 2, sum, pw);
        sum = mul(add(identity(n), pw), sum);
        pw = mul(pw, pw);
    }
}

bool is_prime(int p){
    for (int d = 2; d * d <= p; d++){
        if (p % d == 0) return false;
    }
    return p > 1;
}

long long power(long long b, long long e, long long m){
    long long r = 1;
    for (b %= m; e; e >>= 1, b = b * b % m){
        if (e & 1) r = r * b % m;
    }
    return r;
}

/// The count is exact and huge, up to (X - 1)^(N - 1) * N^(N - 2) < 2^2930, so it is computed
/// modulo 100 primes just below 1e9 (product > 2^2989) and rebuilt with Garner's algorithm
int main(){
    int n, m;
    long long k;
    scanf("%d %d %lld %lld", &n, &m, &X, &k);
    matrix g(n, vector<long long>(n, 0));
    while (m--){
        int u, v;
        scanf("%d %d", &u, &v);
        u--, v--;
        g[u][v] = g[v][u] = 1 % X;
    }

    matrix walks, pw;
    geometric(g, k + 1, walks, pw);

    /// Walk counts are edge multiplicities of up to ~7.7e6, too many for add_edge, so build the matrix tree minor directly
    vector<vector<long long>> minor(n - 1, vector<long long>(n - 1, 0));
    for (int i = 1; i < n; i++){
        for (int j = 0; j < n; j++){
            if (i == j) continue;
            minor[i - 1][i - 1] += walks[i][j];
            if (j > 0) minor[i - 1][j - 1] -= walks[i][j];
        }
    }

    vector<long long> primes, coef;
    for (int p = 999999999; (int)primes.size() < 100; p -= 2){
        if (!is_prime(p)) continue;
        long long x = determinant_mod(minor, p), prod = 1;
        for (int i = 0; i < (int)primes.size(); i++){
            x = ((x - coef[i] * prod) % p + p) % p;
            prod = prod * primes[i] % p;
        }
        primes.push_back(p);
        coef.push_back(x * power(prod, p - 2, p) % p);
    }

    const long long BASE = 1000000000;
    vector<long long> big(1, 0);
    for (int i = (int)primes.size() - 1; i >= 0; i--){
        long long carry = coef[i];
        for (auto& limb : big){
            long long cur = limb * primes[i] + carry;
            limb = cur % BASE;
            carry = cur / BASE;
        }
        for (; carry; carry /= BASE) big.push_back(carry % BASE);
    }
    while (big.size() > 1 && big.back() == 0) big.pop_back();

    printf("%lld", big.back());
    for (int i = (int)big.size() - 2; i >= 0; i--) printf("%09lld", big[i]);
    puts("");
    return 0;
}
