#include "../common.h"

#define main library_main
#include "../../code_library/linear_algebra/gauss_bitset.cpp"
#undef main

bool satisfies(const vector<bitset<MAX>>& equations, int m, const bitset<MAX>& x){
    for (const auto& e : equations){
        int v = 0;
        for (int j = 0; j < m; j++) v ^= e[j] & x[j];
        if (v != e[m]) return false;
    }
    return true;
}

/// Rank of the coefficient columns by a leading-bit XOR basis, an elimination order independent of the library's
int rank_of(const vector<bitset<MAX>>& equations, int m){
    vector<bitset<MAX>> pivot(m);
    vector<bool> used(m, false);
    int rank = 0;

    for (auto e : equations){
        e[m] = 0;
        for (int j = m - 1; j >= 0; j--){
            if (!e[j]) continue;
            if (!used[j]){
                pivot[j] = e, used[j] = true, rank++;
                break;
            }
            e ^= pivot[j];
        }
    }

    return rank;
}

/// Row i of x * y is the XOR of the rows of y selected by row i of x
bool is_identity_product(const vector<bitset<MAX>>& x, const vector<bitset<MAX>>& y){
    int n = x.size();
    for (int i = 0; i < n; i++){
        bitset<MAX> row;
        for (int k = 0; k < n; k++) if (x[i][k]) row ^= y[k];
        bitset<MAX> unit;
        unit[i] = 1;
        if (row != unit) return false;
    }
    return true;
}

/// A full rank n x n matrix: a unit upper triangular one with shuffled rows mixed by random row additions
vector<bitset<MAX>> random_invertible(int n){
    vector<bitset<MAX>> a(n);
    for (int i = 0; i < n; i++){
        a[i][i] = 1;
        for (int j = i + 1; j < n; j++) a[i][j] = stress::rand_int(0, 1);
    }
    shuffle(a.begin(), a.end(), stress::rng());
    for (int t = 0; n >= 2 && t < 3 * n; t++){
        int k = stress::rand_int(0, n - 1), l = stress::rand_int(0, n - 2);
        a[k] ^= a[l + (l >= k)];
    }
    return a;
}

int main(){
    /// Small matrices: rank from the kernel size 2^(m - rank), invertible exactly when the kernel is trivial
    for (long long it = 0; it < stress::scaled(3000); it++){
        int m = stress::rand_int(0, 10), n = stress::rand_int(0, 10), density = stress::rand_int(1, 3);
        vector<bitset<MAX>> a(n);
        for (auto& row : a) for (int j = 0; j < m; j++) row[j] = stress::rand_int(0, density) == 0;
        if (n >= 2 && stress::rand_int(0, 1)) a[stress::rand_int(0, n - 1)] = a[0] ^ a[n - 1];

        long long kernel = 0;
        for (int mask = 0; mask < (1 << m); mask++){
            bitset<MAX> x(mask);
            bool zero = true;
            for (const auto& row : a) zero &= (row & x).count() % 2 == 0;
            kernel += zero;
        }
        int rank = matrix_rank(m, a);
        assert(rank >= 0 && kernel == 1LL << (m - rank));

        if (n != m) continue;
        vector<bitset<MAX>> inv(1);
        bool ok = matrix_inverse(n, a, inv);
        assert(ok == (kernel == 1));
        if (ok){
            assert((int)inv.size() == n && is_identity_product(a, inv) && is_identity_product(inv, a));
            for (const auto& row : inv) assert((row >> n).none());
        }
    }

    /// Up to the MAX boundary: built invertible matrices, then made singular by replacing a row with a sum of two others
    for (long long it = 0; it < stress::scaled(12); it++){
        int n = it < 2 ? MAX : stress::rand_int(1, MAX);
        vector<bitset<MAX>> a = random_invertible(n), inv;
        assert(matrix_inverse(n, a, inv) && is_identity_product(a, inv) && is_identity_product(inv, a));
        assert(matrix_rank(n, a) == n);

        if (n >= 3){
            int r = stress::rand_int(0, n - 1);
            a[r] = a[(r + 1) % n] ^ a[(r + 2) % n];
            assert(!matrix_inverse(n, a, inv) && matrix_rank(n, a) == n - 1);
            if (n < MAX) assert(rank_of(a, n) == n - 1);
        }
    }

    for (long long it = 0; it < stress::scaled(3000); it++){
        int m = stress::rand_int(0, 12), n = stress::rand_int(0, 14), density = stress::rand_int(1, 3);
        vector<bitset<MAX>> equations(n);
        for (auto& e : equations){
            for (int j = 0; j <= m; j++) e[j] = stress::rand_int(0, density) == 0;
        }
        /// Duplicated and combined rows make dependent systems common
        if (n >= 2 && stress::rand_int(0, 1)) equations[stress::rand_int(0, n - 1)] = equations[0] ^ equations[n - 1];

        long long solutions = 0;
        for (int mask = 0; mask < (1 << m); mask++) solutions += satisfies(equations, m, bitset<MAX>(mask));

        bitset<MAX> res;
        res.set();  /// stale bits must be cleared by gauss
        int f_var = gauss(m, equations, res);
        if (!solutions) assert(f_var == -1);
        else{
            assert(f_var >= 0 && solutions == 1LL << f_var);
            assert(satisfies(equations, m, res));
            assert((res >> m).none());
        }
    }

    /// Large systems: a planted solution must be recovered as a valid one with exactly m - rank free variables
    for (long long it = 0; it < stress::scaled(30); it++){
        int m = stress::rand_int(1, MAX - 1), n = stress::rand_int(1, 300), rank = stress::rand_int(0, min(n, m));
        bitset<MAX> planted;
        for (int j = 0; j < m; j++) planted[j] = stress::rand_int(0, 1);

        vector<bitset<MAX>> basis(rank), equations(n);
        for (int r = 0; r < rank; r++){
            basis[r][stress::rand_int(0, m - 1)] = 1;  /// may still be dependent, the true rank is computed below
            for (int j = 0; j < m; j++) if (stress::rand_int(0, 1)) basis[r][j] = 1;
        }
        for (auto& e : equations){
            for (int r = 0; r < rank; r++) if (stress::rand_int(0, 1)) e ^= basis[r];
            e[m] = (e & planted).count() & 1;
        }

        bitset<MAX> res;
        int f_var = gauss(m, equations, res);
        assert(f_var == m - rank_of(equations, m) && satisfies(equations, m, res));
    }

    /// The m = MAX - 1 boundary, the right-hand side in the last bit: an invertible system recovers its planted solution exactly
    for (long long it = 0; it < stress::scaled(3); it++){
        int m = MAX - 1;
        vector<bitset<MAX>> equations = random_invertible(m);
        bitset<MAX> planted;
        for (int j = 0; j < m; j++) planted[j] = stress::rand_int(0, 1);
        for (auto& e : equations) e[m] = (e & planted).count() & 1;

        bitset<MAX> res;
        assert(gauss(m, equations, res) == 0 && res == planted);
    }

    return 0;
}
