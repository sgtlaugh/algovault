#include "../common.h"
#include <sys/wait.h>
#include <unistd.h>

#define main library_main
#include "../../code_library/linear_algebra/gauss_bitset.cpp"
#undef main

template <size_t W>
bool satisfies(const vector<bitset<W>>& equations, int m, const bitset<W>& x){
    for (const auto& e : equations){
        int v = 0;
        for (int j = 0; j < m; j++) v ^= e[j] & x[j];
        if (v != e[m]) return false;
    }
    return true;
}

/// Rank of the first m columns by a leading-bit XOR basis, an elimination order independent of the library's
/// Only bits below m are read, so a right-hand side in bit m is ignored and m = W works
template <size_t W>
int rank_of(const vector<bitset<W>>& rows, int m){
    vector<bitset<W>> pivot(m);
    vector<bool> used(m, false);
    int rank = 0;

    for (auto e : rows){
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
template <size_t W>
bool is_identity_product(const vector<bitset<W>>& x, const vector<bitset<W>>& y){
    int n = x.size();
    for (int i = 0; i < n; i++){
        bitset<W> row;
        for (int k = 0; k < n; k++) if (x[i][k]) row ^= y[k];
        bitset<W> unit;
        unit[i] = 1;
        if (row != unit) return false;
    }
    return true;
}

/// A full rank n x n matrix: a unit upper triangular one with shuffled rows mixed by random row additions
template <size_t W>
vector<bitset<W>> random_invertible(int n){
    vector<bitset<W>> a(n);
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

/// A built invertible n x n matrix, then made singular by replacing a row with the sum of two others, which leaves rank exactly n - 1
/// One product direction is enough here: a square a * inv = I implies inv * a = I, and each product costs as much as the inverse
template <size_t W>
void check_invertible(int n){
    vector<bitset<W>> a = random_invertible<W>(n), inv;
    assert(matrix_inverse(n, a, inv) && is_identity_product(a, inv));
    assert(matrix_rank(n, a) == n);

    if (n >= 3){
        int r = stress::rand_int(0, n - 1);
        a[r] = a[(r + 1) % n] ^ a[(r + 2) % n];
        assert(!matrix_inverse(n, a, inv) && matrix_rank(n, a) == n - 1);
    }
}

/// An invertible system in m = W - 1 variables, the right-hand side in the last bit, recovers its planted solution exactly
template <size_t W>
void check_full_width_system(){
    int m = W - 1;
    vector<bitset<W>> equations = random_invertible<W>(m);
    bitset<W> planted;
    for (int j = 0; j < m; j++) planted[j] = stress::rand_int(0, 1);
    for (auto& e : equations) e[m] = (e & planted).count() & 1;

    bitset<W> res;
    assert(gauss(m, equations, res) == 0 && res == planted);
}

/// Reduced row echelon form: pivots at rows 0, 1, ..., rank - 1 in column order, each pivot column a unit column,
/// no bit left of a pivot, zero rows below the pivot rows, and the same row space as a
template <size_t W>
void check_reduced(const vector<bitset<W>>& a, int m, int rank){
    vector<bitset<W>> r = a;
    vector<int> pos = eliminate_gf2(r, m);
    int n = a.size(), pivots = 0;
    assert((int)pos.size() == m);

    for (int j = 0; j < m; j++){
        if (pos[j] == -1) continue;
        assert(pos[j] == pivots++);
        for (int k = 0; k < n; k++) assert(r[k][j] == (k == pos[j]));
        for (int l = 0; l < j; l++) assert(!r[pos[j]][l]);
    }
    assert(pivots == rank);
    for (int k = rank; k < n; k++) for (int j = 0; j < m; j++) assert(!r[k][j]);

    vector<bitset<W>> both = a;
    both.insert(both.end(), r.begin(), r.end());
    assert(rank_of(both, m) == rank);
}

/// matrix_inverse must abort on n > W before it writes the identity past the end of the rows
/// A single-word bitset maps every index to its one word, so only W > 64 lets the late check be seen as an overflow
void check_inverse_rejects_too_wide(){
    pid_t pid = fork();
    assert(pid >= 0);
    if (pid == 0){
        assert(freopen("/dev/null", "w", stderr));
        vector<bitset<128>> a(400), inv;
        matrix_inverse(400, a, inv);
        _exit(0);
    }

    int status;
    assert(waitpid(pid, &status, 0) == pid);
    assert(WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT);
}

int main(){
    check_inverse_rejects_too_wide();

    /// Small matrices at width W = 10, so m = W and n = W are reached: rank from the kernel size 2^(m - rank), invertible exactly when the kernel is trivial
    const size_t SMALL = 10;
    for (long long it = 0; it < stress::scaled(3000); it++){
        int m = stress::rand_int(0, SMALL), n = stress::rand_int(0, SMALL), density = stress::rand_int(1, 3);
        vector<bitset<SMALL>> a(n);
        for (auto& row : a) for (int j = 0; j < m; j++) row[j] = stress::rand_int(0, density) == 0;
        if (n >= 2 && stress::rand_int(0, 1)) a[stress::rand_int(0, n - 1)] = a[0] ^ a[n - 1];

        long long kernel = 0;
        for (int mask = 0; mask < (1 << m); mask++){
            bitset<SMALL> x(mask);
            bool zero = true;
            for (const auto& row : a) zero &= (row & x).count() % 2 == 0;
            kernel += zero;
        }
        int rank = matrix_rank(m, a);
        assert(rank >= 0 && kernel == 1LL << (m - rank));
        check_reduced(a, m, rank);

        if (n != m) continue;
        vector<bitset<SMALL>> inv(1);
        bool ok = matrix_inverse(n, a, inv);
        assert(ok == (kernel == 1));
        if (ok){
            assert((int)inv.size() == n && is_identity_product(a, inv) && is_identity_product(inv, a));
            for (const auto& row : inv) assert((row >> n).none());
        }
    }

    /// Widths above 1024: random sizes at W = 2000 (not a multiple of 64), then the full default width n = MAX
    for (long long it = 0; it < stress::scaled(8); it++) check_invertible<2000>(it == 0 ? 2000 : stress::rand_int(1, 2000));
    check_invertible<MAX>(MAX);

    /// Small systems at width W = 13, so m = W - 1 is reached, checked against all 2^m assignments
    const size_t NARROW = 13;
    for (long long it = 0; it < stress::scaled(3000); it++){
        int m = stress::rand_int(0, NARROW - 1), n = stress::rand_int(0, 14), density = stress::rand_int(1, 3);
        vector<bitset<NARROW>> equations(n);
        for (auto& e : equations){
            for (int j = 0; j <= m; j++) e[j] = stress::rand_int(0, density) == 0;
        }
        /// Duplicated and combined rows make dependent systems common
        if (n >= 2 && stress::rand_int(0, 1)) equations[stress::rand_int(0, n - 1)] = equations[0] ^ equations[n - 1];

        long long solutions = 0;
        for (int mask = 0; mask < (1 << m); mask++) solutions += satisfies(equations, m, bitset<NARROW>(mask));

        bitset<NARROW> res;
        res.set();  /// stale bits must be cleared by gauss
        int f_var = gauss(m, equations, res);
        if (!solutions) assert(f_var == -1);
        else{
            assert(f_var >= 0 && solutions == 1LL << f_var);
            assert(satisfies(equations, m, res));
            assert((res >> m).none());
        }
    }

    /// Large systems at W = 2000: a planted solution must be recovered as a valid one with exactly m - rank free variables
    const size_t WIDE = 2000;
    for (long long it = 0; it < stress::scaled(30); it++){
        int m = stress::rand_int(1, WIDE - 1), n = stress::rand_int(1, 300), rank = stress::rand_int(0, min(n, m));
        bitset<WIDE> planted;
        for (int j = 0; j < m; j++) planted[j] = stress::rand_int(0, 1);

        vector<bitset<WIDE>> basis(rank), equations(n);
        for (int r = 0; r < rank; r++){
            basis[r][stress::rand_int(0, m - 1)] = 1;  /// may still be dependent, the true rank is computed below
            for (int j = 0; j < m; j++) if (stress::rand_int(0, 1)) basis[r][j] = 1;
        }
        for (auto& e : equations){
            for (int r = 0; r < rank; r++) if (stress::rand_int(0, 1)) e ^= basis[r];
            e[m] = (e & planted).count() & 1;
        }

        bitset<WIDE> res;
        int f_var = gauss(m, equations, res), true_rank = rank_of(equations, m);
        assert(f_var == m - true_rank && satisfies(equations, m, res));
        check_reduced(equations, m, true_rank);
    }

    /// The m = W - 1 boundary above 1024 and at the default width
    for (long long it = 0; it < stress::scaled(2); it++) check_full_width_system<WIDE>();
    check_full_width_system<MAX>();

    return 0;
}
