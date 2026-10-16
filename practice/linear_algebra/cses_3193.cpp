/***
 *
 * Gauss-Jordan Elimination over GF(2)
 * Linear systems, rank and inverse of matrices over GF(2), with bitset rows
 *
 * Complexity: O(min(n, m) * n * W / 64) for gauss and matrix_rank, O(n^2 * W / 64) for matrix_inverse
 * Rows are bitset<W>, W is deduced from the arguments and defaults to GF2_WIDTH = 4096, several widths can be used in one program
 * Every row operation costs W / 64 whatever m is, so pick the smallest W that fits: bitset<128> rows for 100 columns XOR 2 words per row operation instead of 64
 *
 * gauss(m, equations, res): equations in m variables, bit j of equations[i] is the coefficient of x_j, bit m the right-hand side
 *     needs m < W, returns -1 when there is no solution
 *     otherwise returns the number of free variables, 0 for a unique solution, and res holds one solution
 * matrix_rank(m, a): rank of the n x m matrix in bits 0 to m - 1 of the rows of a, needs m <= W
 * matrix_inverse(n, a, inv): for the n x n matrix in a, needs n <= W
 *     returns false when it is singular, otherwise inv = a^-1
 * eliminate_gf2(a, cols, side): reduces the first cols columns of a to reduced row echelon form in place, needs cols <= W
 *     returns pos, pos[j] = row of the pivot in column j or -1, pivot rows are 0, 1, ..., rank - 1 in column order
 *     every row operation is applied to *side as well when side is not null
 * rows_of<W>(strings): rows from strings of '0' and '1', character j is bit j
 *
 * For instance, the system of linear equations modulo 2
 *
 *     (2x + y - z)   % 2 = 0  ----->  equations[0] = bits {0, 1, 1, 0}
 *     (-3x - y + 2z) % 2 = 1  ----->  equations[1] = bits {1, 1, 0, 1}
 *     (-2x + y + 2z) % 2 = 1  ----->  equations[2] = bits {0, 1, 0, 1}
 *
 * has the unique solution x = 0, y = 1, z = 1
 *
 * To make it run even faster in practice, use custom bitset class with unsigned long long
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

const int GF2_WIDTH = 4096;

template <size_t W>
vector<int> eliminate_gf2(vector<bitset<W>>& a, int cols, vector<bitset<W>>* side = nullptr){
    assert(cols <= (int)W);
    int n = a.size(), i = 0;
    vector<int> pos(cols, -1);

    for (int j = 0; j < cols && i < n; j++){
        int p = i;
        while (p < n && !a[p][j]) p++;
        if (p == n) continue;

        swap(a[p], a[i]);
        if (side) swap((*side)[p], (*side)[i]);
        for (int k = 0; k < n; k++){
            if (k == i || !a[k][j]) continue;
            a[k] ^= a[i];
            if (side) (*side)[k] ^= (*side)[i];
        }

        pos[j] = i++;
    }

    return pos;
}

template <size_t W = GF2_WIDTH>
int gauss(int m, vector<bitset<W>> equations, bitset<W>& res){
    assert(m < (int)W);
    int n = equations.size(), f_var = 0;
    vector<int> pos = eliminate_gf2(equations, m);

    res.reset();
    int rank = m - count(pos.begin(), pos.end(), -1);

    /// Only rows without a pivot can be inconsistent, as 0 = 1
    for (int k = rank; k < n; k++){
        if (equations[k][m]) return -1;
    }

    for (int j = 0; j < m; j++){
        if (pos[j] == -1) f_var++;
        else res[j] = equations[pos[j]][m];
    }

    return f_var;
}

template <size_t W = GF2_WIDTH>
bool matrix_inverse(int n, vector<bitset<W>> a, vector<bitset<W>>& inv){
    /// Checked before the identity is written, side[i][i] with i >= W is past the end of the row
    assert((int)a.size() == n && n <= (int)W);
    vector<bitset<W>> side(n);
    for (int i = 0; i < n; i++) side[i][i] = 1;

    vector<int> pos = eliminate_gf2(a, n, &side);
    if (n && pos[n - 1] != n - 1) return false;

    inv = side;
    return true;
}

/// W is not deduced from an empty braced list, matrix_rank(0, {}), so it falls back to GF2_WIDTH
template <size_t W = GF2_WIDTH>
int matrix_rank(int m, vector<bitset<W>> a){
    vector<int> pos = eliminate_gf2(a, m);
    return m - count(pos.begin(), pos.end(), -1);
}

template <size_t W = GF2_WIDTH>
vector<bitset<W>> rows_of(const vector<string>& bits){
    vector<bitset<W>> res;
    for (const string& s : bits){
        bitset<W> row;
        for (int j = 0; j < (int)s.size(); j++) row[j] = s[j] == '1';
        res.push_back(row);
    }
    return res;
}

int main(){
    const int MAXV = 5000, MOD = 1000000007;
    vector<int> primes;
    for (int p = 2; p <= MAXV; p++){
        bool prime = true;
        for (int q = 2; q * q <= p && prime; q++) prime = p % q != 0;
        if (prime) primes.push_back(p);
    }

    int n;
    scanf("%d", &n);

    /// A product is a square iff its prime exponent parity vector is zero, so the answer is 2^(n - rank)
    vector<bitset<672>> rows(n);
    for (int i = 0; i < n; i++){
        int x;
        scanf("%d", &x);
        for (int j = 0; j < (int)primes.size(); j++){
            while (x % primes[j] == 0){
                x /= primes[j];
                rows[i].flip(j);
            }
        }
    }

    int free_count = n - matrix_rank(primes.size(), rows);
    long long ans = 1;
    for (int i = 0; i < free_count; i++) ans = ans * 2 % MOD;
    printf("%lld\n", ans);
    return 0;
}
