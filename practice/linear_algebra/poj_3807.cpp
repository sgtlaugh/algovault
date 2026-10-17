/***
 *
 * Gauss-Jordan Elimination
 * Linear systems, rank and inverse of matrices over the reals
 *
 * Complexity: O(min(n, m) n m) for gauss and matrix_rank, O(n^3) for matrix_inverse
 *
 * gauss(equations, eps = 1e-9): {solvable, x, free_variables} for n >= 1 equations in m variables, equations[i][m] is the right-hand side
 * matrix_rank(a, eps = 1e-9): rank of any n x m matrix
 * matrix_inverse(a, eps = 1e-9): a^-1 for a square a, nullopt when a is singular
 *
 * x is one solution, empty when solvable is false, and free_variables is m minus the rank either way
 * eps is relative, scaled by max(1, largest |coefficient| of the input), a pivot within it counts as 0
 * The right-hand side of gauss is left out of that scale, so x = 1e9 with coefficients of 1 still solves
 * A row of gauss that reduces to 0 = r is consistent when |r| is within eps scaled by the largest |entry| including it
 * T may be an exact type such as a fraction with eps = 0, and long double or __float128 give more precision
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <class T>
struct GaussSolution{
    bool solvable;
    vector<T> x;
    int free_variables;
};

/// Relative, an absolute eps rejects solvable systems with large coefficients
/// Over the first cols columns, a large right-hand side in the scale would push ordinary pivots under it
template <class T>
T tolerance(const vector<vector<T>>& a, int cols, T eps){
    T tol = 1;
    for (auto& row : a) for (int j = 0; j < cols; j++) tol = max(tol, (T)abs(row[j]));
    return tol * eps;
}

/// Brings the first cols columns of a to row echelon form with partial pivoting, pivots are not scaled to 1
/// Returns pos, pos[j] = row of the pivot in column j or -1, pivot rows are 0, 1, ..., rank - 1 in column order
/// Only rows below a pivot are cleared, about 2/3 of the work of clearing above too
template <class T>
vector<int> row_reduce(vector<vector<T>>& a, int cols, T tol){
    int n = a.size(), i = 0;
    vector<int> pos(cols, -1);

    for (int j = 0; j < cols && i < n; j++){
        int p = i;
        for (int k = i; k < n; k++){
            if (abs(a[k][j]) > abs(a[p][j])) p = k;
        }
        if (abs(a[p][j]) <= tol) continue;

        swap(a[p], a[i]);
        int width = a[i].size();
        for (int k = i + 1; k < n; k++){
            if (a[k][j] == 0) continue;
            T factor = a[k][j] / a[i][j];
            for (int l = j; l < width; l++) a[k][l] -= a[i][l] * factor;
        }

        pos[j] = i++;
    }

    return pos;
}

template <class T>
GaussSolution<T> gauss(vector<vector<T>> equations, const T eps=1e-9){
    int n = equations.size(), m = equations[0].size() - 1;
    T tol = tolerance(equations, m + 1, eps);
    vector<int> pos = row_reduce(equations, m, tolerance(equations, m, eps));
    int rank = m - count(pos.begin(), pos.end(), -1);

    /// Only rows without a pivot can be inconsistent, as 0 = nonzero
    for (int k = rank; k < n; k++){
        if (abs(equations[k][m]) > tol) return {false, {}, m - rank};
    }

    /// Back substitution from the last pivot up, free variables stay 0
    vector<T> x(m, 0);
    for (int j = m - 1; j >= 0; j--){
        if (pos[j] == -1) continue;

        vector<T>& row = equations[pos[j]];
        T value = row[m];
        for (int l = j + 1; l < m; l++) value -= row[l] * x[l];
        x[j] = value / row[j];
    }

    return {true, x, m - rank};
}

template <class T>
int matrix_rank(vector<vector<T>> a, const T eps=1e-9){
    int cols = a.empty() ? 0 : a[0].size();
    vector<int> pos = row_reduce(a, cols, tolerance(a, cols, eps));
    return cols - count(pos.begin(), pos.end(), -1);
}

template <class T>
optional<vector<vector<T>>> matrix_inverse(const vector<vector<T>>& a, const T eps=1e-9){
    int n = a.size();
    vector<vector<T>> aug(n, vector<T>(2 * n, 0));
    for (int i = 0; i < n; i++){
        assert((int)a[i].size() == n);
        for (int j = 0; j < n; j++) aug[i][j] = a[i][j];
        aug[i][n + i] = 1;
    }

    vector<int> pos = row_reduce(aug, n, tolerance(a, n, eps));
    if (n && pos[n - 1] != n - 1) return nullopt;

    /// Clears above each pivot from the last up, the left half of a pivot row is then 0 past its pivot
    for (int i = n - 1; i >= 0; i--){
        for (int k = 0; k < i; k++){
            if (aug[k][i] == 0) continue;
            T factor = aug[k][i] / aug[i][i];
            for (int l = n; l < 2 * n; l++) aug[k][l] -= aug[i][l] * factor;
        }
    }

    vector<vector<T>> inv(n, vector<T>(n));
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++) inv[i][j] = aug[i][n + j] / aug[i][i];
    }

    return inv;
}

/// Exact rationals so gauss runs with eps = 0, 128-bit parts leave headroom for the products before reducing
struct Fraction{
    __int128 num, den;

    Fraction(long long x = 0) : num(x), den(1) {}

    Fraction(__int128 x, __int128 y) : num(x), den(y){
        if (den < 0) num = -num, den = -den;
        __int128 g = gcd128(num < 0 ? -num : num, den);
        if (g > 1) num /= g, den /= g;
    }

    static __int128 gcd128(__int128 a, __int128 b){
        while (b) a %= b, swap(a, b);
        return a;
    }

    Fraction operator - (const Fraction& o) const{
        __int128 g = gcd128(den, o.den);
        return Fraction(num * (o.den / g) - o.num * (den / g), den / g * o.den);
    }

    Fraction operator * (const Fraction& o) const{
        __int128 g1 = gcd128(num < 0 ? -num : num, o.den), g2 = gcd128(o.num < 0 ? -o.num : o.num, den);
        return Fraction((num / g1) * (o.num / g2), (den / g2) * (o.den / g1));
    }

    Fraction operator / (const Fraction& o) const{
        return *this * Fraction(o.den, o.num);
    }

    Fraction& operator -= (const Fraction& o){
        return *this = *this - o;
    }

    bool operator < (const Fraction& o) const{
        return num * o.den < o.num * den;
    }

    bool operator > (const Fraction& o) const{
        return o < *this;
    }

    bool operator <= (const Fraction& o) const{
        return !(o < *this);
    }

    bool operator == (const Fraction& o) const{
        return num == o.num && den == o.den;
    }
};

Fraction abs(const Fraction& f){
    return Fraction(f.num < 0 ? -f.num : f.num, f.den);
}

/// Adds mult times the atoms of molecule s to counts, parentheses are never nested
void add_atoms(const string& s, int mult, map<string, int>& counts){
    for (size_t i = 0; i < s.size();){
        size_t end = i + 1;
        if (s[i] == '(') end = s.find(')', i) + 1;
        else if (end < s.size() && islower(s[end])) end++;

        int k = 0;
        size_t digits = end;
        while (digits < s.size() && isdigit(s[digits])) k = k * 10 + (s[digits++] - '0');
        if (digits == end) k = 1;

        if (s[i] == '(') add_atoms(s.substr(i + 1, end - i - 2), mult * k, counts);
        else counts[s.substr(i, end - i)] += mult * k;
        i = digits;
    }
}

vector<string> split(const string& s, char separator){
    vector<string> parts(1);
    for (char c: s){
        if (c == separator) parts.emplace_back();
        else parts.back() += c;
    }
    return parts;
}

int main(){
    string line;
    while (cin >> line && line != "."){
        line.pop_back();
        size_t arrow = line.find("->");
        vector<string> molecules = split(line.substr(0, arrow), '+');
        int left = molecules.size();
        for (auto& s: split(line.substr(arrow + 2), '+')) molecules.push_back(s);
        int m = molecules.size();

        /// Products count negative so every element balances to 0
        map<string, vector<long long>> element_rows;
        for (int j = 0; j < m; j++){
            map<string, int> counts;
            add_atoms(molecules[j], j < left ? 1 : -1, counts);
            for (auto& [element, count]: counts){
                auto& row = element_rows[element];
                row.resize(m + 1, 0);
                row[j] = count;
            }
        }

        /// The balanced coefficients are one ray, so pinning the first to 1 leaves a unique rational solution
        vector<vector<Fraction>> equations;
        for (auto& [element, row]: element_rows) equations.emplace_back(row.begin(), row.end());
        equations.emplace_back(m + 1, Fraction(0));
        equations.back()[0] = equations.back()[m] = 1;

        vector<Fraction> x = gauss(equations, Fraction(0)).x;
        __int128 lcm = 1, g = 0;
        for (auto& f: x) lcm = lcm / Fraction::gcd128(lcm, f.den) * f.den;
        vector<long long> coefficients;
        for (auto& f: x){
            coefficients.push_back((long long)(f.num * (lcm / f.den)));
            g = Fraction::gcd128(coefficients.back(), g);
        }
        for (int j = 0; j < m; j++) printf("%lld%c", coefficients[j] / (long long)g, j + 1 < m ? ' ' : '\n');
    }
    return 0;
}
