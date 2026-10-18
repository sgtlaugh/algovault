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

/// Sign and magnitude in base 1e9 limbs, least significant first, zero has no limbs and is never negative
struct BigInt{
    static const int BASE = 1000000000;

    bool negative = false;
    vector<int> limbs;

    BigInt(long long x = 0){
        negative = x < 0;
        unsigned long long mag = negative ? -(unsigned long long)x : x;
        for (; mag; mag /= BASE) limbs.push_back(mag % BASE);
    }

    static int compare_magnitude(const vector<int>& a, const vector<int>& b){
        if (a.size() != b.size()) return a.size() < b.size() ? -1 : 1;
        for (int i = (int)a.size() - 1; i >= 0; i--){
            if (a[i] != b[i]) return a[i] < b[i] ? -1 : 1;
        }
        return 0;
    }

    static void trim(vector<int>& a){
        while (!a.empty() && a.back() == 0) a.pop_back();
    }

    static vector<int> add_magnitude(const vector<int>& a, const vector<int>& b){
        vector<int> res;
        long long carry = 0;
        for (size_t i = 0; i < max(a.size(), b.size()) || carry; i++){
            if (i < a.size()) carry += a[i];
            if (i < b.size()) carry += b[i];
            res.push_back(carry % BASE);
            carry /= BASE;
        }
        return res;
    }

    /// Requires |a| >= |b|
    static vector<int> subtract_magnitude(const vector<int>& a, const vector<int>& b){
        vector<int> res = a;
        long long borrow = 0;
        for (size_t i = 0; i < res.size(); i++){
            long long cur = res[i] - borrow - (i < b.size() ? b[i] : 0);
            borrow = cur < 0;
            res[i] = cur + borrow * BASE;
        }
        trim(res);
        return res;
    }

    /// a -= b k in place, requires a >= b k
    static void subtract_multiple(vector<int>& a, const vector<int>& b, long long k){
        long long borrow = 0;
        for (size_t i = 0; i < a.size(); i++){
            long long cur = a[i] - borrow - (i < b.size() ? b[i] * k : 0);
            borrow = cur < 0 ? (-cur + BASE - 1) / BASE : 0;
            a[i] = cur + borrow * BASE;
        }
        trim(a);
    }

    static BigInt from(bool negative, vector<int> limbs){
        BigInt res;
        trim(limbs);
        res.negative = negative && !limbs.empty();
        res.limbs = std::move(limbs);
        return res;
    }

    BigInt operator - () const{
        return from(!negative, limbs);
    }

    BigInt operator + (const BigInt& o) const{
        if (negative == o.negative) return from(negative, add_magnitude(limbs, o.limbs));
        if (compare_magnitude(limbs, o.limbs) >= 0) return from(negative, subtract_magnitude(limbs, o.limbs));
        return from(o.negative, subtract_magnitude(o.limbs, limbs));
    }

    BigInt operator - (const BigInt& o) const{
        return *this + (-o);
    }

    BigInt operator * (const BigInt& o) const{
        vector<int> res(limbs.size() + o.limbs.size(), 0);
        for (size_t i = 0; i < limbs.size(); i++){
            long long carry = 0;
            for (size_t j = 0; j < o.limbs.size(); j++){
                long long cur = res[i + j] + carry + (long long)limbs[i] * o.limbs[j];
                res[i + j] = cur % BASE;
                carry = cur / BASE;
            }
            res[i + o.limbs.size()] = carry;
        }
        return from(negative != o.negative, std::move(res));
    }

    /// a / BASE^shift from the top three limbs, enough to guess a quotient limb within 1
    static long double leading(const vector<int>& a, int shift){
        long double res = 0;
        for (int i = (int)a.size() - 1; i >= max(0, (int)a.size() - 3); i--) res += a[i] * powl(BASE, i - shift);
        return res;
    }

    /// Truncates toward zero like int, schoolbook long division with each quotient limb guessed then corrected
    static pair<BigInt, BigInt> divmod(const BigInt& a, const BigInt& b){
        if (compare_magnitude(a.limbs, b.limbs) < 0) return {0, a};

        /// Most gcds are small, and dividing by one limb needs no quotient guessing
        if (b.limbs.size() == 1){
            vector<int> quotient(a.limbs.size());
            long long rem = 0;
            for (int i = (int)a.limbs.size() - 1; i >= 0; i--){
                long long cur = rem * BASE + a.limbs[i];
                quotient[i] = cur / b.limbs[0];
                rem = cur % b.limbs[0];
            }
            return {from(a.negative != b.negative, quotient), from(a.negative, vector<int>{(int)rem})};
        }

        /// The top b.size() - 1 limbs of a are below b, so the quotient starts at limb a.size() - b.size()
        int shift = b.limbs.size() - 1, start = a.limbs.size() - b.limbs.size();
        vector<int> quotient(start + 1), rem(a.limbs.begin() + start + 1, a.limbs.end());
        for (int i = start; i >= 0; i--){
            rem.insert(rem.begin(), a.limbs[i]);
            trim(rem);

            /// The guess is at most 1 too high, so one less never overshoots and a few corrections finish it
            long long q = max(0LL, min<long long>(BASE - 1, leading(rem, shift) / leading(b.limbs, shift)) - 1);
            subtract_multiple(rem, b.limbs, q);
            while (compare_magnitude(rem, b.limbs) >= 0){
                q++;
                subtract_multiple(rem, b.limbs, 1);
            }
            quotient[i] = q;
        }
        return {from(a.negative != b.negative, quotient), from(a.negative, rem)};
    }

    BigInt operator / (const BigInt& o) const{
        return divmod(*this, o).first;
    }

    BigInt operator % (const BigInt& o) const{
        return divmod(*this, o).second;
    }

    bool operator < (const BigInt& o) const{
        if (negative != o.negative) return negative;
        int cmp = compare_magnitude(limbs, o.limbs);
        return negative ? cmp > 0 : cmp < 0;
    }

    bool operator == (const BigInt& o) const{
        return negative == o.negative && limbs == o.limbs;
    }

    string str() const{
        if (limbs.empty()) return "0";
        string res = (negative ? "-" : "") + to_string(limbs.back());
        for (int i = (int)limbs.size() - 2; i >= 0; i--){
            string part = to_string(limbs[i]);
            res += string(9 - part.size(), '0') + part;
        }
        return res;
    }
};

BigInt abs(const BigInt& x){
    return x.negative ? -x : x;
}

long long floor_div(long long a, long long b){
    return a / b - ((a % b != 0) && ((a < 0) != (b < 0)));
}

/// Lehmer's gcd replays a run of Euclid quotients found on the top two limbs, plain Euclid took 90% of the run time
BigInt gcd(BigInt a, BigInt b){
    a = abs(a), b = abs(b);
    if (a < b) swap(a, b);

    while (b.limbs.size() >= 2){
        if (b.limbs.size() + 1 < a.limbs.size()){
            a = a % b;
            swap(a, b);
            continue;
        }

        int top = a.limbs.size() - 1;
        long long x = (long long)a.limbs[top] * BigInt::BASE + a.limbs[top - 1];
        long long y = (top < (int)b.limbs.size() ? (long long)b.limbs[top] * BigInt::BASE : 0) + b.limbs[top - 1];

        /// The quotient is exact while both ends of the cofactor interval agree on it
        long long A = 1, B = 0, C = 0, D = 1;
        while (y + C != 0 && y + D != 0){
            long long q = floor_div(x + A, y + C);
            if (q != floor_div(x + B, y + D)) break;
            long long t = A - q * C;
            A = C, C = t;
            t = B - q * D;
            B = D, D = t;
            t = x - q * y;
            x = y, y = t;
        }

        if (B == 0){
            a = a % b;
            swap(a, b);
        }
        else{
            BigInt next_a = a * A + b * B, next_b = a * C + b * D;
            a = abs(next_a), b = abs(next_b);
            if (a < b) swap(a, b);
        }
    }

    if (b == 0) return a;
    long long small = b.limbs[0];
    BigInt rest = a % b;
    return std::gcd(small, rest.limbs.empty() ? 0LL : (long long)rest.limbs[0]);
}

/// Exact rationals so matrix_inverse runs with eps = 0, big integers because answers on dense networks pass 60 digits
struct Fraction{
    BigInt num, den;

    Fraction(long long x = 0) : num(x), den(1) {}

    Fraction(BigInt x, BigInt y) : num(x), den(y){
        if (den.negative) num = -num, den = -den;
        BigInt g = gcd(num, den);
        if (!(g == 1)) num = num / g, den = den / g;
    }

    /// Knuth's form: two gcds on denominators instead of one on the double-length cross product
    Fraction operator - (const Fraction& o) const{
        if (o.num == 0) return *this;
        BigInt g = gcd(den, o.den), diff = num * (o.den / g) - o.num * (den / g);
        if (diff == 0) return Fraction();

        BigInt g2 = gcd(diff, g);
        Fraction res;
        res.num = diff / g2;
        res.den = (den / g) * (o.den / g2);
        return res;
    }

    /// Cross gcds run Euclid on the smaller factors instead of on their products
    Fraction operator * (const Fraction& o) const{
        if (num == 0 || o.num == 0) return Fraction();
        BigInt g1 = gcd(num, o.den), g2 = gcd(o.num, den);
        Fraction res;
        res.num = (num / g1) * (o.num / g2);
        res.den = (den / g2) * (o.den / g1);
        return res;
    }

    Fraction operator / (const Fraction& o) const{
        return Fraction(num * o.den, den * o.num);
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
    Fraction res = f;
    res.num = abs(res.num);
    return res;
}

int main(){
    int cases;
    scanf("%d", &cases);
    for (int tc = 1; tc <= cases; tc++){
        int n, m;
        scanf("%d %d", &n, &m);

        /// Conductances scaled by lcm(1..9) = 2520 stay integers
        const int scale = 2520;
        vector<vector<long long>> conductance(n, vector<long long>(n, 0));
        while (m--){
            int u, v, r;
            scanf("%d %d %d", &u, &v, &r);
            if (u == v) continue;
            conductance[u][v] += scale / r;
            conductance[v][u] += scale / r;
        }

        vector<int> component(n, -1);
        vector<bool> ground(n, false);
        for (int root = 0; root < n; root++){
            if (component[root] != -1) continue;
            ground[root] = true;
            component[root] = root;
            vector<int> stack = {root};
            while (!stack.empty()){
                int u = stack.back();
                stack.pop_back();
                for (int v = 0; v < n; v++){
                    if (conductance[u][v] && component[v] == -1) component[v] = root, stack.push_back(v);
                }
            }
        }

        /// Kirchhoff's current law with one grounded node per component, one inverse serves every query
        vector<vector<Fraction>> laplacian(n, vector<Fraction>(n));
        for (int u = 0; u < n; u++){
            if (ground[u]){
                laplacian[u][u] = 1;
                continue;
            }

            long long total = 0;
            for (int v = 0; v < n; v++){
                total += conductance[u][v];
                if (conductance[u][v] && !ground[v]) laplacian[u][v] = -conductance[u][v];
            }
            laplacian[u][u] = total;
        }
        vector<vector<Fraction>> inv = *matrix_inverse(laplacian, Fraction(0));

        int queries;
        scanf("%d", &queries);
        printf("Case #%d:\n", tc);
        while (queries--){
            int source, sink;
            scanf("%d %d", &source, &sink);
            if (component[source] != component[sink]){
                printf("Resistance between %d and %d is 1/0\n", source, sink);
                continue;
            }

            /// Potential at x for a unit current injected at y, a grounded y absorbs it
            auto potential = [&](int x, int y){
                return ground[y] ? Fraction(0) : inv[x][y];
            };
            Fraction voltage = (potential(source, source) - potential(sink, source)) - (potential(source, sink) - potential(sink, sink));
            voltage = voltage * Fraction(scale);
            printf("Resistance between %d and %d is %s/%s\n", source, sink, voltage.num.str().c_str(), voltage.den.str().c_str());
        }
        puts("");
    }
    return 0;
}
