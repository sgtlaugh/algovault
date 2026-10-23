/***
 *
 * Linear Recurrence
 * Finds the linear recurrence of a sequence from its first terms and computes far terms, modulo any mod
 *
 * Complexity, for a recurrence of degree k:
 *   - finding the recurrence from m terms: O(m^2) for a prime mod, O(m^2 log mod + sqrt(mod)) for a composite one
 *   - nth_term: O(k log k log n)
 *   - nth_terms(n, c): O(k log k log n + k c)
 *
 * LinearRecurrence lr(sequence, mod, recurrence): f(i) = recurrence[0] f(i - 1) + ... + recurrence[k - 1] f(i - k), sequence holds f(0), f(1), ..., at least k terms
 * LinearRecurrence lr(sequence, mod): finds a recurrence generating sequence and stores it in lr.recurrence, in the same order
 * lr.nth_term(n): f(n), 0 <= n <= 1e18
 * lr.nth_terms(n, c): f(n), ..., f(n + c - 1), 1 <= c <= k
 *
 * 1 <= mod and 2 * mod < INT_MAX, values of any sign are reduced modulo mod
 * 2k terms are enough to find a recurrence of degree k
 * With a prime mod the found recurrence is the shortest one generating every given term, from any number of terms
 * With a composite mod fewer than 2k terms may give a recurrence missing given terms, and it may end in zero coefficients
 *
 * Thanks to anta(https://codeforces.com/profile/anta) for first showing me this cool trick and recommending the Black Box Linear Algebra book
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

namespace nt{
    map<int, int> get_factors(int mod){
        long long d, x = mod;
        map<int, int> factors;
        for (d = 2; d * d <= x; d += (d & 1) + 1){
            while (x % d == 0){
                factors[d]++;
                x /= d;
            }
        }
        if (x > 1) factors[x]++;

        return factors;
    }

    int extended_gcd(int a, int b, int& x, int& y){
        if (!b){
            y = 0, x = 1;
            return a;
        }

        int g = extended_gcd(b, a % b, y, x);
        y -= ((a / b) * x);
        return g;
    }

    int mod_inverse(int a, int mod){
        int x, y;
        extended_gcd(a, mod, x, y);
        return (x + mod) % mod;
    }

    int chinese_remainder_theorem(vector<int> ar, vector<int> mods){
        int i, res = 0, mod = 1;
        for (auto m: mods) mod *= m;

        for (i = 0; i < (int)ar.size(); i++){
            long long x = mod / mods[i];
            long long y = mod_inverse(x, mods[i]);
            res = (res + (x * ar[i] % mod) * y) % mod;
        }

        return res;
    }

    int convolution(int* A, int* B, int n, int mod){
        unsigned long long res = 0;
        int i, j, lim = min(ULONG_LONG_MAX / ((unsigned long long)mod * mod), (unsigned long long)n);

        for (i = 0; (i + lim) <= n; res %= mod){
            for (j = 0; j < lim; j++, i++) res += (unsigned long long)A[i] * B[i];
        }

        for (; i < n; i++) res += (unsigned long long)A[i] * B[i];

        return res % mod;
    }
}

namespace fft{
    struct ComplexNum{
        double real, img;

        inline ComplexNum(double real=0.0, double img=0.0) : real(real), img(img) {}

        inline ComplexNum conjugate(){
            return ComplexNum(real, -img);
        }

        inline ComplexNum operator + (ComplexNum other){
            return ComplexNum(real + other.real, img + other.img);
        }

        inline ComplexNum operator - (ComplexNum other){
            return ComplexNum(real - other.real, img - other.img);
        }

        inline ComplexNum operator * (ComplexNum other){
            return ComplexNum((real * other.real) - (img * other.img), (real * other.img) + (img * other.real));
        }
    };

    /// roots[k + j] = e^(i * pi * j / k) for every power of two k, grown on demand and never changed afterwards
    vector<ComplexNum> roots = {ComplexNum(), ComplexNum(1)};

    void extend_roots(int len){
        for (int k = roots.size(); k < len; k <<= 1){
            roots.resize(2 * k);
            for (int j = 0; j < k; j++){
                long double theta = acosl(-1.0L) * j / k;  /// computed directly, a chain of products would accumulate rounding error
                roots[k + j] = ComplexNum(cosl(theta), sinl(theta));
            }
        }
    }

    void transform(vector<ComplexNum>& a){
        int i, j, k, len = a.size();
        for (i = 1, j = 0; i < len; i++){
            for (k = len >> 1; j & k; k >>= 1) j ^= k;
            j ^= k;
            if (i < j) swap(a[i], a[j]);
        }

        for (k = 1; k < len; k <<= 1){
            for (i = 0; i < len; i += (k << 1)){
                for (j = 0; j < k; j++){
                    auto z = a[i + j + k] * roots[j + k];
                    a[i + j + k] = a[i + j] - z;
                    a[i + j] = a[i + j] + z;
                }
            }
        }
    }

    /// Splits the transforms at i and j (= -i) into the low and high 15 bit halves and writes the products of the halves to u and v
    inline void combine(ComplexNum fi, ComplexNum fj, ComplexNum gi, ComplexNum gj, int len, ComplexNum& u, ComplexNum& v){
        auto a1 = (fi + fj.conjugate()) * ComplexNum(0.5, 0);
        auto a2 = (fi - fj.conjugate()) * ComplexNum(0, -0.5);
        auto b1 = (gi + gj.conjugate()) * ComplexNum(0.5 / len, 0);
        auto b2 = (gi - gj.conjugate()) * ComplexNum(0, -0.5 / len);
        v = a1 * b2 + a2 * b1;
        u = a1 * b1 + a2 * b2 * ComplexNum(0, 1);
    }

    /// parity 0 or 1 returns only the coefficients of p1 * p2 at indices of that parity, -1 returns all of them
    vector<int> mod_multiply(const vector<int>& p1, const vector<int>& p2, int mod, int parity = -1){
        int i, j, a = p1.size(), b = p2.size(), p_len = a + b - 1, step = parity < 0 ? 1 : 2;
        while (a > 1 && p1[a - 1] == 0) a--;
        while (b > 1 && p2[b - 1] == 0) b--;
        bool is_equal = (a == b && equal(p1.begin(), p1.begin() + a, p2.begin()));

        int len = 2;
        while (len < a + b) len <<= 1;
        extend_roots(len);

        vector<ComplexNum> f(len), g(len);
        for (i = 0; i < a; i++) f[i] = ComplexNum((p1[i] % mod) & 32767, (p1[i] % mod) >> 15);
        for (i = 0; i < b; i++) g[i] = ComplexNum((p2[i] % mod) & 32767, (p2[i] % mod) >> 15);

        transform(f);
        if (is_equal) g = f;
        else transform(g);

        for (i = 0; i <= len / 2; i++){  /// in place, so index i and its mirror j are rewritten together
            j = (len - i) & (len - 1);
            auto fi = f[i], fj = f[j], gi = g[i], gj = g[j];
            combine(fi, fj, gi, gj, len, f[j], g[j]);
            if (i != j) combine(fj, fi, gj, gi, len, f[i], g[i]);
        }

        if (parity >= 0){  /// folding the halves of the spectrum leaves the transform of one parity at half the length
            len >>= 1;
            for (i = 0; i < len; i++){
                if (parity) f[i] = (f[i] - f[i + len]) * roots[len + i], g[i] = (g[i] - g[i + len]) * roots[len + i];
                else f[i] = f[i] + f[i + len], g[i] = g[i] + g[i + len];
            }
            f.resize(len), g.resize(len);
        }

        transform(f);
        transform(g);

        vector<int> res((p_len - max(parity, 0) + step - 1) / step, 0);
        for (i = 0; i < min(len, (int)res.size()); i++){
            long long x = f[i].real + 0.5, y = g[i].real + 0.5, z = f[i].img + 0.5;
            res[i] = (x + ((y % mod) << 15) + ((z % mod) << 30)) % mod;
        }

        return res;
    }
}

struct LinearRecurrence{
    long long mod;
    vector<int> sequence, recurrence;

    LinearRecurrence() {}

    LinearRecurrence(vector<int> sequence, int mod, vector<int> recurrence={}) : mod(mod), sequence(sequence), recurrence(recurrence){
        int n = sequence.size(), k = recurrence.size();
        assert(n >= k && mod >= 1 && 2LL * mod < INT_MAX);

        for (auto& x : this->sequence) x = ((long long)x % mod + mod) % mod;
        for (auto& c : this->recurrence) c = ((long long)c % mod + mod) % mod;
        if (recurrence.empty() && n > 0 && mod > 1){
            auto C = find_recurrence(this->sequence, mod);  /// C = [1, C1, ..., Cd] with sum C[j] * f(x - j) = 0
            for (int j = 1; j < (int)C.size(); j++) this->recurrence.push_back((mod - C[j]) % mod);
        }
    }

    int nth_term(long long n){
        return nth_terms(n, 1)[0];
    }

    vector<int> nth_terms(long long n, int k){
        int len = recurrence.size();
        if (n + k <= (long long)sequence.size()) return vector<int>(sequence.begin() + n, sequence.begin() + n + k);
        if (mod == 1 || len == 0) return vector<int>(k, 0);
        assert(1 <= k && k <= len);
        if (k == 1) return {bostan_mori(n)};

        auto p = power_of_x(n);

        vector<int> terms = sequence;
        while ((int)terms.size() < len + k - 1){  /// a given recurrence may come with only its first len terms
            long long val = 0;
            int s = terms.size();
            for (int i = 0; i < len; i++) val = (val + (long long)recurrence[i] * terms[s - 1 - i]) % mod;
            terms.push_back(val);
        }

        vector<int> res(k);
        for (int j = 0; j < k; j++){
            long long val = 0;
            for (int i = 0; i < (int)p.size(); i++) val = (val + (long long)p[i] * terms[i + j]) % mod;
            res[j] = val;
        }

        return res;
    }

    vector<int> find_recurrence(vector<int> sequence, int mod){
        int i, j, e, n, p, m, k = 0;
        auto factors = nt::get_factors(mod);
        n = factors.size();

        if (n == 1 && factors.begin()->second == 1){
            return berlekamp_massey(sequence, mod);
        }

        vector<int> v, mods, vals(n);
        vector<vector<int>> coeffs;

        for (auto it: factors){
            tie(p, e) = it;
            m = powl(p, e) + 0.5;
            mods.push_back(m);

            v = reeds_sloane(sequence, p, e);
            k = max(k, (int)v.size());
            coeffs.push_back(v);
        }

        vector<int> res(k);
        for (i = 0; i < k; i++){
            for (j = 0; j < n; j++){
                vals[j] = i < (int)coeffs[j].size() ? coeffs[j][i] : 0;
            }
            res[i] = nt::chinese_remainder_theorem(vals, mods);
        }

        res.resize(sequence.size() / 2 + 1, 0);
        return res;
    }

    /// Bostan-Mori: f(n) = [x^n] P(x) / Q(x), each step multiplies both by Q(-x) and keeps every second coefficient to halve n
    int bostan_mori(long long n){
        int k = recurrence.size();
        vector<int> q(k + 1), v(k + 1);
        q[0] = 1;
        for (int i = 1; i <= k; i++) q[i] = (mod - recurrence[i - 1]) % mod;

        auto p = fft::mod_multiply(vector<int>(sequence.begin(), sequence.begin() + k), q, mod);
        p.resize(k);

        for (; n; n >>= 1){
            for (int i = 0; i <= k; i++) v[i] = (i & 1) && q[i] ? mod - q[i] : q[i];
            p = fft::mod_multiply(p, v, mod, n & 1);
            q = fft::mod_multiply(q, v, mod, 0);
        }

        return p[0];
    }

    /// x^n modulo the characteristic polynomial x^k - c_1 x^(k - 1) - ... - c_k, lowest coefficient first, for n >= 1 and k >= 2
    vector<int> power_of_x(long long n){
        int k = recurrence.size();
        vector<int> charpoly(k + 1, 1);
        for (int j = 0; j < k; j++) charpoly[j] = (mod - recurrence[k - 1 - j]) % mod;
        vector<int> reversed(charpoly.rbegin(), charpoly.rend());

        /// Newton's step inv = inv (2 - reversed inv) doubles the correct prefix of 1 / reversed, k - 1 terms cover every quotient
        vector<int> inv = {1};
        while ((int)inv.size() < k - 1){
            int len = min(2 * (int)inv.size(), k - 1);
            auto t = fft::mod_multiply(inv, vector<int>(reversed.begin(), reversed.begin() + len), mod);
            t.resize(len);
            for (auto& x : t) x = x ? mod - x : 0;
            t[0] = (t[0] + 2) % mod;
            inv = fft::mod_multiply(inv, t, mod);
            inv.resize(len);
        }

        vector<int> res = {1};
        for (int bit = __lg(n); bit >= 0; bit--){
            res = fft::mod_multiply(res, res, mod);

            if ((int)res.size() > k){  /// the top m coefficients reversed, times 1 / reversed charpoly, give the quotient reversed
                int m = res.size() - k;
                auto quot = fft::mod_multiply(vector<int>(res.rbegin(), res.rbegin() + m), vector<int>(inv.begin(), inv.begin() + m), mod);
                quot.resize(m);
                reverse(quot.begin(), quot.end());

                auto sub = fft::mod_multiply(vector<int>(charpoly.begin(), charpoly.begin() + k), quot, mod);
                res.resize(k);
                for (int j = 0; j < k; j++) res[j] = (res[j] - sub[j] + mod) % mod;
            }

            if (n >> bit & 1){
                res.insert(res.begin(), 0);
                if ((int)res.size() > k){  /// x^k reduces to x^k - charpoly
                    long long top = res.back();
                    res.pop_back();
                    for (int j = 0; j < k; j++) res[j] = (res[j] + (mod - charpoly[j]) * top) % mod;
                }
            }
        }

        return res;
    }

    vector<int> berlekamp_massey(vector<int> sequence, long long mod){
        int n = sequence.size();
        vector<int> T, B(n + 1, 0), C(n + 1, 0);
        reverse(sequence.begin(), sequence.end());

        C[0] = 1, B[0] = 1;
        int i, j, d, l = 0, m = 1, deg = 0;
        unsigned long long x, b = 1;

        for (i = 0; i < n; i++){
            d = sequence[n - i - 1];
            if (l > 0) d = (d + nt::convolution(&C[1], &sequence[n - i], l, mod)) % mod;
            if (d == 0) m++;
            else{
                if (l * 2 <= i) T = C, T.resize(l + 1);
                x = b * (mod - d) % mod;

                for (j = 0, m++; j <= deg; j++){
                    C[m + j - 1] = (C[m + j - 1] + x * B[j]) % mod;
                }

                if (l * 2 <= i){
                    B.swap(T);
                    deg = B.size() - 1;
                    b = nt::mod_inverse(d, mod), m = 1, l = i - l + 1;
                }
            }
        }

        C.resize(l + 1);
        return C;
    }

    vector<int> reeds_sloane(vector<int> sequence, int p, int e){
        vector<int> T(e), U(e), R(e), TO(e, 1), UO(e), POW(e + 1, 1);
        int c, d, i, j, k, l, m, len = sequence.size(), mod = powl(p, e) + 0.5;

        vector<vector<int>> A[3], B[3];
        for (i = 0; i < 3; i++){
            A[i] = vector<vector<int>>(e);
            B[i] = vector<vector<int>>(e);
        }

        for (i = 0; i < len; i++) sequence[i] %= mod;
        for (i = 1; i <= e; i++) POW[i] = (long long)POW[i - 1] * p % mod;
        for (i = 0; i < e; i++){
            T[i] = (long long)sequence[0] * POW[i] % mod;
            A[0][i] = {POW[i]}, A[1][i] = {POW[i]};
            B[0][i] = {0}, B[1][i] = {(int)((long long)sequence[0] * POW[i] % mod)};
            for (U[i] = 0; T[i] && T[i] % p == 0; T[i] /= p, U[i]++) {}
            if (!T[i]) T[i] = 1, U[i] = e;
        }

        for (k = 1; k < len; k++){
            for (j = 0; j < e; j++){
                if (compare(A[0][j], B[0][j]) < compare(A[1][j], B[1][j])){
                    R[j] = k - 1;
                    A[2][j] = A[0][e - 1 - U[j]], B[2][j] = B[0][e - 1 - U[j]];
                    UO[j] = U[e - 1 - U[j]], TO[j] = T[e - 1 - U[j]];
                }
            }

            A[0] = A[1], B[0] = B[1];
            for (l = 0; l < e; l++){
                for (d = 0, i = min(k, (int)A[0][l].size() - 1); i >= 0; i--){
                    d = (d + A[0][l][i] * (long long)sequence[k - i]) % mod;
                }

                if (!d) T[l] = 1, U[l] = e;
                else{
                    for (U[l] = 0, T[l] = d; T[l] % p == 0; T[l] /= p, U[l]++) {}
                    j = e - 1 - U[l];

                    if (!compare(A[0][j], B[0][j])){
                        if ((k + 1) > (int)B[1][l].size()) B[1][l].resize(k + 1, 0);
                        B[1][l][k] = (B[1][l][k] + d) % mod;
                    }
                    else{
                        m = k - R[j];
                        c = (long long)T[l] * nt::mod_inverse(TO[j], mod) % mod;
                        c = (long long)c * (long long)POW[U[l] - UO[j]] % mod;
                        update(A[1][l], A[2][j], m, c);
                        update(B[1][l], B[2][j], m, c);
                    }
                }
            }
        }

        return A[1][0];
    }

    long long compare(const vector<int> &a, const vector<int> &b){
        return max(find_rank(a), find_rank(b) + 1);
    }

    void update(vector<int>& u, vector<int>& v, int m, long long c){
        int i, k = v.size();
        if (k + m > (int)u.size()) u.resize(k + m, 0);

        for (i = 0; i < k; i++){
            u[i + m] -= c * v[i] % mod;
            if (u[i + m] < 0) u[i + m] += mod;
        }

        while (u.size() && !u.back()) u.pop_back();
    }

    long long find_rank(const vector<int> & v){
        return (v.size() > 1 || (v.size() == 1 && v[0])) ? v.size() - 1 : (long long)INT_MIN;
    }
};

int main(){
    /// Fibonacci from its recurrence f(i) = f(i - 1) + f(i - 2) and first two terms
    LinearRecurrence fib({0, 1}, 1000000007, {1, 1});
    assert(fib.nth_term(10) == 55);
    assert(fib.nth_terms(9, 2) == vector<int>({34, 55}));

    /// Fibonacci mod 10 repeats every 60 terms and 1e18 % 60 == 40, so f(1e18) ends like f(40) = 102334155
    LinearRecurrence fib10({0, 1}, 10, {1, 1});
    assert(fib10.nth_term(1000000000000000000LL) == 5);

    /// The recurrence found from the first terms
    LinearRecurrence powers({1, 2, 4, 8, 16}, 1000000007);
    assert(powers.recurrence == vector<int>({2}));
    assert(powers.nth_term(30) == 73741817);  /// 2^30 - 1000000007

    LinearRecurrence found({0, 1, 1, 2}, 1000000007);
    assert(found.recurrence == vector<int>({1, 1}));

    /// A composite mod, 55 % 6 == 1
    LinearRecurrence fib6({0, 1, 1, 2}, 6);
    assert(fib6.nth_term(10) == 1);

    return 0;
}
