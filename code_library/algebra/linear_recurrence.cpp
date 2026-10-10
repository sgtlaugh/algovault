/***
 * Generalized linear recurrence solver for recurrence with high degrees
 * Algorithm uses a combination of Berlekamp-Massey and Reeds-Sloane to derive the recurrence
 * Then it evaluates the n'th term with Bostan-Mori, and a window of terms with x^n mod the characteristic polynomial, using fast fourier transform
 *
 * Overall complexity is roughly O(k^2) + O(k * log k * log n) with high constant factors
 * O(k^2) is for initialization only and can be skipped if the recurrence is known
 * O(k * log k * log n) accounts for each query to calculate the n'th term, nth_terms(n, c) adds O(k * c)
 * Any modulus with 2 * mod < INT_MAX, FFT lengths grow as needed: nth_term with k = 1e5 and n = 1e18 takes about 3 s at -O2
 *
 * Thanks to anta(https://codeforces.com/profile/anta) for first showing me this cool trick and recommending the Black Box Linear Algebra book
 *
***/

#include <stdio.h>
#include <bits/stdtr1c++.h>

using namespace std;

namespace nt{
    map<int, map<int, int>> factors_map;

    map<int, int> get_factors(int mod){
        if (factors_map.count(mod)) return factors_map[mod];

        long long d, x = mod;
        map<int, int> factors;
        for (d = 2; d * d <= x; d += (d & 1) + 1){
            while (x % d == 0){
                factors[d]++;
                x /= d;
            }
        }
        if (x > 1) factors[x]++;

        return factors_map[mod] = factors;
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

        for (j = 0; i < n; i++) res += (unsigned long long)A[i] * B[i];

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

int _POLYNOMIAL_MOD;

struct Polynomial{
    vector<int> coefficient;

    Polynomial(vector<int> coefficient={}) : coefficient(coefficient) {}

    inline int* data(){
        return coefficient.empty() ? 0 : &coefficient[0];
    }

    inline void copy_poly(int* to, int* from, int n){
        if (n > 0) memcpy(to, from, n * sizeof(int));  /// an empty polynomial's data is null, memcpy from null is UB even for 0 bytes
    }

    inline void subtract(int* res, int* P, int pn, int* Q, int qn){
        for (int i = 0; i < qn; i++){
            res[i] = P[i] - Q[i];
            if (res[i] < 0) res[i] += _POLYNOMIAL_MOD;
        }

        copy_poly(res + qn, P + qn, pn - qn);
    }

    inline void shift(int* res, int* P, int n, int k){
        for (int i = 0; i < n; i++){
            res[i] = ((long long)P[i] << k) % _POLYNOMIAL_MOD;
        }
    }

    inline void reverse_poly(int* res, int* P, int n){
        if (&res[0] == &P[0]) return reverse(res, res + n);
        for (int i = 0; i < n; i++) res[n - i - 1] = P[i];
    }

    inline int size(){
        return coefficient.size();
    }

    inline void normalize(){
        while (size() && !coefficient.back()) coefficient.pop_back();
    }

    void multiply(Polynomial &res, Polynomial& p, Polynomial& q){
        if (&res == &p || &res == &q){
            Polynomial temp;
            multiply(temp, p, q);
            res = temp;
            return;
        }

        res.coefficient.clear();
        if (p.size() && q.size()){
            res.coefficient = fft::mod_multiply(p.coefficient, q.coefficient, _POLYNOMIAL_MOD);
        }
    }

    Polynomial inverse(int n){
        Polynomial res({n});
        res.coefficient.resize(n);
        find_inverse(res.data(), n, data(), size());
        return res;
    }

    void divide(Polynomial &quot, Polynomial &rem, Polynomial &p, Polynomial &q, Polynomial &inv){
        int pn = p.size(), qn = q.size();
        quot.coefficient.resize(max((int)0, pn - qn + 1)), rem.coefficient.resize(qn - 1);
        divide_remainder_inverse(quot.data(), rem.data(), p.data(), pn, q.data(), qn, inv.data());
        quot.normalize(), rem.normalize();
    }

    Polynomial remainder(Polynomial &q, Polynomial &inv){
        Polynomial quot, rem;
        divide(quot, rem, *this, q, inv);
        return rem;
    }

    Polynomial power(Polynomial &q, long long k){
        int qn = q.size();
        if (qn == 1) return Polynomial();
        if (k == 0) return Polynomial({1});

        Polynomial inv = q.inverse(max(size() - qn + 1, qn));
        Polynomial p = this->remainder(q, inv);

        Polynomial res = p;
        for (int l = 62 - __builtin_clzll(k); l >= 0; l--){
            multiply(res, res, res);
            res = res.remainder(q, inv);

            if (k & (1LL << l)){
                multiply(res, res, p);
                res = res.remainder(q, inv);
            }
        }

        return res;
    }

    void multiply(int* res, int* P, int pn, int* Q, int qn){
        auto P1 = Polynomial(vector<int>(P, P + pn));
        auto P2 = Polynomial(vector<int>(Q, Q + qn));

        Polynomial P_res;
        P_res.multiply(P_res, P1, P2);
        for (int i = 0; i < (int)P_res.coefficient.size(); i++){
            res[i] = P_res.coefficient[i];
        }
    }

    void inverse_power_series(int* res, int res_n, int* P, int pn){
        if (res_n == 0) return;
        unique_ptr <int[]> ptr(new int[4 * res_n]);
        int* u = ptr.get(), *v = u + res_n * 2, cur = 1, nxt = 1;

        for (int i = 0; i < res_n; i++) res[i] = 0;
        res[0] = P[0];

        while (cur < res_n){
            nxt = min(res_n, cur * 2);
            multiply(u, res, cur, res, cur);
            multiply(v, u, min(nxt, cur * 2 - 1), P, min(nxt, pn));
            shift(res, res, cur, 1);
            subtract(res, res, nxt, v, nxt);
            cur = nxt;
        }
    }

    void find_inverse(int* res, int res_n, int* P, int pn){
        unique_ptr <int[]> ptr(new int[pn]);
        int* tmp = ptr.get();
        reverse_poly(tmp, P, pn);
        inverse_power_series(res, res_n, tmp, pn);
    }

    void divide_inverse(int* res, int res_n, int* revp, int pn, int* inv){
        unique_ptr <int[]> ptr(new int[pn + res_n]);
        int* tmp = ptr.get();
        multiply(tmp, revp, min(pn, res_n), inv, res_n);  /// only the first res_n coefficients are used
        reverse_poly(res, tmp, res_n);
    }

    void divide_remainder_inverse(int* quot, int* rem, int* P, int pn, int* Q, int qn, int* inv){
        if (pn < qn){
            copy_poly(rem, P, pn);
            for (int i = 0; i < qn - pn - 1; i++) rem[i + pn] = 0;
            return;
        }

        if (qn == 1) return;

        int quot_n = pn - qn + 1;
        int rn = qn - 1, tn = min(quot_n, rn), un = tn + rn;
        unique_ptr <int[]> ptr(new int[pn + un + (quot != 0 ? 0 : quot_n)]);

        int* revp = ptr.get(), *qmul = revp + pn;
        if (quot == 0) quot = qmul + un;

        reverse_poly(revp, P, pn);
        divide_inverse(quot, quot_n, revp, pn, inv);
        multiply(qmul, Q, rn, quot, tn);
        subtract(rem, P, rn, qmul, rn);
    }
};

struct LinearRecurrence{
    long long mod;
    vector<int> raw_base_sequence, base_sequence, recurrence;

    LinearRecurrence() {}

    /***
     * Construct linear recurrence from its first terms
     * sequence: without a recurrence, the first 2*n terms where n >= k and the recurrence has degree k
     *           with a recurrence of degree k, at least its first k terms
     * mod: all values are considered modulo this number
     * recurrence: an optional recurrence vector, calculated if not given
    ***/

    LinearRecurrence(vector<int> sequence, int mod, vector<int> recurrence={}) : mod(mod), recurrence(recurrence){
        int n = sequence.size(), m = recurrence.empty() ? n >> 1 : n;

        if (recurrence.empty() && (n == 0 || (n % 2) != 0)){
            throw std::invalid_argument("base sequence must be non-empty and of even length");
        }

        if (n < (int)recurrence.size()){
            throw std::invalid_argument("base sequence must have at least as many terms as the recurrence degree");
        }

        if ((long long)mod * 2 >= INT_MAX){
            throw std::out_of_range("mod value out of range, consider replacing long long with __int128 and int with long long if need to support higher ranges and improve fft multiplication precision");
        }

        for (int i = 0; i < n; i++){
            sequence[i] = ((long long)sequence[i] % mod + mod) % mod;
        }

        raw_base_sequence = sequence;
        base_sequence = {sequence.begin(), sequence.begin() + m};

        if (this->recurrence.size()){
            for (int i = 0; i < (int)this->recurrence.size(); i++){
                this->recurrence[i] %= mod;
                this->recurrence[i] = ((long long)mod - this->recurrence[i] + mod) % mod;
            }
        }
        else if (mod > 1){
            this->recurrence = decompose(sequence, mod);
            reverse(this->recurrence.begin(), this->recurrence.end());
            this->recurrence.pop_back();
        }
    }

    long long find_rank(const vector<int> & v){
        return (v.size() > 1 || (v.size() == 1 && v[0])) ? v.size() - 1 : (long long)INT_MIN;
    }

    long long compare(const vector<int> &a, const vector<int> &b){
        return max(find_rank(a), find_rank(b) + 1);
    };

    void update(vector<int>& u, vector<int>& v, int m, long long c){
        int i, k = v.size();
        if (k + m > (int)u.size()) u.resize(k + m, 0);

        for (i = 0; i < k; i++){
            u[i + m] -= c * v[i] % mod;
            if (u[i + m] < 0) u[i + m] += mod;
        }

        while (u.size() && !u.back()) u.pop_back();
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

    vector<int> berlekamp_massey(vector<int> sequence, long long mod){
        int n = sequence.size();
        vector<int> T, B(n + 1, 0), C(n + 1, 0);
        reverse(sequence.begin(), sequence.end());

        C[0] = 1, B[0] = 1;
        int i, j, d, l = 0, m = 1, deg = 0;
        unsigned long long x, b = nt::mod_inverse(1, mod);

        for (i = 0; i < n; i++){
            d = sequence[n - i - 1];
            if (l > 0) d = (d + nt::convolution(&C[1], &sequence[n - i], l, mod)) % mod;
            if (d == 0) m++;
            else{
                if (l * 2 <= i) T = C, T.resize(l + 1);
                x = (b * (mod - d) % mod + mod) % mod;

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

    vector<int> decompose(vector<int> sequence, int mod){
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
        for (int i = 1; i <= k; i++) q[i] = recurrence[k - i];

        auto p = fft::mod_multiply(vector<int>(raw_base_sequence.begin(), raw_base_sequence.begin() + k), q, mod);
        p.resize(k);

        for (; n; n >>= 1){
            for (int i = 0; i <= k; i++) v[i] = (i & 1) && q[i] ? mod - q[i] : q[i];
            p = fft::mod_multiply(p, v, mod, n & 1);
            q = fft::mod_multiply(q, v, mod, 0);
        }

        return p[0];
    }

    /***
     * Calculates the n'th, n+1'th, ... , n+k-1'th term of the recurrence
     * k must be at most the recurrence degree
    ***/

    vector<int> nth_terms(long long n, int k){
        vector<int> u = recurrence, v, res;

        if (mod == 1) return vector<int> (k, 0);
        if ((n + k - 1) < (int)base_sequence.size()){
            for (int i = 0; i < k; i++) res.push_back(base_sequence[(int)n + i]);
            return res;
        }

        int len = min(recurrence.size(), base_sequence.size());
        if (!len) return vector<int> (k, 0);
        if (k == 1) return {bostan_mori(n)};

        _POLYNOMIAL_MOD = mod;  /// set per call so recurrences with different mods can be interleaved
        u.push_back(1);
        auto p = Polynomial(u);
        p = Polynomial({0, 1}).power(p, n);

        /// p drops zero high coefficients, pad them back
        v.assign(recurrence.size(), 0);
        for (int i = 0; i < min((int)p.size(), (int)recurrence.size()); i++){
            v[i] = p.coefficient[i];
        }

        assert(len >= k);
        vector<int> terms = raw_base_sequence;
        while ((int)terms.size() < len + k - 1){  /// a given recurrence may come with only its first len terms
            long long val = 0;
            int s = terms.size();
            for (int i = 0; i < len; i++) val = (val + (mod - recurrence[i]) * terms[s - len + i]) % mod;
            terms.push_back(val);
        }

        for (int j = 0; j < k; j++){
            long long val = 0;
            for (int i = 0; i < len; i++){
                assert((i + j) < (int)terms.size());
                val = (val + (long long)v[i] * terms[i + j]) % mod;
            }
            res.push_back(val);
        }

        return res;
    }

    int nth_term(long long n){
        auto v = nth_terms(n, 1);
        return v[0];
    }
};

int main(){
    int i, j, k, v, mod;
    LinearRecurrence lr;
    vector<int> sequence, recurrence;

    /// recurrence given, f(x) = f(x - 3) + 2*f(x - 2) + 4*f(x-1), x > 2
    recurrence = {1, 2, 4};
    sequence = {0, 1, 1, 6, 27, 121}, mod = 100;

    lr = LinearRecurrence(sequence, mod, recurrence);
    assert(lr.nth_term(6) == 44);
    assert(lr.nth_term(9) == 90);
    assert(lr.nth_term(10) == 83);
    assert(lr.nth_term(1e18) == 19);

    /// recurrence given, f(x) = f(x - 1) + f(x - 2), x > 1, fibonacci series
    recurrence = {1, 1};
    sequence = {0, 1, 1, 2}, mod = 1000000007;

    lr = LinearRecurrence(sequence, mod, recurrence);
    assert(lr.nth_term(8) == 21);
    assert(lr.nth_term(9) == 34);
    assert(lr.nth_term(10) == 55);
    assert(lr.nth_term(1e18) == 209783453);

    /// alternatively, the recurrence of degree k can be derived from the first 2k terms
    lr = LinearRecurrence(sequence, mod);
    assert(lr.nth_term(8) == 21);
    assert(lr.nth_term(9) == 34);
    assert(lr.nth_term(10) == 55);
    assert(lr.nth_term(1e18) == 209783453);

    /// consecutive terms come from a single polynomial power
    assert(lr.nth_terms(8, 2) == vector<int>({21, 34}));

    /// with the recurrence given, its first k terms are enough
    lr = LinearRecurrence({0, 1}, mod, recurrence);
    assert(lr.nth_term(20) == 6765);
    assert(lr.nth_terms(9, 2) == vector<int>({34, 55}));
    assert(lr.raw_base_sequence == vector<int>({0, 1}));  /// queries leave the given terms untouched
    assert(lr.nth_term(1e18) == 209783453);

    /***
     *
     * Test performance on a recurrence with large degree
     * f(n) = n for n < k
     * f(n) = f(n-1) + 2*f(n-2) + 3*(f-3) + ... + k*f(n-k) for n >= k
     *
    ***/

    sequence.clear();
    k = 10000, mod = 1000000007;

    for (i = 0; i < 2 * k + 10; i++){
        if (i < k) sequence.push_back(i);
        else{
            for (v = 0, j = 1; j <= k; j++){
                v = (v + (long long)sequence[i - j] * j) % mod;
            }
            sequence.push_back(v);
        }
    }

    clock_t start = clock();
    lr = LinearRecurrence(sequence, mod);
    assert((int)lr.recurrence.size() == k);
    assert(lr.nth_term(1000000000000000000LL) == 255380209);

    fprintf(stderr, "\nTime taken = %0.5f\n", (clock() - start) / (double)CLOCKS_PER_SEC);  /// Time taken = 0.44215

    return 0;
}
