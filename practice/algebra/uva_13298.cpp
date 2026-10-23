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
 * nth_term with k = 1e5 and n = 1e18 takes about 3 s at -O2
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

    long long find_rank(const vector<int> & v){
        return (v.size() > 1 || (v.size() == 1 && v[0])) ? v.size() - 1 : (long long)INT_MIN;
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
};

/***
 *
 * Fast I/O
 * Buffered stdin and stdout through fread and fwrite, about 8x faster than scanf for integers and strings
 *
 * Complexity: O(1) amortized per character, one fread or fwrite call per 8192 bytes
 *
 * fio::read(a, b, ...): reads whitespace separated values into a, b, ... and returns how many it read before the input ended
 * read accepts integers of any width including __int128, float, double, long double, strings, chars and vectors of these
 * fio::read_line(s): reads one line into s without its newline or a trailing \r, false only at the end of input
 * fio::write(a, b, ...): writes a, b, ... separated by single spaces and ends the line, a vector argument on its own line
 * write accepts integers, strings, chars, vectors of these and fio::fixed values
 * fio::fixed(x, digits): a floating x printed with exactly digits digits after the point, rounded like printf("%.*f")
 * fio::flush(): writes the pending output now, which also happens at normal exit
 *
 * Once fio handles a stream, use nothing else on it: scanf, printf, cin or cout would reorder the input or output
 * A char reads one non-space character
 * A vector reads exactly its current size, so read a length first: read(n); a.resize(n); read(a);
 * Writing a float without fixed() does not compile, so its precision is always chosen
 * Floats print at printf's speed, since the formatting itself dominates; a floating token keeps its first 127 characters
 * Output is lost on abnormal termination such as abort() or _exit()
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// Because is_integral is not always true for __int128
#ifdef __SIZEOF_INT128__
    namespace std{
        template<> struct is_integral<__int128_t>: true_type {};
    }
#endif

/// Global buffers on purpose: stdin and stdout are one stream each per process, so a second reader instance
/// would split the input between two buffers and a second writer would reorder the output
namespace fio{
    const int BUF_SIZE = 8192;

    int buf_len = 0, inptr = 0, outptr = 0;
    char inbuf[BUF_SIZE], outbuf[BUF_SIZE], tmpbuf[128];

    /// Returns int so that EOF never collides with a 0xFF byte or with an unsigned char platform
    inline int read_char(){
        if (inptr >= buf_len){
            inptr = 0, buf_len = fread(inbuf, 1, BUF_SIZE, stdin);
            if (buf_len == 0) return EOF;
        }

        return (unsigned char)inbuf[inptr++];
    }

    template <typename T, typename=typename enable_if<is_integral<T>::value, T>::type>
    bool read_one(T &x){
        int c = ' ', neg = 0;
        while (c != '-' && !isdigit(c) && c != EOF) c = read_char();
        if (c == '-') neg = 1, c = read_char();
        if (c == EOF) return false;

        /// Accumulate as negative so that the minimum value (e.g. INT_MIN) does not overflow
        for (x = 0; isdigit(c); c = read_char()){
            x = x * 10 - (c - '0');
        }
        if (!neg) x = -x;

        return true;
    }

    /// enable_if in the return type, since a second default template argument would redeclare the integral overload
    template <typename T>
    typename enable_if<is_floating_point<T>::value, bool>::type read_one(T& x){
        int c = ' ', len = 0;
        while (isspace(c) && c != EOF) c = read_char();
        if (c == EOF) return false;

        for (; !isspace(c) && c != EOF; c = read_char()){
            if (len + 1 < (int)sizeof(tmpbuf)) tmpbuf[len++] = c;
        }
        tmpbuf[len] = 0;

        /// Parsing at the target's own width, since strtold then a cast to double rounds twice and can miss by one ulp
        if (is_same<T, float>::value) x = strtof(tmpbuf, nullptr);
        else if (is_same<T, double>::value) x = strtod(tmpbuf, nullptr);
        else x = strtold(tmpbuf, nullptr);
        return true;
    }

    bool read_one(string &s){
        int c = ' ';
        while (isspace(c) && c != EOF) c = read_char();
        if (c == EOF) return false;

        for (s.clear(); !isspace(c) && c != EOF; c = read_char()){
            s.push_back(c);
        }
        return true;
    }

    bool read_one(char& ch){
        int c;
        do c = read_char(); while (isspace(c));
        if (c == EOF) return false;

        ch = c;
        return true;
    }

    template <typename T>
    bool read_one(vector<T>& v){
        for (auto &&x: v){
            if (!read_one(x)) return false;
        }
        return true;
    }

    /// Returns false only at the end of input, a last line without a trailing newline is still returned
    bool read_line(string &s){
        s.clear();
        int c = read_char();
        if (c == EOF) return false;

        for (; c != '\n' && c != EOF; c = read_char()) s.push_back(c);
        if (!s.empty() && s.back() == '\r') s.pop_back();  /// CRLF line endings
        return true;
    }

    int read(){
        return 0;
    }

    template <typename T, typename ...Args>
    int read(T &x, Args& ...args){
        if (!read_one(x)) return 0;
        return read(args...) + 1;
    }

    void flush(){
        fwrite(outbuf, 1, outptr, stdout);
        outptr = 0;
    }

    struct AutoFlush{ ~AutoFlush(){ flush(); } } auto_flush;

    inline void write_char(const char& c){
        if (outptr == BUF_SIZE) flush();
        outbuf[outptr++] = c;
    }

    void write_one(char c){
        write_char(c);
    }

    void write_one(const char* s){
        for (int j = 0; s[j]; j++) write_char(s[j]);
    }

    void write_one(const string& s){
        for (auto &&c: s) write_char(c);
    }

    template <typename T, typename=typename enable_if<is_integral<T>::value, T>::type>
    void write_one(T x){
        int l = 0;
        if (x < 0){
            /// Peel one digit while still negative so that -x cannot overflow for the minimum value
            write_char('-');
            tmpbuf[l++] = '0' - x % 10;
            x = -(x / 10);
        }

        while (x || !l){
            tmpbuf[l++] = (x % 10) + '0';
            x /= 10;
        }
        while (l) write_char(tmpbuf[--l]);
    }

    /// write(x) for a floating x does not compile on purpose: the digits after the point must be chosen
    struct Fixed{
        long double value;
        int digits;
    };

    Fixed fixed(long double value, int digits){
        return {value, digits};
    }

    void write_one(const Fixed& f){
        int len = snprintf(tmpbuf, sizeof(tmpbuf), "%.*Lf", f.digits, f.value);
        if (len < (int)sizeof(tmpbuf)){
            for (int j = 0; j < len; j++) write_char(tmpbuf[j]);
            return;
        }

        string s(len + 1, '\0');
        snprintf(&s[0], len + 1, "%.*Lf", f.digits, f.value);
        s.resize(len);
        write_one(s);
    }

    template <typename T>
    void write_one(const vector<T>& v){
        for (int i = 0; i < (int)v.size(); i++){
            if (i) write_char(' ');
            write_one(v[i]);
        }
    }

    void write() {}

    template <typename T> struct is_vector : false_type {};
    template <typename T, typename A> struct is_vector<vector<T, A>> : true_type {};

    template <typename T, typename ...Args>
    void write(const T& x, const Args& ...args){
        write_one(x);
        write_char(sizeof...(args) && !is_vector<T>::value ? ' ' : '\n');
        write(args...);
    }
}

using namespace fio;

#define MOD 1000000009

int main(){
    int k, i, j;
    long long n;
    map<int, LinearRecurrence> family;

    while (read(k, n) == 2 && (k || n)){
        if (!family.count(k)){
            vector<int> first(k, 0);
            for (i = 0; i < k; i++){
                first[i] = (i == 0);
                for (j = max(0, i - k); j < i; j++) first[i] = (first[i] + first[j]) % MOD;
            }
            family.emplace(k, LinearRecurrence(first, MOD, vector<int>(k, 1)));
        }
        write(family.at(k).nth_term(n));
    }
    return 0;
}
