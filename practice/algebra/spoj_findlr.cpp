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
 * LinearRecurrence lr(sequence, mod, recurrence): f(i) = recurrence[0] f(i - 1) + ... + recurrence[k - 1] f(i - k)
 * LinearRecurrence lr(sequence, mod): finds a recurrence generating sequence, stored in lr.recurrence in the same order
 * lr.nth_term(n): f(n), 0 <= n <= 1e18
 * lr.nth_terms(n, c): f(n), ..., f(n + c - 1), 1 <= c <= k
 *
 * Notes:
 *   - sequence holds f(0), f(1), ..., at least k terms when the recurrence is given
 *   - 1 <= mod and 2 * mod < LLONG_MAX, values of any sign are reduced modulo mod
 *   - nth_term and nth_terms split values into 15 bit halves for the FFT, so need mod < 2^30
 *   - 2k terms are enough to find a recurrence of degree k
 *   - With a prime mod the found recurrence is the shortest one generating every given term, from any number of terms
 *   - With a composite mod the recurrence may end in zero coefficients, and may miss given terms when fewer than 2k
 *   - Requires __int128
 *   - Thanks to anta (codeforces.com/profile/anta) for showing me this trick and the Black Box Linear Algebra book
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

namespace nt{
    map<long long, long long> get_factors(long long mod){
        __int128 d, x = mod;
        map<long long, long long> factors;

        for (d = 2; d * d <= x; d += (d & 1) + 1){
            while (x % d == 0){
                factors[d]++;
                x /= d;
            }
        }
        if (x > 1) factors[x]++;

        return factors;
    }

    long long extended_gcd(long long a, long long b, long long& x, long long& y){
        if (!b){
            y = 0, x = 1;
            return a;
        }

        long long g = extended_gcd(b, a % b, y, x);
        y -= ((a / b) * x);
        return g;
    }

    long long mod_inverse(long long a, long long mod){
        long long x, y;
        extended_gcd(a, mod, x, y);
        return (x + mod) % mod;
    }

    long long chinese_remainder_theorem(vector<long long> ar, vector<long long> mods){
        long long i, res = 0, mod = 1;

        for (auto m: mods) mod *= m;

        for (i = 0; i < (long long)ar.size(); i++){
            __int128 x = mod / mods[i];
            __int128 y = mod_inverse(x, mods[i]);
            res = (res + (x * ar[i] % mod) * y) % mod;
        }

        return res;
    }

    long long convolution(long long* A, long long* B, long long n, long long mod){
        unsigned __int128 res = 0;
        long long i, j, lim = min(ULONG_LONG_MAX / ((unsigned __int128)mod * mod), (unsigned __int128)n);

        for (i = 0; (i + lim) <= n; res %= mod){
            for (j = 0; j < lim; j++, i++) res += (unsigned __int128)A[i] * B[i];
        }

        for (; i < n; i++) res += (unsigned __int128)A[i] * B[i];

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

    void extend_roots(long long len){
        for (long long k = roots.size(); k < len; k <<= 1){
            roots.resize(2 * k);
            for (long long j = 0; j < k; j++){
                long double theta = acosl(-1.0L) * j / k;  /// direct, a chain of products accumulates rounding error
                roots[k + j] = ComplexNum(cosl(theta), sinl(theta));
            }
        }
    }

    void transform(vector<ComplexNum>& a){
        long long i, j, k, len = a.size();

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

    /// Splits the transforms at i and j (= -i) into low and high 15 bit halves, products of the halves go to u and v
    inline void combine(ComplexNum fi, ComplexNum fj, ComplexNum gi, ComplexNum gj, long long len, ComplexNum& u, ComplexNum& v){
        auto a1 = (fi + fj.conjugate()) * ComplexNum(0.5, 0);
        auto a2 = (fi - fj.conjugate()) * ComplexNum(0, -0.5);
        auto b1 = (gi + gj.conjugate()) * ComplexNum(0.5 / len, 0);
        auto b2 = (gi - gj.conjugate()) * ComplexNum(0, -0.5 / len);
        v = a1 * b2 + a2 * b1;
        u = a1 * b1 + a2 * b2 * ComplexNum(0, 1);
    }

    /// parity 0 or 1 returns only the coefficients of p1 * p2 at indices of that parity, -1 returns all of them
    vector<long long> mod_multiply(const vector<long long>& p1, const vector<long long>& p2, long long mod, long long parity = -1){
        long long i, j, a = p1.size(), b = p2.size(), p_len = a + b - 1, step = parity < 0 ? 1 : 2;

        while (a > 1 && p1[a - 1] == 0) a--;
        while (b > 1 && p2[b - 1] == 0) b--;
        bool is_equal = (a == b && equal(p1.begin(), p1.begin() + a, p2.begin()));

        long long len = 2;
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

        vector<long long> res((p_len - max(parity, 0LL) + step - 1) / step, 0);
        for (i = 0; i < min(len, (long long)res.size()); i++){
            __int128 x = f[i].real + 0.5, y = g[i].real + 0.5, z = f[i].img + 0.5;
            res[i] = (x + ((y % mod) << 15) + ((z % mod) << 30)) % mod;
        }

        return res;
    }
}

struct LinearRecurrence{
    __int128 mod;
    vector<long long> sequence, recurrence;

    LinearRecurrence() {}

    LinearRecurrence(vector<long long> sequence, long long mod, vector<long long> recurrence={}) : mod(mod), sequence(sequence), recurrence(recurrence){
        long long n = sequence.size(), k = recurrence.size();

        assert(n >= k && mod >= 1 && 2 * mod < LLONG_MAX);

        for (auto& x : this->sequence) x = ((__int128)x % mod + mod) % mod;
        for (auto& c : this->recurrence) c = ((__int128)c % mod + mod) % mod;
        if (recurrence.empty() && n > 0 && mod > 1){
            auto C = find_recurrence(this->sequence, mod);  /// C = [1, C1, ..., Cd] with sum C[j] * f(x - j) = 0
            for (long long j = 1; j < (long long)C.size(); j++) this->recurrence.push_back((mod - C[j]) % mod);
        }
    }

    long long nth_term(__int128 n){
        return nth_terms(n, 1)[0];
    }

    vector<long long> nth_terms(__int128 n, long long k){
        long long len = recurrence.size();

        if (n + k <= (__int128)sequence.size()) return vector<long long>(sequence.begin() + n, sequence.begin() + n + k);
        if (mod == 1 || len == 0) return vector<long long>(k, 0);
        assert(1 <= k && k <= len);
        if (k == 1) return {bostan_mori(n)};

        auto p = power_of_x(n);

        vector<long long> terms = sequence;
        while ((long long)terms.size() < len + k - 1){  /// a given recurrence may come with only its first len terms
            __int128 val = 0;
            long long s = terms.size();
            for (long long i = 0; i < len; i++) val = (val + (__int128)recurrence[i] * terms[s - 1 - i]) % mod;
            terms.push_back(val);
        }

        vector<long long> res(k);
        for (long long j = 0; j < k; j++){
            __int128 val = 0;
            for (long long i = 0; i < (long long)p.size(); i++) val = (val + (__int128)p[i] * terms[i + j]) % mod;
            res[j] = val;
        }

        return res;
    }

    vector<long long> find_recurrence(vector<long long> sequence, long long mod){
        long long i, j, e, n, p, m, k = 0;
        auto factors = nt::get_factors(mod);

        n = factors.size();

        if (n == 1 && factors.begin()->second == 1){
            return berlekamp_massey(sequence, mod);
        }

        vector<long long> v, mods, vals(n);
        vector<vector<long long>> coeffs;

        for (auto it: factors){
            tie(p, e) = it;
            m = powl(p, e) + 0.5;
            mods.push_back(m);

            v = reeds_sloane(sequence, p, e);
            k = max(k, (long long)v.size());
            coeffs.push_back(v);
        }

        vector<long long> res(k);
        for (i = 0; i < k; i++){
            for (j = 0; j < n; j++){
                vals[j] = i < (long long)coeffs[j].size() ? coeffs[j][i] : 0;
            }
            res[i] = nt::chinese_remainder_theorem(vals, mods);
        }

        res.resize(sequence.size() / 2 + 1, 0);
        return res;
    }

    /// Bostan-Mori: f(n) = [x^n] P(x) / Q(x), each step multiplies both by Q(-x) and keeps every second coefficient
    long long bostan_mori(__int128 n){
        long long k = recurrence.size();
        vector<long long> q(k + 1), v(k + 1);

        q[0] = 1;
        for (long long i = 1; i <= k; i++) q[i] = (mod - recurrence[i - 1]) % mod;

        auto p = fft::mod_multiply(vector<long long>(sequence.begin(), sequence.begin() + k), q, mod);
        p.resize(k);

        for (; n; n >>= 1){
            for (long long i = 0; i <= k; i++) v[i] = (i & 1) && q[i] ? mod - q[i] : q[i];
            p = fft::mod_multiply(p, v, mod, n & 1);
            q = fft::mod_multiply(q, v, mod, 0);
        }

        return p[0];
    }

    /// x^n modulo x^k - c_1 x^(k - 1) - ... - c_k, lowest coefficient first, for n >= 1 and k >= 2
    vector<long long> power_of_x(__int128 n){
        long long k = recurrence.size();
        vector<long long> charpoly(k + 1, 1);

        for (long long j = 0; j < k; j++) charpoly[j] = (mod - recurrence[k - 1 - j]) % mod;
        vector<long long> reversed(charpoly.rbegin(), charpoly.rend());

        /// Newton step inv = inv (2 - reversed inv) doubles the correct prefix, k - 1 terms cover every quotient
        vector<long long> inv = {1};
        while ((long long)inv.size() < k - 1){
            long long len = min(2 * (long long)inv.size(), k - 1);
            auto t = fft::mod_multiply(inv, vector<long long>(reversed.begin(), reversed.begin() + len), mod);
            t.resize(len);
            for (auto& x : t) x = x ? mod - x : 0;
            t[0] = (t[0] + 2) % mod;
            inv = fft::mod_multiply(inv, t, mod);
            inv.resize(len);
        }

        vector<long long> res = {1};
        for (long long bit = __lg((long long)n); bit >= 0; bit--){
            res = fft::mod_multiply(res, res, mod);

            if ((long long)res.size() > k){  /// reversed quotient = reversed top m coefficients * (1 / reversed charpoly)
                long long m = res.size() - k;
                auto quot = fft::mod_multiply(vector<long long>(res.rbegin(), res.rbegin() + m), vector<long long>(inv.begin(), inv.begin() + m), mod);
                quot.resize(m);
                reverse(quot.begin(), quot.end());

                auto sub = fft::mod_multiply(vector<long long>(charpoly.begin(), charpoly.begin() + k), quot, mod);
                res.resize(k);
                for (long long j = 0; j < k; j++) res[j] = (res[j] - sub[j] + mod) % mod;
            }

            if (n >> bit & 1){
                res.insert(res.begin(), 0);
                if ((long long)res.size() > k){  /// x^k reduces to x^k - charpoly
                    __int128 top = res.back();
                    res.pop_back();
                    for (long long j = 0; j < k; j++) res[j] = (res[j] + (mod - charpoly[j]) * top) % mod;
                }
            }
        }

        return res;
    }

    vector<long long> berlekamp_massey(vector<long long> sequence, __int128 mod){
        long long n = sequence.size();
        vector<long long> T, B(n + 1, 0), C(n + 1, 0);

        reverse(sequence.begin(), sequence.end());

        C[0] = 1, B[0] = 1;
        long long i, j, d, l = 0, m = 1, deg = 0;
        unsigned __int128 x, b = 1;

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

    vector<long long> reeds_sloane(vector<long long> sequence, long long p, long long e){
        vector<long long> T(e), U(e), R(e), TO(e, 1), TO_INV(e, 1), UO(e), POW(e + 1, 1);
        long long c, d, i, j, k, l, m, len = sequence.size(), mod = powl(p, e) + 0.5;
        unsigned long long inv = ~0ULL / mod;

        vector<vector<long long>> A[3], B[3];
        for (i = 0; i < 3; i++){
            A[i] = vector<vector<long long>>(e);
            B[i] = vector<vector<long long>>(e);
        }

        for (i = 0; i < len; i++) sequence[i] %= mod;
        for (i = 1; i <= e; i++) POW[i] = (__int128)POW[i - 1] * p % mod;
        for (i = 0; i < e; i++){
            T[i] = (__int128)sequence[0] * POW[i] % mod;
            A[0][i] = {POW[i]}, A[1][i] = {POW[i]};
            B[0][i] = {0}, B[1][i] = {(long long)((__int128)sequence[0] * POW[i] % mod)};
            for (U[i] = 0; T[i] && T[i] % p == 0; T[i] /= p, U[i]++) {}
            if (!T[i]) T[i] = 1, U[i] = e;
        }

        for (k = 1; k < len; k++){
            for (j = 0; j < e; j++){
                if (compare(A[0][j], B[0][j]) < compare(A[1][j], B[1][j])){
                    R[j] = k - 1;
                    A[2][j] = A[0][e - 1 - U[j]], B[2][j] = B[0][e - 1 - U[j]];
                    UO[j] = U[e - 1 - U[j]], TO[j] = T[e - 1 - U[j]];
                    TO_INV[j] = nt::mod_inverse(TO[j], mod);
                }
            }

            A[0] = A[1], B[0] = B[1];
            for (l = 0; l < e; l++){
                for (d = 0, i = min(k, (long long)A[0][l].size() - 1); i >= 0; i--){
                    d = reduce(d + (unsigned __int128)A[0][l][i] * sequence[k - i], mod, inv);
                }

                if (!d) T[l] = 1, U[l] = e;
                else{
                    for (U[l] = 0, T[l] = d; T[l] % p == 0; T[l] /= p, U[l]++) {}
                    j = e - 1 - U[l];

                    if (!compare(A[0][j], B[0][j])){
                        if ((k + 1) > (long long)B[1][l].size()) B[1][l].resize(k + 1, 0);
                        B[1][l][k] = (B[1][l][k] + d) % mod;
                    }
                    else{
                        m = k - R[j];
                        c = (__int128)T[l] * TO_INV[j] % mod;
                        c = (__int128)c * (__int128)POW[U[l] - UO[j]] % mod;
                        update(A[1][l], A[2][j], m, c);
                        update(B[1][l], B[2][j], m, c);
                    }
                }
            }
        }

        return A[1][0];
    }

    __int128 compare(const vector<long long> &a, const vector<long long> &b){
        return max(find_rank(a), find_rank(b) + 1);
    }

    /// Barrett reduction of x < 2^62 modulo m < 2^31 with inv = (2^64 - 1) / m, about 2x faster than % in Reeds-Sloane
    unsigned long long reduce(unsigned long long x, unsigned long long m, unsigned long long inv){
        unsigned long long r = x - (unsigned long long)((unsigned __int128)x * inv >> 64) * m;
        return r >= m ? r - m : r;
    }

    void update(vector<long long>& u, vector<long long>& v, long long m, __int128 c){
        long long i, k = v.size();
        unsigned long long inv = ~0ULL / mod;

        if (k + m > (long long)u.size()) u.resize(k + m, 0);

        for (i = 0; i < k; i++){
            u[i + m] -= (long long)reduce(c * v[i], mod, inv);
            if (u[i + m] < 0) u[i + m] += mod;
        }

        while (u.size() && !u.back()) u.pop_back();
    }

    __int128 find_rank(const vector<long long> & v){
        return (v.size() > 1 || (v.size() == 1 && v[0])) ? v.size() - 1 : (__int128)INT_MIN;
    }
};

/***
 *
 * Fast I/O
 * Buffered stdin and stdout, about 5x faster than scanf for integers, interactive problems included
 *
 * Complexity: O(1) amortized per character, one system call per buffer of up to 8192 bytes
 *
 * fio::flush(): writes the pending output now, which also happens at normal exit
 * fio::read(a, b, ...): reads whitespace separated values into a, b, ..., returns the count read before the input ended
 * fio::read_char(): the next input byte, whitespace included, EOF at the end of input
 * fio::read_line(s): reads one line into s without its newline or a trailing \r, false only at the end of input
 * fio::write(a, b, ...): writes a, b, ... space separated and ends the line, a vector argument on its own line
 * fio::fixed(x, digits): a floating x with exactly digits digits after the point, the same text as printf("%.*f")
 *
 * Notes:
 *   - read takes any integer including __int128, float, double, long double, strings, chars and vectors of these
 *   - Reading into a char array does not compile: read a string, or a grid as vector<string>
 *   - write takes integers, strings, chars, vectors of these and fio::fixed values
 *   - Once fio handles a stream, use nothing else on it: scanf, printf, cin or cout would reorder the input or output
 *   - An integer read skips everything up to the next digit or -, words included, where scanf fails and reads nothing
 *   - A char reads one non-space character
 *   - A number or string read eats the byte after it: read_line next gets the next line after "9\n", "" after "9\r\n"
 *   - A text line after a number line, under any line ending: read_line(s), k = stoi(s), then read_line(text)
 *   - A loop over read_char() parses characters 2-3x faster than read(char) and about 4x faster than read(string)
 *   - A vector reads exactly its current size, so read a length first: read(n); a.resize(n); read(a);
 *   - read(u, ar[u], d) binds ar[u] before u is read, so read an index on its own first: read(u); read(ar[u], d);
 *   - Writing a float without fixed() does not compile, so its precision is always chosen
 *   - fixed() is about 8x faster than printf while |x| * 10^digits < 2^52 and digits <= 18, as fast as printf otherwise
 *   - A floating token read keeps its first 127 characters
 *   - A float, double or long double read gives strtof's, strtod's or strtold's value, 2-3x faster on short decimals
 *   - Interactive problems need no flush(): pending output is written each time the input buffer runs empty
 *   - Output is lost on abnormal termination such as abort() or _exit()
 *
***/

#include <bits/stdtr1c++.h>
#include <unistd.h>  /// read() and write(), which MinGW's unistd.h also declares, so Windows judges compile them too

using namespace std;

/// Global buffers on purpose: a second reader would split stdin between two buffers, a second writer reorder stdout
namespace fio{
    const int BUF_SIZE = 8192;

    int buf_len = 0, inptr = 0, outptr = 0;
    char inbuf[BUF_SIZE], outbuf[BUF_SIZE], tmpbuf[128];

    /// write, not fwrite: fwrite strands bytes in stdout's buffer, and fflush after it costs a second call, 7% slower
    void flush(){
        int done = 0;
        while (done < outptr){
            int len = ::write(1, outbuf + done, outptr - done);
            if (len <= 0) break;
            done += len;
        }
        outptr = 0;
    }

    /// Returns int so that EOF never collides with a 0xFF byte or with an unsigned char platform
    inline int read_char(){
        if (inptr >= buf_len){
            /// An interactive judge answers only after it has read our query, so the query leaves before we wait
            flush();

            /// read returns the bytes available, fread waits for all 8192 or EOF: a deadlock on an interactive pipe
            inptr = 0, buf_len = ::read(0, inbuf, BUF_SIZE);
            if (buf_len <= 0) return EOF;
        }

        return (unsigned char)inbuf[inptr++];
    }

    /// Catches every integer the overloads below miss, so __int128 needs no is_integral specialization (breaks C++20)
    template <typename T>
    bool read_one(T &x){
        int c = ' ', neg = 0;

        while (c != '-' && !isdigit(c) && c != EOF) c = read_char();
        if (c == '-') neg = 1, c = read_char();
        if (c == EOF) return false;

        /// A local, since an int x may alias inptr or buf_len and writing it per digit forces their reload from memory
        T v = 0;

        /// Accumulate as negative so that the minimum value (e.g. INT_MIN) does not overflow
        for (; isdigit(c); c = read_char()){
            v = v * 10 - (c - '0');
        }
        x = neg ? v : -v;

        return true;
    }

    /// Returns the length kept in tmpbuf, 0 only at the end of input
    int read_token(){
        int c = ' ', len = 0;

        while (isspace(c) && c != EOF) c = read_char();
        if (c == EOF) return 0;

        for (; !isspace(c) && c != EOF; c = read_char()){
            if (len + 1 < (int)sizeof(tmpbuf)) tmpbuf[len++] = c;
        }
        tmpbuf[len] = 0;
        return len;
    }

    /// Digits m and 10^|e| both exact in F round once, as parse does: 2-3x faster than parse on short tokens
    template <typename F, F (*parse)(const char *, char **)>
    bool read_floating(F &x, int max_len, unsigned long long max_mantissa, int max_exp){
        int len = read_token();
        if (!len) return false;

        /// x87 arithmetic (FLT_EVAL_METHOD != 0) rounds float and double twice, long double once
        bool own_width = !FLT_EVAL_METHOD || is_same<F, long double>::value;

        /// A longer token skips the digit scan, which would make its parse about 10% slower
        if (len <= max_len && own_width){
            static const F p10[] = {1e0L, 1e1L, 1e2L, 1e3L, 1e4L, 1e5L, 1e6L, 1e7L, 1e8L, 1e9L,
                                    1e10L, 1e11L, 1e12L, 1e13L, 1e14L, 1e15L, 1e16L, 1e17L, 1e18L, 1e19L,
                                    1e20L, 1e21L, 1e22L, 1e23L, 1e24L, 1e25L, 1e26L, 1e27L};
            const unsigned char *s = (const unsigned char *)tmpbuf + (*tmpbuf == '-' || *tmpbuf == '+');
            unsigned long long m = 0;
            int digits = 0, e = 0;
            for (; isdigit(*s); s++) m = m * 10 + (*s - '0'), digits++;
            if (*s == '.'){
                for (s++; isdigit(*s); s++) m = m * 10 + (*s - '0'), digits++, e--;
            }

            /// No exponent digits after the e equals strtod, which then ignores the e
            if (*s == 'e' || *s == 'E'){
                int sign = s[1] == '-' ? -1 : 1, ex = 0;
                s += 1 + (s[1] == '-' || s[1] == '+');
                for (; isdigit(*s) && ex < 1000; s++) ex = ex * 10 + (*s - '0');
                e += sign * ex;
            }

            if (!*s && digits && digits <= 19 && m <= max_mantissa && e >= -max_exp && e <= max_exp){
                F v = e < 0 ? m / p10[-e] : m * p10[e];
                x = *tmpbuf == '-' ? -v : v;
                return true;
            }
        }

        x = parse(tmpbuf, nullptr);
        return true;
    }

    /// Each float parses at its own width, since strtold then a cast to double rounds twice and can miss by one ulp
    bool read_one(float &x){
        return read_floating<float, strtof>(x, 10, 1 << 24, 10);
    }

    /// 22 characters keep 16-digit scientific tokens like -1.234567890123456e-05 on the fast path
    bool read_one(double &x){
        return read_floating<double, strtod>(x, 22, 1ULL << 53, 22);
    }

    /// Only an x87 64-bit mantissa takes the fast path: a double-double long double does not round a division once
    bool read_one(long double &x){
        return read_floating<long double, strtold>(x, 21, ~0ULL, LDBL_MANT_DIG == 64 ? 27 : -1);
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

    bool read_one(char &ch){
        int c;
        do c = read_char(); while (isspace(c));
        if (c == EOF) return false;

        ch = c;
        return true;
    }

    template <typename T>
    bool read_one(vector<T> &v){
        for (auto &&x: v){
            if (!read_one(x)) return false;
        }

        return true;
    }

    template <typename T>
    int read(T &x){
        return read_one(x);
    }

    template <typename T, typename ...Args>
    int read(T &x, Args &...args){
        if (!read_one(x)) return 0;
        return read(args...) + 1;
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

    struct AutoFlush{ ~AutoFlush(){ flush(); } } auto_flush;

    inline void write_char(char c){
        if (outptr == BUF_SIZE) flush();
        outbuf[outptr++] = c;
    }

    /// Integers, like read_one: the plain overloads below win for every other type
    template <typename T>
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

    void write_one(char c){
        write_char(c);
    }

    void write_one(const char *s){
        for (int j = 0; s[j]; j++) write_char(s[j]);
    }

    /// Without it a char * picks the integer template, an exact match that beats the const char * overload
    void write_one(char *s){
        write_one((const char *)s);
    }

    void write_one(const string &s){
        for (auto &&c: s) write_char(c);
    }

    /// write(x) for a floating x does not compile on purpose: the digits after the point must be chosen
    struct Fixed{
        long double value;
        int digits;
    };

    void write_one(const Fixed &f){
        static const long long p10[] = {1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000,
            10000000000, 100000000000, 1000000000000, 10000000000000, 100000000000000, 1000000000000000,
            10000000000000000, 100000000000000000, 1000000000000000000};

        if (f.digits >= 0 && f.digits <= 18){
            long double scaled = fabsl(f.value) * p10[f.digits];

            /// Rounding never moves the product across k + 0.5, exact below 2^52, so only an exact k + 0.5 needs printf
            if (scaled < 4503599627370496.0L && scaled - (long long)scaled != 0.5L){
                long long q = scaled + 0.5L;
                if (signbit(f.value)) write_char('-');
                write_one(q / p10[f.digits]);
                if (!f.digits) return;

                long long r = q % p10[f.digits];
                for (int j = f.digits - 1; j >= 0; j--) tmpbuf[j] = '0' + r % 10, r /= 10;
                write_char('.');
                for (int j = 0; j < f.digits; j++) write_char(tmpbuf[j]);
                return;
            }
        }

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
    void write_one(const vector<T> &v){
        for (int i = 0; i < (int)v.size(); i++){
            if (i) write_char(' ');
            write_one(v[i]);
        }
    }

    template <typename T>
    char separator_after(const T &){
        return ' ';
    }

    template <typename T>
    char separator_after(const vector<T> &){
        return '\n';
    }

    template <typename T>
    void write(const T &x){
        write_one(x);
        write_char('\n');
    }

    template <typename T, typename ...Args>
    void write(const T &x, const Args &...args){
        write_one(x);
        write_char(separator_after(x));
        write(args...);
    }

    Fixed fixed(long double value, int digits){
        return {value, digits};
    }
}

using namespace fio;

int main(){
    int t, k, m, i;

    read(t);
    while (t--){
        read(k, m);
        vector<long long> seq(2 * k);
        read(seq);

        LinearRecurrence lr(seq, m);

        long long res = 0;
        for (i = 0; i < (int)lr.recurrence.size(); i++){
            res = (res + (__int128)lr.recurrence[i] * seq[2 * k - 1 - i]) % m;
        }
        write(res);
    }

    return 0;
}
