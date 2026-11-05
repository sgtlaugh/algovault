/***
 *
 * Wildcard Pattern Matching with NTT
 * Finds every position where pattern occurs in text, a wildcard in either string matches any character
 *
 * Complexity: O(n log n) for text length n, three convolutions modulo 998244353
 *
 * wildcard_match(text, pattern, wildcard = '?'): the 0-indexed start positions in increasing order
 *
 * Notes:
 *   - Works on strings and string literals mixed freely, and on vectors of any ordered type given the wildcard value
 *   - An empty pattern or one longer than the text returns no positions, as kmp_search does
 *   - Requires n <= 2^23, the largest NTT length modulo 998244353
 *   - A real match is always reported, a mismatching position is reported with probability at most 4 / 998244352
 *   - A call reports some false position with probability at most (n - m + 1) * 4 / 998244352, about 3% at n = 2^23
 *   - Small alphabets err far less: with two letters a call errs only if both random weights collide
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

namespace wildcard_ntt{
    const unsigned MOD = 998244353;
    const unsigned ROOT = 3;
    const int max_log = 23;

    /// BEGIN COPY ntt_core from code_library/algebra/ntt.cpp
    constexpr unsigned long long power(unsigned long long b, unsigned long long e, unsigned long long mod){
        unsigned long long res = 1;
        b %= mod;
        while (e){
            if (e & 1) res = res * b % mod;
            b = b * b % mod;
            e >>= 1;
        }

        return res;
    }

    /// Newton's step y (2 - p y) doubles the correct low bits of 1 / p modulo 2^32, y = p starts with 3 of them
    constexpr unsigned inverse_mod_2_32(unsigned p, unsigned y, int steps){
        return steps ? inverse_mod_2_32(p, y * (2 - p * y), steps - 1) : y;
    }

    /// Values live in [0, BOUND), 2 MOD below 2^30 to skip the final reduction, else MOD, so u - v + BOUND < 2^32
    template<unsigned MOD>
    constexpr unsigned bound(){
        return MOD < (1u << 30) ? 2 * MOD : MOD;
    }

    /// Montgomery reduction, x / 2^32 modulo MOD for x < MOD 2^32, in [0, bound<MOD>()), replacing the slower x % MOD
    template<unsigned MOD>
    unsigned redc(unsigned long long x){
        constexpr unsigned NEG_INV = 0u - inverse_mod_2_32(MOD, MOD, 4);
        unsigned r = (x + (unsigned long long)((unsigned)x * NEG_INV) * MOD) >> 32;
        return MOD < (1u << 30) || r < MOD ? r : r - MOD;
    }

    /// rt[k + j] = w^j R modulo MOD for each power of two k < n, w a primitive 2k-th root of unity or its inverse
    template<unsigned MOD, unsigned ROOT>
    vector<unsigned> roots(int n, bool invert){
        vector<unsigned> rt(max(n, 2), ((unsigned long long)1 << 32) % MOD);
        for (int k = 2; k < n; k <<= 1){
            unsigned long long z = power(ROOT, (MOD - 1) / (2 * k), MOD);
            if (invert) z = power(z, MOD - 2, MOD);
            z = (z << 32) % MOD;
            for (int i = k; i < 2 * k; i++){
                unsigned w = redc<MOD>(rt[i >> 1] * z);
                rt[i] = i & 1 ? (w >= MOD ? w - MOD : w) : rt[i >> 1];
            }
        }

        return rt;
    }

    /// Decimation in frequency, natural order in, bit-reversed order out, a.size() a power of two at most rt.size()
    template<unsigned MOD>
    void forward_transform(vector<unsigned>& a, const vector<unsigned>& rt){
        const unsigned BOUND = bound<MOD>();
        int n = a.size();
        for (int k = n / 2; k >= 1; k >>= 1){
            for (int i = 0; i < n; i += 2 * k){
                for (int j = 0; j < k; j++){
                    unsigned u = a[i + j], v = a[i + j + k];
                    a[i + j] = u + v >= BOUND ? u + v - BOUND : u + v;
                    a[i + j + k] = redc<MOD>((unsigned long long)(u - v + BOUND) * rt[j + k]);
                }
            }
        }
    }

    /// Decimation in time with the inverse roots, bit-reversed order in, natural order out, scaled by n
    template<unsigned MOD>
    void inverse_transform(vector<unsigned>& a, const vector<unsigned>& rt){
        const unsigned BOUND = bound<MOD>();
        int n = a.size();
        for (int k = 1; k < n; k <<= 1){
            for (int i = 0; i < n; i += 2 * k){
                for (int j = 0; j < k; j++){
                    unsigned u = a[i + j], v = redc<MOD>((unsigned long long)a[i + j + k] * rt[j + k]);
                    a[i + j] = u + v >= BOUND ? u + v - BOUND : u + v;
                    a[i + j + k] = u >= v ? u - v : u + BOUND - v;
                }
            }
        }
    }
    /// END COPY ntt_core

    /// Transform of w^k placed at the front of a zero-padded array of length len
    vector<unsigned> transformed_powers(const vector<unsigned>& weights, int k, int len, const vector<unsigned>& rt){
        vector<unsigned> a(len);
        for (int i = 0; i < (int)weights.size(); i++) a[i] = power(weights[i], k, MOD);
        forward_transform<MOD>(a, rt);
        return a;
    }
}

template<typename Container>
vector<int> wildcard_match(const Container& text, const Container& pattern, typename Container::value_type wildcard = '?'){
    using namespace wildcard_ntt;
    using T = typename Container::value_type;

    int n = text.size(), m = pattern.size();
    vector<int> positions;
    if (m == 0 || m > n) return positions;
    assert(n <= (1 << max_log));

    mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());
    auto random_weight = [&]{ return uniform_int_distribution<unsigned>(1, MOD - 1)(rng); };
    vector<unsigned> wt(n), wp(m);

    /// Equal values share one weight, from a byte table or runs of a sort: a map lookup per element cost 4x the NTTs
    if constexpr (is_integral_v<T> && sizeof(T) == 1){
        array<unsigned, 256> weight_of;
        for (auto& w : weight_of) w = random_weight();
        weight_of[(unsigned char)wildcard] = 0;
        for (int i = 0; i < n; i++) wt[i] = weight_of[(unsigned char)text[i]];
        for (int j = 0; j < m; j++) wp[j] = weight_of[(unsigned char)pattern[m - 1 - j]];
    }
    else{
        vector<pair<T, int>> items;
        items.reserve(n + m);
        for (int i = 0; i < n; i++){
            if (!(text[i] == wildcard)) items.emplace_back(text[i], i);
        }
        for (int j = 0; j < m; j++){
            if (!(pattern[m - 1 - j] == wildcard)) items.emplace_back(pattern[m - 1 - j], n + j);
        }
        sort(items.begin(), items.end(), [](const pair<T, int>& a, const pair<T, int>& b){ return a.first < b.first; });

        unsigned w = 0;
        for (int k = 0; k < (int)items.size(); k++){
            if (k == 0 || items[k - 1].first < items[k].first) w = random_weight();
            int at = items[k].second;
            (at < n ? wt[at] : wp[at - n]) = w;
        }
    }

    /// A cyclic length >= n only wraps the convolution onto indices below m - 1, which are never read
    int len = 1;
    while (len < n) len <<= 1;
    auto rt = roots<MOD, ROOT>(len, false);
    uint64_t len_inv = power(len, MOD - 2, MOD);

    vector<unsigned> score(len);
    for (int k = 1; k <= 3; k++){
        auto a = transformed_powers(wt, k, len, rt), b = transformed_powers(wp, 4 - k, len, rt);
        uint64_t coef = (k == 2 ? MOD - 2 : 1) * len_inv % MOD;
        for (int i = 0; i < len; i++){
            score[i] = (score[i] + (uint64_t)a[i] * b[i] % MOD * coef) % MOD;
        }
    }
    inverse_transform<MOD>(score, roots<MOD, ROOT>(len, true));

    for (int i = m - 1; i < n; i++){
        if (score[i] % MOD == 0) positions.push_back(i - m + 1);
    }

    return positions;
}

vector<int> wildcard_match(const string& text, const string& pattern, char wildcard = '?'){
    return wildcard_match<string>(text, pattern, wildcard);
}

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

    bool read_token(){
        int c = ' ', len = 0;

        while (isspace(c) && c != EOF) c = read_char();
        if (c == EOF) return false;

        for (; !isspace(c) && c != EOF; c = read_char()){
            if (len + 1 < (int)sizeof(tmpbuf)) tmpbuf[len++] = c;
        }
        tmpbuf[len] = 0;
        return true;
    }

    /// Each float parses at its own width, since strtold then a cast to double rounds twice and can miss by one ulp
    bool read_one(float &x){
        if (!read_token()) return false;
        x = strtof(tmpbuf, nullptr);
        return true;
    }

    /// Digits m <= 2^53 times or over an exact 10^k, k <= 22, round once as strtod: 2-3x faster, 6% slower at 17 digits
    bool read_one(double &x){
        if (!read_token()) return false;

        static const double p10[] = {1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11,
                                     1e12, 1e13, 1e14, 1e15, 1e16, 1e17, 1e18, 1e19, 1e20, 1e21, 1e22};
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

        /// x87 arithmetic (FLT_EVAL_METHOD != 0) rounds twice, so it always takes strtod
        if (*s || !digits || digits > 19 || m > 1ULL << 53 || e < -22 || e > 22 || FLT_EVAL_METHOD){
            x = strtod(tmpbuf, nullptr);
            return true;
        }

        double v = e < 0 ? m / p10[-e] : m * p10[e];
        x = *tmpbuf == '-' ? -v : v;
        return true;
    }

    bool read_one(long double &x){
        if (!read_token()) return false;
        x = strtold(tmpbuf, nullptr);
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
    int n, m, r, c, w, i, j;

    read(n, m);
    vector<string> table(n);
    read(table);
    read(r, c);
    vector<string> pat(r);
    read(pat);

    w = m + c - 1;
    string text, pattern;
    for (i = 0; (i + 1) < (n + r); i++){
        for (j = 0; j < w; j++) text += table[i % n][j % m];
    }
    for (i = 0; i < r; i++){
        pattern += pat[i];
        if ((i + 1) < r) pattern += string(w - c, '?');
    }

    vector<string> res(n, string(m, '0'));
    for (int p : wildcard_match(text, pattern)){
        if ((p % w) < m && (p / w) < n) res[p / w][p % w] = '1';
    }

    for (i = 0; i < n; i++) write(res[i]);
    return 0;
}
