/***
 *
 * Anti Polynomial Hash
 * Builds two different strings with the same polynomial hash, to hack solutions or break weak tests
 *
 * Complexity, with chunk = 15 and lim = 1 + log2(n / 60) rounded down, at least 1:
 *   - AntiHash: O(n log n)
 *   - solve: O(2^chunk chunk + n + chunk 2^lim) per randomized trial, trials repeat until a collision, with no bound
 *   - 512 KB in a fixed array inside the struct, plus 16 n lim bytes, about 600 KB at n = 1000
 *
 * AntiHash hack(base, mod, n): the attack on h = (h * base + s[i]) % mod, n >= 30, base + mod < 2^63, both positive
 * hack.solve(): {s, t}, two different n-character strings of '0' and '1' with equal hashes
 *
 * Notes:
 *   - solve averaged under 10 ms at n = 1000 for mod <= 1e14 and at n = 2000 for mod < 2^61
 *   - Too small an n for mod slows solve down fast: n = 1000 with mod near 2^61 took 2 s on average and up to 9 s
 *   - The strings differ by a +-1 combination of base powers, so they collide whatever values '0' and '1' hash to
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct AntiHash{
    /// Magic constant! Each trial meets in the middle over 2^chunk subset sums

    const static int chunk = 15;

    int n, lim, len, found, sample[2 * chunk];
    long long base, mod;
    string res[2];
    pair<long long, int> ar[1 << chunk];
    vector<vector<pair<long long, int>>> dp;

    AntiHash(long long base, long long mod, int n) : n(n), base(base), mod(mod){
        assert(n >= 2 * chunk);
        generate();
    }

    inline long long fast_modmul(long long a, long long b, long long m){
        if (a >= m) a %= m;
        if (b >= m) b %= m;
        if (m < (long long)UINT_MAX) return (uint64_t)a * b % m;

        #ifdef __SIZEOF_INT128__
            return __int128(a) * b % m;
        #endif

        long double x = (long double)a * b;
        long long c = x / m;

        a = (uint64_t)a * b - (uint64_t)c * m;
        if (a >= m) a -= m;
        if (a < 0) a += m;

        return a;
    }

    inline long long get_hash(const string& s){
        long long h = 0;
        for (int i = 0; i < (int)s.length(); i++){
            h = fast_modmul(h, base, mod);
            if (s[i] == '0') h = (h + 17) % mod;
            if (s[i] == '1') h = (h + 1000000009) % mod;
        }

        return h;
    }

    inline void generate(){
        res[0] = res[1] = string(n, '0');

        /// Halving stops before the last level drops under 4 chunk values, as the tuned n = 2000, lim = 6 left it
        for (lim = 1, len = n; len / 2 >= 4 * chunk; lim++) len /= 2;
        dp.assign(lim, vector<pair<long long, int>>(n));

        len = 0, found = 0;
        dp[0][len++] = make_pair(1, n - 1);
        for (int i = 1; i < n; i++){
            dp[0][len++] = make_pair(fast_modmul(dp[0][i - 1].first, base, mod), n - i - 1);
        }
        sort(dp[0].begin(), dp[0].begin() + len);

        for (int i = 1; i < lim; i++){
            for (int j = 0, k = 0; k + 1 < len; j++, k += 2){
                long long x = (dp[i - 1][k + 1].first - dp[i - 1][k].first + mod) % mod;
                dp[i][j] = make_pair(x, j);
            }
            len >>= 1;
            sort(dp[i].begin(), dp[i].begin() + len);
        }
    }

    inline void dfs(unsigned int l, int i, int pos){
        if (l == 0) res[pos][dp[0][i].second] = '1';
        else{
            dfs(l - 1, dp[l][i].second * 2 + 1, pos);
            dfs(l - 1, dp[l][i].second * 2, pos ^ 1);
        }
    }

    inline void insert_hash(int i, long long sum, long long mask){
        if (i == chunk) ar[mask] = make_pair(sum, mask);
        else{
            insert_hash(i + 1, sum, mask);
            insert_hash(i + 1, (sum + dp[lim - 1][sample[i]].first) % mod, mask | (1LL << i));
        }
    }

    inline void find_hash(int i, long long sum, long long mask){
        if (found) return;
        if (i == chunk){
            int idx = lower_bound(ar + 1, ar + (1 << chunk), pair<long long, int>(sum, 0)) - ar;

            if (idx < (1 << chunk) && ar[idx].first == sum){
                mask = (mask << (long long)chunk) | ar[idx].second;
                for (i = 0; i < 2 * chunk; i++){
                    if (mask & (1LL << i)){
                        dfs(lim - 1, sample[i], (i < chunk));
                    }
                }

                if (get_hash(res[0]) == get_hash(res[1]) && res[0] != res[1]) found = 1;
                else {
                    for (i = 0; i < n; i++) res[0][i] = res[1][i] = '0';
                }

                return;
            }
        }
        else{
            find_hash(i + 1, sum, mask);
            find_hash(i + 1, (sum + dp[lim - 1][sample[i + chunk]].first) % mod, mask | (1LL << i));
        }
    }

    pair<string, string> solve(){
        while (!found){
            vector<int> v;
            for (int i = 0; i < len; i++){
                v.push_back(i < (2 * chunk));
            }
            shuffle(v.begin(), v.end(), mt19937(rand()));
            for (int i = 0, k = 0; i < (int)v.size(); i++){
                if (v[i] == 1) sample[k++] = i;
            }
            shuffle(sample, sample + chunk * 2, mt19937(rand()));

            insert_hash(0, 0, 0);
            sort(ar, ar + (1 << chunk));
            find_hash(0, 0, 0);
        }

        return {res[0], res[1]};
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
    long long b, m;

    freopen("hash.in", "r", stdin);
    freopen("hash.out", "w", stdout);

    read(b, m);
    AntiHash hack(b, m, 1000);
    pair<string, string> res = hack.solve();

    write(res.first);
    write(res.second);
    return 0;
}
