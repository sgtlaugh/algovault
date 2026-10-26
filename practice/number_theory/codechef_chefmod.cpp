/***
 *
 * Discrete Logarithm
 * Generalized discrete logarithm with Shanks' baby step giant step algorithm
 *
 * Complexity: O(sqrt(mod)) expected
 *
 * discrete_log(a, b, mod): the smallest x >= 0 with a^x = b (modulo mod), -1 if no solution exists
 *
 * Notes:
 *   - a and b are any int (negative ones too), 1 <= mod <= INT_MAX, a, b and mod need not be pairwise co-prime
 *   - a^0 = 1 (including 0^0), so x = 0 whenever b = 1 (modulo mod), which is every b when mod = 1
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// The table masks the low bits, and residues like b * 3^i mod 2^30 all share them
struct mix_hash{
    size_t operator()(int x) const{
        return (unsigned long long)x * 0x9E3779B97F4A7C15ULL >> 32;
    }
};

int discrete_log(int a, int b, int mod){
    int n = sqrt(mod + 0.5) + 1;

    auto expo = [mod](long long base, int k){
        long long res = 1 % mod;
        for (; k; k >>= 1, base = base * base % mod){
            if (k & 1) res = res * base % mod;
        }
        return (int)res;
    };

    a %= mod, b %= mod;
    if (a < 0) a += mod;
    if (b < 0) b += mod;

    /// Linear probing at most half full, -1 marks an empty slot: 51% faster than gp_hash_table on Library Checker data
    int cap = 2 << __lg(2 * n + 5);
    vector<int> keys(cap, -1), vals(cap);
    auto slot = [&](int key){
        size_t h = mix_hash()(key) & (cap - 1);
        while (keys[h] != -1 && keys[h] != key) h = (h + 1) & (cap - 1);
        return h;
    };

    int cur = 1 % mod;
    for (int i = 0; i < n + 3; i++){
        if (cur == b) return i;
        int key = (long long)b * cur % mod;
        size_t h = slot(key);
        keys[h] = key, vals[h] = i;
        cur = (long long)cur * a % mod;
    }

    int giant = expo(a, n);
    cur = giant;
    for (int i = 2; i < n + 3; i++){
        cur = (long long)cur * giant % mod;
        size_t h = slot(cur);
        if (keys[h] == cur){
            int x = (((long long)n * i - vals[h]) % mod + mod) % mod;
            return expo(a, x) == b ? x : -1;
        }
    }

    return -1;
}

/***
 *
 * Extended GCD
 * Extended gcd, modular inverse, and linear diophantine equations (solve and count solutions in a box)
 *
 * Complexity: O(log m) for mod_inverse, O(log min(|a|, |b|)) for the others
 *
 * count_diophantine(a, b, c, x1, x2, y1, y2): number of (x, y) with a * x + b * y = c in [x1, x2] x [y1, y2]
 * diophantine(a, b, c): some (x, y) with a * x + b * y = c, nullopt if there is none
 * extended_gcd(a, b): {g, x, y} with g = gcd(a, b) and a * x + b * y = g, for a, b >= 0
 * mod_inverse(a, m): inverse of a modulo m in [0, m), -1 if gcd(a, m) != 1, for any a and 1 <= m < 2^63
 *
 * Notes:
 *   - count_diophantine: arguments in [-1e9, 1e9], an empty range gives 0, a = b = c = 0 counts the whole box
 *   - diophantine: a, b, c may be negative, |a|, |b|, |c| <= 1e9 keeps x and y within long long
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// BEGIN SHARED extended_gcd
array<long long, 3> extended_gcd(long long a, long long b){
    assert(a >= 0 && b >= 0);

    long long x = 1, y = 0, x1 = 0, y1 = 1;
    while (b){
        long long q = a / b;
        tie(x, x1) = make_pair(x1, x - q * x1);
        tie(y, y1) = make_pair(y1, y - q * y1);
        tie(a, b) = make_pair(b, a - q * b);
    }
    return {a, x, y};
}

long long mod_inverse(long long a, long long m){
    assert(m >= 1);

    a %= m;
    if (a < 0) a += m;

    auto [g, x, y] = extended_gcd(a, m);
    if (g != 1) return -1;
    return x < 0 ? x + m : x;
}
/// END SHARED extended_gcd

optional<pair<long long, long long>> diophantine(long long a, long long b, long long c){
    if (a == 0 && b == 0){
        if (c) return nullopt;
        return pair(0LL, 0LL);
    }

    auto [g, x, y] = extended_gcd(llabs(a), llabs(b));
    if (c % g) return nullopt;

    x *= c / g, y *= c / g;
    if (a < 0) x = -x;
    if (b < 0) y = -y;
    return pair(x, y);
}

long long count_diophantine(long long a, long long b, long long c, long long x1, long long x2, long long y1, long long y2){
    auto sol = diophantine(a, b, c);
    if (!sol) return 0;
    if (a == 0 && b == 0) return max(0LL, x2 - x1 + 1) * max(0LL, y2 - y1 + 1);

    /// shrinks [lo, hi] to the k with lo_val <= base + k * step <= hi_val
    long long g = gcd(llabs(a), llabs(b)), lo = LLONG_MIN, hi = LLONG_MAX;
    auto restrict_steps = [&](long long base, long long step, long long lo_val, long long hi_val){
        if (step == 0){
            if (base < lo_val || base > hi_val) lo = 1, hi = 0;
            return;
        }

        if (step < 0) step = -step, base = -base, tie(lo_val, hi_val) = make_pair(-hi_val, -lo_val);
        auto floor_div = [&](long long num){ return num / step - (num % step < 0); };
        lo = max(lo, -floor_div(base - lo_val));
        hi = min(hi, floor_div(hi_val - base));
    };
    restrict_steps(sol->first, b / g, x1, x2);
    restrict_steps(sol->second, -a / g, y1, y2);
    return max(0LL, hi - lo + 1);
}

/***
 *
 * Fast I/O
 * Buffered stdin and stdout through fread and fwrite, about 8x faster than scanf for integers and strings
 *
 * Complexity: O(1) amortized per character, one fread or fwrite call per 8192 bytes
 *
 * fio::flush(): writes the pending output now, which also happens at normal exit
 * fio::read(a, b, ...): reads whitespace separated values into a, b, ..., returns the count read before the input ended
 * fio::read_line(s): reads one line into s without its newline or a trailing \r, false only at the end of input
 * fio::write(a, b, ...): writes a, b, ... space separated and ends the line, a vector argument on its own line
 * fio::fixed(x, digits): a floating x printed with exactly digits digits after the point, rounded like printf("%.*f")
 *
 * Notes:
 *   - read takes any integer including __int128, float, double, long double, strings, chars and vectors of these
 *   - write takes integers, strings, chars, vectors of these and fio::fixed values
 *   - Once fio handles a stream, use nothing else on it: scanf, printf, cin or cout would reorder the input or output
 *   - An integer read skips everything up to the next digit or -, words included, where scanf fails and reads nothing
 *   - A char reads one non-space character
 *   - A vector reads exactly its current size, so read a length first: read(n); a.resize(n); read(a);
 *   - Writing a float without fixed() does not compile, so its precision is always chosen
 *   - Floats print at printf's speed, a floating token read keeps its first 127 characters
 *   - Output is lost on abnormal termination such as abort() or _exit()
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// Global buffers on purpose: a second reader would split stdin between two buffers, a second writer reorder stdout
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

    bool read_one(double &x){
        if (!read_token()) return false;
        x = strtod(tmpbuf, nullptr);
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

    void flush(){
        fwrite(outbuf, 1, outptr, stdout);
        outptr = 0;
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

    void write_one(const string &s){
        for (auto &&c: s) write_char(c);
    }

    /// write(x) for a floating x does not compile on purpose: the digits after the point must be chosen
    struct Fixed{
        long double value;
        int digits;
    };

    void write_one(const Fixed &f){
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

#define LIM 3000

const int P = 100000007;

int g, spf[LIM], lg[LIM];

long long expo(long long x, long long k){
    long long res = 1;
    for (; k; k >>= 1, x = x * x % P){
        if (k & 1) res = res * x % P;
    }
    return res;
}

int primitive_root(){
    int i, x, m = P - 1;
    vector<int> factors;

    for (i = 2; (i * i) <= m; i++){
        if (m % i) continue;
        factors.push_back(i);
        while (m % i == 0) m /= i;
    }
    if (m > 1) factors.push_back(m);

    for (x = 2; ; x++){
        bool flag = true;
        for (int f: factors) flag &= (expo(x, (P - 1) / f) != 1);
        if (flag) return x;
    }
}

void gen(){
    int i, j;

    g = primitive_root();
    for (i = 2; i < LIM; i++){
        if (spf[i]) continue;
        for (j = i; j < LIM; j += i){
            if (!spf[j]) spf[j] = i;
        }
    }

    for (i = 2; i < LIM; i++){
        if (spf[i] == i) lg[i] = discrete_log(g, i, P);
        else lg[i] = (lg[spf[i]] + lg[i / spf[i]]) % (P - 1);
    }
}

long long solve(int A){
    int v, a, d, n;
    long long i, inv, res = 0, cur = A, low = P, best = LLONG_MAX;

    for (; low >= LIM; cur = cur * A % P){
        if (cur < low) res += cur, low = cur;
    }

    a = discrete_log(g, A, P), d = gcd(a, P - 1), n = (P - 1) / d;
    inv = mod_inverse(a / d, n);
    for (v = 1; v < low; v++){
        if (lg[v] % d) continue;
        i = (lg[v] / d) * inv % n;
        if (!i) i = n;
        if (i < best) best = i, res += v;
    }

    return res;
}

int main(){
    int t, A;

    gen();
    read(t);
    while (t--){
        read(A);
        write(solve(A));
    }
    return 0;
}
