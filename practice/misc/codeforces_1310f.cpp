/***
 *
 * Nimber Arithmetic
 * Nim product and inverse of 64-bit nimbers, the field GF(2^64) of Conway's nim arithmetic
 *
 * Complexity: O(1) per nim_product (46 table lookups), 126 products per nim_inverse, 78 KB of tables
 *
 * nim_product(a, b): a times b in nim arithmetic, commutative, associative and distributive over the nim sum a ^ b
 * nim_inverse(a): the x with nim_product(a, x) == 1, a != 0
 *
 * Notes:
 *   - Nimbers below 2^(2^k) are closed under both operations, a subfield: GF(4) below 4, GF(16) below 16
 *   - Products of 2^(2^k) for distinct k are ordinary products, and 2^(2^k) squared is 3 2^(2^k - 1): 4 * 4 = 6
 *   - The tables are built on the first call, about 0.3 ms
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// With x = 2^k, x * x is the nim sum of x and 2^(k - 1), so a product of 2k-bit values needs 3 k-bit products, not 4
struct NimberTables{
    unsigned char small[256][256] = {};
    unsigned long long half[6][256]; /// rows 0..1 for k = 16 and 2..5 for k = 32: byte i of c, times 2^(k - 1)

    NimberTables(){
        small[1][1] = 1;
        for (int h = 2; h <= 8; h *= 2){
            for (int a = 0; a < (1 << h); a++){
                for (int b = 0; b < (1 << h); b++){
                    if (max(a, b) >= (1 << (h / 2))) small[a][b] = split(a, b, h);
                }
            }
        }

        for (int k = 16; k <= 32; k *= 2){
            for (int i = 0; i < k / 8; i++){
                for (int v = 0; v < 256; v++){
                    half[k / 8 - 2 + i][v] = split((unsigned long long)v << (8 * i), 1ULL << (k - 1), k);
                }
            }
        }
    }

    unsigned long long product(unsigned long long a, unsigned long long b, int h) const{
        return h <= 8 ? small[a][b] : split(a, b, h);
    }

    /// a, b < 2^h, and every table entry it reads for width h / 2 is already built
    unsigned long long split(unsigned long long a, unsigned long long b, int h) const{
        int k = h / 2;
        unsigned long long mask = (1ULL << k) - 1;
        unsigned long long ah = a >> k, al = a & mask, bh = b >> k, bl = b & mask;

        unsigned long long low = product(al, bl, k), high = product(ah, bh, k), mid = product(ah ^ al, bh ^ bl, k);
        return ((mid ^ low) << k) ^ low ^ times_half(high, k);
    }

    /// c * 2^(k - 1) for c < 2^k, linear in c, so one lookup per byte of c
    unsigned long long times_half(unsigned long long c, int k) const{
        if (k <= 8) return small[c][1 << (k - 1)];

        unsigned long long res = 0;
        for (int i = 0; i < k / 8; i++) res ^= half[k / 8 - 2 + i][c >> (8 * i) & 255];
        return res;
    }
};

unsigned long long nim_product(unsigned long long a, unsigned long long b){
    static const NimberTables tables;
    return tables.product(a, b, 64);
}

/// The multiplicative group has order 2^64 - 1, so a^(2^64 - 2) is the inverse
unsigned long long nim_inverse(unsigned long long a){
    assert(a != 0);

    unsigned long long res = 1;
    for (int i = 1; i < 64; i++){
        a = nim_product(a, a);
        res = nim_product(res, a);
    }

    return res;
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

const unsigned long long ORDER = ULLONG_MAX;
const unsigned long long primes[] = {3, 5, 17, 257, 641, 65537, 6700417};

unsigned long long power(unsigned long long a, unsigned long long e){
    unsigned long long res = 1;

    while (e){
        if (e & 1) res = nim_product(res, a);
        a = nim_product(a, a);
        e >>= 1;
    }

    return res;
}

long long bsgs(unsigned long long a, unsigned long long b, unsigned long long p){
    unsigned long long i, m = sqrt(p) + 1, cur = 1;
    vector<pair<unsigned long long, unsigned long long>> baby;

    for (i = 0; i < m; i++){
        baby.push_back({cur, i});
        cur = nim_product(cur, a);
    }
    sort(baby.begin(), baby.end());

    unsigned long long step = nim_inverse(cur);
    cur = b;
    for (i = 0; i < m; i++){
        auto it = lower_bound(baby.begin(), baby.end(), make_pair(cur, 0ULL));
        if (it != baby.end() && it->first == cur) return (i * m) + it->second;
        cur = nim_product(cur, step);
    }

    return -1;
}

unsigned long long mod_power(unsigned long long a, unsigned long long e, unsigned long long m){
    unsigned long long res = 1;

    a %= m;
    while (e){
        if (e & 1) res = (res * a) % m;
        a = (a * a) % m;
        e >>= 1;
    }

    return res;
}

bool solve(unsigned long long a, unsigned long long b, unsigned long long &res){
    int i;
    unsigned long long p, q, r;

    res = 0;
    for (i = 0; i < 7; i++){
        p = primes[i], q = ORDER / p;
        unsigned long long ap = power(a, q), bp = power(b, q);
        if (ap == 1 && bp != 1) return false;
        if (ap == 1) continue;

        long long x = bsgs(ap, bp, p);
        if (x == -1) return false;

        r = (x * mod_power(q % p, p - 2, p)) % p;
        res = (res + ((unsigned __int128)q * r)) % ORDER;
    }

    return power(a, res) == b;
}

int main(){
    int t;
    unsigned long long a, b, res;

    read(t);
    while (t--){
        read(a, b);
        if (solve(a, b, res)) write(res);
        else write(-1);
    }

    return 0;
}
