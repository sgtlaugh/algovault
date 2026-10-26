/***
 *
 * Continued Fractions
 * Continued fraction expansion, convergents, best rational approximation and fraction binary search
 *
 * Complexity: O(log q) for the expansion, convergents and best approximation, O(log N) predicate calls for the search
 *
 * best_approximation(p, q, N): the fraction closest to p / q among all fractions with denominator <= N
 *     Ties go to the smaller denominator, then to the smaller value
 *     Requires p > LLONG_MIN, q >= 1, N >= 1 and __int128
 *
 * continued_fraction(p, q): [a0; a1, ..., ak] with p / q = a0 + 1 / (a1 + 1 / (... + 1 / ak))
 *     a0 = floor(p / q) can be negative, a1..ak >= 1 and ak >= 2 when k >= 1, so the expansion is unique
 *     p / q need not be reduced, requires q >= 1, any p works
 *
 * convergents(a): the fractions [a0; a1, ..., ai] for every i, reduced, the last one is p / q reduced
 *     Requires a non-empty expansion, never overflows on continued_fraction(p, q) for p > LLONG_MIN
 *
 * fraction_search(f, N): the smallest fraction p / q with 0 <= p <= N and 1 <= q <= N for which f(p, q) is true
 *     f must be monotone over the non-negative rationals: false up to some point, true from there on
 *     Returns {1, 0} if f is false for every such fraction
 *     f is only called with reduced fractions inside the bounds, at most 4 log2(N) + 10 times, N <= 1e18
 *
 * Notes:
 *   - Fractions are {numerator, denominator} pairs, the denominator is always positive
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

vector<long long> continued_fraction(long long p, long long q){
    assert(q > 0);
    vector<long long> a;

    while (q){
        long long d = p / q, r = p % q;
        if (r < 0) d--, r += q;
        a.push_back(d);
        p = q, q = r;
    }

    return a;
}

pair<long long, long long> best_approximation(long long p, long long q, long long n){
    assert(q > 0 && n > 0);
    auto a = continued_fraction(p, q);
    long long pp = 1, pq = 0, cp = a[0], cq = 1;

    for (size_t i = 1; i < a.size(); i++){
        if (a[i] <= (n - pq) / cq){
            tie(pp, pq, cp, cq) = make_tuple(cp, cq, a[i] * cp + pp, a[i] * cq + pq);
            continue;
        }

        /// The convergent and the largest semiconvergent are the neighbours of p / q among denominators <= n
        long long t = (n - pq) / cq, sp = pp + t * cp, sq = pq + t * cq;
        __int128 dc = (__int128)p * cq - (__int128)cp * q, ds = (__int128)p * sq - (__int128)sp * q;
        if (dc < 0) dc = -dc;
        if (ds < 0) ds = -ds;

        __int128 lhs = dc * sq, rhs = ds * cq;
        if (lhs != rhs) return lhs < rhs ? make_pair(cp, cq) : make_pair(sp, sq);
        if (cq != sq) return cq < sq ? make_pair(cp, cq) : make_pair(sp, sq);
        return (__int128)cp * sq < (__int128)sp * cq ? make_pair(cp, cq) : make_pair(sp, sq);
    }

    return {cp, cq};
}

vector<pair<long long, long long>> convergents(const vector<long long>& a){
    assert(!a.empty());
    vector<pair<long long, long long>> res;
    long long pp = 1, pq = 0, cp = a[0], cq = 1;
    res.push_back({cp, cq});

    for (size_t i = 1; i < a.size(); i++){
        tie(pp, pq, cp, cq) = make_tuple(cp, cq, a[i] * cp + pp, a[i] * cq + pq);
        res.push_back({cp, cq});
    }

    return res;
}

template<typename F>
pair<long long, long long> fraction_search(F f, long long n){
    assert(n > 0);
    if (f(0LL, 1LL)) return {0, 1};

    /// Invariant: f(lp / lq) is false, f(rp / rq) is true, with 1 / 0 standing for infinity
    long long lp = 0, lq = 1, rp = 1, rq = 0;
    while (lp + rp <= n && lq + rq <= n){
        bool side = f(lp + rp, lq + rq);
        long long& bp = side ? rp : lp;
        long long& bq = side ? rq : lq;
        long long dp = side ? lp : rp, dq = side ? lq : rq;

        auto valid = [&](long long k){
            if (dp && k > (n - bp) / dp) return false;
            if (dq && k > (n - bq) / dq) return false;
            return f(bp + k * dp, bq + k * dq) == side;
        };

        long long k = 1, step = 1;
        while (valid(k + step)) k += step, step *= 2;
        for (step /= 2; step; step /= 2){
            if (valid(k + step)) k += step;
        }

        bp += k * dp, bq += k * dq;
    }

    return {rp, rq};
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

const long long INF = LLONG_MAX;

char buf[64];
long long lp, lq, rp, rq;

bool inside(pair<long long, long long> c){
    return (__int128)lp * c.second < (__int128)c.first * lq && (__int128)c.first * rq < (__int128)rp * c.second;
}

vector<long long> other(vector<long long> a){
    a.back()--;
    a.push_back(1);
    return a;
}

pair<long long, long long> candidate(const vector<long long>& x, const vector<long long>& y){
    vector<long long> a;

    for (size_t i = 0; ; i++){
        long long u = i < x.size() ? x[i] : INF, v = i < y.size() ? y[i] : INF;
        if (u != v){
            a.push_back(min(u, v) + 1);
            break;
        }
        a.push_back(u);
    }

    return convergents(a).back();
}

pair<long long, long long> simplest(){
    vector<long long> x = continued_fraction(lp, lq), y = continued_fraction(rp, rq);
    pair<long long, long long> res = candidate(x, y);

    if (inside(res)) return res;

    vector<pair<long long, long long>> c = {candidate(other(x), y), candidate(x, other(y)), candidate(other(x), other(y))};
    res = {-1, -1};
    for (int i = 0; i < 3; i++){
        if (inside(c[i]) && (res.second == -1 || c[i].second < res.second)) res = c[i];
    }

    return res;
}

int main(){
    int t;
    long long p, x;

    while (read(t)){
        while (t--){
            read(p, x);
            lp = x - 1, lq = p, rp = x, rq = p;

            pair<long long, long long> res = simplest();
            long long a = (__int128)res.second * x - (__int128)res.first * p;
            snprintf(buf, sizeof(buf), "%lld/%lld", a, res.second);
            write(string(buf));
        }
    }
    return 0;
}
