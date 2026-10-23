/***
 *
 * Walsh Hadamard Transform
 * xor, or and and convolutions of two vectors, exact or modulo mod
 *
 * Complexity: O(n log n)
 *
 * WalshHadamard fwht: a convolver for any number of calls of any length
 * fwht.xor_convolution(A, B, mod = 0): C where C[k] is the sum of A[i] * B[j] over all i ^ j = k, reduced modulo mod when mod > 0
 * fwht.or_convolution(A, B, mod = 0): the same with i | j = k
 * fwht.and_convolution(A, B, mod = 0): the same with i & j = k
 * WalshHadamard::walsh_transform(ar, n, conv_type, mod = 0): transforms ar[0..n) in place, conv_type WalshHadamard::OR, AND or XOR
 * WalshHadamard::inverse_walsh_transform(ar, n, conv_type, mod = 0): undoes walsh_transform
 * With a mod, the two transforms need every ar[i] already in [0, mod)
 *
 * A and B must have the same length n, and n must be a power of 2
 * Exact mode (mod = 0): no overflow while n^2 * max|A[i]| * max|B[j]| < 2^62
 * Modular mode: 1 <= mod < 2^31, any input values (negative or >= mod are reduced), results in [0, mod)
 * xor needs an odd mod, or and and take any mod
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct WalshHadamard{
    enum ConvType{OR, AND, XOR};

    vector<long long> buffer;  /// reused across calls: allocating a fresh 8 MB vector per call at n = 2^20 costs ~10% in page faults

    vector<long long> xor_convolution(const vector<long long>& A, const vector<long long>& B, long long mod = 0){
        return convolution(A, B, XOR, mod);
    }

    vector<long long> or_convolution(const vector<long long>& A, const vector<long long>& B, long long mod = 0){
        return convolution(A, B, OR, mod);
    }

    vector<long long> and_convolution(const vector<long long>& A, const vector<long long>& B, long long mod = 0){
        return convolution(A, B, AND, mod);
    }

    static void walsh_transform(long long* ar, int n, ConvType conv_type, long long mod = 0){
        dispatch<false>(ar, n, conv_type, mod);
    }

    static void inverse_walsh_transform(long long* ar, int n, ConvType conv_type, long long mod = 0){
        dispatch<true>(ar, n, conv_type, mod);
    }

    vector<long long> convolution(const vector<long long>& A, const vector<long long>& B, ConvType conv_type, long long mod = 0){
        int n = A.size();
        assert(A.size() == B.size() && __builtin_popcount(n) == 1);
        vector<long long> res = A;
        buffer.assign(B.begin(), B.end());
        if (mod){
            for (auto& x : res) if ((x %= mod) < 0) x += mod;
            for (auto& x : buffer) if ((x %= mod) < 0) x += mod;
        }

        walsh_transform(res.data(), n, conv_type, mod);
        walsh_transform(buffer.data(), n, conv_type, mod);
        if (mod) for (int i = 0; i < n; i++) res[i] = res[i] * buffer[i] % mod;
        else for (int i = 0; i < n; i++) res[i] = res[i] * buffer[i];
        inverse_walsh_transform(res.data(), n, conv_type, mod);

        return res;
    }

    /// conv_type and mod != 0 become template arguments so the butterfly loops carry no per-element branches
    template<bool inverse>
    static void dispatch(long long* ar, int n, ConvType conv_type, long long mod){
        assert(mod >= 0 && mod < (1LL << 31) && (conv_type != XOR || !mod || mod % 2 == 1));
        if (conv_type == OR) mod ? transform<OR, inverse, true>(ar, n, mod) : transform<OR, inverse, false>(ar, n, mod);
        if (conv_type == AND) mod ? transform<AND, inverse, true>(ar, n, mod) : transform<AND, inverse, false>(ar, n, mod);
        if (conv_type == XOR) mod ? transform<XOR, inverse, true>(ar, n, mod) : transform<XOR, inverse, false>(ar, n, mod);
    }

    template<ConvType conv_type, bool inverse, bool modular>
    static void transform(long long* ar, int n, long long mod){
        if (n < 2) return;

        int m = n >> 1;
        transform<conv_type, inverse, modular>(ar, m, mod);
        transform<conv_type, inverse, modular>(ar + m, m, mod);

        for (int i = 0; i < m; i++){
            long long x = ar[i], y = ar[i + m];
            if constexpr (conv_type == OR) ar[i + m] = inverse ? sub<modular>(y, x, mod) : add<modular>(x, y, mod);
            if constexpr (conv_type == AND) ar[i] = inverse ? sub<modular>(x, y, mod) : add<modular>(x, y, mod);
            if constexpr (conv_type == XOR){
                ar[i] = add<modular>(x, y, mod), ar[i + m] = sub<modular>(x, y, mod);
                if constexpr (inverse) ar[i] = half<modular>(ar[i], mod), ar[i + m] = half<modular>(ar[i + m], mod);
            }
        }
    }

    /// modular operands are in [0, mod); (s >> 63) & mod adds mod back exactly when s went negative
    template<bool modular>
    static long long add(long long x, long long y, long long mod){
        long long s = x + y;
        if constexpr (modular) s -= mod, s += (s >> 63) & mod;
        return s;
    }

    template<bool modular>
    static long long sub(long long x, long long y, long long mod){
        long long d = x - y;
        if constexpr (modular) d += (d >> 63) & mod;
        return d;
    }

    /// exact mode: v is even, so >> 1 divides; modular mode: mod is odd, so adding it to an odd v makes v + mod even
    template<bool modular>
    static long long half(long long v, long long mod){
        if constexpr (modular) v += (v & 1) * mod;
        return v >> 1;
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

#define MAX 8192
#define LOG 6

int n, mod, a[MAX], b[MAX];
long long hat[MAX][1 << LOG], S[MAX];

void add_frequency(int k){
    int j, d, q = 0;
    bitset<MAX> zero, odd, valid;

    for (j = n + 1; j < MAX; j++) zero[j] = 1;
    for (j = 1; j <= n; j++){
        if (!hat[j][k]) zero[j] = 1;
        else if (hat[j][k] < 0) q ^= 1;
        odd[j] = q;
    }

    valid.set();
    valid >>= (MAX - n);
    for (d = 1; d <= n; d++){
        valid &= ~(zero >> d);
        if (valid.none()) break;
        S[d] += (long long)valid.count() - 2 * (valid & (odd ^ (odd >> d))).count();
    }
}

int main(){
    using namespace fio;
    int i, k;

    read(n, mod);
    for (i = 1; i <= n; i++) read(a[i]);
    for (i = 1; i <= n; i++) read(b[i]);

    for (i = 1; i <= n; i++){
        hat[i][a[i]]++, hat[i][b[i]]++;
        WalshHadamard::walsh_transform(hat[i], 1 << LOG, WalshHadamard::XOR);
    }
    for (k = 0; k < (1 << LOG); k++) add_frequency(k);

    long long pw = 1;
    for (i = 1; i <= n; i++){
        if (i < LOG) write(((S[i] << i) >> LOG) % mod);
        else{
            write((S[i] % mod + mod) % mod * pw % mod);
            pw = pw * 2 % mod;
        }
    }
    return 0;
}
