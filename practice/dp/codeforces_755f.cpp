/***
 *
 * Bounded Subset Sum
 * Which sums up to W are reachable by a sub-multiset of the items, with item ids for a chosen sum
 *
 * Complexity:
 *   - bounded_subset_sums: O(W K / 64 + W) time, O(n + W) memory, K <= 2 min(W, D (1 + log2 W), sqrt(2 S)) copies
 *   - SubsetSum: O(W (D + 1) + n log n) time, O(n + W) memory
 *
 * bounded_subset_sums(items, W): can[s] for s in [0, W], items any range of (value, count) pairs like map<int, int>
 * SubsetSum ss(a, W): sub-multisets of the values a[i] up to W
 * ss.reachable(s): whether some sub-multiset sums to s, for any int s
 * ss.subset(s): distinct ids i with the a[i] summing to s, needs reachable(s), empty for s = 0
 *
 * Notes:
 *   - D is the number of distinct values in [1, W], S the sum of all values with their counts
 *   - D <= sqrt(2 S), so both are O(W sqrt W) when the values sum to O(W), like component sizes summing to n
 *   - Values and counts are non-negative, counts up to 1e18, W >= 0, and equal values add their counts
 *   - Values 0 and values above W are never used
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename Items = vector<pair<int, long long>>>
vector<char> bounded_subset_sums(const Items& items, int W){
    assert(W >= 0);
    vector<long long> cnt(W + 1, 0);
    for (const auto& [v, c] : items){
        assert(v >= 0 && c >= 0);
        if (v > 0 && v <= W) cnt[v] = min<long long>(W, cnt[v] + min<long long>(c, W));
    }

    int words = W / 64 + 1, top = 0;  /// words above top are still all zero
    vector<unsigned long long> bits(words, 0);
    bits[0] = 1;
    for (int v = 1; v <= W; v++){
        if (cnt[v] > 2){
            long long pairs = (cnt[v] - 1) / 2;
            cnt[v] -= 2 * pairs;
            if (v <= W / 2) cnt[2 * v] = min<long long>(W, cnt[2 * v] + pairs);
        }

        int q = v / 64, r = v % 64;
        for (int k = 0; k < cnt[v]; k++){
            int hi = min(words - 1, top + q + 1);
            if (r == 0){
                for (int i = hi; i >= q; i--) bits[i] |= bits[i - q];
            } else{
                for (int i = hi; i > q; i--) bits[i] |= bits[i - q] << r | bits[i - q - 1] >> (64 - r);
                bits[q] |= bits[0] << r;
            }
            top = hi;
        }
    }

    vector<char> can(W + 1);
    for (int s = 0; s <= W; s++) can[s] = bits[s / 64] >> (s % 64) & 1;
    return can;
}

struct SubsetSum{
    vector<int> a, last;
    vector<char> can;

    SubsetSum(const vector<int>& a, int W) : a(a){
        assert(W >= 0);
        last.assign(W + 1, -1);
        can.assign(W + 1, 0);

        vector<int> order;
        for (int i = 0; i < (int)a.size(); i++){
            assert(a[i] >= 0);
            if (a[i] > 0 && a[i] <= W) order.push_back(i);
        }
        sort(order.begin(), order.end(), [&](int i, int j){ return a[i] < a[j]; });

        can[0] = 1;
        for (int lo = 0, hi; lo < (int)order.size(); lo = hi){
            int v = a[order[lo]];
            for (hi = lo; hi < (int)order.size() && a[order[hi]] == v; hi++) {}

            for (int r = 0; r < v; r++){
                int left = 0;
                for (int s = r; s <= W; s += v){
                    if (can[s]) left = hi - lo;
                    else if (left > 0) can[s] = 1, last[s] = order[lo + --left];
                }
            }
        }
    }

    bool reachable(int s) const{
        return s >= 0 && s < (int)can.size() && can[s];
    }

    vector<int> subset(int s) const{
        assert(reachable(s));
        vector<int> ids;
        for (; s > 0; s -= a[ids.back()]) ids.push_back(last[s]);
        return ids;
    }
};

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

#define MAX 1000010

int n, k, p[MAX];
bool visited[MAX];

int main(){
    int i, j, len, pairs = 0;

    read(n, k);
    for (i = 1; i <= n; i++) read(p[i]);

    map<int, int> cnt;
    for (i = 1; i <= n; i++){
        if (visited[i]) continue;
        for (len = 0, j = i; !visited[j]; j = p[j], len++) visited[j] = true;
        cnt[len]++;
        pairs += len >> 1;
    }

    vector<char> can = bounded_subset_sums(cnt, k);
    int low = can[k] ? k : k + 1;
    int high = k <= pairs ? 2 * k : min(n, pairs + k);

    write(low, high);
    return 0;
}
