/***
 *
 * Dancing Links
 * Solves exact cover with Knuth's Algorithm X: select a subset of rows covering every column exactly once
 *
 * Complexity: exponential worst case (exact cover is NP-complete), O(columns + total row length) memory
 *
 * DancingLinks dlx(ncolumns, max_nodes = 0): columns 0 to ncolumns - 1, max_nodes = total row length presizes arrays
 * dlx.add_row(r, columns): adds a row with any int id r over distinct columns, an empty row is ignored
 * dlx.count(limit = LLONG_MAX): number of exact covers, or limit if there are more, searching no further, limit >= 0
 * dlx.exact_cover(): the selected row ids, nullopt if no cover exists, an empty cover for 0 columns
 *
 * Notes:
 *   - Rows beyond max_nodes still fit after a reallocation
 *   - Deterministic: the same rows added in the same order give the same cover
 *   - count(2) == 1 tells a unique cover, like a sudoku with one solution, without searching for all of them
 *   - count and exact_cover restore the links, so they can be called again or rows can still be added
 *   - Recursion depth is at most ncolumns, an 8 MB stack overflows near 5e4 (4e4 under sanitizers)
 *   - sudoku::solution(grid) below reduces an n x n sudoku (n a perfect square) to exact cover with 4 n^2 columns
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct DancingLinks{
    int n, idx, cover_size = 0;
    vector<int> row, col, L, R, U, D, selected_rows, column_count;

    DancingLinks(int ncolumns, int max_nodes = 0) : n(ncolumns), idx(ncolumns + 1), selected_rows(ncolumns + 1), column_count(ncolumns + 1){
        for (auto* v : {&row, &col, &L, &R, &U, &D}) v->resize(n + 1 + max_nodes);

        for (int i = 0; i <= n; i++) U[i] = D[i] = i, L[i] = i - 1, R[i] = i + 1;
        L[0] = n, R[n] = 0;
    }

    void add_row(int r, const vector<int>& columns){
        int first = idx, l = columns.size();

        if (!l) return;  /// an empty row covers nothing, and linking it would corrupt the previous row
        if (idx + l > (int)row.size()){
            for (auto* v : {&row, &col, &L, &R, &U, &D}) v->resize(2 * (idx + l));
        }

        for (int column : columns){
            assert(0 <= column && column < n);
            int c = column + 1;  /// node 0 is the header of the column list
            L[idx] = idx - 1, R[idx] = idx + 1, D[idx] = c, U[idx] = U[c];
            D[U[c]] = idx, U[c] = idx, row[idx] = r, col[idx] = c;
            column_count[c]++, idx++;
        }
        R[idx - 1] = first, L[first] = idx - 1;
    }

    long long count(long long limit = LLONG_MAX){
        assert(limit >= 0);
        return limit ? algorithm_x(0, limit) : 0;
    }

    optional<vector<int>> exact_cover(){
        if (!algorithm_x(0, 1)) return nullopt;
        return vector<int>(selected_rows.begin(), selected_rows.begin() + cover_size);
    }

    /// Counts covers up to limit >= 1, the last one found stays in selected_rows[0, cover_size)
    long long algorithm_x(int depth, long long limit){
        if (R[0] == 0){
            cover_size = depth;
            return 1;
        }

        int c = R[0];
        /// Stops at a count <= 1: the choice is forced, and a 0 count column left behind still fails the branch
        for (int i = R[0]; i != 0 && column_count[c] > 1; i = R[i]){
            if (column_count[i] < column_count[c]) c = i;
        }

        remove(c);
        long long found = 0;
        for (int i = D[c]; i != c && found < limit; i = D[i]){
            selected_rows[depth] = row[i];
            for (int j = R[i]; j != i; j = R[j]) remove(col[j]);
            found += algorithm_x(depth + 1, limit - found);
            for (int j = L[i]; j != i; j = L[j]) restore(col[j]);
        }

        restore(c);
        return found;
    }

    void remove(int c){
        int l = L[c], r = R[c];

        L[r] = l, R[l] = r;

        for (int i = D[c]; i != c; i = D[i]){
            for (int j = R[i]; j != i; j = R[j]){
                int u = U[j], d = D[j];
                column_count[col[j]]--;
                U[d] = u, D[u] = d;
            }
        }
    }

    void restore(int c){
        for (int i = U[c]; i != c; i = U[i]){
            for (int j = L[i]; j != i; j = L[j]){
                int u = U[j], d = D[j];
                column_count[col[j]]++;
                U[d] = j, D[u] = j;
            }
        }

        int l = L[c], r = R[c];
        L[r] = c, R[l] = c;
    }
};

namespace sudoku{
    int encode(int n, int a, int b, int c){
        return a * n * n + b * n + c;
    }

    void decode(int n, int v, int& a, int& b, int& c){
        c = v % n, v /= n;
        b = v % n, a = (v / n) % n;
    }

    /***
     * Returns the solved grid, nullopt if the puzzle is unsolvable; the puzzle itself is not changed
     * Zeros are blanks, the grid must be n x n with n a perfect square
     *
     * The exact cover search is exponential in the worst case
     * 9 x 9 grids, well constrained 16 x 16 grids and even empty 49 x 49 grids solve fast,
     * but sparse 16 x 16 grids (around 75% blank) can take over a minute
     *
    ***/
    optional<vector<vector<int>>> solution(const vector<vector<int>>& grid){
        int i, j, k, n = grid.size(), m = sqrt(n + 0.5);

        assert(m * m == n);
        for (i = 0; i < n; i++) assert((int)grid[i].size() == n);

        int nrows = 0;
        for (auto& r : grid) for (int x : r) nrows += x ? 1 : n;

        DancingLinks dlx(4 * n * n, 4 * nrows);  /// n * n columns each for cells, rows, columns and boxes
        for (i = 0; i < n; i++){
            for (j = 0; j < n; j++){
                for (k = 0; k < n; k++){
                    if (grid[i][j] == 0 || grid[i][j] == (k + 1)){
                        dlx.add_row(encode(n, i, j, k), {encode(n, 0, i, j), encode(n, 1, i, k), encode(n, 2, j, k), encode(n, 3, (i / m) * m + j / m, k)});
                    }
                }
            }
        }

        auto res = dlx.exact_cover();
        if (!res) return nullopt;

        vector<vector<int>> solved = grid;
        for (int v : *res){
            decode(n, v, i, j, k);
            solved[i][j] = k + 1;
        }

        return solved;
    }
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

#define MAXA 36

typedef vector<pair<int, int>> shape;

const int dx[] = {0, 0, -1, 1};
const int dy[] = {-1, 1, 0, 0};

vector<shape> ominoes[7];
vector<int> kind[7];
map<array<int, 3>, bool> memo;
char buf[64];

shape normalize(shape s){
    int low_x = INT_MAX, low_y = INT_MAX;

    for (auto& p : s) low_x = min(low_x, p.first), low_y = min(low_y, p.second);
    for (auto& p : s) p.first -= low_x, p.second -= low_y;
    sort(s.begin(), s.end());
    return s;
}

shape canonical(shape s){
    int i;
    shape res = normalize(s);

    for (i = 0; i < 8; i++){
        for (auto& p : s) p = make_pair(p.second, -p.first);
        if (i == 3){
            for (auto& p : s) p.first = -p.first;
        }
        res = min(res, normalize(s));
    }

    return res;
}

void grow(const shape& s, set<shape>& next){
    int d;

    for (auto& p : s){
        for (d = 0; d < 4; d++){
            shape t = s;
            t.push_back(make_pair(p.first + dx[d], p.second + dy[d]));
            if (find(s.begin(), s.end(), t.back()) == s.end()) next.insert(normalize(t));
        }
    }
}

void gen(){
    int x;
    set<shape> cur = {shape(1, make_pair(0, 0))};

    for (x = 1; x <= 6; x++){
        map<shape, int> id;
        for (auto& s : cur){
            ominoes[x].push_back(s);
            shape c = canonical(s);
            if (!id.count(c)){
                int k = id.size();
                id[c] = k;
            }
            kind[x].push_back(id[c]);
        }

        set<shape> next;
        for (auto& s : cur) grow(s, next);
        cur = next;
    }
}

void add(DancingLinks& dlx, int c, const shape& s, int a, int b, int mark){
    vector<int> columns;

    for (auto& p : s) columns.push_back((a + p.first) * c + b + p.second);
    dlx.add_row(0, columns);
    if (mark < 0) return;

    columns.push_back(mark);
    dlx.add_row(0, columns);
}

bool fits(int x, int r, int c, int k){
    int i, a, b, h, w;
    DancingLinks dlx(r * c + (k >= 0));

    for (i = 0; i < (int)ominoes[x].size(); i++){
        h = w = 0;
        for (auto& p : ominoes[x][i]) h = max(h, p.first + 1), w = max(w, p.second + 1);

        for (a = 0; (a + h) <= r; a++){
            for (b = 0; (b + w) <= c; b++) add(dlx, c, ominoes[x][i], a, b, kind[x][i] == k ? r * c : -1);
        }
    }

    return dlx.count(1) > 0;
}

bool tileable(int x, int r, int c){
    if (!r || !c) return true;
    if ((r * c) % x) return false;

    return fits(x, r, c, -1);
}

bool beats(int x, int r, int c, int k){
    int a, b;

    for (a = 1; a <= r; a++){
        for (b = 1; b <= c && (a * b) <= MAXA; b++){
            if ((a * b) % x || !fits(x, a, b, k)) continue;
            if (tileable(x, a, c - b) && tileable(x, r - a, c)) return true;
            if (tileable(x, r, c - b) && tileable(x, r - a, b)) return true;
        }
    }

    return fits(x, r, c, k);
}

bool gabriel(int x, int r, int c){
    if (x >= 7 || (r * c) % x) return false;
    if (r > c) swap(r, c);
    if (memo.count({x, r, c})) return memo[{x, r, c}];

    int k, kinds = *max_element(kind[x].begin(), kind[x].end()) + 1;
    bool res = true;
    for (k = 0; k < kinds && res; k++) res = beats(x, r, c, k);

    return (memo[{x, r, c}] = res);
}

int main(){
    int t, x, r, c, tc = 0;

    gen();
    read(t);
    while (t--){
        read(x, r, c);
        snprintf(buf, sizeof(buf), "Case #%d: %s", ++tc, gabriel(x, r, c) ? "GABRIEL" : "RICHARD");
        write(string(buf));
    }

    return 0;
}
