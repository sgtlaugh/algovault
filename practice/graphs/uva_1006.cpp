/***
 *
 * Hungarian Algorithm (Kuhn-Munkres)
 * Matches rows to distinct columns of a cost matrix so the sum of the chosen cells is minimum or maximum
 *
 * Complexity: O(min(n, m)^2 max(n, m)) for n rows and m columns
 *
 * hungarian(mat, minimize = true): {total, pairs}, an optimal assignment of the n x m vector<vector<T>> mat
 *
 * Notes:
 *   - pairs has min(n, m) (row, col) pairs sorted by row, no row or column twice, total sums their cells
 *   - minimize = false maximizes the total instead
 *   - T is signed (int, long long or double) and 2 max(n, m) max |mat[i][j]| must fit in T
 *   - An empty mat returns {0, {}}
 *   - With ties, a tall matrix may get a different optimal matching than its transpose, the total is the same
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <class T>
pair<T, vector<pair<int, int>>> hungarian(vector<vector<T>> mat, bool minimize=true){
    if (!mat.size()) return {(T)0, {}};

    int n = mat.size(), m = mat[0].size();
    for (int i = 0; i < n; i++){
        assert((int)mat[i].size() == m);
    }

    if (n > m){
        vector<vector<T>> t(m, vector<T>(n));
        for (int i = 0; i < n; i++){
            for (int j = 0; j < m; j++) t[j][i] = mat[i][j];
        }
        auto res = hungarian(t, minimize);
        for (auto& p: res.second) swap(p.first, p.second);
        sort(res.second.begin(), res.second.end());
        return res;
    }

    if (!minimize){
        for (auto& r: mat){
            for (auto& x: r) x = -x;
        }
    }

    vector<int> way(m + 1, 0), match(n + 1, 0), visited(m + 1, 0), row_of_col(m + 1, 0);
    vector<T> row_potential(n + 1, 0), col_potential(m + 1, 0), minv(m + 1, 0);

    for (int i = 1; i <= n; i++){
        int col = 0;
        row_of_col[0] = i;
        for (int j = 0; j <= m; j++){
            minv[j] = numeric_limits<T>::max(), visited[j] = 0;
        }
        do{
            T w = numeric_limits<T>::max();
            int row = row_of_col[col], next = 0;
            visited[col] = 1;

            for (int j = 1; j <= m; j++){
                if (!visited[j]){
                    T r = mat[row - 1][j - 1] - row_potential[row] - col_potential[j];

                    if (r < minv[j]) minv[j] = r, way[j] = col;
                    if (minv[j] < w) w = minv[j], next = j;
                }
            }

            for (int j = 0; j <= m; j++){
                if (visited[j]) row_potential[row_of_col[j]] += w, col_potential[j] -= w;
                else minv[j] -= w;
            }
            col = next;
        } while (row_of_col[col]);

        do{
            int prev = way[col];
            row_of_col[col] = row_of_col[prev], col = prev;
        } while (col);
    }

    vector<pair<int, int>> matches;

    for (int j = 1; j <= m; j++) match[row_of_col[j]] = j;
    for (int i = 1; i <= n; i++) matches.push_back({i - 1, match[i] - 1});

    return {minimize ? -col_potential[0] : col_potential[0], matches};
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

        /// Accumulate as negative so that the minimum value (e.g. INT_MIN) does not overflow
        for (x = 0; isdigit(c); c = read_char()){
            x = x * 10 - (c - '0');
        }
        if (!neg) x = -x;

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

const long long INF = 1000000000;

int m, n, region[12], cnt[55], S[55][12], Tm[55][12], run_time[55][12], where[55], pos[55], start_at[55];

int main(){
    int i, j, k, T = 0;
    char buf[128];

    while (read(m, n) && (m || n)){
        for (j = 0; j < m; j++) read(region[j]);
        for (i = 0; i < n; i++){
            read(cnt[i]);
            for (k = 0; k < cnt[i]; k++) read(S[i][k], Tm[i][k]);
            for (j = 0; j < m; j++){
                run_time[i][j] = -1;
                for (k = 0; k < cnt[i]; k++){
                    if (S[i][k] <= region[j]) run_time[i][j] = Tm[i][k];
                }
            }
        }

        vector<vector<long long>> mat(n, vector<long long>(m * n));
        for (i = 0; i < n; i++){
            for (j = 0; j < m; j++){
                for (k = 0; k < n; k++) mat[i][j * n + k] = (run_time[i][j] == -1) ? INF : (long long)(k + 1) * run_time[i][j];
            }
        }
        pair<long long, vector<pair<int, int>>> res = hungarian(mat);

        for (i = 0; i < n; i++) where[i] = res.second[i].second / n, pos[i] = res.second[i].second % n;
        for (j = 0; j < m; j++){
            vector<pair<int, int>> order;
            for (i = 0; i < n; i++){
                if (where[i] == j) order.push_back(make_pair(-pos[i], i));
            }
            sort(order.begin(), order.end());

            int now = 0;
            for (k = 0; k < (int)order.size(); k++){
                i = order[k].second;
                start_at[i] = now, now += run_time[i][j];
            }
        }

        snprintf(buf, sizeof(buf), "Case %d", ++T);
        write(string(buf));
        snprintf(buf, sizeof(buf), "Average turnaround time = %.2f", (double)res.first / n);
        write(string(buf));
        for (i = 0; i < n; i++){
            snprintf(buf, sizeof(buf), "Program %d runs in region %d from %d to %d", i + 1, where[i] + 1, start_at[i], start_at[i] + run_time[i][where[i]]);
            write(string(buf));
        }
        write("");
    }

    return 0;
}
