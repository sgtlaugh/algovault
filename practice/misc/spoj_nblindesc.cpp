/***
 *
 * Local Minimization (Nelder-Mead)
 * A local minimum of a function of a few real variables, smooth or not, with no derivatives
 *
 * Complexity: at most 2500 d^2 calls of f for d variables on the stress test's functions, O(d^2) work per call
 *
 * local_minimum(x, f, step, eps = 1e-12): a local minimum of f found from x, step about the distance to it
 *
 * Notes:
 *   - f(const vector<double>&) returns a double, finite at x; x has at least one coordinate and step > 0
 *   - Ends once 3 passes in a row gain below eps * max(1, |f|), each refining moves to eps * max(1, max |x_i|)
 *   - A larger eps ends sooner, as needed at an infimum f only approaches on a domain edge
 *   - A domain: f returns a huge value like 1e100 outside it, and x starts inside it
 *   - Probabilities summing to 1 without walls: p_i = x_i^2 / sum x_j^2, starting from x = (1, ..., 1)
 *   - Within 1e-9 relative on quadratics to condition 1e6, Rosenbrock, sums of distances, convex functions on a simplex
 *   - Max of distances, the enclosing ball: within 1e-9 in 2 dims, up to 1e-6 off in 3 dims
 *   - Local only: restart from several points when f has several minima
 *   - Rounding blurs a slanted wall like x_0 + x_1 <= 1: a minimum on it can end 1e-2 off, change variables instead
 *   - A value f takes only at exact points, like an expected time once a probability is 0, is missed: search each apart
 *   - Tested up to 16 variables
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template<typename F>
vector<double> nelder_mead(vector<double> x, F& f, double step, double eps){
    int d = x.size();
    vector<vector<double>> p(d + 1, x);
    vector<double> v(d + 1);
    for (int i = 0; i < d; i++) p[i + 1][i] += step;
    for (int i = 0; i <= d; i++) v[i] = f(p[i]);

    /// Capped: wedged in a corner of walls it can shrink for 40000 steps, the next run restarts it fresh
    for (int iter = 0; iter < 100 * d * d; iter++){
        vector<int> order(d + 1);
        iota(order.begin(), order.end(), 0);
        sort(order.begin(), order.end(), [&](int a, int b){ return v[a] < v[b]; });
        int lo = order[0], hi = order[d], second = order[d - 1];

        /// Relative: vertices an ulp apart average to one of themselves, so an absolute size might never be reached
        double size = 0, scale = 1;
        for (int i = 0; i <= d; i++){
            for (int j = 0; j < d; j++) size = max(size, abs(p[i][j] - p[lo][j])), scale = max(scale, abs(p[i][j]));
        }
        if (size < eps * scale) return p[lo];

        vector<double> centroid(d, 0);
        for (int i = 0; i <= d; i++){
            for (int j = 0; i != hi && j < d; j++) centroid[j] += p[i][j] / d;
        }
        auto along = [&](double t){
            vector<double> y(d);
            for (int j = 0; j < d; j++) y[j] = centroid[j] + t * (p[hi][j] - centroid[j]);
            return y;
        };

        vector<double> r = along(-1);
        double fr = f(r);
        if (fr < v[lo]){
            vector<double> e = along(-2);
            double fe = f(e);
            if (fe < fr) p[hi] = e, v[hi] = fe;
            else p[hi] = r, v[hi] = fr;
        }
        else if (fr < v[second]) p[hi] = r, v[hi] = fr;
        else{
            vector<double> c = along(fr < v[hi] ? -0.5 : 0.5);
            double fc = f(c);
            if (fc < min(fr, v[hi])) p[hi] = c, v[hi] = fc;
            else{
                for (int i = 0; i <= d; i++){
                    if (i == lo) continue;
                    for (int j = 0; j < d; j++) p[i][j] = (p[i][j] + p[lo][j]) / 2;
                    v[i] = f(p[i]);
                }
            }
        }
    }

    return p[min_element(v.begin(), v.end()) - v.begin()];
}

/// A simplex squeezed flat against a wall stalls, and only moves between coordinate pairs slide along x_0 + x_1 = c
template<typename F>
vector<double> compass_search(vector<double> x, F& f, double step, double eps){
    int d = x.size();
    double fx = f(x);

    /// Capped: zigzagging down a narrow valley could take millions of rounds, the next Nelder-Mead run follows it
    for (int round = 0; round < 100; round++){
        double scale = 1;
        for (double xi : x) scale = max(scale, abs(xi));
        if (step < eps * scale) break;

        bool moved = false;
        for (int i = 0; i < d; i++){
            for (int j = i; j < d; j++){
                for (double s : {step, -step}){
                    vector<double> y = x;
                    y[i] += s;
                    if (j != i) y[j] -= s;
                    double fy = f(y);
                    if (fy < fx) x = y, fx = fy, moved = true;
                }
            }
        }
        step = moved ? step * 2 : step / 2;
    }

    return x;
}

template<typename F>
vector<double> local_minimum(vector<double> x, F f, double step, double eps = 1e-12){
    assert(!x.empty() && step > 0);
    double fx = f(x);

    /// A simplex stalled on a ridge of a max of functions often moves again when restarted flipped and smaller
    for (int pass = 0, idle = 0; idle < 3; pass++){
        double size = (pass % 2 ? -step : step) / (1 << min(pass, 30));
        vector<double> y = compass_search(nelder_mead(x, f, size, eps), f, step, eps);
        double fy = f(y);
        idle = fy < fx - eps * max(1.0, abs(fx)) ? 0 : idle + 1;
        x = y, fx = fy;
    }

    return x;
}

/***
 *
 * Gaussian Elimination on a Band Matrix
 * Solves grid equations where each cell's equation only uses nearby cells, like expected values of a random walk
 *
 * Complexity: O(n m r^2) time and O(n m r + t) memory for t terms, against O(n^3 m^3) for dense gauss on the grid
 *
 * GaussBand<T> grid(n, m): one unknown x(i, j) per cell of an n x m grid, every equation starts as x(i, j) = 0
 * grid.add_constant(i, j, c): adds c to the right side of the equation of x(i, j)
 * grid.add_term(i, j, k, l, w): adds w * x(k, l) to the right side of the equation of x(i, j), (k, l) may be (i, j)
 * grid.solve(): the n x m grid of x(i, j), once per instance and after every add call
 * grid.probabilities: after solve(), true if the weights were probabilities as in the notes
 *
 * Notes:
 *   - Every equation reads x(i, j) = constant + sum of w * x(k, l), so a cell without calls solves to 0
 *   - r is how far apart a cell and its farthest term are in the cell order, r <= min(n, m) + 1 for 8-neighbour moves
 *   - Probabilities: every w passed is >= 0 and each cell's weights sum to at most 1, the rest its chance to stop
 *   - A chance to stop below 1e-14 counts as 0, so a sum that rounding pushed off 1 never stops the walk
 *   - Probabilities solve within 1e-13 relative in double on walks of up to 120 cells with chances down to 1e-15
 *   - A chance s to stop adds about 1e-16 / s relative error in double, as it is read as 1 - (sum of weights)
 *   - Otherwise no row swaps: safe when each cell's |1 - w on itself| exceeds the sum of |w| on the other cells
 *   - A cell that cannot reach a stop makes the system singular and the answer inf or nan
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <class T>  /// If more precision is required, use __float128 from quadmath.h if supported :)
struct GaussBand{
    static constexpr double STOP_ROUNDING = 1e-14;

    struct Term{
        int i, j, k, l;
        T weight;
    };

    int n, m;
    vector<vector<T>> constant;
    vector<Term> terms;
    bool solved = false;

    bool transposed = false, probabilities = false;
    int reach = 0;
    vector<T> band, rhs, stop;

    GaussBand(int n, int m) : n(n), m(m), constant(n, vector<T>(m, 0)){
        assert(n >= 1 && m >= 1);
    }

    void add_constant(int i, int j, T c){
        assert(!solved && inside(i, j));
        constant[i][j] += c;
    }

    void add_term(int i, int j, int k, int l, T weight){
        assert(!solved && inside(i, j) && inside(k, l));
        terms.push_back({i, j, k, l, weight});
    }

    vector<vector<T>> solve(){
        assert(!solved);
        solved = true;

        transposed = false;
        int row_reach = farthest_term();
        transposed = true;
        int column_reach = farthest_term();

        /// Numbering the cells along the shorter side keeps every term near the diagonal, so the band is narrower
        transposed = column_reach < row_reach;
        reach = min(row_reach, column_reach);

        /// The equation of cell u, x(u) - sum of w * x(v) = constant, becomes row u of the band
        int cells = n * m;
        band.assign(cells * (2 * reach + 1), 0);
        rhs.assign(cells, 0);
        for (int i = 0; i < n; i++){
            for (int j = 0; j < m; j++){
                int u = cell_index(i, j);
                at(u, u) = 1;
                rhs[u] = constant[i][j];
            }
        }

        /// The chance the walk stops at u, kept apart: 1 - (weights left) cancels to noise once it is below rounding
        probabilities = true;
        stop.assign(cells, 1);
        for (auto& term : terms){
            int u = cell_index(term.i, term.j);
            at(u, cell_index(term.k, term.l)) -= term.weight;
            stop[u] -= term.weight;
            if (term.weight < 0) probabilities = false;
        }
        for (int u = 0; u < cells; u++){
            if (stop[u] < -STOP_ROUNDING) probabilities = false;
            if (stop[u] < STOP_ROUNDING) stop[u] = 0;
        }

        /// Forward elimination, a row only reaches the next reach rows and columns
        for (int row = 0; row < cells; row++){
            int last = min(cells - 1, row + reach);

            /// GTH: the pivot sums the chances of leaving the cell, never 1 - (chance of staying), so it cannot cancel
            T pivot = at(row, row);
            if (probabilities){
                pivot = stop[row];
                for (int col = row + 1; col <= last; col++) pivot -= at(row, col);
            }
            for (int col = row + 1; col <= last; col++) at(row, col) /= pivot;
            rhs[row] /= pivot, stop[row] /= pivot;

            for (int below = row + 1; below <= last; below++){
                /// += of -factor is bit for bit -= of factor, and 8% to 14% faster on walks of 40000 to 80000 cells
                T factor = -at(below, row);
                if (factor == 0) continue;  /// walls and cells without calls leave many: 18% faster on SPOJ NBLINDESC
                for (int col = row + 1; col <= last; col++) at(below, col) += factor * at(row, col);
                rhs[below] += factor * rhs[row], stop[below] += factor * stop[row];
            }
        }

        /// Back substitution, every row is already divided by its pivot
        for (int row = cells - 1; row >= 0; row--){
            int last = min(cells - 1, row + reach);
            for (int col = row + 1; col <= last; col++) rhs[row] -= at(row, col) * rhs[col];
        }

        vector<vector<T>> x(n, vector<T>(m));
        for (int i = 0; i < n; i++){
            for (int j = 0; j < m; j++) x[i][j] = rhs[cell_index(i, j)];
        }

        return x;
    }

    bool inside(int i, int j){
        return 0 <= i && i < n && 0 <= j && j < m;
    }

    int farthest_term(){
        int farthest = 0;
        for (auto& term : terms){
            farthest = max(farthest, abs(cell_index(term.i, term.j) - cell_index(term.k, term.l)));
        }
        return farthest;
    }

    /// Cells are numbered row by row, or column by column when transposed
    int cell_index(int i, int j){
        return transposed ? j * n + i : i * m + j;
    }

    /// Row and column of the full matrix, each row stores only the columns within reach of its diagonal
    T& at(int row, int col){
        return band[row * (2 * reach + 1) + col - row + reach];
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

const int dx[] = {-1, 0, 1, 0};
const int dy[] = {0, 1, 0, -1};

int n, m, cells, src, dst, id[8][8], ci[64], cj[64], nxt[64][4], order[64];
bool reached[64], escapes[64];
string grid[8];

bool can_escape(const double* p){
    int i, d, u, v, len = 1;

    for (i = 0; i < cells; i++) reached[i] = escapes[i] = false;
    order[0] = src, reached[src] = true;
    for (i = 0; i < len; i++){
        for (d = 0; d < 4 && order[i] != dst; d++){
            v = nxt[order[i]][d];
            if (p[d] > 0 && !reached[v]) reached[v] = true, order[len++] = v;
        }
    }
    if (!reached[dst]) return false;

    order[0] = dst, escapes[dst] = true, len = 1;
    for (i = 0; i < len; i++){
        u = order[i];
        for (d = 0; d < 4; d++){
            v = nxt[u][d];
            if (reached[v] && !escapes[v] && p[d ^ 2] > 0) escapes[v] = true, order[len++] = v;
        }
    }

    for (i = 0; i < cells; i++){
        if (reached[i] && !escapes[i]) return false;
    }

    return true;
}

double expected_time(const double* p){
    int u, d;

    if (!can_escape(p)) return 1e100;

    GaussBand<double> walk(n, m);
    for (u = 0; u < cells; u++){
        if (!reached[u] || u == dst) continue;

        double out = 0;
        for (d = 0; d < 4; d++){
            if (nxt[u][d] != u) out += p[d];
        }

        walk.add_constant(ci[u], cj[u], 1 / out);
        for (d = 0; d < 4; d++){
            if (p[d] > 0 && nxt[u][d] != u) walk.add_term(ci[u], cj[u], ci[nxt[u][d]], cj[nxt[u][d]], p[d] / out);
        }
    }

    return walk.solve()[ci[src]][cj[src]];
}

double solve(int mask){
    vector<int> dirs;

    for (int d = 0; d < 4; d++){
        if ((mask >> d) & 1) dirs.push_back(d);
    }
    int k = dirs.size();

    auto f = [&](const vector<double>& y){
        double p[4] = {0, 0, 0, 0}, last = 1;
        for (int i = 0; (i + 1) < k; i++){
            if (y[i] <= 0) return 1e100;
            p[dirs[i]] = y[i], last -= y[i];
        }
        if (last <= 0) return 1e100;
        p[dirs[k - 1]] = last;
        return expected_time(p);
    };

    vector<double> y(k - 1, 1.0 / k);
    if (k == 1 || f(y) >= 1e100) return f(y);
    return f(local_minimum(y, f, 0.5 / k, 1e-8));
}

int main(){
    int t, i, j, d, k, l, mask;

    read(t);
    while (t--){
        read(n, m);
        for (i = 0; i < n; i++) read(grid[i]);

        cells = 0;
        for (i = 0; i < n; i++){
            for (j = 0; j < m; j++){
                id[i][j] = -1;
                if (grid[i][j] != '#') ci[cells] = i, cj[cells] = j, id[i][j] = cells++;
            }
        }

        for (i = 0; i < n; i++){
            for (j = 0; j < m; j++){
                if (id[i][j] == -1) continue;
                if (grid[i][j] == '@') src = id[i][j];
                if (grid[i][j] == '*') dst = id[i][j];
                for (d = 0; d < 4; d++){
                    k = i + dx[d], l = j + dy[d];
                    bool open = k >= 0 && k < n && l >= 0 && l < m && id[k][l] != -1;
                    nxt[id[i][j]][d] = open ? id[k][l] : id[i][j];
                }
            }
        }

        double res = 1e100;
        for (mask = 1; mask < 16; mask++) res = min(res, solve(mask));
        write(fixed(res, 10));
    }
    return 0;
}
