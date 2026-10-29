/***
 *
 * Simulated Annealing
 * Minimizes by a random walk that accepts worse moves with a cooling probability, as a loop driver or over whole states
 *
 * Complexity: O(s c) for s steps and c the cost of one move, plus one steady_clock read per step in time mode
 *
 * Annealer::by_steps(t_start, t_end, steps, seed = 17): exactly `steps` moves, the same seed gives the same result
 * Annealer::by_time(t_start, t_end, seconds, seed = 17): moves until the time is up, timed from the first running()
 * sa.running(): cools to the next step, false once the steps or the time are used up
 * sa.accept(delta): true when delta <= 0, otherwise with probability exp(-delta / T)
 * simulated_annealing(initial, energy, neighbour, annealer): the lowest energy state it visited
 * sa.rng: the mt19937_64 for the caller's moves
 *
 * Notes:
 *   - T = t_start * (t_end / t_start)^progress, progress = step / steps or elapsed / seconds, t_start >= t_end > 0
 *   - A good t_start passes a typical uphill step about a third of the time, a good t_end the smallest one almost never
 *   - energy(const State&): signed integer or floating value to minimize (negate to maximize), any range is safe
 *   - neighbour(const State&, mt19937_64&): a random nearby state, all randomness drawn from the given rng
 *   - simulated_annealing copies and rescores each state from scratch, use the loop driver for moves with a cheap delta
 *   - Nothing guarantees the optimum, more steps and a move that changes little help
 *   - The loop driver, for TSP with in-place 2-opt moves on a tour p of length cur:
 *         auto sa = Annealer::by_time(2, 0.01, 1.9);
 *         while (sa.running()){
 *             int i = sa.rng() % n, j = sa.rng() % n;
 *             long long d = two_opt_delta(i, j);
 *             if (sa.accept(d)) apply_two_opt(i, j), cur += d, best = min(best, cur);
 *         }
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct Annealer{
    double t_start, t_end, seconds;
    long long steps, step = 0;
    double temperature = 0;
    mt19937_64 rng;
    chrono::steady_clock::time_point start;

    static Annealer by_steps(double t_start, double t_end, long long steps, unsigned long long seed = 17){
        return {t_start, t_end, 0, steps, 0, 0, mt19937_64(seed), {}};
    }

    static Annealer by_time(double t_start, double t_end, double seconds, unsigned long long seed = 17){
        return {t_start, t_end, seconds, 0, 0, 0, mt19937_64(seed), {}};
    }

    bool running(){
        if (step == 0) start = chrono::steady_clock::now();
        double progress = seconds > 0 ? chrono::duration<double>(chrono::steady_clock::now() - start).count() / seconds
                                      : step >= steps ? 1 : (double)step / steps;
        if (progress >= 1) return false;

        temperature = t_start * pow(t_end / t_start, progress), step++;
        return true;
    }

    bool accept(double delta){
        return delta <= 0 || uniform_real_distribution<double>(0, 1)(rng) < exp(-delta / temperature);
    }
};

/// Unsigned subtraction is exact as hi > lo, a double cast first merges values past 2^53, freeing small uphill moves
template<typename E>
double uphill_gap(E lo, E hi){
    if constexpr (is_integral_v<E>) return (double)((make_unsigned_t<E>)hi - (make_unsigned_t<E>)lo);
    else return (double)(hi - lo);
}

template<typename State, typename Energy, typename Neighbour>
State simulated_annealing(State cur, Energy energy, Neighbour neighbour, Annealer sa){
    auto cur_e = energy(cur), best_e = cur_e;
    State best = cur;

    while (sa.running()){
        State next = neighbour(cur, sa.rng);
        auto next_e = energy(next);
        if (next_e <= cur_e || sa.accept(uphill_gap(cur_e, next_e))){
            cur = move(next), cur_e = next_e;
            if (cur_e < best_e) best = cur, best_e = cur_e;
        }
    }

    return best;
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

int n, m;
double binom[10][10];

double wins(const vector<double> &p, int k){
    int i, j, c;
    double dp[10] = {1}, nd[10];

    for (j = 0; j < m; j++){
        for (i = 0; i < n; i++) nd[i] = 0;
        for (i = 0; i < n; i++){
            double pw = 1;
            for (c = 0; (i + c) < n; c++, pw *= p[j]){
                if ((j == k && c > 0) || (j < k && c == 1)) continue;
                nd[i + c] += dp[i] * binom[n - 1 - i][c] * pw;
            }
        }
        for (i = 0; i < n; i++) dp[i] = nd[i];
    }

    return dp[n - 1];
}

double regret(const vector<double> &p){
    double best = 0, worst = 1;

    for (int k = 0; k < m; k++){
        double w = wins(p, k);
        best = max(best, w);
        if (p[k] > 0) worst = min(worst, w);
    }

    return best - worst;
}

int main(){
    int i, j;

    read(n, m);
    for (i = 0; i < 10; i++){
        binom[i][0] = 1;
        for (j = 1; j <= i; j++) binom[i][j] = binom[i - 1][j - 1] + (j < i ? binom[i - 1][j] : 0);
    }

    vector<double> p(m, 1.0 / m), res = p;
    double cur = regret(p), low = cur;

    auto sa = Annealer::by_steps(1e-2, 1e-9, 100000);
    while (m > 1 && sa.running()){
        i = sa.rng() % m, j = sa.rng() % m;
        if (i == j) continue;

        double step = sqrt(sa.temperature);
        vector<double> q = p;
        double amount = min(q[i], step * uniform_real_distribution<double>(0, 1)(sa.rng));
        q[i] -= amount, q[j] += amount;

        double next = regret(q);
        if (sa.accept(next - cur)){
            p = q, cur = next;
            if (cur < low) low = cur, res = p;
        }
    }

    for (i = 0; i < m; i++) write(fixed(res[i], 5));
    return 0;
}
