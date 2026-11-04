/***
 *
 * Link-Cut Tree
 * Dynamic rooted forest with link, cut, re-rooting, LCA and path aggregates
 *
 * Complexity: amortized O(log n) per operation
 *
 * LinkCutTree<T> lct(n): nodes 0..n - 1 and no edges, values start at T() and paths fold with +
 * LinkCutTree<T, Combine> lct(n, identity, combine): values start at identity and paths fold with combine
 * lct.cut(u, v): removes edge u-v, false if there is no such edge
 * lct.link(u, v): adds edge u-v with u as a child of v, false and no change if u and v are already connected
 * lct.make_root(u): re-roots u's tree at u
 * lct.set(u, value): sets the value of u
 * lct.connected(u, v): whether u and v are in the same tree
 * lct.find_root(u): root of u's tree
 * lct.lca(u, v): lowest common ancestor under the current root, -1 if u and v are not connected
 * lct.parent(u): parent of u under the current root, -1 for the root
 * lct.query(u, v): combine of the values on the path u..v, both ends included, u and v must be connected
 *
 * Notes:
 *   - combine must be associative and commutative with identity as its neutral element, since re-rooting reverses paths
 *   - link re-roots u's tree at u and v's tree keeps its root, so link(child, parent) builds a rooted forest for lca
 *   - After cut the piece holding the old root keeps it, the other piece is rooted at its endpoint of the edge
 *   - query leaves the root where it was
 *   - Everything is iterative, deep paths are safe
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T, typename Combine = plus<T>>
struct LinkCutTree{
    struct Node{
        /// par of an auxiliary splay tree's root is the path-parent pointer into the splay tree above it
        int ch[2], par;
        bool flip;
        T value, agg;
    };

    T identity;
    Combine combine;
    vector<Node> nodes;

    LinkCutTree(int n, T identity = T(), Combine combine = Combine())
        : identity(identity), combine(combine), nodes(n, Node{{-1, -1}, -1, false, identity, identity}) {}

    bool cut(int u, int v){
        if (parent(u) != v) swap(u, v);
        if (parent(u) != v) return false;

        access(u);
        nodes[nodes[u].ch[0]].par = -1;
        nodes[u].ch[0] = -1;
        pull(u);
        return true;
    }

    bool link(int u, int v){
        if (connected(u, v)) return false;
        make_root(u);
        nodes[u].par = v;
        return true;
    }

    void make_root(int u){
        access(u);
        nodes[u].flip ^= 1;
    }

    void set(int u, T value){
        splay(u);
        nodes[u].value = value;
        pull(u);
    }

    bool connected(int u, int v){
        return find_root(u) == find_root(v);
    }

    int find_root(int u){
        access(u);
        int x = u;
        for (push(x); nodes[x].ch[0] != -1; push(x)) x = nodes[x].ch[0];

        splay(x);
        return x;
    }

    int lca(int u, int v){
        if (!connected(u, v)) return -1;
        access(u);
        return access(v);
    }

    int parent(int u){
        access(u);
        int x = nodes[u].ch[0];
        if (x == -1) return -1;

        for (push(x); nodes[x].ch[1] != -1; push(x)) x = nodes[x].ch[1];
        splay(x);
        return x;
    }

    T query(int u, int v){
        assert(connected(u, v));
        int root = find_root(u);

        make_root(u);
        access(v);
        T res = nodes[v].agg;

        make_root(root);
        return res;
    }

    /// Makes root..x the preferred path, x ends as the root of its splay tree; returns the last node joined
    int access(int x){
        int last = -1;

        for (int y = x; y != -1; y = nodes[y].par){
            splay(y);
            nodes[y].ch[1] = last;
            pull(y);
            last = y;
        }

        splay(x);
        return last;
    }

    void pull(int x){
        nodes[x].agg = combine(combine(agg_of(nodes[x].ch[0]), nodes[x].value), agg_of(nodes[x].ch[1]));
    }

    /// Pushing only g, p, x is enough: a flip pending higher up reverses the whole subtree and commutes with rotations
    void splay(int x){
        push(x);
        while (!is_root(x)){
            int p = nodes[x].par, g = nodes[p].par;
            if (!is_root(p)) push(g);
            push(p);
            push(x);
            if (!is_root(p)) rotate((nodes[g].ch[1] == p) == (nodes[p].ch[1] == x) ? p : x);
            rotate(x);
        }
    }

    void push(int x){
        if (!nodes[x].flip) return;
        swap(nodes[x].ch[0], nodes[x].ch[1]);
        for (int c : nodes[x].ch){
            if (c != -1) nodes[c].flip ^= 1;
        }
        nodes[x].flip = false;
    }

    T agg_of(int x) const{
        return x == -1 ? identity : nodes[x].agg;
    }

    bool is_root(int x) const{
        int p = nodes[x].par;
        return p == -1 || (nodes[p].ch[0] != x && nodes[p].ch[1] != x);
    }

    void rotate(int x){
        int p = nodes[x].par, g = nodes[p].par, d = nodes[p].ch[1] == x, child = nodes[x].ch[!d];

        if (!is_root(p)) nodes[g].ch[nodes[g].ch[1] == p] = x;
        nodes[x].par = g;
        nodes[p].ch[d] = child;
        if (child != -1) nodes[child].par = p;
        nodes[x].ch[!d] = p;
        nodes[p].par = x;

        pull(p);
        pull(x);
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

    /// Digits m <= 2^53 times or over an exact 10^k, k <= 22, round once as strtod: 2-3x faster, 6% slower at 17 digits
    bool read_one(double &x){
        if (!read_token()) return false;

        static const double p10[] = {1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11,
                                     1e12, 1e13, 1e14, 1e15, 1e16, 1e17, 1e18, 1e19, 1e20, 1e21, 1e22};
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

        /// x87 arithmetic (FLT_EVAL_METHOD != 0) rounds twice, so it always takes strtod
        if (*s || !digits || digits > 19 || m > 1ULL << 53 || e < -22 || e > 22 || FLT_EVAL_METHOD){
            x = strtod(tmpbuf, nullptr);
            return true;
        }

        double v = e < 0 ? m / p10[-e] : m * p10[e];
        x = *tmpbuf == '-' ? -v : v;
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
    int n, m, a, b;
    string op;

    read(n, m);
    LinkCutTree<int> lct(n + 1);
    while (m--){
        read(op, a);
        if (op == "cut") lct.cut(a, lct.parent(a));
        else{
            read(b);
            if (op == "link") lct.link(a, b);
            else write(lct.lca(a, b));
        }
    }

    return 0;
}
