/***
 *
 * Treap
 * Randomized balanced binary search tree, as an ordered multiset and as an implicit-key sequence
 *
 * Complexity: expected O(log n) per operation, plus O(k) to build from or insert k values at once
 *
 * Treap<T> set: an empty multiset of keys
 * set.erase(key): removes one copy of key, false if there is none
 * set.insert(key): adds one copy of key
 * set.count(key): how many copies of key
 * set.count_less(key): how many keys are < key
 * set.kth(k): the k-th smallest key, 0-based, 0 <= k < size()
 * set.size(): number of keys, copies included
 *
 * ImplicitTreap<T = long long, Ops = SumAdd<T>, L = T> seq(values): the sequence of values, empty when omitted
 * seq.erase(pos): removes the value at pos, 0 <= pos < size()
 * seq.erase(l, r): removes positions l..r
 * seq.insert(pos, value): inserts value before position pos, 0 <= pos <= size(), pos = size() appends
 * seq.insert(pos, values): inserts the whole vector values before position pos, in order
 * seq.reverse(l, r): reverses the order of positions l..r
 * seq.update(l, r, lz): applies the lazy update lz to positions l..r, adds lz by default
 * seq.get(pos): the value at pos, 0 <= pos < size()
 * seq.max_right(l, pred): largest r with pred(seq.query(l, r)) true, l - 1 if none, 0 <= l <= size()
 * seq.query(l, r): combine over positions l..r, their sum by default
 * seq.size(): length of the sequence
 * seq.to_vector(): the whole sequence in order
 *
 * Notes:
 *   - Positions are 0-based, ranges need 0 <= l <= r < size(), and every aggregate must fit T
 *   - pred must stay false once false as the range grows, it is never called on an empty range
 *   - Ops holds four static functions, like SegmentTree's, L is the lazy type when it differs from T:
 *       - combine(a, b): fold of two adjacent ranges, a on the left
 *       - apply(val, lz, len): value of a range of len elements after lz, len = 1 for a single value
 *       - compose(old, new): one lazy meaning old then new, also called on a stale old whose result is dropped
 *       - reversed(a): the aggregate of the range read backwards, a itself when combine is commutative
 *   - Best subarray sum with reversals: T = {sum, prefix, suffix, best}, reversed swaps prefix and suffix
 *   - Treap is a portable, augmentable alternative to ordered_set.cpp's OrderedMultiset for judges without pb_ds
 *   - Long runs of equal values, up to 4.6e18 positions: run_length_implicit_treap.cpp
 *   - Priorities are seeded from the clock: a fixed seed lets an adversary force an O(n) deep chain
 *   - Erased nodes stay in nodes, so memory grows with every insert: reserve nodes ahead for millions of them
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T>
struct Treap{
    struct Node{
        T key;
        unsigned int priority;
        int size, l, r;
    };

    vector<Node> nodes;
    int root = -1;
    mt19937 rng{(unsigned int)chrono::steady_clock::now().time_since_epoch().count()};

    bool erase(const T& key){
        return erase(root, key);
    }

    bool erase(int& t, const T& key){
        if (t == -1) return false;
        if (!(nodes[t].key < key) && !(key < nodes[t].key)){
            t = merge(nodes[t].l, nodes[t].r);
            return true;
        }
        bool found = erase(key < nodes[t].key ? nodes[t].l : nodes[t].r, key);
        if (found) pull(t);
        return found;
    }

    void insert(const T& key){
        nodes.push_back({key, (unsigned int)rng(), 1, -1, -1});
        int a, b;
        split(root, key, a, b);
        root = merge(merge(a, (int)nodes.size() - 1), b);
    }

    int count(const T& key) const{
        int res = 0;
        for (int t = root; t != -1; ){
            if (key < nodes[t].key) t = nodes[t].l;
            else res += size(nodes[t].l) + 1, t = nodes[t].r;
        }

        return res - count_less(key);
    }

    int count_less(const T& key) const{
        int res = 0;
        for (int t = root; t != -1; ){
            if (nodes[t].key < key) res += size(nodes[t].l) + 1, t = nodes[t].r;
            else t = nodes[t].l;
        }

        return res;
    }

    T kth(int k) const{
        assert(0 <= k && k < size());
        for (int t = root; ; ){
            int left = size(nodes[t].l);
            if (k < left) t = nodes[t].l;
            else if (k == left) return nodes[t].key;
            else k -= left + 1, t = nodes[t].r;
        }
    }

    int size() const{
        return size(root);
    }

    int size(int t) const{
        return t == -1 ? 0 : nodes[t].size;
    }

    int merge(int a, int b){
        if (a == -1 || b == -1) return a == -1 ? b : a;
        if (nodes[a].priority > nodes[b].priority){
            nodes[a].r = merge(nodes[a].r, b);
            pull(a);
            return a;
        }
        nodes[b].l = merge(a, nodes[b].l);
        pull(b);
        return b;
    }

    void pull(int t){
        nodes[t].size = 1 + size(nodes[t].l) + size(nodes[t].r);
    }

    /// Splits t into keys < key and keys >= key
    void split(int t, const T& key, int& a, int& b){
        if (t == -1){
            a = b = -1;
            return;
        }

        if (nodes[t].key < key){
            split(nodes[t].r, key, nodes[t].r, b);
            a = t;
        }
        else{
            split(nodes[t].l, key, a, nodes[t].l);
            b = t;
        }
        pull(t);
    }
};

/// Range add and range sum, the default Ops of ImplicitTreap
template <typename T>
struct SumAdd{
    static T combine(T a, T b){ return a + b; }
    static T apply(T val, T lz, int len){ return val + lz * len; }
    static T compose(T old_lz, T new_lz){ return old_lz + new_lz; }
    static T reversed(T a){ return a; }
};

template <typename T = long long, typename Ops = SumAdd<T>, typename L = T>
struct ImplicitTreap{
    struct Node{
        T value, agg;
        L lazy;
        unsigned int priority;
        int size, l, r;
        bool has_lazy, flip;
    };

    vector<Node> nodes;
    int root = -1;
    mt19937 rng{(unsigned int)chrono::steady_clock::now().time_since_epoch().count()};

    ImplicitTreap(const vector<T>& values = {}){
        root = build(values);
    }

    void erase(int pos){
        erase(pos, pos);
    }

    void erase(int l, int r){
        assert(0 <= l && l <= r && r < size());
        int a, b, c;
        split(root, l, a, b);
        split(b, r - l + 1, b, c);
        root = merge(a, c);
    }

    void insert(int pos, const T& value){
        assert(0 <= pos && pos <= size());
        int a, b;
        split(root, pos, a, b);
        root = merge(merge(a, new_node(value)), b);
    }

    void insert(int pos, const vector<T>& values){
        assert(0 <= pos && pos <= size());
        int a, b;
        split(root, pos, a, b);
        root = merge(merge(a, build(values)), b);
    }

    void reverse(int l, int r){
        on_range(l, r, [&](int t){ apply_flip(t); });
    }

    void update(int l, int r, const L& lz){
        on_range(l, r, [&](int t){ apply_lazy(t, lz); });
    }

    T get(int pos){
        return query(pos, pos);
    }

    template <typename F>
    int max_right(int l, F pred){
        assert(0 <= l && l <= size());
        int a, b, res = l - 1;
        split(root, l, a, b);

        bool empty = true;
        T acc{};
        for (int t = b; t != -1; ){
            push(t);
            int c = nodes[t].l;
            if (c != -1){
                T nxt = empty ? nodes[c].agg : Ops::combine(acc, nodes[c].agg);
                if (!pred(nxt)){
                    t = c;
                    continue;
                }
                acc = nxt, empty = false, res += nodes[c].size;
            }

            T nxt = empty ? nodes[t].value : Ops::combine(acc, nodes[t].value);
            if (!pred(nxt)) break;
            acc = nxt, empty = false, res++;
            t = nodes[t].r;
        }

        root = merge(a, b);
        return res;
    }

    T query(int l, int r){
        T res{};
        on_range(l, r, [&](int t){ res = nodes[t].agg; });
        return res;
    }

    int size() const{
        return size(root);
    }

    int size(int t) const{
        return t == -1 ? 0 : nodes[t].size;
    }

    vector<T> to_vector(){
        vector<T> res;
        vector<int> stack;

        for (int t = root; t != -1 || !stack.empty(); ){
            for (; t != -1; t = nodes[t].l) push(t), stack.push_back(t);
            t = stack.back();
            stack.pop_back();
            res.push_back(nodes[t].value);
            t = nodes[t].r;
        }

        return res;
    }

    /// Cartesian tree of values on fresh random priorities, built left to right with a stack of the right spine
    int build(const vector<T>& values){
        vector<int> stack;

        for (const T& v : values){
            int t = new_node(v), last = -1;
            while (!stack.empty() && nodes[stack.back()].priority < nodes[t].priority){
                last = stack.back();
                stack.pop_back();
                pull(last);
            }
            nodes[t].l = last;
            if (!stack.empty()) nodes[stack.back()].r = t;
            stack.push_back(t);
        }

        int res = stack.empty() ? -1 : stack[0];
        while (!stack.empty()){
            pull(stack.back());
            stack.pop_back();
        }

        return res;
    }

    int new_node(const T& v){
        nodes.push_back({v, v, L{}, (unsigned int)rng(), 1, -1, -1, false, false});
        return nodes.size() - 1;
    }

    void pull(int t){
        Node& node = nodes[t];
        T agg = node.value;
        node.size = 1;
        if (node.l != -1) agg = Ops::combine(nodes[node.l].agg, agg), node.size += nodes[node.l].size;
        if (node.r != -1) agg = Ops::combine(agg, nodes[node.r].agg), node.size += nodes[node.r].size;
        node.agg = agg;
    }

    template <typename F>
    void on_range(int l, int r, F f){
        assert(0 <= l && l <= r && r < size());
        int a, b, c;
        split(root, l, a, b);
        split(b, r - l + 1, b, c);
        f(b);
        root = merge(merge(a, b), c);
    }

    void apply_lazy(int t, const L& lz){
        if (t == -1) return;
        Node& node = nodes[t];
        node.value = Ops::apply(node.value, lz, 1);
        node.agg = Ops::apply(node.agg, lz, node.size);
        /// compose also runs on a stale lazy and the result is dropped: this select beat a branch by 12%
        L composed = Ops::compose(node.lazy, lz);
        node.lazy = node.has_lazy ? composed : lz;
        node.has_lazy = true;
    }

    /// The children swap at once and agg is read backwards, so a flipped node's own fields are already correct
    void apply_flip(int t){
        if (t == -1) return;
        swap(nodes[t].l, nodes[t].r);
        nodes[t].agg = Ops::reversed(nodes[t].agg);
        nodes[t].flip ^= 1;
    }

    /// Splits t into the first k elements and the rest
    void split(int t, int k, int& a, int& b){
        if (t == -1){
            a = b = -1;
            return;
        }

        push(t);
        if (size(nodes[t].l) < k){
            split(nodes[t].r, k - size(nodes[t].l) - 1, nodes[t].r, b);
            a = t;
        }
        else{
            split(nodes[t].l, k, a, nodes[t].l);
            b = t;
        }
        pull(t);
    }

    int merge(int a, int b){
        if (a == -1 || b == -1) return a == -1 ? b : a;
        if (nodes[a].priority > nodes[b].priority){
            push(a);
            nodes[a].r = merge(nodes[a].r, b);
            pull(a);
            return a;
        }
        push(b);
        nodes[b].l = merge(a, nodes[b].l);
        pull(b);
        return b;
    }

    void push(int t){
        Node& node = nodes[t];
        if (!(node.has_lazy | node.flip)) return;

        int l = node.l, r = node.r;
        if (node.has_lazy){
            apply_lazy(l, node.lazy);
            apply_lazy(r, node.lazy);
            node.has_lazy = false;
        }
        if (node.flip){
            apply_flip(l);
            apply_flip(r);
            node.flip = false;
        }
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

struct Node{
    long long sum, pre, suf, best;
};

struct MaxSubarray{
    static Node combine(Node a, Node b){
        Node res;
        res.sum = a.sum + b.sum;
        res.pre = max(a.pre, a.sum + b.pre);
        res.suf = max(b.suf, b.sum + a.suf);
        res.best = max(max(a.best, b.best), a.suf + b.pre);
        return res;
    }
    static Node apply(Node val, int, int){ return val; }
    static int compose(int, int new_lz){ return new_lz; }
    static Node reversed(Node a){ return a; }
};

Node leaf(long long x){
    return {x, x, x, x};
}

int main(){
    int n, q, i, x, y;
    char op;

    read(n);
    vector<Node> ar(n);
    for (i = 0; i < n; i++){
        read(x);
        ar[i] = leaf(x);
    }

    ImplicitTreap<Node, MaxSubarray, int> seq(ar);
    read(q);
    while (q--){
        read(op, x);
        x--;
        if (op == 'D'){
            seq.erase(x);
            continue;
        }

        read(y);
        if (op == 'I') seq.insert(x, leaf(y));
        else if (op == 'R') seq.erase(x), seq.insert(x, leaf(y));
        else write(seq.query(x, y - 1).best);
    }

    return 0;
}
