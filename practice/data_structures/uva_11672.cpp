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
 * ImplicitTreap seq(values): the sequence of long long values, empty when values is omitted
 * seq.add(l, r, delta): adds delta to positions l..r
 * seq.erase(pos): removes the value at pos, 0 <= pos < size()
 * seq.erase(l, r): removes positions l..r
 * seq.insert(pos, value): inserts value before position pos, 0 <= pos <= size(), pos = size() appends
 * seq.insert(pos, values): inserts the whole vector values before position pos, in order
 * seq.reverse(l, r): reverses the order of positions l..r
 * seq.get(pos): the value at pos, 0 <= pos < size()
 * seq.size(): length of the sequence
 * seq.sum(l, r): sum of positions l..r
 * seq.to_vector(): the whole sequence in order
 *
 * Notes:
 *   - Positions are 0-based, ranges need 0 <= l <= r < size() and sums must fit in long long
 *   - Treap is a portable, augmentable alternative to ordered_set.cpp's OrderedMultiset for judges without pb_ds
 *   - For other values, aggregates or lazy updates use run_length_implicit_treap.cpp and its Runs struct
 *   - Priorities are seeded from the clock: a fixed seed lets an adversary force an O(n) deep chain
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

struct ImplicitTreap{
    struct Node{
        long long value, total, lazy;
        unsigned int priority;
        int size, l, r;
        bool flip;
    };

    vector<Node> nodes;
    int root = -1;
    mt19937 rng{(unsigned int)chrono::steady_clock::now().time_since_epoch().count()};

    ImplicitTreap(const vector<long long>& values = {}){
        root = build(values);
    }

    void add(int l, int r, long long delta){
        on_range(l, r, [&](int t){ apply(t, delta, false); });
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

    void insert(int pos, long long value){
        assert(0 <= pos && pos <= size());
        int a, b;
        split(root, pos, a, b);
        root = merge(merge(a, new_node(value)), b);
    }

    void insert(int pos, const vector<long long>& values){
        assert(0 <= pos && pos <= size());
        int a, b;
        split(root, pos, a, b);
        root = merge(merge(a, build(values)), b);
    }

    void reverse(int l, int r){
        on_range(l, r, [&](int t){ apply(t, 0, true); });
    }

    long long get(int pos){
        return sum(pos, pos);
    }

    int size() const{
        return size(root);
    }

    int size(int t) const{
        return t == -1 ? 0 : nodes[t].size;
    }

    long long sum(int l, int r){
        long long res = 0;
        on_range(l, r, [&](int t){ res = nodes[t].total; });
        return res;
    }

    vector<long long> to_vector(){
        vector<long long> res;
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
    int build(const vector<long long>& values){
        vector<int> stack;
        for (long long v : values){
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

    int new_node(long long v){
        nodes.push_back({v, v, 0, (unsigned int)rng(), 1, -1, -1, false});
        return nodes.size() - 1;
    }

    void pull(int t){
        int l = nodes[t].l, r = nodes[t].r;
        nodes[t].size = 1 + size(l) + size(r);
        nodes[t].total = nodes[t].value + (l == -1 ? 0 : nodes[l].total) + (r == -1 ? 0 : nodes[r].total);
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

    void apply(int t, long long delta, bool flip){
        if (t == -1) return;
        nodes[t].value += delta, nodes[t].total += delta * nodes[t].size, nodes[t].lazy += delta;
        if (flip) swap(nodes[t].l, nodes[t].r), nodes[t].flip ^= 1;
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
        if (nodes[t].lazy || nodes[t].flip){
            apply(nodes[t].l, nodes[t].lazy, nodes[t].flip);
            apply(nodes[t].r, nodes[t].lazy, nodes[t].flip);
            nodes[t].lazy = 0, nodes[t].flip = false;
        }
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

int main(){
    int q, scenario = 0, i, j, k;
    long long n, r, m, a, c;
    char buf[32];
    string cmd;

    while (read(q) && q){
        snprintf(buf, sizeof(buf), "Scenario %d:", ++scenario);
        write(string(buf));
        ImplicitTreap seq;
        seq.nodes.reserve(2000000);
        while (q--){
            read(cmd, i);
            if (cmd[0] == 'i'){
                read(n, r, m, a, c);
                vector<long long> b(n);
                for (k = 0; k < r; k++) read(b[k]);
                for (k = r; k < n; k++){
                    long long x = b[k - r] * a + c;
                    b[k] = x % m, c = x / m % m;
                }
                seq.insert(i, b);
            }
            else{
                read(j);
                if (cmd[0] == 's') write(seq.sum(i, j));
                else seq.erase(i, j);
            }
        }
        write(string());
    }

    return 0;
}
