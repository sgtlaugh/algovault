/***
 *
 * Longest Increasing Subsequence
 * LIS and LDS lengths ending at every index, plus one longest subsequence as indices
 *
 * Complexity: O(n log n)
 *
 * lis_vector(a, allow_equal, cmp)[i]: length of the longest increasing subsequence ending at index i
 * lis_indices(a, allow_equal, cmp): increasing indices of one longest increasing subsequence, empty for empty input
 * lis_length(a, allow_equal): length of the longest increasing subsequence
 * lds_vector(a, allow_equal)[i]: length of the longest decreasing subsequence ending at index i
 * lds_length(a, allow_equal): length of the longest decreasing subsequence
 *
 * Notes:
 *   - cmp is a strict weak ordering that defines "increasing", less<T> by default
 *   - allow_equal = true gives the non-decreasing and non-increasing versions
 *   - lis_indices(a, false, greater<T>()) gives one longest decreasing subsequence
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T, typename Compare = less<T>>
vector<int> lis_vector(const vector<T>& ar, bool allow_equal=false, Compare cmp=Compare()){
    vector<T> tails;
    vector<int> res;
    auto before = [&](const T& tail, const T& x){ return allow_equal ? !cmp(x, tail) : cmp(tail, x); };

    for (const T &x : ar){
        if (tails.empty() || before(tails.back(), x)){
            tails.push_back(x);
            res.push_back(tails.size());
            continue;
        }

        /// Branchless search in [0, size - 1], the back is not before x; 2.5x faster than lower_bound at random n = 1e6
        const T* first = tails.data();
        size_t n = tails.size();
        while (n > 1){
            size_t half = n / 2;
            first += half * before(first[half - 1], x);
            n -= half;
        }

        res.push_back(first - tails.data() + 1);
        tails[first - tails.data()] = x;
    }

    return res;
}

/// The rightmost i before j with len[i] == len[j] - 1 always precedes ar[j], so no comparator check is needed
template <typename T, typename Compare = less<T>>
vector<int> lis_indices(const vector<T>& ar, bool allow_equal=false, Compare cmp=Compare()){
    if (ar.empty()) return {};
    vector<int> len = lis_vector(ar, allow_equal, cmp);
    int last = max_element(len.begin(), len.end()) - len.begin();
    vector<int> res(len[last]);

    res.back() = last;
    for (int i = last - 1, k = len[last] - 2; k >= 0; i--){
        if (len[i] == k + 1) res[k--] = i;
    }

    return res;
}

template <typename T>
vector<int> lds_vector(const vector<T>& ar, bool allow_equal=false){
    return lis_vector(ar, allow_equal, [](const T& a, const T& b){ return b < a; });
}

template <typename T>
int lis_length(const vector<T>& ar, bool allow_equal=false){
    auto lis = lis_vector(ar, allow_equal);
    return lis.empty() ? 0 : *max_element(lis.begin(), lis.end());
}

template <typename T>
int lds_length(const vector<T>& ar, bool allow_equal=false){
    auto lds = lds_vector(ar, allow_equal);
    return lds.empty() ? 0 : *max_element(lds.begin(), lds.end());
}

/***
 *
 * Garsia-Wachs
 * Minimum total cost of merging adjacent piles into one, a merge costing the merged size
 *
 * Complexity: expected O(n log n)
 *
 * garsia_wachs(weights): minimum total merge cost, 0 for at most one pile
 *
 * Notes:
 *   - The same value is the optimal alphabetic tree: min sum weights[i] * depth(i) over trees with the leaves in order
 *   - Weights are non-negative integers with total * max(1, ceil(log2 n)) < LLONG_MAX, which bounds the answer
 *   - Merge costs that are not the merged size: knuth_optimization.cpp
 *   - https://en.wikipedia.org/wiki/Garsia%E2%80%93Wachs_algorithm
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct GarsiaWachs{
    struct Node{
        long long val, max_val;
        unsigned int priority;
        int l, r;
    };

    /// A wall above every merged pile stops the leftward walk, and end above every weight merges the last pairs
    static constexpr long long wall = LLONG_MAX, end = LLONG_MAX - 1;

    vector<Node> nodes;
    mt19937 rng{(unsigned int)chrono::steady_clock::now().time_since_epoch().count()};
    long long cost = 0;

    template <typename T>
    GarsiaWachs(const vector<T>& weights){
        nodes.reserve(weights.size() + 3);
        vector<int> rest = {node(end)};
        for (int i = (int)weights.size() - 1; i >= 0; i--) rest.push_back(node(weights[i]));

        /// The sequence is root then a, b; it never holds s[i] <= s[i + 2], which would merge s[i] and s[i + 1]
        int root = -1, a = node(wall), b = node(wall);
        while (!rest.empty()){
            int v = rest.back();
            rest.pop_back();

            /// Children mark a run of smaller piles put back whole; a run of one is compared as a pile and never merges
            if (nodes[v].l != -1 || nodes[v].r != -1){
                int last = pop_last(v), before_last = pop_last(v);
                root = merge(merge(merge(root, a), b), v);
                a = before_last, b = last;
            }
            else if (nodes[a].val <= nodes[v].val){
                nodes[a].val = nodes[a].max_val = nodes[a].val + nodes[b].val;
                cost += nodes[a].val;

                /// The merged pile moves left past the smaller piles, which need no new comparison except with v
                int smaller = -1;
                split_after_last_at_least(root, nodes[a].val, root, smaller);
                rest.push_back(v);
                if (smaller != -1) rest.push_back(smaller);
                rest.push_back(a);
                b = pop_last(root), a = pop_last(root);
            }
            else{
                root = merge(root, a);
                a = b, b = v;
            }
        }
    }

    int node(long long x){
        nodes.push_back({x, x, (unsigned int)rng(), -1, -1});
        return nodes.size() - 1;
    }

    int merge(int x, int y){
        if (x == -1 || y == -1) return x == -1 ? y : x;
        if (nodes[x].priority > nodes[y].priority){
            nodes[x].r = merge(nodes[x].r, y);
            pull(x);
            return x;
        }
        nodes[y].l = merge(x, nodes[y].l);
        pull(y);
        return y;
    }

    /// The returned node has no children, or the main loop would read it as a run of piles
    int pop_last(int& t){
        if (nodes[t].r == -1){
            int last = t;
            t = nodes[t].l;
            nodes[last].l = -1, nodes[last].max_val = nodes[last].val;
            return last;
        }
        int last = pop_last(nodes[t].r);
        pull(t);
        return last;
    }

    void pull(int t){
        nodes[t].max_val = nodes[t].val;
        if (nodes[t].l != -1) nodes[t].max_val = max(nodes[t].max_val, nodes[nodes[t].l].max_val);
        if (nodes[t].r != -1) nodes[t].max_val = max(nodes[t].max_val, nodes[nodes[t].r].max_val);
    }

    void split_after_last_at_least(int t, long long at_least, int& x, int& y){
        if (t == -1){
            x = y = -1;
            return;
        }

        if (nodes[t].r != -1 && nodes[nodes[t].r].max_val >= at_least){
            split_after_last_at_least(nodes[t].r, at_least, nodes[t].r, y);
            x = t;
        }
        else if (nodes[t].val >= at_least){
            y = nodes[t].r, nodes[t].r = -1;
            x = t;
        }
        else{
            split_after_last_at_least(nodes[t].l, at_least, x, nodes[t].l);
            y = t;
        }
        pull(t);
    }
};

template <typename T>
long long garsia_wachs(const vector<T>& weights){
    return GarsiaWachs(weights).cost;
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

int main(){
    int t, n, i;

    read(t);
    for (int tc = 1; tc <= t; tc++){
        read(n);
        vector<int> a(n);
        read(a);

        vector<int> left = lis_vector(a), right = lis_vector(vector<int>(a.rbegin(), a.rend()));
        vector<long long> power(n);
        for (i = 0; i < n; i++) power[i] = left[i] + right[n - 1 - i] - 1;

        write("Case " + to_string(tc) + ": " + to_string(garsia_wachs(power)));
    }

    return 0;
}
