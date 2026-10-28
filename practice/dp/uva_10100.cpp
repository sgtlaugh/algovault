/***
 *
 * Hirschberg LCS
 * One longest common subsequence of two sequences, in linear memory
 *
 * Complexity: O(n m) time, O(n + m) memory
 *
 * lcs_sequence(a, b): one longest common subsequence of a and b, in their Sequence type (a string for strings)
 *
 * Notes:
 *   - Elements need only ==, a and b any random-access containers with push_back
 *   - Only the length: bit_string_lcs.cpp in O(n m / 64), or hunt_szymanski.cpp when few positions match
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename Sequence>
struct Hirschberg{
    const Sequence &a, &b;
    Sequence res;

    Hirschberg(const Sequence& a, const Sequence& b) : a(a), b(b){
        solve(0, a.size(), 0, b.size());
    }

    void solve(int a0, int a1, int b0, int b1){
        if (a0 >= a1 || b0 >= b1) return;
        if (a1 - a0 == 1){
            for (int j = b0; j < b1; j++){
                if (a[a0] == b[j]){
                    res.push_back(a[a0]);
                    return;
                }
            }

            return;
        }

        /// An optimal subsequence splits b where the first half of a ends; the rows are freed before recursing
        int mid = (a0 + a1) / 2, split = b0;
        {
            vector<int> left = forward_row(a0, mid, b0, b1), right = backward_row(mid, a1, b0, b1);
            for (int j = 0, best = -1; j <= b1 - b0; j++){
                if (left[j] + right[j] > best) best = left[j] + right[j], split = b0 + j;
            }
        }

        solve(a0, mid, b0, split);
        solve(mid, a1, split, b1);
    }

    /// row[j] = LCS length of a[a0..a1) and b[b0..b0 + j)
    vector<int> forward_row(int a0, int a1, int b0, int b1) const{
        vector<int> row(b1 - b0 + 1, 0);
        for (int i = a0; i < a1; i++){
            int diag = 0;
            for (int j = b0; j < b1; j++){
                int up = row[j - b0 + 1];
                row[j - b0 + 1] = a[i] == b[j] ? diag + 1 : max(up, row[j - b0]);
                diag = up;
            }
        }

        return row;
    }

    /// row[j] = LCS length of a[a0..a1) and b[b0 + j..b1)
    vector<int> backward_row(int a0, int a1, int b0, int b1) const{
        vector<int> row(b1 - b0 + 1, 0);
        for (int i = a1 - 1; i >= a0; i--){
            int diag = 0;
            for (int j = b1 - 1; j >= b0; j--){
                int down = row[j - b0];
                row[j - b0] = a[i] == b[j] ? diag + 1 : max(down, row[j - b0 + 1]);
                diag = down;
            }
        }

        return row;
    }
};

template <typename Sequence>
Sequence lcs_sequence(const Sequence& a, const Sequence& b){
    /// The rows run over the second sequence, so make it the shorter one
    if (a.size() < b.size()) return Hirschberg<Sequence>(b, a).res;
    return Hirschberg<Sequence>(a, b).res;
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

vector<string> split_words(const string& line){
    vector<string> words;
    string cur;

    for (char c: line){
        if (isalnum((unsigned char)c)) cur.push_back(c);
        else if (!cur.empty()) words.push_back(cur), cur.clear();
    }
    if (!cur.empty()) words.push_back(cur);

    return words;
}

int main(){
    int T = 0;
    char buf[64];
    string a, b;

    while (read_line(a) && read_line(b)){
        if (a.empty() || b.empty()) snprintf(buf, sizeof(buf), "%2d. Blank!", ++T);
        else snprintf(buf, sizeof(buf), "%2d. Length of longest match: %d", ++T, (int)lcs_sequence(split_words(a), split_words(b)).size());
        write(string(buf));
    }

    return 0;
}
