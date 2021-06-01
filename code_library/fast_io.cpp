/***
 *
 * Fast i/o with fread and fwrite
 *
 * Function read() can take multiple arguments and returns the number of objects read
 * Works for strings or any integral types (int, long long, __int128 etc)
 * Also works for vector, but the type has to be string or integral
 *
 * Function read_line() reads one line at a time
 * It does not skip empty lines
 *
 * Similarly, write() can take multiple arguments, which can be string or integral
 * The write functions prints them separated by a single space, ending with a newline
 *
 * The write function also works for vectors just like read
 * Vector elements are separated by a single space and ends with a newline
 *
 * If multiple vector arguments are passed however, write will print them on separate lines
 * In other words we can do write(v1, v2, v3) and they will be printed on separate lines
 * It is discouraged to mix vectors with other types in the same write call because of formatting difference
 *
 * Don't forget to flush the buffer at the end after writing
 *
 * Example usage - Read four numbers till EOF and print their sum
 *
 * int a, c;
 * long long b;
 * __int128 d;
 *
 * while (read(a, b, c, d)){
 *     write(a + b + c + d);
 * }
 *
 * flush();
 *
***/

#include <bits/stdc++.h>

using namespace std;

/// Because is_integral is not always true for __int128
#ifdef __SIZEOF_INT128__
    namespace std {
	    template<> struct is_integral<__int128_t>: true_type {};
	}
#endif

namespace fio{
    const int BUF_SIZE = 8192;

    int buf_len = 0, inptr = 0, outptr = 0;
    char inbuf[BUF_SIZE], outbuf[BUF_SIZE], tmpbuf[128];

    inline char read_char(){
        if (inptr >= buf_len){
            inptr = 0, buf_len = fread(inbuf, 1, BUF_SIZE, stdin);
            if (buf_len == 0) return EOF;
        }

        return inbuf[inptr++];
    }

    template <typename T, typename=typename enable_if<is_integral<T>::value, T>::type>
    bool read_one(T &x){
        int c = ' ', neg = 0;
        while (c != '-' && !isdigit(c) && c != EOF) c = read_char();
        if (c == '-') neg = 1, c = read_char();
        if (c == EOF) return false;

        for (x = 0; isdigit(c); c = read_char()){
            x = x * 10 + c - '0';
        }
        if (neg) x = -x;

        return true;
    }

    bool read_one(string &s){
        int c = ' ';
        while (isspace(c) && c != EOF) c = read_char();
        if (c == EOF) return false;

        for (s.clear(); !isspace(c); c = read_char()){
            s.push_back(c);
        }
        return true;
    }

    template <typename T>
    bool read_one(vector<T>& v, bool read_length=true){
        if (read_length){
            int n;
            if (!read_one(n)) return false;
            v.resize(n);
        }

        for (auto &&x: v){
            if (!read_one(x)) return false;
        }
        return true;
    }

    bool read_line(string &s){
        s.clear();
        int c = '\n';
        while (c == '\n' || c == '\r') c = read_char();

        while (c != '\n' && c != '\r' && c != EOF){
            s.push_back(c);
            c = read_char();
        }
        return c != EOF;
    }

    int read(){
        return 0;
    }

    template <typename T, typename ...Args>
    int read(T &x, Args& ...args){
        if (!read_one(x)) return 0;
        return read(args...) + 1;
    }

    /* End of read methods, write methods start below */

    void flush(){
        fwrite(outbuf, 1, outptr, stdout);
        outptr = 0;
    }

    inline void write_char(const char& c){
        if (outptr == BUF_SIZE) flush();
        outbuf[outptr++] = c;
    }

    void write_one(const char* s){
        for (int j = 0; s[j]; j++) write_char(s[j]);
    }

    void write_one(const string& s){
        for (auto &&c: s) write_char(c);
    }

    template <typename T, typename=typename enable_if<is_integral<T>::value, T>::type>
    void write_one(T x){
        if (x < 0) x = -x, write_char('-');

        int l = 0;
        while (x || !l){
            tmpbuf[l++] = (x % 10) + '0';
            x /= 10;
        }
        while (l) write_char(tmpbuf[--l]);
    }

    template <typename T>
    void write_one(const vector<T>& v){
        for (int i = 0; i < (int)v.size(); i++){
            if (i) write_char(' ');
            write_one(v[i]);
        }
    }

    void write(){}

    template <typename T, typename ...Args>
    void write(const T x, const Args ...args){
        write_one(x);
        write_char(sizeof...(args) && is_trivial<T>::value ? ' ' : '\n');
        write(args...);
    }
}

int main(){
    /// SPOJ Enormous Input and Output Test (https://www.spoj.com/problems/INOUTEST/)

    using namespace fio;
    int n, a, b;

    read(n);
    while (n--){
        read(a, b);
        write(a * b);
    }
    flush();

    return 0;
}
