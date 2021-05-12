/***
 *
 * Fast i/o with fread and fwrite
 *
 * read() can take multiple arguments and returns the number of objects read
 * read() works for strings or any integral types (int, long long, etc)
 * read() also works for vector, but the type has to be string or integral
 *
 * read_line() reads one line at a time, it does not skip empty lines
 *
 * Similarly, write() can take multiple arguments, arguments can be string or integral
 * The write functions prints them separated by a single space, ending with a newline
 * Write also works for vectors just like read
 * Vector elements are separated by a single space and ends with a newline
 *
 * Don't forget to flush the buffer at the end after writing
 *
 * Example usage - Read three numbers till EOF and print their sum
 *
 * int a, c;
 * long long b;
 *
 * while (read(a, b, c)){
 *     write(a + b + c);
 * }
 * flush();
 *
***/

#include <bits/stdc++.h>

using namespace std;

namespace fstio{
    const int BUF_SIZE = 8192;

    int len = 0, inptr = 0, outptr = 0;
    char in[BUF_SIZE], out[BUF_SIZE], tmp[128];

    char read_char(){
        if (inptr >= len){
            inptr = 0, len = fread(in, 1, BUF_SIZE, stdin);
            if (len == 0) return EOF;
        }

        return in[inptr++];
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

    void flush(){
        fwrite(out, 1, outptr, stdout);
        outptr = 0;
    }

    void write_char(const char& c){
        if (outptr == BUF_SIZE) flush();
        out[outptr++] = c;
    }

    void write_one(const string& s){
        for (auto &&c: s) write_char(c);
    }

    template <typename T>
    void write_one(T x){
        if (x < 0) x = -x, write_char('-');

        int l = 0;
        while (x || !l){
            tmp[l++] = (x % 10) + '0';
            x /= 10;
        }
        while (l) write_char(tmp[--l]);
    }

    template <typename T>
    void write_one(const vector<T>& v, char delimiter=' ', char last_char=0){
        for (int i = 0; i < (int)v.size(); i++){
            if (i) write_char(delimiter);
            write_one(v[i]);
        }
        if (last_char) write_char(last_char);
    }

    void write(){}

    template <typename T, typename ...Args>
    void write(const T x, const Args ...args){
        write_one(x);
        write_char(sizeof...(args) ? ' ' : '\n');
        write(args...);
    }
}

int main(){
    /// SPOJ Enormous Input and Output Test (https://www.spoj.com/problems/INOUTEST/)

    using namespace fstio;

    int n, a, b;

    read(n);
    while (n--){
        read(a, b);
        write(a * b);
    }
    flush();

    return 0;
}
