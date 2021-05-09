/***
 *
 * Fast i/o with fread and fwrite
 * readint(), readstr() can take multiple arguments and returns the number of objects read
 * readline() reads one line at a time, it does not skip empty lines
 *
 * Similarly, writeint() and writestr() takes multiple arguments
 * The write functions prints them separated by a single space and ending with a newline
 * Don't forget to flush the buffer at the end
 *
 * Example usage - Read two numbers till EOF and print their sum
 *
 * int a, b;
 * while (readint(a, b)){
 *     writeint(a + b);
 * }
 * flush();
 *
***/

#include <bits/stdc++.h>

using namespace std;

namespace fstio{
    const int BUFSIZE = 65536;

    int len = 0, inptr = 0, outptr = 0;
    char in[BUFSIZE], out[BUFSIZE], tmp[256];

    char readchar(){
        if (inptr >= len){
            inptr = 0, len = fread(in, 1, BUFSIZE, stdin);
            if (len == 0) return EOF;
        }

        return in[inptr++];
    }

    int readint(){
        return 0;
    }

    template <typename T, typename ...Args>
    int readint(T &x, Args& ...args){
        int c = ' ', neg = 0;
        while (c != '-' && !isdigit(c) && c != EOF) c = readchar();
        if (c == '-') neg = 1, c = readchar();
        if (c == EOF) return 0;

        for (x = 0; isdigit(c); c = readchar()){
            x = x * 10 + c - '0';
        }
        if (neg) x = -x;

        return readint(args...) + 1;
    }

    int readstr(){
        return 0;
    }

    template <typename T, typename ...Args>
    int readstr(T &s, Args& ...args){
        s.clear();
        int c = ' ';
        while (isspace(c) && c != EOF) c = readchar();
        if (c == EOF) return 0;

        for (; !isspace(c); c = readchar()){
            s.push_back(c);
        }
        return readstr(args...) + 1;
    }

    int readline(string &s){
        s.clear();
        int c = ' ';
        while (c == '\n' || c == '\r') c = readchar();

        while (c != '\n' && c != '\r'){
            if (c == EOF) return 0;
            s.push_back(c);
            c = readchar();
        }
        return 1;
    }

    void flush(){
        fwrite(out, 1, outptr, stdout);
        outptr = 0;
    }

    void writechar(char c){
        if (outptr == BUFSIZE) flush();
        out[outptr++] = c;
    }

    void writestr(){}

    template <typename T, typename ...Args>
    void writestr(const T& s, const Args ...args){
        for (int j = 0; s[j]; j++) writechar(s[j]);

        const int n = sizeof...(args);
        writechar(n == 0 ? '\n' : ' ');
        writestr(args...);
    }

    void writeint(){}

    template <typename T, typename ...Args>
    void writeint(T x, Args ...args){
        if (x < 0) x = -x, writechar('-');

        int l = 0;
        while (x || !l){
            tmp[l++] = (x % 10) + '0';
            x /= 10;
        }
        while (l) writechar(tmp[--l]);

        const int n = sizeof...(args);
        writechar(n == 0 ? '\n' : ' ');
        writeint(args...);
    }
}

int main(){
    /// SPOJ Enormous Input and Output Test (https://www.spoj.com/problems/INOUTEST/)

    using namespace fstio;

    int n, a, b;

    readint(n);
    while (n--){
        readint(a, b);
        writeint(a * b);
    }
    flush();

    return 0;
}
