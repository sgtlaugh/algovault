// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/matrix_product
// competitive-verifier: TLE 10
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/linear_algebra/matrix.cpp"
#undef main

long long read_int(){
    int c = getchar();
    while (c < '0' || c > '9') c = getchar();
    long long x = 0;
    while (c >= '0' && c <= '9'){
        x = x * 10 + c - '0';
        c = getchar();
    }
    return x;
}

void write_int(long long x, char end){
    char s[24];
    int len = 0;
    do{
        s[len++] = '0' + x % 10;
        x /= 10;
    } while (x);
    while (len) putchar(s[--len]);
    putchar(end);
}

int main(){
    const long long MOD = 998244353;
    int n = read_int(), m = read_int(), k = read_int();

    Matrix a(n, m, MOD), b(m, k, MOD);
    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++) a[i][j] = read_int();
    }
    for (int i = 0; i < m; i++){
        for (int j = 0; j < k; j++) b[i][j] = read_int();
    }

    Matrix c = a * b;

    for (int i = 0; i < n; i++){
        for (int j = 0; j < k; j++) write_int(c[i][j], j + 1 == k ? '\n' : ' ');
    }
    return 0;
}
