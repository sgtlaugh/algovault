// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/find_linear_recurrence
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/algebra/linear_recurrence.cpp"
#undef main

long long read_int(){
    int c = getchar();
    while (c != '-' && !isdigit(c)) c = getchar();

    bool neg = c == '-';
    if (neg) c = getchar();

    long long x = 0;
    for (; isdigit(c); c = getchar()) x = x * 10 + (c - '0');
    return neg ? -x : x;
}

template<typename T>
void print_line(const vector<T>& v){
    string out;
    for (size_t i = 0; i < v.size(); i++){
        if (i) out += ' ';
        out += to_string(v[i]);
    }
    out += '\n';
    fwrite(out.data(), 1, out.size(), stdout);
}

int main(){
    const int mod = 998244353;
    int n = read_int();
    vector<int> a(n);
    for (auto& x : a) x = read_int();

    /// berlekamp_massey returns C = [1, C1, ..., Cd] with sum C[j] * a[i - j] = 0, so c_j = -C[j]
    auto C = LinearRecurrence().berlekamp_massey(a, mod);
    vector<int> c;
    for (size_t j = 1; j < C.size(); j++) c.push_back((mod - C[j]) % mod);

    printf("%d\n", (int)c.size());
    print_line(c);
    return 0;
}
