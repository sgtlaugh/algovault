// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/division_of_polynomials
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/algebra/polynomial.cpp"
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
    int n = read_int(), m = read_int();
    vector<int> f(n), g(m);
    for (auto& x : f) x = read_int();
    for (auto& x : g) x = read_int();

    auto [q, r] = Polynomial<>(f).divmod(Polynomial<>(g));
    printf("%d %d\n", q.size(), r.size());
    print_line(q.a);
    print_line(r.a);
    return 0;
}
