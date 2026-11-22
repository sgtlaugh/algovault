// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/0350
// competitive-verifier: TLE 1
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/algebra/fraction.cpp"
#undef main

long long power_of_ten(int k){
    long long res = 1;
    while (k--) res *= 10;
    return res;
}

int main(){
    char buf[16];
    if (scanf("%15s", buf) != 1) return 0;
    string s = buf;

    size_t dot = s.find('.'), open = s.find('(');
    string fixed = s.substr(dot + 1, open == string::npos ? string::npos : open - dot - 1);
    Fraction res(stoll(s.substr(0, dot)));
    if (!fixed.empty()) res += Fraction(stoll(fixed), power_of_ten(fixed.size()));

    /// 0.0(r) with r of k digits is r / (10^k - 1), shifted right past the fixed digits
    if (open != string::npos){
        string cycle = s.substr(open + 1, s.size() - open - 2);
        res += Fraction(stoll(cycle), power_of_ten(fixed.size()) * (power_of_ten(cycle.size()) - 1));
    }

    printf("%lld/%lld\n", res.num, res.den);
    return 0;
}
