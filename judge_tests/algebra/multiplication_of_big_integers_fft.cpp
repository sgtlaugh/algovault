// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/multiplication_of_big_integers
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/algebra/fft.cpp"
#undef main

/// Base 1000 keeps 2 * 10^6 digits within MAX = 2^21 transform points and every coefficient below 7 * 10^11
vector<long long> to_groups(const string& s, size_t from){
    vector<long long> res;
    for (size_t j = s.size(); j > from; j -= min<size_t>(3, j - from)){
        long long group = 0;
        for (size_t k = j - min<size_t>(3, j - from); k < j; k++) group = group * 10 + (s[k] - '0');
        res.push_back(group);
    }
    return res;
}

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int t;
    cin >> t;

    string out;
    while (t--){
        string x, y;
        cin >> x >> y;

        bool neg = (x[0] == '-') != (y[0] == '-');
        auto res = fft::multiply(to_groups(x, x[0] == '-'), to_groups(y, y[0] == '-'));

        long long carry = 0;
        for (auto& v : res){
            v += carry;
            carry = v / 1000, v %= 1000;
        }
        for (; carry; carry /= 1000) res.push_back(carry % 1000);
        while (res.size() > 1 && res.back() == 0) res.pop_back();

        if (neg && res.back() != 0) out += '-';
        out += to_string(res.back());
        for (int i = (int)res.size() - 2; i >= 0; i--){
            char digits[4];
            snprintf(digits, sizeof(digits), "%03lld", res[i]);
            out += digits;
        }
        out += '\n';
    }

    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
