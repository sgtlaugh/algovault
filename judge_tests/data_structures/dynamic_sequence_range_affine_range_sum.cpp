// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/dynamic_sequence_range_affine_range_sum
// competitive-verifier: TLE 10
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/run_length_implicit_treap.cpp"
#undef main

const long long MOD = 998244353;

struct AffineRuns{
    using Value = long long;
    using Agg = long long;
    using Lazy = pair<long long, long long>;

    static Agg run_agg(Value x, long long count){ return x * count % MOD; }
    static Agg combine(Agg a, Agg b){ return (a + b) % MOD; }
    static Value apply(Value x, Lazy lz){ return (lz.first * x + lz.second) % MOD; }
    static Agg apply(Agg a, Lazy lz, long long count){ return (lz.first * a + lz.second * (count % MOD)) % MOD; }
    static Lazy compose(Lazy f, Lazy g){ return {g.first * f.first % MOD, (g.first * f.second + g.second) % MOD}; }
    static Value reversed(Value x, long long){ return x; }
    static Agg reversed(Agg a){ return a; }
    static pair<Value, Value> cut(Value x, long long){ return {x, x}; }
};

int read_int(){
    int c = getchar();
    while (c < '0' || c > '9') c = getchar();

    int x = 0;
    for (; c >= '0' && c <= '9'; c = getchar()) x = x * 10 + c - '0';
    return x;
}

int main(){
    int n = read_int(), q = read_int();
    vector<pair<long long, long long>> runs(n);
    for (auto& [x, count] : runs) x = read_int(), count = 1;

    RunLengthImplicitTreap<AffineRuns> seq(runs);
    string out;
    while (q--){
        int type = read_int();
        if (type == 0){
            int i = read_int(), x = read_int();
            seq.insert(i, x, 1);
        }
        else if (type == 1){
            int i = read_int();
            seq.erase(i, i);
        }
        else if (type == 2){
            int l = read_int(), r = read_int();
            seq.reverse(l, r - 1);
        }
        else if (type == 3){
            int l = read_int(), r = read_int(), b = read_int(), c = read_int();
            seq.update(l, r - 1, {b, c});
        }
        else{
            int l = read_int(), r = read_int();
            out += to_string(seq.query(l, r - 1)) + '\n';
        }
    }

    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
