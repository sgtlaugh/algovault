// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/stern_brocot_tree
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/number_theory/continued_fractions.cpp"
#undef main

/// Run lengths of the path from 1/1, alternating R, L, R, ... starting with R, so runs[0] may be 0
/// The path to [a0; a1, ..., ak] is R^a0 L^a1 ... with the last term one short
vector<long long> encode(long long a, long long b){
    auto runs = continued_fraction(a, b);
    runs.back()--;
    return runs;
}

pair<long long, long long> decode(vector<long long> runs){
    runs.back()++;
    return convergents(runs).back();
}

vector<long long> truncate(const vector<long long>& runs, long long depth){
    vector<long long> res;
    for (long long run : runs){
        res.push_back(min(run, depth));
        depth -= res.back();
        if (!depth) break;
    }
    return res;
}

int main(){
    int t;
    if (scanf("%d", &t) != 1) return 0;

    char op[16];
    while (t--){
        if (scanf("%s", op) != 1) return 0;
        string s = op;

        if (s == "ENCODE_PATH"){
            long long a, b;
            if (scanf("%lld %lld", &a, &b) != 2) return 0;
            auto runs = encode(a, b);

            vector<pair<char, long long>> out;
            for (size_t i = 0; i < runs.size(); i++){
                if (runs[i]) out.push_back({i % 2 ? 'L' : 'R', runs[i]});
            }
            printf("%d", (int)out.size());
            for (auto [c, n] : out) printf(" %c %lld", c, n);
            puts("");
        }
        else if (s == "DECODE_PATH"){
            int k;
            if (scanf("%d", &k) != 1) return 0;

            vector<long long> runs;
            for (int i = 0; i < k; i++){
                char c;
                long long n;
                if (scanf(" %c %lld", &c, &n) != 2) return 0;
                if (i == 0 && c == 'L') runs.push_back(0);
                runs.push_back(n);
            }
            if (runs.empty()) runs.push_back(0);

            auto [a, b] = decode(runs);
            printf("%lld %lld\n", a, b);
        }
        else if (s == "LCA"){
            long long a, b, c, d;
            if (scanf("%lld %lld %lld %lld", &a, &b, &c, &d) != 4) return 0;
            auto x = encode(a, b), y = encode(c, d);

            vector<long long> common;
            for (size_t i = 0; i < min(x.size(), y.size()); i++){
                common.push_back(min(x[i], y[i]));
                if (x[i] != y[i]) break;
            }

            auto [f, g] = decode(common);
            printf("%lld %lld\n", f, g);
        }
        else if (s == "ANCESTOR"){
            long long k, a, b;
            if (scanf("%lld %lld %lld", &k, &a, &b) != 3) return 0;
            auto runs = encode(a, b);

            if (accumulate(runs.begin(), runs.end(), 0LL) < k){
                puts("-1");
                continue;
            }
            auto [f, g] = decode(truncate(runs, k));
            printf("%lld %lld\n", f, g);
        }
        else{
            long long a, b;
            if (scanf("%lld %lld", &a, &b) != 2) return 0;
            auto conv = convergents(continued_fraction(a, b));

            /// a / b is the mediant of its two bounds, one of which is the previous convergent (1/0 before the first)
            pair<long long, long long> lo = conv.size() > 1 ? conv[conv.size() - 2] : make_pair(1LL, 0LL);
            pair<long long, long long> hi = {a - lo.first, b - lo.second};
            if (lo.first * hi.second > hi.first * lo.second) swap(lo, hi);
            printf("%lld %lld %lld %lld\n", lo.first, lo.second, hi.first, hi.second);
        }
    }
    return 0;
}
