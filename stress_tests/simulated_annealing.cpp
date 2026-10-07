#include "common.h"

#define main library_main
#include "../code_library/simulated_annealing.c"
#undef main

const string IN = "/tmp/algovault_annealing_in.txt", OUT = "/tmp/algovault_annealing_out.txt";

int lnds(const vector<int>& v){
    vector<int> tails;  /// longest non-decreasing subsequence by patience sorting
    for (int x : v){
        auto it = upper_bound(tails.begin(), tails.end(), x);
        if (it == tails.end()) tails.push_back(x);
        else *it = x;
    }
    return tails.size();
}

/// Every subset of positions reversed in place, the exact answer the annealing searches for
int best_reversal(const vector<int>& v){
    int n = v.size(), best = 0;
    for (int mask = 0; mask < (1 << n); mask++){
        vector<int> idx, w = v;
        for (int i = 0; i < n; i++) if (mask >> i & 1) idx.push_back(i);
        for (size_t i = 0, j = idx.size(); i + 1 < j; i++, j--) swap(w[idx[i]], w[idx[j - 1]]);
        best = max(best, lnds(w));
    }
    return best;
}

vector<int> run(const string& input){
    ofstream(IN) << input;
    assert(freopen(IN.c_str(), "r", stdin) && freopen(OUT.c_str(), "w", stdout));
    library_main();
    fflush(stdout);
    vector<int> res;
    ifstream out(OUT);
    for (int x; out >> x; ) res.push_back(x);
    return res;
}

int main(){
    for (long long it = 0; it < stress::scaled(3000); it++){
        n = stress::rand_int(1, 50);
        vector<int> v(n);
        for (int i = 0; i < n; i++) v[i] = ar[i] = stress::rand_int(-3, 3) * (it % 2 ? 1 : 1000000);
        assert(solve() == lnds(v));
    }

    vector<vector<int>> cases;
    string input;
    for (long long k = 0; k < stress::scaled(6); k++){
        vector<int> v(stress::rand_int(1, 11));
        int range = stress::rand_int(1, 10);
        for (auto& x : v) x = stress::rand_int(1, range);
        cases.push_back(v);
        input += to_string(v.size()) + "\n";
        for (int x : v) input += to_string(x) + " ";
        input += "\n";
    }

    auto res = run(input);
    assert(equal(cases.back().begin(), cases.back().end(), ar));  /// each candidate reversal is applied twice, which must restore the array
    assert(res.size() == cases.size() && run(input) == res);  /// srand(17) per run, so a rerun must repeat itself
    for (size_t k = 0; k < cases.size(); k++){
        assert(lnds(cases[k]) <= res[k] && res[k] <= best_reversal(cases[k]));
    }

    /// Only reversing positions 0 and 2 reaches 3, with 8 states and 180000 steps the search cannot miss it
    assert(run("3\n9 1 1\n") == vector<int>{3});
    return 0;
}
