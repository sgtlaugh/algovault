#include "../common.h"

#define main library_main
#include "../../code_library/strings/hirschberg_lcs.cpp"
#undef main

template <typename Sequence>
int lcs_length(const Sequence& a, const Sequence& b){
    vector<vector<int>> dp(a.size() + 1, vector<int>(b.size() + 1, 0));
    for (size_t i = 1; i <= a.size(); i++){
        for (size_t j = 1; j <= b.size(); j++){
            dp[i][j] = a[i - 1] == b[j - 1] ? dp[i - 1][j - 1] + 1 : max(dp[i - 1][j], dp[i][j - 1]);
        }
    }

    return dp[a.size()][b.size()];
}

template <typename Sequence>
bool is_subsequence(const Sequence& sub, const Sequence& of){
    size_t k = 0;
    for (size_t i = 0; i < of.size() && k < sub.size(); i++){
        if (of[i] == sub[k]) k++;
    }

    return k == sub.size();
}

/// Many LCSs may tie, so the oracle is the optimal length plus being a common subsequence
template <typename Sequence>
void check(const Sequence& a, const Sequence& b){
    Sequence res = lcs_sequence(a, b);
    assert((int)res.size() == lcs_length(a, b));
    assert(is_subsequence(res, a));
    assert(is_subsequence(res, b));
}

string random_string(int n, int sigma){
    string s(n, 'a');
    for (auto& c : s) c = 'a' + stress::rand_int(0, sigma - 1);
    return s;
}

int main(){
    /// Tiny alphabets make long ties, both length orders exercise the swap to the shorter row
    for (long long it = 0; it < stress::scaled(20000); it++){
        int sigma = stress::rand_int(1, 4);
        check(random_string(stress::rand_int(0, 12), sigma), random_string(stress::rand_int(0, 12), sigma));
    }

    for (long long it = 0; it < stress::scaled(300); it++){
        int sigma = stress::rand_int(1, 26);
        check(random_string(stress::rand_int(0, 200), sigma), random_string(stress::rand_int(0, 200), sigma));

        vector<int> a(stress::rand_int(0, 60)), b(stress::rand_int(0, 60));
        for (auto& x : a) x = stress::rand_int(-3, 3);
        for (auto& x : b) x = stress::rand_int(-3, 3);
        check(a, b);
    }

    /// Larger inputs, where the recursion splits many levels deep
    for (long long it = 0; it < stress::scaled(2); it++){
        check(random_string(2000, 4), random_string(2000, 4));
        check(random_string(2000, 26), random_string(1500, 26));
    }

    return 0;
}
