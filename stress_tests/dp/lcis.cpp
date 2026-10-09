#include "../common.h"

#define main library_main
#include "../../code_library/dp/lcis.cpp"
#undef main

/// Every strictly increasing subsequence of A, checked for being a subsequence of B
int brute(const vector<int>& A, const vector<int>& B){
    int n = A.size(), best = 0;
    for (int mask = 1; mask < (1 << n); mask++){
        vector<int> seq;
        for (int i = 0; i < n; i++){
            if (mask >> i & 1) seq.push_back(A[i]);
        }
        bool increasing = true;
        for (size_t i = 1; i < seq.size(); i++) increasing &= seq[i - 1] < seq[i];
        if (!increasing || (int)seq.size() <= best) continue;
        size_t k = 0;
        for (int b : B){
            if (k < seq.size() && b == seq[k]) k++;
        }
        if (k == seq.size()) best = seq.size();
    }
    return best;
}

int main(){
    for (long long it = 0; it < stress::scaled(6000); it++){
        int n = stress::rand_int(0, 12), m = stress::rand_int(0, 14), range = stress::rand_int(1, it % 2 ? 4 : 30);
        vector<int> A(n), B(m);
        for (auto& x : A) x = stress::rand_int(-range, range);
        for (auto& x : B) x = stress::rand_int(-range, range);
        assert(lcis(A, B) == brute(A, B));
    }
    vector<int> A(3000), B(3000);
    iota(A.begin(), A.end(), 0), iota(B.begin(), B.end(), 0);
    assert(lcis(A, B) == 3000);
    return 0;
}
