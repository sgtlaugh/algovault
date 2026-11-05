#include "../common.h"

#define main library_main
#include "../../code_library/strings/hunt_szymanski.cpp"
#undef main

template <typename Seq>
int dp_lcs(const Seq& a, const Seq& b){
    vector<int> prev(b.size() + 1), cur(b.size() + 1);
    for (auto c : a){
        for (size_t j = 0; j < b.size(); j++) cur[j + 1] = c == b[j] ? prev[j] + 1 : max(prev[j + 1], cur[j]);
        swap(prev, cur);
    }

    return prev[b.size()];
}

string random_string(int len, const string& letters){
    string s(len, 0);
    for (auto& c : s) c = letters[stress::rand_int(0, letters.size() - 1)];
    return s;
}

int main(){
    const string alphabets[] = {"a", "ab", "acgt", "abcdefghijklmnopqrstuvwxyz", "a\x80\xff\x01 Z~"};

    for (long long it = 0; it < stress::scaled(3000); it++){
        const string& letters = alphabets[stress::rand_int(0, 4)];
        int len = it % 10 ? 40 : 400;
        string a = random_string(stress::rand_int(0, len), letters);
        string b = random_string(stress::rand_int(0, len), letters);
        assert(lcs(a, b) == dp_lcs(a, b));
    }

    /// Values no byte holds, 0, negatives and LLONG_MIN, from few values for dense matches or many for sparse ones
    for (long long it = 0; it < stress::scaled(3000); it++){
        long long range = stress::rand_int(0, 1) ? 3 : 1000000;
        vector<long long> a(stress::rand_int(0, 60)), b(stress::rand_int(0, 60));
        for (auto& x : a) x = stress::rand_int(0, 9) ? stress::rand_int(-range, range) * 1000003 : LLONG_MIN;
        for (auto& x : b) x = stress::rand_int(0, 9) ? stress::rand_int(-range, range) * 1000003 : LLONG_MIN;
        assert(lcs(a, b) == dp_lcs(a, b));
    }

    /// Long enough to overflow a fixed-size buffer, a wide alphabet keeps the match count R small
    string letters;
    for (int c = 1; c < 256; c++) letters += (char)c;
    string a = random_string(52000, letters), b = a;
    b.erase(stress::rand_int(0, b.size() - 1), 1);
    assert(lcs(a, b) == (int)b.size());

    /// Permutations of 1..n reduce LCS to a longest increasing subsequence of positions, checked by patience sorting
    int n = 200000;
    vector<int> p(n), q(n), where(n + 1), lis;
    iota(p.begin(), p.end(), 1), iota(q.begin(), q.end(), 1);
    shuffle(p.begin(), p.end(), mt19937(stress::rand_int(0, 1000000)));
    shuffle(q.begin(), q.end(), mt19937(stress::rand_int(0, 1000000)));
    for (int i = 0; i < n; i++) where[q[i]] = i;
    for (int x : p){
        auto it = lower_bound(lis.begin(), lis.end(), where[x]);
        if (it == lis.end()) lis.push_back(where[x]);
        else *it = where[x];
    }
    assert(lcs(p, q) == (int)lis.size());

    return 0;
}
