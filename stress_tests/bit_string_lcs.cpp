#include "common.h"

#define main library_main
#include "../code_library/bit_string_lcs.cpp"
#undef main

int dp_lcs(const string& a, const string& b){
    vector<int> prev(b.size() + 1), cur(b.size() + 1);
    for (char c : a){
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
    for (long long it = 0; it < stress::scaled(2000); it++){
        const string& letters = alphabets[stress::rand_int(0, 4)];
        int n = stress::rand_int(0, it % 10 ? 70 : 600), m = stress::rand_int(0, it % 10 ? 200 : 600);  /// m crosses several 64-bit blocks
        string a = random_string(n, letters), b = random_string(m, letters);
        assert(lcs(a.c_str(), b.c_str()) == dp_lcs(a, b));
    }

    /// A longer than the old fixed carry buffer
    string a = random_string(150000, "ab"), b = random_string(stress::rand_int(1, 130), "ab");
    assert(lcs(a.c_str(), b.c_str()) == dp_lcs(a, b));
    return 0;
}
