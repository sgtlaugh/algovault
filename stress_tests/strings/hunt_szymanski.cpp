#include "../common.h"

#define main library_main
#include "../../code_library/strings/hunt_szymanski.cpp"
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
    for (long long it = 0; it < stress::scaled(3000); it++){
        const string& letters = alphabets[stress::rand_int(0, 4)];
        string a = random_string(stress::rand_int(0, it % 10 ? 40 : 400), letters), b = random_string(stress::rand_int(0, it % 10 ? 40 : 400), letters);
        assert(lcs(a.c_str(), b.c_str()) == dp_lcs(a, b));
    }

    /// Longer than the old fixed MAX buffer, a wide alphabet keeps the match count R small
    string letters;
    for (int c = 1; c < 256; c++) letters += (char)c;
    string a = random_string(52000, letters), b = a;
    b.erase(stress::rand_int(0, b.size() - 1), 1);
    assert(lcs(a.c_str(), b.c_str()) == (int)b.size());
    return 0;
}
