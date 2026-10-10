#include "../common.h"

#define main library_main
#include "../../code_library/strings/dynamic_aho_corasick.cpp"
#undef main

string random_word(int len, int alphabet){
    string s(len, 'a');
    for (auto& c : s) c = 'a' + stress::rand_int(0, alphabet - 1);
    return s;
}

long long brute(const vector<string>& patterns, const string& text){
    long long res = 0;
    for (const auto& p : patterns){
        for (size_t i = 0; i + p.size() <= text.size(); i++) res += text.compare(i, p.size(), p) == 0;
    }

    return res;
}

/// Queried after every insert so that each merge of the binary decomposition is exercised
int main(){
    for (long long it = 0; it < stress::scaled(200); it++){
        int alphabet = stress::rand_int(1, it % 5 ? 3 : 26);
        vector<string> patterns;
        DynamicAhoCorasick ac;

        for (int ins = stress::rand_int(0, it % 10 ? 20 : 150); ins; ins--){
            patterns.push_back(random_word(stress::rand_int(1, 6), alphabet));
            ac.insert(patterns.back());
            assert((int)ac.ar.size() == __lg((int)patterns.size()) + 1);

            string text = random_word(stress::rand_int(0, 60), alphabet);
            for (auto& c : text) if (stress::rand_int(0, 15) == 0) c = "A{ .\xe9\x80"[stress::rand_int(0, 5)];  /// outside the alphabet, including bytes above 0x7f
            assert(ac.count(text) == brute(patterns, text));
        }
    }
    return 0;
}
