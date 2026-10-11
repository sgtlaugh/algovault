#include "../common.h"

#define main library_main
#include "../../code_library/strings/aho_corasick.cpp"
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

vector<int> brute_each(const vector<string>& patterns, const string& text){
    vector<int> res;
    for (const auto& p : patterns){
        int c = 0;
        for (size_t i = 0; i + p.size() <= text.size(); i++) c += text.compare(i, p.size(), p) == 0;
        res.push_back(c);
    }

    return res;
}

int main(){
    for (long long it = 0; it < stress::scaled(1000); it++){
        int alphabet = stress::rand_int(1, it % 5 ? 3 : 26);
        vector<string> patterns(stress::rand_int(0, it % 10 ? 8 : 60));
        for (auto& p : patterns) p = random_word(stress::rand_int(1, 6), alphabet);

        AhoCorasick ac;
        for (const auto& p : patterns){
            int v = ac.insert(p), cur = 0;
            for (char c : p) cur = ac.go[cur][c - 'a'];  /// still the plain trie before build()
            assert(cur == v);
        }
        ac.build();
        assert(ac.size() == (int)patterns.size());

        for (int q = 0; q < 10; q++){
            string text = random_word(stress::rand_int(0, it % 10 ? 40 : 2000), alphabet);
            for (auto& c : text) if (stress::rand_int(0, 15) == 0) c = "A{ .\xe9\x80"[stress::rand_int(0, 5)];  /// outside the alphabet, including bytes above 0x7f
            assert(ac.count(text) == brute(patterns, text));
            assert(ac.occurrences(text) == brute_each(patterns, text));
        }
    }

    /// A shifted alphabet: digits, with letters and other bytes in the text acting as separators
    for (long long it = 0; it < stress::scaled(300); it++){
        vector<string> patterns(stress::rand_int(1, 10));
        for (auto& p : patterns){
            p = random_word(stress::rand_int(1, 4), 3);
            for (auto& c : p) c = c - 'a' + '0';
        }

        AhoCorasick<10, '0'> ac;
        for (const auto& p : patterns) ac.insert(p);
        ac.build();
        string text = random_word(stress::rand_int(0, 60), 3);
        for (auto& c : text) c = stress::rand_int(0, 7) ? c - 'a' + '0' : "a/:9\xff"[stress::rand_int(0, 4)];
        assert(ac.count(text) == brute(patterns, text));
        assert(ac.occurrences(text) == brute_each(patterns, text));
    }

    return 0;
}
