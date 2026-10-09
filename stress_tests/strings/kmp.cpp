#include "../common.h"

#define main library_main
#include "../../code_library/strings/kmp.cpp"
#undef main

template <typename Container>
void check(const Container& text, const Container& pattern){
    int n = text.size(), m = pattern.size();
    vector<int> expected;
    for (int i = 0; m && i + m <= n; i++){
        if (equal(pattern.begin(), pattern.end(), text.begin() + i)) expected.push_back(i);
    }
    assert(kmp_search(text, pattern) == expected);

    auto fail = kmp_failure(pattern);
    assert((int)fail.size() == m);
    for (int i = 0; i < m; i++){
        int border = i;
        while (border > 0 && !equal(pattern.begin(), pattern.begin() + border, pattern.begin() + i + 1 - border)) border--;
        assert(fail[i] == border - 1);
    }
}

/// one naive KMP step: longest prefix of pattern (at most m) that is a suffix of pattern[0..state) + c
int naive_step(const vector<int>& pattern, int state, int c){
    int m = pattern.size();
    vector<int> read(pattern.begin(), pattern.begin() + state);
    read.push_back(c);

    for (int k = min(m, (int)read.size()); k > 0; k--){
        if (equal(pattern.begin(), pattern.begin() + k, read.end() - k)) return k;
    }

    return 0;
}

/// automaton table against naive_step for every state and symbol, then a walk over text against brute-force matches
void check_automaton(const vector<int>& pattern, const vector<int>& text, int alphabet, int first){
    int m = pattern.size();
    vector<int> shifted(pattern);
    for (auto& x : shifted) x += first;
    auto aut = kmp_automaton(shifted, alphabet, first);
    assert((int)aut.size() == m + 1);

    for (int i = 0; i <= m; i++){
        assert((int)aut[i].size() == alphabet);
        for (int c = 0; c < alphabet; c++) assert(aut[i][c] == naive_step(pattern, i, c));
    }

    int state = 0;
    for (int i = 0; i < (int)text.size(); i++){
        state = aut[state][text[i]];
        int expected = 0;
        for (int k = min(m, i + 1); k > 0; k--){
            if (equal(pattern.begin(), pattern.begin() + k, text.begin() + i + 1 - k)){
                expected = k;
                break;
            }
        }
        assert(state == expected);
    }
}

int main(){
    for (long long it = 0; it < stress::scaled(5000); it++){
        int alphabet = stress::rand_int(1, 3);
        string text(stress::rand_int(0, it % 10 ? 40 : 1000), 'a'), pattern(stress::rand_int(0, 8), 'a');
        for (auto& c : text) c = 'a' + stress::rand_int(0, alphabet - 1);
        for (auto& c : pattern) c = 'a' + stress::rand_int(0, alphabet - 1);
        if (it % 4 == 0 && !text.empty()){
            int l = stress::rand_int(0, text.size() - 1);
            pattern = text.substr(l, stress::rand_int(1, text.size() - l));  /// guaranteed to occur
        }
        check(text, pattern);

        vector<long long> vt(text.begin(), text.end()), vp(pattern.begin(), pattern.end());
        for (auto& x : vt) x = x * 1000000007LL - 5;
        for (auto& x : vp) x = x * 1000000007LL - 5;
        check(vt, vp);
    }

    for (long long it = 0; it < stress::scaled(3000); it++){
        int alphabet = stress::rand_int(1, it % 3 ? 3 : 26);
        int first = stress::rand_int(-5, 100);
        vector<int> pattern(stress::rand_int(0, it % 10 ? 10 : 40)), text(stress::rand_int(0, 60));
        for (auto& x : pattern) x = stress::rand_int(0, alphabet - 1);
        for (auto& x : text) x = stress::rand_int(0, alphabet - 1);
        check_automaton(pattern, text, alphabet, first);
    }

    for (int len = 0; len <= 7; len++){
        for (int mask = 0; mask < (1 << len); mask++){
            vector<int> pattern(len);
            for (int i = 0; i < len; i++) pattern[i] = mask >> i & 1;
            check_automaton(pattern, {}, 2, 'a');
        }
    }

    string word = "abracadabra";
    auto aut = kmp_automaton(word);
    int state = 0, matches = 0;
    string text = "abracadabracadabrabracadabra";
    for (char ch : text) matches += (state = aut[state][ch - 'a']) == (int)word.size();
    assert(matches == (int)kmp_search(text, word).size());

    return 0;
}
