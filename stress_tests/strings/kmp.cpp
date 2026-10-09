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
    return 0;
}
