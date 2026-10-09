#include "../common.h"

#define main library_main
#include "../../code_library/strings/wildcard_matching.cpp"
#undef main

template<typename Container>
vector<int> naive_match(const Container& text, const Container& pattern, typename Container::value_type wildcard){
    int n = text.size(), m = pattern.size();
    vector<int> positions;
    for (int i = 0; m && i + m <= n; i++){
        bool ok = true;
        for (int j = 0; j < m && ok; j++) ok = text[i + j] == wildcard || pattern[j] == wildcard || text[i + j] == pattern[j];
        if (ok) positions.push_back(i);
    }
    return positions;
}

string random_string(int len, int alphabet, int wildcard_percent){
    string s(len, '?');
    for (auto& c : s){
        if (stress::rand_int(1, 100) > wildcard_percent) c = 'a' + stress::rand_int(0, alphabet - 1);
    }
    return s;
}

/// Plants copies of pattern (with its wildcards replaced by random letters) so long inputs have many matches
void plant(string& text, const string& pattern, int copies){
    int n = text.size(), m = pattern.size();
    for (int k = 0; k < copies && m <= n; k++){
        int at = stress::rand_int(0, n - m);
        for (int j = 0; j < m; j++) text[at + j] = pattern[j] == '?' ? 'a' + stress::rand_int(0, 1) : pattern[j];
    }
}

/// Fastest of three runs, so one slow scheduling slice cannot fail the timing check
template<typename F>
double best_seconds(F run){
    double best = 1e18;
    for (int rep = 0; rep < 3; rep++){
        auto start = chrono::steady_clock::now();
        run();
        best = min(best, chrono::duration<double>(chrono::steady_clock::now() - start).count());
    }
    return best;
}

/// Compares against the naive O(n * m) scan on random strings and int vectors, then checks the n = 2^23 boundary and the cost of many distinct values
int main(){
    for (long long it = 0; it < stress::scaled(4000); it++){
        int n = it < 500 ? it % 12 : stress::rand_int(0, it % 20 ? 60 : 3000);
        int m = stress::rand_int(0, it % 3 ? 6 : n + 1);
        int alphabet = stress::rand_int(1, 3), wildcard_percent = stress::rand_int(0, 4) * 25;
        string text = random_string(n, alphabet, wildcard_percent), pattern = random_string(m, alphabet, wildcard_percent);
        if (it % 2) plant(text, pattern, stress::rand_int(1, 5));
        assert(wildcard_match(text, pattern) == naive_match(text, pattern, '?'));
        auto high_bytes = [](string s){
            for (auto& c : s) c = c == '?' ? '\xff' : (char)(c - 'a' + 0x80);
            return s;
        };
        assert(wildcard_match(high_bytes(text), high_bytes(pattern), '\xff') == naive_match(text, pattern, '?'));

        vector<int> vt(n), vp(m);
        for (auto& x : vt) x = stress::rand_int(0, 7) ? (int)(stress::rand_int(-alphabet, alphabet) * 1000000000 / 3) : INT_MIN;
        for (auto& x : vp) x = stress::rand_int(0, 7) ? (int)(stress::rand_int(-alphabet, alphabet) * 1000000000 / 3) : INT_MIN;
        if (m && m <= n && it % 2) copy(vp.begin(), vp.end(), vt.begin() + stress::rand_int(0, n - m));
        assert(wildcard_match(vt, vp, INT_MIN) == naive_match(vt, vp, INT_MIN));
    }

    for (int it = 0; it < 6; it++){
        int n = stress::rand_int(100000, 300000), m = stress::rand_int(1, it < 3 ? 20 : 300);
        string pattern = random_string(m, 2, 30), text = random_string(n, 2, it % 2 ? 0 : 20);
        plant(text, pattern, 2000);
        assert(wildcard_match(text, pattern) == naive_match(text, pattern, '?'));
    }

    int n = 1 << 23;
    string text = random_string(n, 2, 10), pattern = "ab?ba?b";
    plant(text, pattern, 5000);
    assert(wildcard_match(text, pattern) == naive_match(text, pattern, '?'));

    /// Weighing values by one map lookup per element made this call 3x the string time under the sanitizers, it is now 1.3x
    int big = 1 << 20;
    vector<int> vt(big);
    iota(vt.begin(), vt.end(), 0);
    shuffle(vt.begin(), vt.end(), stress::rng());
    int at = stress::rand_int(0, big - 1000);
    vector<int> vp(vt.begin() + at, vt.begin() + at + 1000);
    for (int k = 0; k < 100; k++) vp[stress::rand_int(0, 999)] = -1;
    string st = random_string(big, 2, 10), sp = random_string(1000, 2, 10);

    vector<int> found;
    double vector_time = best_seconds([&]{ found = wildcard_match(vt, vp, -1); });
    double string_time = best_seconds([&]{ wildcard_match(st, sp); });
    /// Distinct values are the worst case for false positions (< 0.5% per call), so a few extra ones are tolerated
    assert(find(found.begin(), found.end(), at) != found.end() && found.size() <= 4);
    fprintf(stderr, "distinct ints %.2fs, string %.2fs\n", vector_time, string_time);
    assert(vector_time < 2 * string_time);

    return 0;
}
