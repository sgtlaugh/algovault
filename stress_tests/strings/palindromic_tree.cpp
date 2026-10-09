#include "../common.h"

#define main library_main
#include "../../code_library/strings/palindromic_tree.cpp"
#undef main

bool is_palindrome(const string& s, int l, int r){
    while (l < r){
        if (s[l++] != s[r--]) return false;
    }

    return true;
}

/// Distinct palindromes, occurrence counts, longest and freshness per character against enumerating substrings
void check(const string& s){
    PalindromicTree<3, 'a'> tree;
    set<string> seen;
    int longest = 0;
    for (int i = 0; i < (int)s.size(); i++){
        bool fresh = tree.add(s[i]);
        bool expected = false;
        for (int l = 0; l <= i; l++){
            if (is_palindrome(s, l, i)){
                expected |= seen.insert(s.substr(l, i - l + 1)).second;
                longest = max(longest, i - l + 1);
            }
        }
        assert(fresh == expected);
        assert(tree.distinct() == (int)seen.size() && tree.longest == longest);
    }

    map<string, long long> count;
    for (int l = 0; l < (int)s.size(); l++){
        for (int r = l; r < (int)s.size(); r++){
            if (is_palindrome(s, l, r)) count[s.substr(l, r - l + 1)]++;
        }
    }

    auto occ = tree.occurrences();
    map<int, vector<long long>> by_length_tree, by_length_brute;
    for (int v = 2; v < (int)tree.len.size(); v++) by_length_tree[tree.len[v]].push_back(occ[v]);
    for (auto& [p, c] : count) by_length_brute[p.size()].push_back(c);
    for (auto& [len, v] : by_length_tree) sort(v.begin(), v.end());
    for (auto& [len, v] : by_length_brute) sort(v.begin(), v.end());
    assert(by_length_tree == by_length_brute);

    int n = s.size();
    vector<int> best(n + 1, INT_MAX);
    best[0] = 0;
    for (int i = 1; i <= n; i++){
        for (int j = 0; j < i; j++){
            if (is_palindrome(s, j, i - 1)) best[i] = min(best[i], best[j] + 1);
        }
    }

    auto fact = min_palindromic_factorization<3, 'a'>(s);
    for (int i = 0; i < n; i++) assert(fact[i] == best[i + 1]);
}

int main(){
    for (long long it = 0; it < stress::scaled(5000); it++){
        int n = stress::rand_int(0, 40), sigma = stress::rand_int(1, 3);
        string s;
        for (int i = 0; i < n; i++) s += char('a' + stress::rand_int(0, sigma - 1));
        check(s);
    }

    for (int n = 0; n <= 12; n++){
        for (int mask = 0; mask < (1 << n); mask++){
            string s;
            for (int i = 0; i < n; i++) s += mask >> i & 1 ? 'b' : 'a';
            check(s);
        }
    }

    /// Long inputs: Fibonacci words have many palindromes, "aaa..." is the worst case for suffix runs
    string fib_a = "a", fib_b = "ab";
    while (fib_b.size() < 200000) tie(fib_a, fib_b) = make_pair(fib_b, fib_b + fib_a);
    for (const string& s : {fib_b, string(200000, 'a')}){
        auto fact = min_palindromic_factorization<3, 'a'>(s);
        assert((int)fact.size() == (int)s.size());
        if (s[1] == 'a') assert(fact.back() == 1);
    }

    return 0;
}
