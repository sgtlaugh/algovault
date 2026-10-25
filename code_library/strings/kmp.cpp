/***
 *
 * Knuth-Morris-Pratt (KMP) Algorithm
 * Prefix function, every or the first occurrence of a pattern in a text, and the prefix-function automaton
 *
 * Complexity: O(n + m) for text length n and pattern length m, O(m alphabet) for kmp_automaton
 *
 * prefix_function(s)[i]: length of the longest proper prefix of s[0..i] that is also its suffix
 * kmp_find(text, pattern, from = 0): start of the first occurrence at index >= from, -1 if none or the pattern is empty
 * kmp_search(text, pattern): starting indices of every occurrence, overlaps included, none for an empty pattern
 * kmp_automaton(pattern, alphabet = 26, first = 'a'): aut[i][c] = the state after reading symbol first + c in state i
 *
 * Notes:
 *   - State i is the length of the longest pattern prefix that is a suffix of the text read, so m is a full match
 *   - State m continues through the border, so overlapping matches count
 *   - Pattern symbols must lie in [first, first + alphabet)
 *   - kmp_find stops at that match, so a match near from costs far less than kmp_search over the whole text
 *   - The smallest period of s is p = n - prefix_function(s)[n - 1], s repeats s[0..p) iff n % p == 0
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template<typename Container>
vector<int> prefix_function(const Container& c){
    int n = c.size();
    vector<int> pi(n, 0);

    for (int i = 1, k = 0; i < n; i++){
        while (k > 0 && c[k] != c[i]) k = pi[k - 1];
        if (c[k] == c[i]) k++;
        pi[i] = k;
    }

    return pi;
}

/// Its own loop rather than a limit on kmp_search, so kmp_search's hot loop gains no per-match check
template<typename Container>
int kmp_find(const Container& text, const Container& pattern, int from = 0){
    int m = pattern.size();

    if (m == 0) return -1;

    auto pi = prefix_function(pattern);
    int k = 0;

    for (int i = max(from, 0); i < (int)text.size(); i++){
        while (k > 0 && text[i] != pattern[k]) k = pi[k - 1];
        if (pattern[k] == text[i]) k++;
        if (k == m) return i - m + 1;
    }

    return -1;
}

template<typename Container>
vector<int> kmp_search(const Container& text, const Container& pattern){
    vector<int> positions;

    if (pattern.size() == 0) return positions;

    auto pi = prefix_function(pattern);
    int k = 0;
    int m = pattern.size();

    for (int i = 0; i < (int)text.size(); i++){
        while (k > 0 && text[i] != pattern[k]) k = pi[k - 1];
        if (pattern[k] == text[i]) k++;
        if (k == m){
            positions.push_back(i - m + 1);
            k = pi[k - 1];
        }
    }

    return positions;
}

/// aut[i][c] = longest pattern prefix that is a suffix of pattern[0..i) + c, for 0 <= i <= m
template<typename Container>
vector<vector<int>> kmp_automaton(const Container& pattern, int alphabet = 26, int first = 'a'){
    int m = pattern.size();
    auto pi = prefix_function(pattern);
    vector<vector<int>> aut(m + 1, vector<int>(alphabet, 0));

    for (int i = 0; i <= m; i++){
        int border = i ? pi[i - 1] : 0;
        int next = i < m ? int(pattern[i] - first) : -1;
        for (int c = 0; c < alphabet; c++){
            if (c == next) aut[i][c] = i + 1;
            else if (i > 0) aut[i][c] = aut[border][c];
        }
    }

    return aut;
}

vector<int> prefix_function(const char* s){
    return prefix_function(string(s));
}

int kmp_find(const string& text, const string& pattern, int from = 0){
    return kmp_find<string>(text, pattern, from);
}

int kmp_find(const char* text, const char* pattern, int from = 0){
    return kmp_find(string(text), string(pattern), from);
}

/// Lets a std::string be mixed with a literal, kmp_search(s, "ab"), where the template cannot deduce one Container
vector<int> kmp_search(const string& text, const string& pattern){
    return kmp_search<string>(text, pattern);
}

vector<int> kmp_search(const char* text, const char* pattern){
    return kmp_search(string(text), string(pattern));
}

vector<vector<int>> kmp_automaton(const char* pattern, int alphabet = 26, int first = 'a'){
    return kmp_automaton(string(pattern), alphabet, first);
}

int main(){
    auto count_avoiding = [](const vector<vector<int>>& aut, int n){
        int m = aut.size() - 1;
        vector<long long> dp(m + 1, 0);
        dp[0] = 1;

        for (int step = 0; step < n; step++){
            vector<long long> nxt(m + 1, 0);
            for (int i = 0; i < m; i++){
                for (int to : aut[i]) nxt[to] += dp[i];
            }
            dp = nxt;
        }

        return accumulate(dp.begin(), dp.begin() + m, 0LL);
    };

    assert((prefix_function("ababaca") == vector<int>{0, 0, 1, 2, 3, 0, 1}));
    assert((prefix_function("aaa") == vector<int>{0, 1, 2}));
    assert((prefix_function("abcde") == vector<int>{0, 0, 0, 0, 0}));
    assert((prefix_function("") == vector<int>{}));
    assert((prefix_function(vector<int>{7, 7, 3, 7, 7}) == vector<int>{0, 1, 0, 1, 2}));

    assert((kmp_search("ababcababa", "aba") == vector<int>{0, 5, 7}));
    assert((kmp_search("hello world", "o") == vector<int>{4, 7}));
    assert((kmp_search("aaaaaa", "aa") == vector<int>{0, 1, 2, 3, 4}));
    assert((kmp_search("test", "xyz") == vector<int>{}));
    assert((kmp_search("", "a") == vector<int>{}));
    assert((kmp_search(string("abab"), "ab") == vector<int>{0, 2}));
    assert((kmp_search("abab", string("bab")) == vector<int>{1}));
    assert((kmp_search(vector<int>{1, 2, 1, 2, 1}, vector<int>{1, 2, 1}) == vector<int>{0, 2}));

    assert(kmp_find("ababcababa", "aba") == 0);
    assert(kmp_find("ababcababa", "aba", 1) == 5);
    assert(kmp_find("ababcababa", "aba", 8) == -1);
    assert(kmp_find(string("abab"), "bab") == 1);

    assert((kmp_automaton("aba", 2) == vector<vector<int>>{{1, 0}, {1, 2}, {3, 0}, {1, 2}}));
    assert((kmp_automaton("aab", 2) == vector<vector<int>>{{1, 0}, {2, 0}, {2, 3}, {1, 0}}));
    assert((kmp_automaton("", 3) == vector<vector<int>>{{0, 0, 0}}));
    assert((kmp_automaton(vector<int>{1, 1}, 2, 0) == vector<vector<int>>{{0, 1}, {0, 2}, {0, 2}}));
    assert(count_avoiding(kmp_automaton("aa", 2), 4) == 8);
    assert(count_avoiding(kmp_automaton("aa", 2), 10) == 144);
    assert(count_avoiding(kmp_automaton("ab", 26), 3) == 17576 - 2 * 26);

    return 0;
}
