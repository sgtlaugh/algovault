/***
 *
 * Knuth-Morris-Pratt (KMP) Algorithm
 * Fast pattern matching using failure function
 * Finds all occurrences of pattern in text without backtracking
 *
 * Complexity: O(n + m) where n = text length, m = pattern length
 *
 * Failure function: fail[i] = (length of longest proper prefix of pattern[0..i] that is also suffix) - 1, -1 if none
 * Allows skipping redundant comparisons when mismatch occurs
 *
 * Prefix-function automaton: aut[i][c] = KMP state after reading symbol c in state i
 * State i = length of the longest pattern prefix that is a suffix of the text read so far, 0 <= i <= m
 * State m is a full match, and its transitions continue through the border so overlapping matches count
 * Symbols are pattern[i] - first and must lie in [0, alphabet), build is O(m * alphabet) time and space
 * Use it for string DP: dp over (text position, KMP state) with one table lookup per transition
 *
 *     auto aut = kmp_automaton("aba", 2);
 *     int state = 0;
 *     for (char ch : text) state = aut[state][ch - 'a'];  // state == 3 where an "aba" ends
 *
 * Applications:
 *   - String pattern matching in linear time
 *   - Multiple pattern searches by rebuilding failure function
 *   - Period detection in strings
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// fail[i] = (length of longest proper prefix of c[0..i] that is also suffix) - 1, -1 if none
template<typename Container>
vector<int> kmp_failure(const Container& c){
    int n = c.size();
    vector<int> fail(n);
    if (n == 0) return fail;

    int k = fail[0] = -1;
    for (int i = 1; i < n; i++){
        while (k >= 0 && c[k + 1] != c[i]) k = fail[k];
        if (c[i] == c[k + 1]) k++;
        fail[i] = k;
    }

    return fail;
}

/// Find all occurrences of pattern in text, return starting indices
template<typename Container>
vector<int> kmp_search(const Container& text, const Container& pattern){
    vector<int> positions;
    if (pattern.size() == 0) return positions;

    auto fail = kmp_failure(pattern);
    int k = 0;
    int m = pattern.size();

    for (int i = 0; i < (int)text.size(); i++){
        while (k > 0 && text[i] != pattern[k]) k = fail[k - 1] + 1;
        if (pattern[k] == text[i]) k++;
        if (k == m){
            positions.push_back(i - m + 1);
            k = fail[k - 1] + 1;
        }
    }

    return positions;
}

/// aut[i][c] = longest pattern prefix that is a suffix of pattern[0..i) + c, for 0 <= i <= m
template<typename Container>
vector<vector<int>> kmp_automaton(const Container& pattern, int alphabet = 26, int first = 'a'){
    int m = pattern.size();
    auto fail = kmp_failure(pattern);
    vector<vector<int>> aut(m + 1, vector<int>(alphabet, 0));

    for (int i = 0; i <= m; i++){
        int border = i ? fail[i - 1] + 1 : 0;
        int next = i < m ? int(pattern[i] - first) : -1;
        for (int c = 0; c < alphabet; c++){
            if (c == next) aut[i][c] = i + 1;
            else if (i > 0) aut[i][c] = aut[border][c];
        }
    }

    return aut;
}

vector<int> kmp_failure(const char* s){
    return kmp_failure(string(s));
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

    assert((kmp_failure("ababaca") == vector<int>{-1, -1, 0, 1, 2, -1, 0}));
    assert((kmp_failure("aaa") == vector<int>{-1, 0, 1}));
    assert((kmp_failure("abcde") == vector<int>{-1, -1, -1, -1, -1}));

    assert((kmp_search("ababcababa", "aba") == vector<int>{0, 5, 7}));
    assert((kmp_search("hello world", "o") == vector<int>{4, 7}));
    assert((kmp_search("aaaaaa", "aa") == vector<int>{0, 1, 2, 3, 4}));
    assert((kmp_search("test", "xyz") == vector<int>{}));
    assert((kmp_search("", "a") == vector<int>{}));

    assert((kmp_automaton("aba", 2) == vector<vector<int>>{{1, 0}, {1, 2}, {3, 0}, {1, 2}}));
    assert((kmp_automaton("aab", 2) == vector<vector<int>>{{1, 0}, {2, 0}, {2, 3}, {1, 0}}));
    assert((kmp_automaton("", 3) == vector<vector<int>>{{0, 0, 0}}));
    assert((kmp_automaton(vector<int>{1, 1}, 2, 0) == vector<vector<int>>{{0, 1}, {0, 2}, {0, 2}}));
    assert(count_avoiding(kmp_automaton("aa", 2), 4) == 8);
    assert(count_avoiding(kmp_automaton("aa", 2), 10) == 144);
    assert(count_avoiding(kmp_automaton("ab", 26), 3) == 17576 - 2 * 26);

    return 0;
}
