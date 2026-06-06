/***
 *
 * Knuth-Morris-Pratt (KMP) Algorithm
 * Fast pattern matching using failure function
 * Finds all occurrences of pattern in text without backtracking
 *
 * Complexity: O(n + m) where n = text length, m = pattern length
 *
 * Failure function: fail[i] = length of longest proper prefix of pattern[0..i] that is also suffix
 * Allows skipping redundant comparisons when mismatch occurs
 *
 * Applications:
 *   - String pattern matching in linear time
 *   - Multiple pattern searches by rebuilding failure function
 *   - Period detection in strings
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// fail[i] = length of longest proper prefix of c[0..i] that is also suffix, -1 for index 0
template<typename Container>
vector<int> kmp_failure(const Container& c){
    int n = c.size();
    vector<int> fail(n);

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

vector<int> kmp_failure(const char* s){
    return kmp_failure(string(s));
}

vector<int> kmp_search(const char* text, const char* pattern){
    return kmp_search(string(text), string(pattern));
}

int main(){
    assert((kmp_failure("ababaca") == vector<int>{-1, -1, 0, 1, 2, -1, 0}));
    assert((kmp_failure("aaa") == vector<int>{-1, 0, 1}));
    assert((kmp_failure("abcde") == vector<int>{-1, -1, -1, -1, -1}));

    assert((kmp_search("ababcababa", "aba") == vector<int>{0, 5, 7}));
    assert((kmp_search("hello world", "o") == vector<int>{4, 7}));
    assert((kmp_search("aaaaaa", "aa") == vector<int>{0, 1, 2, 3, 4}));
    assert((kmp_search("test", "xyz") == vector<int>{}));
    assert((kmp_search("", "a") == vector<int>{}));

    return 0;
}
