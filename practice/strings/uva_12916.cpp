/***
 *
 * Knuth-Morris-Pratt (KMP) Algorithm
 * Prefix function, all occurrences of a pattern in a text, and the prefix-function automaton
 *
 * Complexity: O(n + m) where n = text length, m = pattern length
 *
 * prefix_function(s)[i] = length of the longest proper prefix of s[0..i] that is also its suffix, 0 if none
 * Smallest period of s: p = n - prefix_function(s)[n - 1], s is a repetition of s[0..p) iff n % p == 0
 * kmp_search(text, pattern): starting indices of every occurrence, overlaps included, none for an empty pattern
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
    int t;
    static char buf[1005];
    scanf("%d", &t);
    while (t--){
        scanf("%s", buf);
        string s = buf;
        int n = s.size(), p = n - prefix_function(s).back();
        printf("%d\n", n % p == 0 ? p : n);
    }
    return 0;
}
