/***
 *
 * Palindromic Tree (Eertree)
 * One node per distinct palindromic substring, built online one character at a time
 *
 * Complexity: O(n * SIGMA) memory, amortized O(1) per add, O(n log n) for min_palindromic_factorization
 *
 * PalindromicTree<SIGMA, BASE> tree; characters must lie in [BASE, BASE + SIGMA), defaults to 'a'..'z'
 * tree.add(c): appends c, returns true if a new distinct palindrome appeared
 * tree.distinct(): number of distinct non-empty palindromic substrings so far
 * tree.longest: length of the longest palindromic substring so far
 * tree.occurrences(): occurrences[v] = how many times node v occurs as a substring, nodes 2.. are real palindromes
 * Node v: len[v], link[v] = longest proper palindromic suffix, next[v][c] = c + palindrome(v) + c
 * Nodes 0 and 1 are the roots: node 0 is the empty palindrome, node 1 the imaginary one of length -1
 *
 * min_palindromic_factorization(s): res[i] = fewest palindromes that the prefix s[0..i] splits into
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <int SIGMA = 26, char BASE = 'a'>
struct PalindromicTree{
    string s;
    vector<array<int, SIGMA>> next;
    vector<int> len, link, diff, series_link, ends_here;
    int last = 0, longest = 0;

    PalindromicTree(){
        new_node(0, 1);
        new_node(-1, 0);
    }

    int new_node(int length, int suffix){
        next.push_back({});
        next.back().fill(-1);
        len.push_back(length), link.push_back(suffix), ends_here.push_back(0);
        int v = len.size() - 1;
        diff.push_back(length - (v < 2 ? length : len[suffix]));
        series_link.push_back(v < 2 || diff[v] != diff[suffix] ? suffix : series_link[suffix]);
        return v;
    }

    int extend_from(int v) const{
        int i = s.size() - 1;
        while (i - len[v] - 1 < 0 || s[i - len[v] - 1] != s[i]) v = link[v];
        return v;
    }

    bool add(char ch){
        int c = ch - BASE;
        assert(0 <= c && c < SIGMA);
        s.push_back(ch);

        int v = extend_from(last);
        bool fresh = next[v][c] == -1;
        if (fresh){
            int suffix = len[v] == -1 ? 0 : next[extend_from(link[v])][c];
            int u = new_node(len[v] + 2, suffix);
            next[v][c] = u;
        }
        last = next[v][c];
        ends_here[last]++;
        longest = max(longest, len[last]);
        return fresh;
    }

    int distinct() const{
        return len.size() - 2;
    }

    /// Suffix links always point to older nodes, so one backwards pass pushes counts up
    vector<long long> occurrences() const{
        vector<long long> res(ends_here.begin(), ends_here.end());
        for (int v = len.size() - 1; v >= 2; v--) res[link[v]] += res[v];
        return res;
    }
};

template <int SIGMA = 26, char BASE = 'a'>
vector<int> min_palindromic_factorization(const string& s){
    int n = s.size();
    PalindromicTree<SIGMA, BASE> tree;
    vector<int> best(n + 1, INT_MAX), series(2, 0), res(n);
    best[0] = 0;

    for (int i = 1; i <= n; i++){
        tree.add(s[i - 1]);
        series.resize(tree.len.size(), 0);
        /// Palindromic suffixes split into O(log n) runs with equal length differences, each run reuses its link's answer
        for (int v = tree.last; tree.len[v] > 0; v = tree.series_link[v]){
            series[v] = best[i - (tree.len[tree.series_link[v]] + tree.diff[v])];
            if (tree.diff[v] == tree.diff[tree.link[v]]) series[v] = min(series[v], series[tree.link[v]]);
            best[i] = min(best[i], series[v] + 1);
        }
        res[i - 1] = best[i];
    }

    return res;
}

int main(){
    PalindromicTree<> tree;
    for (char c : string("abacaba")) tree.add(c);
    assert(tree.distinct() == 7);  /// a b c aba aca bacab abacaba
    assert(tree.longest == 7);

    auto occ = tree.occurrences();
    map<int, long long> by_length;
    for (int v = 2; v < (int)tree.len.size(); v++) by_length[tree.len[v]] += occ[v];
    assert(by_length[1] == 7);
    assert(by_length[3] == 3);  /// aba twice, aca once

    /// the prefix abac needs aba | c, the whole word is one palindrome
    assert((min_palindromic_factorization("abacaba") == vector<int>{1, 2, 1, 2, 3, 2, 1}));
    return 0;
}
