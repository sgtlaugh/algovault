/***
 *
 * De Bruijn Sequence
 * Cyclic sequence over [0, k) in which every length n word appears exactly once
 *
 * Complexity: O(k^n) time and memory
 *
 * de_bruijn(k, n): the lexicographically smallest De Bruijn sequence B(k, n), of length k^n, for k >= 1, n >= 1
 *
 * Notes:
 *   - A window may wrap around the end, append the first n - 1 symbols for a linear sequence holding every word
 *   - Map the symbols to any alphabet afterwards, e.g. "ACGT"[x]
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

vector<int> de_bruijn(int k, int n){
    /// The extension loop below would still grow word to length n, costing O(n) for a length 1 answer
    if (k == 1) return {0};

    long long len = 1;
    for (int i = 0; i < n; i++) len *= k;
    vector<int> seq, word = {-1};
    seq.reserve(len);

    /// FKM: bump the last symbol, extend the Lyndon word of length p periodically to n, drop trailing k - 1 symbols
    while (!word.empty()){
        word.back()++;
        int p = word.size();
        if (n % p == 0) seq.insert(seq.end(), word.begin(), word.end());

        while ((int)word.size() < n) word.push_back(word[word.size() - p]);
        while (!word.empty() && word.back() == k - 1) word.pop_back();
    }

    return seq;
}

int main(){
    assert((de_bruijn(2, 1) == vector<int>{0, 1}));
    assert((de_bruijn(2, 2) == vector<int>{0, 0, 1, 1}));
    assert((de_bruijn(2, 3) == vector<int>{0, 0, 0, 1, 0, 1, 1, 1}));
    assert((de_bruijn(2, 4) == vector<int>{0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1}));
    assert((de_bruijn(3, 2) == vector<int>{0, 0, 1, 0, 2, 1, 1, 2, 2}));
    assert((de_bruijn(4, 2) == vector<int>{0, 0, 1, 0, 2, 0, 3, 1, 1, 2, 1, 3, 2, 2, 3, 3}));
    assert((de_bruijn(5, 1) == vector<int>{0, 1, 2, 3, 4}));

    assert((de_bruijn(1, 1) == vector<int>{0}));
    assert((de_bruijn(1, 6) == vector<int>{0}));
    assert((de_bruijn(1, 1000000000) == vector<int>{0}));

    assert(de_bruijn(2, 20).size() == 1u << 20);
    assert(de_bruijn(10, 6).size() == 1000000u);

    return 0;
}
