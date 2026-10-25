/***
 *
 * De Bruijn Sequence
 * Cyclic sequence over [0, k) in which every length n word appears exactly once
 *
 * Complexity: O(k^n) time and memory
 *
 * de_bruijn(k, n): the lexicographically smallest De Bruijn sequence B(k, n), of length k^n, for k >= 1, n >= 1
 * Built by the FKM algorithm: concatenates in lexicographic order the Lyndon words whose length divides n
 * Read cyclically, so a window may wrap around the end. For a linear sequence holding every word,
 * append the first n - 1 symbols (length k^n + n - 1)
 * Map the symbols to any alphabet afterwards, e.g. "ACGT"[x]
 *
 * de_bruijn(2, 3) = {0, 0, 0, 1, 0, 1, 1, 1}
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

    /// Walks the prenecklaces in lexicographic order: bump the last symbol, then word is a Lyndon word of
    /// length p, extend it periodically to length n and drop trailing k - 1 symbols to reach the next one
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
    int m, n;

    while (scanf("%d %d", &m, &n) == 2){
        string ring;
        for (int x : de_bruijn(n, m)) ring += char('0' + x);
        puts(ring.c_str());
    }
    return 0;
}
