/***
 *
 * Permutation Rank
 * Lexicographic rank of a sequence of distinct values, and the permutation of 0..n-1 with a given rank
 *
 * Complexity: O(n^2) for permutation_rank (n <= 20), O(n + 21^2) for kth_permutation
 *
 * permutation_rank(p): 0-based rank of p among all orderings of its values, in [0, n!)
 *   p holds n <= 20 distinct values of any kind (1..n, 0..n-1, arbitrary), duplicates give no meaningful rank
 * kth_permutation(n, k): 0-based k-th permutation of 0..n-1, any n >= 0 and 0 <= k < n! (any k >= 0 when n > 20)
 *   permutation_rank(kth_permutation(n, k)) == k, add 1 to each element for 1-based problems
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename Permutation>
long long permutation_rank(const Permutation& p){
    int n = p.size();
    assert(n <= 20);

    long long res = 0;
    for (int i = 0; i < n; i++){
        int smaller = 0;
        for (int j = i + 1; j < n; j++) smaller += p[j] < p[i];
        res = res * (n - i) + smaller;
    }
    return res;
}

/// Only the last 21 positions can move for k < 2^63, the prefix stays the identity
vector<int> kth_permutation(int n, long long k){
    long long fact[21] = {1};
    for (int i = 1; i <= 20; i++) fact[i] = fact[i - 1] * i;
    assert(n >= 0 && k >= 0 && (n > 20 || k < fact[n]));

    vector<int> res(n);
    iota(res.begin(), res.end(), 0);
    for (int i = max(0, n - 21); i < n; i++){
        int j = i + k / fact[n - 1 - i];
        k %= fact[n - 1 - i];
        rotate(res.begin() + i, res.begin() + j, res.begin() + j + 1);
    }
    return res;
}

int main(){
    int t;
    scanf("%d", &t);
    for (int tc = 1; tc <= t; tc++){
        char s[32];
        long long k;
        scanf("%s %lld", s, &k);
        int n = strlen(s);
        /// s[i] is the p[i]-th letter of the alphabet for the (k - 1)-th permutation p
        vector<int> p = kth_permutation(n, k - 1);
        string alphabet(n, ' ');
        for (int i = 0; i < n; i++) alphabet[p[i]] = s[i];
        printf("Case %d: %s\n", tc, alphabet.c_str());
    }
    return 0;
}
