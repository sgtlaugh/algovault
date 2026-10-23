/***
 *
 * Permutation Rank
 * Lexicographic rank of a sequence of distinct values, and the permutation of 0..n-1 with a given rank
 *
 * Complexity: O(n^2) for permutation_rank (n <= 20), O(n + 21^2) for kth_permutation
 *
 * kth_permutation(n, k): 0-based k-th permutation of 0..n-1, any n >= 0 and 0 <= k < n! (any k >= 0 when n > 20)
 *   permutation_rank(kth_permutation(n, k)) == k, add 1 to each element for 1-based problems
 * permutation_rank(p): 0-based rank of p among all orderings of its values, in [0, n!)
 *   p holds n <= 20 distinct values of any kind (1..n, 0..n-1, arbitrary), duplicates give no meaningful rank
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

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

int main(){
    /// All 3! permutations in lexicographic order, ranked and unranked
    vector<vector<int>> order = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
    for (int r = 0; r < 6; r++){
        assert(permutation_rank(order[r]) == r);
        assert(kth_permutation(3, r) == order[r]);
    }

    assert(permutation_rank(vector<int>{3, 1, 2}) == 4);
    assert(permutation_rank(string("bca")) == 3);
    assert(permutation_rank(vector<int>{}) == 0);
    assert(kth_permutation(0, 0).empty());

    assert((kth_permutation(15, 1000000000006LL) == vector<int>{11, 6, 8, 10, 1, 2, 14, 4, 5, 13, 3, 12, 9, 0, 7}));
    assert(permutation_rank(vector<int>{11, 6, 8, 10, 1, 2, 14, 4, 5, 13, 3, 12, 9, 0, 7}) == 1000000000006LL);

    vector<int> last_two(25);
    iota(last_two.begin(), last_two.end(), 0);
    swap(last_two[23], last_two[24]);
    assert(kth_permutation(25, 1) == last_two);

    vector<int> big = {0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 25, 26, 21, 10, 29, 20, 13, 22, 14, 18, 28, 11, 24, 17, 19, 9, 16, 15, 27, 23};
    assert(kth_permutation(30, LLONG_MAX) == big);

    return 0;
}
