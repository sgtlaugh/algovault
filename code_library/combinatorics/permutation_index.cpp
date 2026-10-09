/***
 *
 * Permutation Index
 * Lexicographic rank of a permutation of 0..n-1 in O(n) with no allocation, a perfect hash for permutation states
 *
 * Complexity: O(n) per query, no build cost and no tables, 1 <= n <= 20 (20! still fits in long long)
 *
 * PermutationIndex idx(n); idx.index(p): rank of p among all n! permutations in [0, n!), p holds 0..n-1
 * Use it as a dense state id in BFS / DP over permutations, e.g. visited[idx.index(p)]
 * permutation_rank.cpp computes the same rank for any n in O(n^2)
 *
 * rank(p) = sum of c_i * (n - 1 - i)!, c_i = number of values below p[i] not used by p[0..i-1]
 * c_i is one popcount over the bitmask of unused values, the sum is evaluated Horner style
 * index is compiled for the POPCNT instruction (every x86 CPU since 2008): without it GCC calls a software
 * popcount and the old 12 MB table version wins at n <= 10, with it this beats the table 1.5x to 6x at every n
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct PermutationIndex{
    int n;

    PermutationIndex(int n) : n(n){
        assert(1 <= n && n <= 20);
    }

    template <typename Permutation>
    __attribute__((target("popcnt"))) long long index(const Permutation& p) const{
        long long res = 0;
        int rest = (1 << n) - 1;
        for (int i = 0; i < n; i++){
            res = res * (n - i) + __builtin_popcount(rest & ((1 << p[i]) - 1));
            rest ^= 1 << p[i];
        }
        return res;
    }
};

int main(){
    PermutationIndex three(3);
    vector<vector<int>> order = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
    for (int r = 0; r < 6; r++) assert(three.index(order[r]) == r);

    PermutationIndex one(1);
    assert(one.index(vector<int>{0}) == 0);

    PermutationIndex twelve(12);
    vector<int> identity(12), reversed(12);
    iota(identity.begin(), identity.end(), 0);
    iota(reversed.rbegin(), reversed.rend(), 0);
    assert(twelve.index(identity) == 0 && twelve.index(reversed) == 479001599);
    int arr[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 10};
    assert(twelve.index(arr) == 1);

    PermutationIndex twenty(20);
    vector<int> last(20);
    iota(last.rbegin(), last.rend(), 0);
    assert(twenty.index(last) == 2432902008176639999LL);

    return 0;
}
