/***
 *
 * Permutation Index
 * Lexicographic rank of a permutation of 0..n-1 in O(n) with no allocation, a perfect hash for permutation states
 *
 * Complexity: O(n) per query, no build cost and no tables
 *
 * PermutationIndex idx(n): ranks permutations of 0..n-1 for 1 <= n <= 20
 * idx.index(p): rank of p among all n! permutations in [0, n!), p holds 0..n-1
 *
 * Notes:
 *   - index uses the POPCNT instruction, present on every x86 CPU since 2008
 *   - permutation_rank.cpp ranks any distinct values (n <= 20) and unranks for any n
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct PermutationIndex{
    int n;

    PermutationIndex(int n) : n(n){
        assert(1 <= n && n <= 20);
    }

    /// Without target popcnt GCC calls a software popcount and a 12 MB table wins at n <= 10, with it 1.5x-6x faster
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
    /// All 3! permutations in lexicographic order get the indices 0 to 5
    PermutationIndex three(3);
    vector<vector<int>> order = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
    for (int r = 0; r < 6; r++) assert(three.index(order[r]) == r);

    PermutationIndex one(1);
    assert(one.index(vector<int>{0}) == 0);

    PermutationIndex twelve(12);
    vector<int> identity(12), reversed(12);
    iota(identity.begin(), identity.end(), 0);
    iota(reversed.rbegin(), reversed.rend(), 0);
    assert(twelve.index(identity) == 0);
    assert(twelve.index(reversed) == 479001599);  /// 12! - 1
    int arr[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 10};
    assert(twelve.index(arr) == 1);

    PermutationIndex twenty(20);
    vector<int> last(20);
    iota(last.rbegin(), last.rend(), 0);
    assert(twenty.index(last) == 2432902008176639999LL);

    return 0;
}
