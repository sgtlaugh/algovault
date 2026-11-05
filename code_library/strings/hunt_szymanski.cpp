/***
 *
 * Hunt-Szymanski Algorithm for LCS
 * Longest common subsequence of two sequences, fast when few pairs of positions match
 *
 * Complexity: O((R + N + M) log M), R = number of position pairs where the sequences match, up to N M
 *
 * lcs(A, B): the LCS length of A and B, two strings or two vectors of any comparable values
 *
 * Notes:
 *   - Values are compared with <, so any byte, negative numbers and huge values all work, like permutations up to 1e9
 *   - Use bit_string_lcs.cpp for dense inputs (small alphabet, large R), it is O(N M / 64)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename Seq>
int lcs(const Seq& A, const Seq& B){
    int i, j, k, n = A.size(), m = B.size();
    vector<typename Seq::value_type> values(B.begin(), B.end());
    sort(values.begin(), values.end());
    values.erase(unique(values.begin(), values.end()), values.end());

    vector<vector<int>> adj(values.size());
    for (i = 0; i < m; i++) adj[lower_bound(values.begin(), values.end(), B[i]) - values.begin()].push_back(i);

    /// rank[i] = -1 when A[i] never occurs in B
    vector<int> rank(n), ar(1, -1);
    for (i = 0; i < n; i++){
        auto it = lower_bound(values.begin(), values.end(), A[i]);
        rank[i] = it == values.end() || A[i] < *it ? -1 : it - values.begin();
    }

    for (i = 0; i < n; i++){
        if (rank[i] == -1) continue;

        const vector<int>& pos = adj[rank[i]];
        /// Positions come in decreasing order, so each one lands at or before the slot of the previous
        for (j = (int)pos.size() - 1, k = ar.size(); j >= 0; j--){
            int x = pos[j];
            if (x > ar.back()) ar.push_back(x), k = ar.size();
            else{
                k = lower_bound(ar.begin(), ar.begin() + k, x) - ar.begin();
                ar[k++] = x;
            }
        }
    }

    return ar.size() - 1;
}

int main(){
    /// ABCBDAB and BDCABA share BCBA, and no common subsequence is longer
    assert(lcs(string("ABCBDAB"), string("BDCABA")) == 4);

    /// Permutations with values beyond a byte: 1 3 4 is common
    assert(lcs(vector<int>{5, 1, 1000000, 3, 4}, vector<int>{1, 3, 7, 4, 1000000}) == 3);

    return 0;
}
