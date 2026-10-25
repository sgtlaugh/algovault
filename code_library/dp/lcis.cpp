/***
 *
 * Longest Common Increasing Subsequence
 * One longest strictly increasing sequence that is a subsequence of both A and B
 *
 * Complexity: O(n m) time, O(m + r) memory, r = number of dp improvements <= m * answer
 *
 * lcis(A, B): the subsequence in the same Sequence type (a string for strings), empty if either is empty
 *
 * Notes:
 *   - Elements need only < and ==
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename Sequence>
Sequence lcis(const Sequence& A, const Sequence& B){
    int m = B.size();
    vector<int> dp(m, 0), last(m, -1);
    vector<pair<int, int>> node;  /// (index in B, parent node), never modified once pushed

    for (const auto& a : A){
        int best = 0, from = -1;
        for (int j = 0; j < m; j++){
            if (B[j] == a && best + 1 > dp[j]) dp[j] = best + 1, last[j] = node.size(), node.push_back({j, from});
            else if (B[j] < a && dp[j] > best) best = dp[j], from = last[j];
        }
    }

    Sequence res;
    int v = m ? last[max_element(dp.begin(), dp.end()) - dp.begin()] : -1;
    for (; v != -1; v = node[v].second) res.push_back(B[node[v].first]);
    reverse(res.begin(), res.end());
    return res;
}

int main(){
    assert((lcis(vector<int>{3, 4, 9, 1}, vector<int>{5, 3, 8, 9, 10, 2, 1}) == vector<int>{3, 9}));
    assert((lcis(vector<int>{1, 2, 3, 4}, vector<int>{1, 2, 3, 4}) == vector<int>{1, 2, 3, 4}));
    assert(lcis(vector<int>{4, 3, 2, 1}, vector<int>{4, 3, 2, 1}).size() == 1);
    assert((lcis(vector<int>{2, 2, 2}, vector<int>{2, 2}) == vector<int>{2}));
    assert(lcis(vector<int>{}, vector<int>{1, 2}).empty());
    assert(lcis(vector<int>{1, 2}, vector<int>{}).empty());
    assert(lcis(vector<int>{1, 2}, vector<int>{3, 4}).empty());

    string s = lcis(string("acbdf"), string("abcdf"));
    assert(s.size() == 4 && (s == "abdf" || s == "acdf"));
    return 0;
}
