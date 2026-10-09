/***
 *
 * Longest Common Increasing Subsequence
 * Length of the longest strictly increasing sequence that is a subsequence of both A and B
 *
 * Complexity: O(n * m) time, O(m) memory
 *
 * lcis(A, B): works for any sequences whose elements compare with < and ==, 0 if either is empty
 *
 * dp[j] = longest common increasing subsequence ending exactly at B[j]; while scanning B for a fixed A[i],
 * best holds the longest one ending at some B[j'] < A[i], ready to be extended when B[j] == A[i]
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename Sequence>
int lcis(const Sequence& A, const Sequence& B){
    vector<int> dp(B.size(), 0);
    for (const auto& a : A){
        int best = 0;
        for (size_t j = 0; j < B.size(); j++){
            if (B[j] == a) dp[j] = max(dp[j], best + 1);
            else if (B[j] < a) best = max(best, dp[j]);
        }
    }
    return dp.empty() ? 0 : *max_element(dp.begin(), dp.end());
}

int main(){
    assert(lcis(vector<int>{3, 4, 9, 1}, vector<int>{5, 3, 8, 9, 10, 2, 1}) == 2);
    assert(lcis(vector<int>{1, 2, 3, 4}, vector<int>{1, 2, 3, 4}) == 4);
    assert(lcis(vector<int>{4, 3, 2, 1}, vector<int>{4, 3, 2, 1}) == 1);
    assert(lcis(vector<int>{2, 2, 2}, vector<int>{2, 2}) == 1);
    assert(lcis(vector<int>{}, vector<int>{1, 2}) == 0);
    assert(lcis(string("acbdf"), string("abcdf")) == 4);
    return 0;
}
