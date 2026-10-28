/***
 *
 * Longest Common Increasing Subsequence
 * One longest strictly increasing sequence that is a subsequence of both A and B
 *
 * Complexity: O(n * m) time, O(m + r) memory, r = number of dp improvements <= m * answer
 *
 * lcis(A, B): works for any sequences whose elements compare with < and ==, returns the subsequence
 * as the same Sequence type (a string for strings), so its size() is the length; empty if either is empty
 *
 * dp[j] = longest common increasing subsequence ending exactly at B[j]; while scanning B for a fixed A[i],
 * best holds the longest one ending at some B[j'] < A[i], ready to be extended when B[j] == A[i]
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
    int t, n;

    scanf("%d", &t);
    while (t--){
        vector<int> log[2];
        for (auto& v : log){
            scanf("%d", &n);
            v.resize(n);
            for (auto& x : v) scanf("%d", &x);
        }
        printf("%d\n", (int)lcis(log[0], log[1]).size());
    }

    return 0;
}
