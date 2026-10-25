/***
 *
 * Stable Marriage (Gale-Shapley)
 * Stable matching between n proposers and n receivers from complete preference lists
 *
 * Complexity: O(n^2) time and memory
 *
 * stable_marriage(proposer_pref, receiver_pref): match, match[p] = receiver matched to proposer p
 *
 * Notes:
 *   - proposer_pref[p] lists receivers from most to least preferred, receiver_pref[r] lists proposers likewise
 *   - Both must be n lists, each a permutation of 0..n - 1
 *   - Stable means no proposer and receiver both prefer each other to their partners
 *   - Each proposer gets its best partner over all stable matchings, each receiver its worst
 *   - Swap the two arguments for the receiver optimal matching
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

vector<int> stable_marriage(const vector<vector<int>>& proposer_pref, const vector<vector<int>>& receiver_pref){
    int n = proposer_pref.size();
    vector<vector<int>> rank(n, vector<int>(n));

    for (int r = 0; r < n; r++){
        for (int i = 0; i < n; i++) rank[r][receiver_pref[r][i]] = i;
    }

    vector<int> match(n, -1), partner(n, -1), next_choice(n, 0), free_proposers(n);
    iota(free_proposers.rbegin(), free_proposers.rend(), 0);
    while (!free_proposers.empty()){
        int p = free_proposers.back();
        int r = proposer_pref[p][next_choice[p]++], q = partner[r];
        if (q != -1 && rank[r][q] < rank[r][p]) continue;

        free_proposers.pop_back();
        if (q != -1) match[q] = -1, free_proposers.push_back(q);
        match[p] = r, partner[r] = p;
    }

    return match;
}

int main(){
    /// Both proposers want receiver 0, who prefers proposer 1, so proposer 0 settles for receiver 1
    assert((stable_marriage({{0, 1}, {0, 1}}, {{1, 0}, {0, 1}}) == vector<int>{1, 0}));

    /// The first choices conflict, so the side that proposes gets them
    vector<vector<int>> pa = {{0, 1}, {1, 0}};
    vector<vector<int>> pb = {{1, 0}, {0, 1}};
    assert((stable_marriage(pa, pb) == vector<int>{0, 1}));
    assert((stable_marriage(pb, pa) == vector<int>{1, 0}));

    return 0;
}
