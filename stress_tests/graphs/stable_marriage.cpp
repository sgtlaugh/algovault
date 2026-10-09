#include "../common.h"

#define main library_main
#include "../../code_library/graphs/stable_marriage.cpp"
#undef main

/// The original Library implementation (proposers are w, receivers m, 1-indexed), unchanged, as a second reference
namespace shohag{
    const int N = 1010;
    int m[N][N], w[N][N], id[N], marry[N], ans[N], pref[N][N];
    void SM(int n) {
      queue<int>q;
      for(int i = 1; i <= n; i++) {
        q.push(i);
        id[i] = 1;
      }
      memset(marry, -1, sizeof marry);
      for(int i = 1; i <= n; i++) {
        for(int j = 1; j <= n; j++) pref[i][m[i][j]] = j;
      }
      while(!q.empty()) {
        int cw = q.front();
        q.pop();
        int cm = w[cw][id[cw]];
        if(marry[cm] != -1) {
          if(pref[cm][marry[cm]] > pref[cm][cw]) {
            id[marry[cm]]++;
            q.push(marry[cm]);
            marry[cm] = cw;
          } else {
            id[cw]++;
            q.push(cw);
          }
        } else marry[cm] = cw;
      }
      for(int i = 1; i <= n; i++) ans[i] = marry[i];
    }
}

vector<int> run_original(const vector<vector<int>>& proposer_pref, const vector<vector<int>>& receiver_pref){
    int n = proposer_pref.size();
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            shohag::w[i + 1][j + 1] = proposer_pref[i][j] + 1;
            shohag::m[i + 1][j + 1] = receiver_pref[i][j] + 1;
        }
    }
    shohag::SM(n);

    vector<int> match(n, -1);
    for (int r = 1; r <= n; r++) match[shohag::ans[r] - 1] = r - 1;
    return match;
}

vector<vector<int>> random_prefs(int n){
    vector<vector<int>> pref(n, vector<int>(n));
    for (auto& row : pref){
        iota(row.begin(), row.end(), 0);
        shuffle(row.begin(), row.end(), stress::rng());
    }
    return pref;
}

vector<vector<int>> ranks(const vector<vector<int>>& pref){
    int n = pref.size();
    vector<vector<int>> rank(n, vector<int>(n));
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++) rank[i][pref[i][j]] = j;
    }
    return rank;
}

bool is_stable(const vector<int>& match, const vector<vector<int>>& rank_p, const vector<vector<int>>& rank_r){
    int n = match.size();
    vector<int> partner(n, -1);
    for (int p = 0; p < n; p++){
        if (match[p] < 0 || match[p] >= n || partner[match[p]] != -1) return false;
        partner[match[p]] = p;
    }

    for (int p = 0; p < n; p++){
        for (int r = 0; r < n; r++){
            if (rank_p[p][r] < rank_p[p][match[p]] && rank_r[r][p] < rank_r[r][partner[r]]) return false;
        }
    }
    return true;
}

void check_brute(const vector<vector<int>>& proposer_pref, const vector<vector<int>>& receiver_pref){
    int n = proposer_pref.size();
    auto rank_p = ranks(proposer_pref), rank_r = ranks(receiver_pref);
    auto res = stable_marriage(proposer_pref, receiver_pref);
    assert(is_stable(res, rank_p, rank_r));

    vector<int> perm(n), best(n, n);
    iota(perm.begin(), perm.end(), 0);
    do{
        if (!is_stable(perm, rank_p, rank_r)) continue;
        for (int p = 0; p < n; p++) best[p] = min(best[p], rank_p[p][perm[p]]);
    } while (next_permutation(perm.begin(), perm.end()));

    for (int p = 0; p < n; p++) assert(rank_p[p][res[p]] == best[p]);
}

int main(){
    for (int n = 0; n <= 3; n++){
        int k = 1;
        for (int i = 2; i <= n; i++) k *= i;
        int total = 1;
        for (int i = 0; i < 2 * n; i++) total *= k;

        for (int code = 0; code < total; code++){
            vector<vector<int>> pref(2 * n, vector<int>(n));
            int c = code;
            for (auto& row : pref){
                iota(row.begin(), row.end(), 0);
                for (int s = c % k; s > 0; s--) next_permutation(row.begin(), row.end());
                c /= k;
            }
            vector<vector<int>> proposer_pref(pref.begin(), pref.begin() + n), receiver_pref(pref.begin() + n, pref.end());
            check_brute(proposer_pref, receiver_pref);
        }
    }

    for (long long it = 0; it < stress::scaled(1000); it++){
        int n = stress::rand_int(1, 7);
        check_brute(random_prefs(n), random_prefs(n));
    }

    for (long long it = 0; it < stress::scaled(60); it++){
        int n = stress::rand_int(1, 1000);
        auto proposer_pref = random_prefs(n), receiver_pref = random_prefs(n);
        auto res = stable_marriage(proposer_pref, receiver_pref);
        assert(is_stable(res, ranks(proposer_pref), ranks(receiver_pref)));
        assert(res == run_original(proposer_pref, receiver_pref));
    }

    /// Identical proposer lists against reversed receiver lists: n^2 / 2 proposals, the O(n^2) worst case
    int n = 3000;
    vector<vector<int>> same(n, vector<int>(n)), reversed(n, vector<int>(n));
    for (int i = 0; i < n; i++){
        iota(same[i].begin(), same[i].end(), 0);
        iota(reversed[i].rbegin(), reversed[i].rend(), 0);
    }
    auto res = stable_marriage(same, reversed);
    assert(is_stable(res, ranks(same), ranks(reversed)));
    for (int p = 0; p < n; p++) assert(res[p] == n - 1 - p);

    return 0;
}
