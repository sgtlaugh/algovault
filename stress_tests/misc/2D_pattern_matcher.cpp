#include "../common.h"

#define main library_main
#include "../../code_library/misc/2D_pattern_matcher.cpp"
#undef main

vector<string> random_grid(int rows, int cols){
    vector<string> g(rows, string(cols, 'a'));
    for (auto& row : g) for (auto& ch : row) ch = 'a' + stress::rand_int(0, 1);
    return g;
}

/// Every match position in column major order, the order solve reports them in
vector<pair<int, int>> brute(const vector<string>& text, const vector<string>& pattern){
    vector<pair<int, int>> res;
    int n = text.size(), m = text[0].size(), r = pattern.size(), c = pattern[0].size();
    for (int j = 0; j + c <= m; j++){
        for (int i = 0; i + r <= n; i++){
            bool ok = true;
            for (int k = 0; k < r && ok; k++) ok = text[i + k].compare(j, c, pattern[k]) == 0;
            if (ok) res.push_back({i, j});
        }
    }
    return res;
}

/// Consecutive calls on unrelated grids, so any state left over from the previous call shows up
int main(){
    for (long long it = 0; it < stress::scaled(120000); it++){
        int n = stress::rand_int(1, 8), m = stress::rand_int(1, 8);
        int r = stress::rand_int(1, n + 2), c = stress::rand_int(1, m + 2);  /// sometimes larger than the text
        auto text = random_grid(n, m), pattern = random_grid(r, c);
        assert(pm::solve(text, pattern) == brute(text, pattern));
    }
    return 0;
}
