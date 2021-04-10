/***
 * 2D Pattern Matcher, n x m text, r x c pattern
 * Returns the top-left coordinates (i, j) of all matched blocks such that
 * For every (k, l), where k >= 0 and k < r and l >= 0 and l < c, (i + k) < n and (j + l) < m and text[i + k][j + l] == pattern[k][l]
 *
 * For example:
 *
 * text:
 * ababa
 * babab
 * abaaa
 *
 * pattern:
 * ba
 * ab
 *
 * returns top left coordinate of all matched positions
 * {{1, 0}, {0, 1}, {0, 3}}
 *
 *
 * Complexity: O(n * m)
 * Algorithm uses hashing modulo 2^64, so might be susceptible to anti-hash attacks
 * Can be modified to use random modulo if required but keeping it simple for now
 *

***/


#include <bits/stdtr1c++.h>

using namespace std;

namespace pm{
    const unsigned long long base1 = 1995433697;
    const unsigned long long base2 = 2117566807;

    int n, m, r, c;
    vector <unsigned long long> dp;
    unsigned long long base_pow, pattern_hash, pattern_pow;

    void build(const vector<string>& text, const vector<string>& pattern) {
        n = text.size(), r = pattern.size();
        m = text[0].length(), c = pattern[0].length();
        for (int i = 0; i < n; i++) assert((int)text[i].length() == m);
        for (int i = 0; i < r; i++) assert((int)pattern[i].length() == c);

        base_pow = 1, pattern_pow = 1, pattern_hash = 0;
        for (int i = 1; i < r; i++) base_pow = base_pow * base1;
        for (int i = 1; i < c; i++) pattern_pow = pattern_pow * base2;

        for (int i = 0; i < r; i++){
            unsigned long long h = 0;
            for (int j = 0; j < c; j++) h = h * base2 + pattern[i][j];
            pattern_hash = pattern_hash * base1 + h;
        }

        dp.resize(n + 1, 0);
        for (int i = 0; i < n; i++){
            for (int j = 0; j < c; j++) dp[i] = dp[i] * base2 + text[i][j];
        }
    }

    vector<pair<int, int>> solve(const vector<string>& text, const vector<string>& pattern) {
        build(text, pattern);

        int i, j;
        unsigned long long h, x, y;
        vector<pair<int, int>> matched_pos;

        for (j = 0; (j + c) <= m; j++){
            for (x = 0, i = 0; i < r; i++) x = x * base1 + dp[i];
            for (i = 0; (i + r) <= n; i++){
                if (x == pattern_hash) matched_pos.push_back({i, j});
                y = x - base_pow * dp[i];
                x = y * base1 + dp[i + r];
            }

            for (i = 0; i < n; i++){
                h = dp[i] - pattern_pow * text[i][j];
                dp[i] = h * base2 + text[i][j + c];
            }
        }

        return matched_pos;
    }
}

int main(){
    vector <string> text = {"ababa", "babab", "abaaa"};
    vector <string> pattern = {"ba", "ab"};

    auto matched_pos = pm::solve(text, pattern);
    vector<pair<int, int>> expected_matched_pos{{1, 0}, {0, 1}, {0, 3}};
    assert(matched_pos == expected_matched_pos);

    return 0;
}
