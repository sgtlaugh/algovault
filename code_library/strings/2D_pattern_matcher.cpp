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
 * Hashing modulo 2^61 - 1 with random bases, a false match has probability about (n * m) / 2^61
 * Modulo 2^64 is avoided on purpose, Thue-Morse patterns collide there for every odd base
 *

***/


#include <bits/stdtr1c++.h>

using namespace std;

namespace pm{
    const uint64_t mod = (1ULL << 61) - 1;
    mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
    const uint64_t base1 = rng() % (mod / 3) + mod / 3;
    const uint64_t base2 = rng() % (mod / 3) + mod / 3;

    int n, m, r, c;
    vector <uint64_t> dp;
    uint64_t base_pow, pattern_hash, pattern_pow;

    /// Every argument must already be reduced below mod
    inline uint64_t mul(uint64_t a, uint64_t b){
        __uint128_t x = (__uint128_t)a * b;
        uint64_t res = (x & mod) + (x >> 61);
        return res >= mod ? res - mod : res;
    }

    inline uint64_t push(uint64_t h, uint64_t base, uint64_t v){
        h = mul(h, base) + v;
        return h >= mod ? h - mod : h;
    }

    inline uint64_t pop(uint64_t h, uint64_t p, uint64_t v){
        v = mul(p, v);
        return h >= v ? h - v : h + mod - v;
    }

    void build(const vector<string>& text, const vector<string>& pattern) {
        n = text.size(), r = pattern.size();
        m = text[0].length(), c = pattern[0].length();
        for (int i = 0; i < n; i++) assert((int)text[i].length() == m);
        for (int i = 0; i < r; i++) assert((int)pattern[i].length() == c);

        base_pow = 1, pattern_pow = 1, pattern_hash = 0;
        for (int i = 1; i < r; i++) base_pow = mul(base_pow, base1);
        for (int i = 1; i < c; i++) pattern_pow = mul(pattern_pow, base2);

        for (int i = 0; i < r; i++){
            uint64_t h = 0;
            for (int j = 0; j < c; j++) h = push(h, base2, (unsigned char)pattern[i][j]);
            pattern_hash = push(pattern_hash, base1, h);
        }

        dp.assign(n + 1, 0);
        for (int i = 0; i < n; i++){
            for (int j = 0; j < c; j++) dp[i] = push(dp[i], base2, (unsigned char)text[i][j]);
        }
    }

    vector<pair<int, int>> solve(const vector<string>& text, const vector<string>& pattern) {
        assert(!pattern.empty());
        if (pattern.size() > text.size() || pattern[0].size() > text[0].size()) return {};  /// build would read past the text
        build(text, pattern);

        int i, j;
        uint64_t x;
        vector<pair<int, int>> matched_pos;

        for (j = 0; (j + c) <= m; j++){
            for (x = 0, i = 0; i < r; i++) x = push(x, base1, dp[i]);
            for (i = 0; (i + r) <= n; i++){
                if (x == pattern_hash) matched_pos.push_back({i, j});
                x = push(pop(x, base_pow, dp[i]), base1, dp[i + r]);
            }

            for (i = 0; i < n; i++){
                dp[i] = push(pop(dp[i], pattern_pow, (unsigned char)text[i][j]), base2, (unsigned char)text[i][j + c]);
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
