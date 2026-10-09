/***
 *
 * 2D Pattern Matcher
 * Finds every occurrence of an r x c pattern in an n x m text by rolling 2D polynomial hashes
 *
 * Complexity: O(r * c) to build, O(n * m) per find
 *
 * PatternMatcher2D matcher(pattern) hashes the pattern once, find(text) can then be called on any number of texts
 * find returns the top-left corners (i, j) of all matches in column-major order (by j, then by i):
 * text[i + k][j + l] == pattern[k][l] for every 0 <= k < r, 0 <= l < c
 * All rows of the pattern must have the same non-zero length, the same for the rows of a text
 * A pattern larger than the text in either dimension has no match
 *
 * Hashing modulo 2^61 - 1 with two random bases per instance
 * A false match has probability of about 3 * n * m * (r + c) / 2^61 (Schwartz-Zippel over the base range)
 * Modulo 2^64 is avoided on purpose, Thue-Morse patterns collide there for every odd base
 *
 * Example:
 *   text = {"ababa", "babab", "abaaa"}, pattern = {"ba", "ab"}
 *   PatternMatcher2D(pattern).find(text) == {{1, 0}, {0, 1}, {0, 3}}
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct PatternMatcher2D{
    static constexpr uint64_t mod = (1ULL << 61) - 1;

    int r, c;
    uint64_t base1, base2, base_pow, pattern_pow, pattern_hash;

    explicit PatternMatcher2D(const vector<string>& pattern){
        r = pattern.size();
        assert(r > 0);
        c = pattern[0].size();
        assert(c > 0);
        for (int i = 0; i < r; i++) assert((int)pattern[i].size() == c);

        mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
        base1 = rng() % (mod / 3) + mod / 3;
        base2 = rng() % (mod / 3) + mod / 3;

        base_pow = 1, pattern_pow = 1, pattern_hash = 0;
        for (int i = 1; i < r; i++) base_pow = mul(base_pow, base1);
        for (int i = 1; i < c; i++) pattern_pow = mul(pattern_pow, base2);

        for (int i = 0; i < r; i++){
            uint64_t h = 0;
            for (int j = 0; j < c; j++) h = push(h, base2, (unsigned char)pattern[i][j]);
            pattern_hash = push(pattern_hash, base1, h);
        }
    }

    vector<pair<int, int>> find(const vector<string>& text) const{
        int n = text.size();
        int m = n > 0 ? text[0].size() : 0;
        for (int i = 0; i < n; i++) assert((int)text[i].size() == m);
        if (r > n || c > m) return {};

        vector<uint64_t> row_hash(n + 1, 0);
        for (int i = 0; i < n; i++){
            for (int j = 0; j < c; j++) row_hash[i] = push(row_hash[i], base2, (unsigned char)text[i][j]);
        }

        vector<pair<int, int>> matched_pos;
        for (int j = 0; j + c <= m; j++){
            uint64_t x = 0;
            for (int i = 0; i < r; i++) x = push(x, base1, row_hash[i]);

            for (int i = 0; i + r <= n; i++){
                if (x == pattern_hash) matched_pos.push_back({i, j});
                x = push(pop(x, base_pow, row_hash[i]), base1, row_hash[i + r]);
            }

            /// text[i][m] is the string's terminating '\0', read once on the last column and never used
            for (int i = 0; i < n; i++){
                row_hash[i] = push(pop(row_hash[i], pattern_pow, (unsigned char)text[i][j]), base2, (unsigned char)text[i][j + c]);
            }
        }

        return matched_pos;
    }

    /// Every argument must already be reduced below mod
    static uint64_t mul(uint64_t a, uint64_t b){
        __uint128_t x = (__uint128_t)a * b;
        uint64_t res = (x & mod) + (x >> 61);
        return res >= mod ? res - mod : res;
    }

    static uint64_t pop(uint64_t h, uint64_t p, uint64_t v){
        v = mul(p, v);
        return h >= v ? h - v : h + mod - v;
    }

    static uint64_t push(uint64_t h, uint64_t base, uint64_t v){
        h = mul(h, base) + v;
        return h >= mod ? h - mod : h;
    }
};

int main(){
    vector<string> text = {"ababa", "babab", "abaaa"};
    PatternMatcher2D matcher({"ba", "ab"});
    assert((matcher.find(text) == vector<pair<int, int>>{{1, 0}, {0, 1}, {0, 3}}));
    assert((matcher.find({"ab", "ba"}) == vector<pair<int, int>>{}));
    assert((matcher.find({"ba", "ab"}) == vector<pair<int, int>>{{0, 0}}));
    assert((matcher.find({"b"}) == vector<pair<int, int>>{}));
    assert((matcher.find({}) == vector<pair<int, int>>{}));

    PatternMatcher2D single({"a"}), column({"x", "y"});
    assert((single.find({"aba", "baa"}) == vector<pair<int, int>>{{0, 0}, {1, 1}, {0, 2}, {1, 2}}));
    assert((column.find({"xyx", "yxy", "xxy"}) == vector<pair<int, int>>{{0, 0}, {0, 2}}));
    assert((single.find({"aaa"}) == vector<pair<int, int>>{{0, 0}, {0, 1}, {0, 2}}));

    return 0;
}
