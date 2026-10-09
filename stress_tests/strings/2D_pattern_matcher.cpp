#include "../common.h"

#define main library_main
#include "../../code_library/strings/2D_pattern_matcher.cpp"
#undef main

vector<string> random_grid(int rows, int cols, int alphabet){
    vector<string> g(rows, string(cols, 'a'));
    for (auto& row : g) for (auto& ch : row) ch = 'a' + stress::rand_int(0, alphabet - 1);
    return g;
}

/// Every match position in column major order, the order find reports them in
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

/// Hashing mod 2^64 maps a Thue-Morse word and its complement to the same value for every odd base
void check_thue_morse(){
    string word, complement;
    for (int i = 0; i < 2048; i++){
        int bit = __builtin_popcount(i) & 1;
        word += "ab"[bit], complement += "ba"[bit];
    }
    assert(PatternMatcher2D({word}).find({complement}).empty());

    vector<string> column, column_complement;
    for (int i = 0; i < 2048; i++) column.push_back(string(1, word[i])), column_complement.push_back(string(1, complement[i]));
    PatternMatcher2D column_matcher(column);
    assert(column_matcher.find(column_complement).empty());
    assert(column_matcher.find(column).size() == 1);
}

/// Planted copies of the pattern in a large text, so matches exist at every scale
void check_large(){
    for (int it = 0; it < stress::scaled(20); it++){
        int n = stress::rand_int(50, 300), m = stress::rand_int(50, 300);
        int r = stress::rand_int(1, 6), c = stress::rand_int(1, 6);
        auto text = random_grid(n, m, 2), pattern = random_grid(r, c, 2);
        for (int t = 0; t < 30; t++){
            int i = stress::rand_int(0, n - r), j = stress::rand_int(0, m - c);
            for (int k = 0; k < r; k++) text[i + k].replace(j, c, pattern[k]);
        }
        assert(PatternMatcher2D(pattern).find(text) == brute(text, pattern));
    }
}

/// Two live matchers queried in interleaved order, each reused on several texts
int main(){
    for (long long it = 0; it < stress::scaled(60000); it++){
        int n = stress::rand_int(1, 8), m = stress::rand_int(1, 8), alphabet = stress::rand_int(1, 3);
        auto pattern1 = random_grid(stress::rand_int(1, n + 2), stress::rand_int(1, m + 2), alphabet);  /// sometimes larger than the text
        auto pattern2 = random_grid(stress::rand_int(1, n), stress::rand_int(1, m), alphabet);
        PatternMatcher2D matcher1(pattern1), matcher2(pattern2);

        for (int q = 0; q < 2; q++){
            auto text = random_grid(n, m, alphabet);
            assert(matcher1.find(text) == brute(text, pattern1));
            assert(matcher2.find(text) == brute(text, pattern2));
        }
    }

    check_large();
    check_thue_morse();
    return 0;
}
