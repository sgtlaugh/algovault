#include "../common.h"

#define main library_main
#include "../../code_library/dp/cky.cpp"
#undef main

/// Every string a CNF grammar derives with up to max_len characters, by expanding nonterminals bottom up
vector<set<string>> derivable(int r, const vector<pair<int, char>>& terminals, const vector<array<int, 3>>& binaries, int max_len){
    vector<vector<set<string>>> by_len(max_len + 1, vector<set<string>>(r));
    for (auto [a, c] : terminals) by_len[1][a].insert(string(1, c));
    for (int len = 2; len <= max_len; len++){
        for (auto [a, b, c] : binaries){
            for (int k = 1; k < len; k++){
                for (auto& x : by_len[k][b]){
                    for (auto& y : by_len[len - k][c]) by_len[len][a].insert(x + y);
                }
            }
        }
    }
    vector<set<string>> res(r);
    for (int len = 1; len <= max_len; len++){
        for (int a = 0; a < r; a++) res[a].insert(by_len[len][a].begin(), by_len[len][a].end());
    }
    return res;
}

int main(){
    for (long long it = 0; it < stress::scaled(400); it++){
        int r = stress::rand_int(1, 4), alphabet = stress::rand_int(1, 2), max_len = 6;
        CKY parser(r);
        vector<pair<int, char>> terminals;
        vector<array<int, 3>> binaries;
        for (int i = stress::rand_int(1, 4); i; i--){
            terminals.push_back({(int)stress::rand_int(0, r - 1), char('a' + stress::rand_int(0, alphabet - 1))});
            parser.add_terminal(terminals.back().first, terminals.back().second);
        }
        for (int i = stress::rand_int(0, 6); i; i--){
            binaries.push_back({(int)stress::rand_int(0, r - 1), (int)stress::rand_int(0, r - 1), (int)stress::rand_int(0, r - 1)});
            parser.add_binary(binaries.back()[0], binaries.back()[1], binaries.back()[2]);
        }
        auto lang = derivable(r, terminals, binaries, max_len);

        for (int len = 0; len <= max_len; len++){
            for (int mask = 0; mask < (1 << len) * (alphabet == 2 ? 1 : 0) + (alphabet == 1 ? 1 : 0); mask++){
                string s;
                for (int i = 0; i < len; i++) s += alphabet == 2 && (mask >> i & 1) ? 'b' : 'a';
                for (int start = 0; start < r; start++) assert(parser.accepts(s, start) == (lang[start].count(s) > 0));
            }
        }
    }
    return 0;
}
