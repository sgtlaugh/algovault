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

/// The plain O(n^3 B) table, as a reference for strings longer than the brute force reaches: which symbols derive s
vector<char> table_derives(int r, const vector<pair<int, char>>& terminals, const vector<array<int, 3>>& binaries, const string& s){
    int n = s.size();

    vector<vector<vector<char>>> can(n + 1, vector<vector<char>>(n, vector<char>(r, 0)));
    for (int i = 0; i < n; i++){
        for (auto [a, c] : terminals) if (c == s[i]) can[1][i][a] = 1;
    }
    for (int len = 2; len <= n; len++){
        for (int i = 0; i + len <= n; i++){
            for (int k = 1; k < len; k++){
                for (auto [a, b, c] : binaries) if (can[k][i][b] && can[len - k][i + k][c]) can[len][i][a] = 1;
            }
        }
    }

    return can[n][0];
}

bool balanced(const string& s){
    int depth = 0;
    for (char ch : s){
        depth += ch == '(' ? 1 : -1;
        if (depth < 0) return false;
    }

    return depth == 0 && !s.empty();
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

    /// Strings up to 200 characters, so splits cross 64-bit words: random grammars against the plain table
    for (long long it = 0; it < stress::scaled(60); it++){
        int r = stress::rand_int(1, 8), alphabet = stress::rand_int(1, 3), n = stress::rand_int(1, 200);
        CKY parser(r);
        vector<pair<int, char>> terminals;
        vector<array<int, 3>> binaries;
        for (int i = stress::rand_int(1, 6); i; i--){
            terminals.push_back({(int)stress::rand_int(0, r - 1), char('a' + stress::rand_int(0, alphabet - 1))});
            parser.add_terminal(terminals.back().first, terminals.back().second);
        }
        for (int i = stress::rand_int(0, 16); i; i--){
            binaries.push_back({(int)stress::rand_int(0, r - 1), (int)stress::rand_int(0, r - 1), (int)stress::rand_int(0, r - 1)});
            parser.add_binary(binaries.back()[0], binaries.back()[1], binaries.back()[2]);
        }

        string s;
        for (int i = 0; i < n; i++) s += char('a' + stress::rand_int(0, alphabet - 1));
        vector<char> derives = table_derives(r, terminals, binaries, s);
        for (int start = 0; start < r; start++) assert(parser.accepts(s, start) == (bool)derives[start]);
    }

    /// Balanced parentheses at every length to 130, near-misses included, against a depth counter
    CKY parens(4);
    const int S = 0, T = 1, L = 2, R = 3;
    parens.add_terminal(L, '('), parens.add_terminal(R, ')');
    parens.add_binary(S, L, R), parens.add_binary(S, L, T), parens.add_binary(S, S, S), parens.add_binary(T, S, R);
    for (int n = 1; n <= 130; n++){
        for (int rep = 0; rep < 4; rep++){
            string s;
            int open = 0;
            for (int i = 0; i < n; i++){
                bool close = open > 0 && (n - i <= open || stress::rand_int(0, 1));
                s += close ? ')' : '(';
                open += close ? -1 : 1;
            }
            if (rep == 3) s[stress::rand_int(0, n - 1)] ^= '(' ^ ')';
            assert(parens.accepts(s, S) == balanced(s));
        }
    }

    return 0;
}
