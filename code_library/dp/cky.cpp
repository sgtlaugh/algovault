/***
 *
 * CKY (Cocke Younger Kasami)
 * Decides whether a string belongs to a context-free language given in Chomsky normal form
 *
 * Complexity: O(n^3 * |binary rules| + n * |terminal rules|), O(n^2 * r) memory
 *
 * CKY parser(r): nonterminals 0..r-1
 * parser.add_terminal(A, c): rule A -> c
 * parser.add_binary(A, B, C): rule A -> B C
 * parser.accepts(s, start): true if start derives s, the empty string is never accepted (no epsilon rules in CNF)
 *
 * can[len][i][A] = A derives s[i .. i + len - 1], built from shorter spans split at every point
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct CKY{
    int r;
    vector<pair<int, char>> terminals;
    vector<array<int, 3>> binaries;

    CKY(int r) : r(r) {}

    void add_terminal(int a, char c){
        terminals.push_back({a, c});
    }

    void add_binary(int a, int b, int c){
        binaries.push_back({a, b, c});
    }

    bool accepts(const string& s, int start) const{
        int n = s.size();
        if (n == 0) return false;
        vector<vector<vector<char>>> can(n + 1, vector<vector<char>>(n, vector<char>(r, 0)));
        for (int i = 0; i < n; i++){
            for (auto [a, c] : terminals){
                if (c == s[i]) can[1][i][a] = 1;
            }
        }
        for (int len = 2; len <= n; len++){
            for (int i = 0; i + len <= n; i++){
                for (int k = 1; k < len; k++){
                    for (auto [a, b, c] : binaries){
                        if (can[k][i][b] && can[len - k][i + k][c]) can[len][i][a] = 1;
                    }
                }
            }
        }
        return can[n][0][start];
    }
};

int main(){
    /***
     * Balanced parentheses, at least one pair, in CNF:
     * S -> L R | L T | S S,  T -> S R,  L -> '(',  R -> ')'
    ***/
    CKY parens(4);
    const int S = 0, T = 1, L = 2, R = 3;
    parens.add_terminal(L, '('), parens.add_terminal(R, ')');
    parens.add_binary(S, L, R), parens.add_binary(S, L, T), parens.add_binary(S, S, S), parens.add_binary(T, S, R);

    for (string good : {"()", "(())", "()()", "(()())", "((()))()"}) assert(parens.accepts(good, S));
    for (string bad : {"", "(", ")(", "(()", "())(", "((("}) assert(!parens.accepts(bad, S));
    assert(!parens.accepts("()", T));
    return 0;
}
