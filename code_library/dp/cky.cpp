/***
 *
 * CKY (Cocke Younger Kasami)
 * Decides whether a string belongs to a context-free language given in Chomsky normal form
 *
 * Complexity: O(n^3 B + n T) for B binary and T terminal rules, O(n^2 r) memory
 *
 * CKY parser(r): nonterminals 0..r-1
 * parser.add_binary(A, B, C): rule A -> B C
 * parser.add_terminal(A, c): rule A -> c
 * parser.accepts(s, start): true if start derives s, never for the empty string (no epsilon rules in CNF)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct CKY{
    int r;
    vector<pair<int, char>> terminals;
    vector<array<int, 3>> binaries;

    CKY(int r) : r(r) {}

    void add_binary(int a, int b, int c){
        binaries.push_back({a, b, c});
    }

    void add_terminal(int a, char c){
        terminals.push_back({a, c});
    }

    bool accepts(const string& s, int start) const{
        int n = s.size();

        if (n == 0) return false;
        /// One flat table: nested vectors cost n^2 allocations and a pointer chase per rule test
        vector<char> can((size_t)(n + 1) * n * r, 0);
        auto at = [&](int len, int i, int a) -> char&{ return can[((size_t)len * n + i) * r + a]; };

        for (int i = 0; i < n; i++){
            for (auto [a, c] : terminals){
                if (c == s[i]) at(1, i, a) = 1;
            }
        }

        for (int len = 2; len <= n; len++){
            for (int i = 0; i + len <= n; i++){
                for (int k = 1; k < len; k++){
                    for (auto [a, b, c] : binaries){
                        if (at(k, i, b) && at(len - k, i + k, c)) at(len, i, a) = 1;
                    }
                }
            }
        }

        return at(n, 0, start);
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
