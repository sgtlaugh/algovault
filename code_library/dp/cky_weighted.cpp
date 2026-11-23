/***
 *
 * Weighted CKY
 * The best score, the number or the total weight of the parses of a sequence under a grammar in Chomsky normal form
 *
 * Complexity: O(n^3 (r + B)) for B binary rules, far less when few symbols derive each substring
 * Memory: one (symbol, value) pair per symbol deriving each substring, at most r n (n + 1) / 2 pairs
 *
 * WeightedCKY<T, Ops> parser(r, ops = Ops()): nonterminals 0..r-1, values of type T combined by ops
 * parser.add_binary(A, B, C): rule A -> B C, the value of A is mul(value of B, value of C)
 * parser.add_binary(A, B, C, w): rule A -> B C with weight w, the value of A is mul(mul(value of B, value of C), w)
 * parser.parse(tokens): {symbol, value} for every symbol deriving the whole sequence, by increasing symbol
 *
 * Notes:
 *   - tokens[i] lists the {symbol, weight} pairs that position i can be; for a plain string, one symbol per position
 *   - A value is add over all parse trees of mul over the tree's leaf weights and rule weights
 *   - Ops has add(x, y) and mul(x, y), add associative and commutative, mul associative and distributing over add
 *   - Yes/no: add = or, mul = and; counting parses mod p: + and * mod p; best parse: add = max, mul = * or +
 *   - A symbol with no parse is left out, so Ops needs no zero; duplicate rules count as separate parse trees
 *   - Products of many probabilities underflow a double: use logs with add = max and mul = +
 *   - Ties: parse lists symbols in increasing order, so the first best one is the smallest symbol
 *   - An empty sequence gives an empty list; for yes/no alone cky.cpp is 2x faster on sparse grammars, 150x on dense
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T, typename Ops>
struct WeightedCKY{
    struct Rule{
        int a, c;
        bool weighted;
        T w;
    };

    int r;
    Ops ops;
    vector<vector<Rule>> by_left;  /// by_left[b] holds the rules a -> b c

    WeightedCKY(int r, Ops ops = Ops()) : r(r), ops(ops), by_left(r) {}

    void add_binary(int a, int b, int c){
        by_left[b].push_back({a, c, false, T()});
    }

    void add_binary(int a, int b, int c, T w){
        by_left[b].push_back({a, c, true, w});
    }

    vector<pair<int, T>> parse(const vector<vector<pair<int, T>>>& tokens) const{
        int n = tokens.size();
        if (n == 0) return {};

        /// cell[id(i, len)] lists the symbols deriving the len tokens from i, each once, with their values
        vector<vector<pair<int, T>>> cell((size_t)n * (n + 1) / 2);
        auto id = [&](int i, int len){ return (size_t)(len - 1) * n - (size_t)(len - 1) * (len - 2) / 2 + i; };
        vector<T> value(r), right(r);
        vector<int> seen(r, -1);
        vector<long long> right_seen(r, -1);
        int stamp = 0;
        auto add = [&](vector<pair<int, T>>& out, int a, const T& x){
            if (seen[a] != stamp) seen[a] = stamp, value[a] = x, out.push_back({a, T()});
            else value[a] = ops.add(value[a], x);
        };
        auto finish = [&](vector<pair<int, T>>& out){
            for (auto& [a, x] : out) x = value[a];
        };

        for (int i = 0; i < n; i++, stamp++){
            for (auto [a, w] : tokens[i]) add(cell[id(i, 1)], a, w);
            finish(cell[id(i, 1)]);
        }

        /// The right cell is spread into a symbol-indexed array, so each rule finds its right child in O(1)
        long long split = 0;
        for (int len = 2; len <= n; len++){
            for (int i = 0; i + len <= n; i++, stamp++){
                vector<pair<int, T>>& out = cell[id(i, len)];
                for (int k = 1; k < len; k++, split++){
                    const vector<pair<int, T>> &left = cell[id(i, k)], &rest = cell[id(i + k, len - k)];
                    if (left.empty() || rest.empty()) continue;

                    for (auto [c, y] : rest) right[c] = y, right_seen[c] = split;
                    for (auto [b, x] : left){
                        for (const Rule& rule : by_left[b]){
                            if (right_seen[rule.c] != split) continue;

                            T v = ops.mul(x, right[rule.c]);
                            add(out, rule.a, rule.weighted ? ops.mul(v, rule.w) : v);
                        }
                    }
                }
                finish(out);
            }
        }

        vector<pair<int, T>> res = cell[id(0, n)];
        sort(res.begin(), res.end(), [](const pair<int, T>& x, const pair<int, T>& y){ return x.first < y.first; });
        return res;
    }
};

int main(){
    struct CountMod{
        long long mod;

        long long add(long long x, long long y) const{ return (x + y) % mod; }
        long long mul(long long x, long long y) const{ return x * y % mod; }
    };

    struct MaxTimes{
        double add(double x, double y) const{ return max(x, y); }
        double mul(double x, double y) const{ return x * y; }
    };

    /***
     * Balanced parentheses, S -> L R | L T | S S, T -> S R: "()()()" is S S S, which S -> S S splits as
     * (S S) S or S (S S), so S has 2 parse trees and T none
    ***/
    const int S = 0, T = 1, L = 2, R = 3;
    WeightedCKY<long long, CountMod> count(4, CountMod{1000000007});
    count.add_binary(S, L, R), count.add_binary(S, L, T), count.add_binary(S, S, S), count.add_binary(T, S, R);
    vector<vector<pair<int, long long>>> tokens;
    for (char ch : string("()()()")) tokens.push_back({{ch == '(' ? L : R, 1}});
    assert((count.parse(tokens) == vector<pair<int, long long>>{{S, 2}}));

    /***
     * Colours: blue + yellow -> green, yellow + red -> orange, blue + orange -> brown
     * Tokens: blue for sure, then yellow 0.75 or red 0.25, then red for sure
     * brown = blue * (yellow * red) = 0.75; green needs blue next to yellow and leaves red alone, no full parse
    ***/
    const int BLUE = 0, YELLOW = 1, RED = 2, GREEN = 3, ORANGE = 4, BROWN = 5;
    WeightedCKY<double, MaxTimes> colours(6);
    colours.add_binary(GREEN, BLUE, YELLOW), colours.add_binary(ORANGE, YELLOW, RED);
    colours.add_binary(BROWN, BLUE, ORANGE);
    vector<pair<int, double>> top = colours.parse({{{BLUE, 1.0}}, {{YELLOW, 0.75}, {RED, 0.25}}, {{RED, 1.0}}});
    assert((top == vector<pair<int, double>>{{BROWN, 0.75}}));

    /// A rule weight: the same grammar where every orange made costs a factor 0.5
    WeightedCKY<double, MaxTimes> costly(6);
    costly.add_binary(ORANGE, YELLOW, RED, 0.5), costly.add_binary(BROWN, BLUE, ORANGE);
    vector<pair<int, double>> cheaper = costly.parse({{{BLUE, 1.0}}, {{YELLOW, 0.75}}, {{RED, 1.0}}});
    assert((cheaper == vector<pair<int, double>>{{BROWN, 0.375}}));

    return 0;
}
