#include "../common.h"

#define main library_main
#include "../../code_library/dp/cky_weighted.cpp"
#undef main

struct Grammar{
    int r;
    vector<array<long long, 4>> rules;  /// {a, b, c, weight}, weight 0 for an unweighted rule
    vector<vector<pair<int, long long>>> tokens;
};

/// Every parse tree of tokens[i..j] rooted at a, each as the list of its leaf and rule weights
vector<vector<long long>> trees(const Grammar& g, int i, int j, int a){
    vector<vector<long long>> res;
    if (i == j){
        for (auto [s, w] : g.tokens[i]) if (s == a) res.push_back({w});
        return res;
    }

    for (auto [ra, b, c, w] : g.rules){
        if (ra != a) continue;
        for (int k = i; k < j; k++){
            for (auto& x : trees(g, i, k, b)){
                for (auto& y : trees(g, k + 1, j, c)){
                    vector<long long> t = x;
                    t.insert(t.end(), y.begin(), y.end());
                    if (w) t.push_back(w);
                    res.push_back(t);
                }
            }
        }
    }

    return res;
}

struct Bool{
    bool add(bool x, bool y) const{ return x || y; }
    bool mul(bool x, bool y) const{ return x && y; }
};

struct CountMod{
    long long mod;

    long long add(long long x, long long y) const{ return (x + y) % mod; }
    long long mul(long long x, long long y) const{ return x * y % mod; }
};

struct MaxTimes{
    long long add(long long x, long long y) const{ return max(x, y); }
    long long mul(long long x, long long y) const{ return x * y; }
};

struct MinPlus{
    long long add(long long x, long long y) const{ return min(x, y); }
    long long mul(long long x, long long y) const{ return x + y; }
};

/// What parse(tokens) must return, folding each tree with mul and the trees with add; forest[a] is every tree of a
template <typename T, typename Ops, typename Leaf>
vector<pair<int, T>> brute(const vector<vector<vector<long long>>>& forest, Ops ops, Leaf leaf){
    vector<pair<int, T>> res;
    for (int a = 0; a < (int)forest.size(); a++){
        bool any = false;
        T total = T();
        for (auto& t : forest[a]){
            T v = leaf(t[0]);
            for (size_t k = 1; k < t.size(); k++) v = ops.mul(v, leaf(t[k]));
            total = any ? ops.add(total, v) : v;
            any = true;
        }
        if (any) res.push_back({a, total});
    }

    return res;
}

template <typename T, typename Ops, typename Leaf>
vector<pair<int, T>> run(const Grammar& g, Ops ops, Leaf leaf){
    WeightedCKY<T, Ops> parser(g.r, ops);
    for (auto [a, b, c, w] : g.rules){
        if (w) parser.add_binary(a, b, c, leaf(w));
        else parser.add_binary(a, b, c);
    }

    vector<vector<pair<int, T>>> tokens;
    for (auto& options : g.tokens){
        tokens.push_back({});
        for (auto [s, w] : options) tokens.back().push_back({s, leaf(w)});
    }

    return parser.parse(tokens);
}

/// The plain O(n^3 B) table with an optional value per cell and symbol, a reference for longer sequences
template <typename T, typename Ops>
vector<pair<int, T>> table(const Grammar& g, Ops ops){
    int n = g.tokens.size(), r = g.r;
    vector<vector<vector<optional<T>>>> best(n + 1, vector<vector<optional<T>>>(n, vector<optional<T>>(r)));
    auto add = [&](optional<T>& cell, T v){ cell = cell ? ops.add(*cell, v) : v; };
    for (int i = 0; i < n; i++){
        for (auto [s, w] : g.tokens[i]) add(best[1][i][s], w);
    }
    for (int len = 2; len <= n; len++){
        for (int i = 0; i + len <= n; i++){
            for (int k = 1; k < len; k++){
                for (auto [a, b, c, w] : g.rules){
                    auto &x = best[k][i][b], &y = best[len - k][i + k][c];
                    if (!x || !y) continue;
                    T v = ops.mul(*x, *y);
                    add(best[len][i][a], w ? ops.mul(v, w) : v);
                }
            }
        }
    }

    vector<pair<int, T>> res;
    for (int a = 0; a < r && n > 0; a++) if (best[n][0][a]) res.push_back({a, *best[n][0][a]});
    return res;
}

Grammar random_grammar(int n, int r, int rules, int options, int max_weight){
    Grammar g{r, {}, {}};
    for (int k = 0; k < rules; k++){
        long long w = stress::rand_int(0, 2) ? 0 : stress::rand_int(1, max_weight);
        g.rules.push_back({stress::rand_int(0, r - 1), stress::rand_int(0, r - 1), stress::rand_int(0, r - 1), w});
    }
    for (int i = 0; i < n; i++){
        g.tokens.push_back({});
        for (int k = stress::rand_int(0, options); k; k--) g.tokens.back().push_back({(int)stress::rand_int(0, r - 1), stress::rand_int(1, max_weight)});
    }

    return g;
}

int main(){
    auto same = [](long long x){ return x; };
    auto truthy = [](long long){ return true; };

    /// Short sequences against every parse tree, in four semirings, duplicate rules and tokens included
    for (long long it = 0; it < stress::scaled(1500); it++){
        Grammar g = random_grammar(stress::rand_int(0, 6), stress::rand_int(1, 3), stress::rand_int(0, 6), 2, 3);
        long long mod = stress::rand_int(0, 1) ? 7 : 1000000007;
        int n = g.tokens.size();
        vector<vector<vector<long long>>> forest(g.r);
        for (int a = 0; a < g.r && n > 0; a++) forest[a] = trees(g, 0, n - 1, a);

        assert((run<bool>(g, Bool{}, truthy) == brute<bool>(forest, Bool{}, truthy)));
        assert((run<long long>(g, CountMod{mod}, same) == brute<long long>(forest, CountMod{mod}, same)));
        assert((run<long long>(g, MaxTimes{}, same) == brute<long long>(forest, MaxTimes{}, same)));
        assert((run<long long>(g, MinPlus{}, same) == brute<long long>(forest, MinPlus{}, same)));
    }

    /// Longer sequences against the plain table, sparse and dense grammars, so stamps and splits run many rounds
    for (long long it = 0; it < stress::scaled(40); it++){
        int r = stress::rand_int(1, 12);
        Grammar g = random_grammar(stress::rand_int(1, 70), r, stress::rand_int(0, 40), stress::rand_int(1, 3), 5);
        assert((run<long long>(g, CountMod{998244353}, same) == table<long long>(g, CountMod{998244353})));
        assert((run<long long>(g, MinPlus{}, same) == table<long long>(g, MinPlus{})));
    }

    return 0;
}
