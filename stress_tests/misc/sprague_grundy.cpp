#include "../common.h"

#define main library_main
#include "../../code_library/misc/sprague_grundy.cpp"
#undef main

int brute_mex(const set<int>& values){
    int res = 0;
    while (values.count(res)) res++;
    return res;
}

/// Octal game: digit[k] bit 0 takes a whole heap of k, bit 1 takes k leaving one heap, bit 2 leaving two heaps
struct OctalBrute{
    vector<int> digit;
    map<vector<int>, int> memo;

    /// The whole position searched as one game, so the xor of heap values is checked, never assumed
    int grundy(vector<int> heaps){
        sort(heaps.begin(), heaps.end());
        auto it = memo.find(heaps);
        if (it != memo.end()) return it->second;

        set<int> reach;
        for (int h = 0; h < (int)heaps.size(); h++){
            for (int k = 1; k < (int)digit.size(); k++){
                int rest = heaps[h] - k;
                vector<int> next = heaps;
                next.erase(next.begin() + h);
                if ((digit[k] & 1) && rest == 0) reach.insert(grundy(next));
                if ((digit[k] & 2) && rest > 0){
                    next.push_back(rest);
                    reach.insert(grundy(next));
                    next.pop_back();
                }
                for (int a = 1; (digit[k] & 4) && 2 * a <= rest; a++){
                    next.push_back(a), next.push_back(rest - a);
                    reach.insert(grundy(next));
                    next.pop_back(), next.pop_back();
                }
            }
        }

        return memo[heaps] = brute_mex(reach);
    }
};

vector<int> octal_table(int n, const vector<int>& digit){
    return grundy_table(n, [&](int i, const vector<int>& g, Mex& mex){
        for (int k = 1; k < (int)digit.size() && k <= i; k++){
            int rest = i - k;
            if ((digit[k] & 1) && rest == 0) mex.add(0);
            if ((digit[k] & 2) && rest > 0) mex.add(g[rest]);
            for (int a = 1; (digit[k] & 4) && 2 * a <= rest; a++) mex.add(g[a] ^ g[rest - a]);
        }
    });
}

void partitions(int left, int largest, vector<int>& cur, vector<vector<int>>& out){
    out.push_back(cur);
    for (int h = 1; h <= min(left, largest); h++){
        cur.push_back(h);
        partitions(left - h, h, cur, out);
        cur.pop_back();
    }
}

void test_octal(){
    vector<vector<int>> positions;
    vector<int> cur;
    partitions(11, 11, cur, positions);

    for (int it = 0; it < stress::scaled(60); it++){
        OctalBrute brute;
        brute.digit.assign(stress::rand_int(2, 5), 0);
        for (int k = 1; k < (int)brute.digit.size(); k++) brute.digit[k] = stress::rand_int(0, 7);
        if (it == 0) brute.digit = {0, 7, 7};  /// Kayles: 0.77 in octal notation, the 0 for k = 0 unused

        vector<int> g = octal_table(11, brute.digit);
        for (const vector<int>& heaps : positions){
            int x = 0;
            for (int h : heaps) x ^= g[h];
            assert(x == brute.grundy(heaps));
        }
    }
}

/// The 11840 Treblecross recurrence as written by hand in an accepted solution, first value above 63 at n = 1318
void test_treblecross(){
    int n = 2000;
    vector<int> want(n + 1), seen(n + 2, -1);
    for (int i = 1; i <= n; i++){
        for (int j = 0; j < i; j++) seen[want[max(0, j - 2)] ^ want[max(0, i - j - 3)]] = i;
        while (seen[want[i]] == i) want[i]++;
    }

    vector<int> g = grundy_table(n, [](int i, const vector<int>& g, Mex& mex){
        for (int j = 0; j < i; j++) mex.add(g[max(0, j - 2)] ^ g[max(0, i - j - 3)]);
    });
    assert(g == want);
    assert(*max_element(g.begin(), g.end()) > 63);
}

/// Random moves to smaller sizes plus moves into fixed side games, so even size 0 can have moves
void test_random_table(){
    for (int it = 0; it < stress::scaled(300); it++){
        int n = stress::rand_int(0, 40);
        vector<vector<int>> to(n + 1), side(n + 1);
        for (int i = 0; i <= n; i++){
            for (int e = stress::rand_int(0, 6); i > 0 && e > 0; e--) to[i].push_back(stress::rand_int(0, i - 1));
            for (int e = stress::rand_int(0, 2); e > 0; e--) side[i].push_back(stress::rand_int(0, 3));
        }

        vector<int> want(n + 1);
        for (int i = 0; i <= n; i++){
            set<int> reach(side[i].begin(), side[i].end());
            for (int j : to[i]) reach.insert(want[j]);
            want[i] = brute_mex(reach);
        }

        vector<int> g = grundy_table(n, [&](int i, const vector<int>& g, Mex& mex){
            for (int j : to[i]) mex.add(g[j]);
            for (int x : side[i]) mex.add(x);
        });
        assert(g == want);
    }
}

/// Nim heaps, take any positive amount: g[i] == i, far above a 64-bit mask
void test_nim(){
    vector<int> g = grundy_table(500, [](int i, const vector<int>& g, Mex& mex){
        for (int j = 0; j < i; j++) mex.add(g[j]);
    });
    for (int i = 0; i <= 500; i++) assert(g[i] == i);
}

void test_mex(){
    for (int it = 0; it < stress::scaled(2000); it++){
        Mex mex;
        set<int> values;
        int ops = stress::rand_int(1, 200), top = stress::rand_int(1, 300);
        for (int op = 0; op < ops; op++){
            int kind = stress::rand_int(0, 20);
            if (kind == 0){
                mex.clear();
                values.clear();
            } else if (kind <= 3){
                int x = stress::rand_int(0, 5000);
                mex.add(x), values.insert(x);
            } else {
                int x = stress::rand_int(0, top);
                mex.add(x), values.insert(x);
            }
            assert(mex.value() == brute_mex(values));
        }
    }
}

vector<vector<int>> random_dag(int n, int edges){
    vector<int> perm(n);
    iota(perm.begin(), perm.end(), 0);
    shuffle(perm.begin(), perm.end(), stress::rng());

    vector<vector<int>> adj(n);
    for (int e = 0; e < edges && n > 1; e++){
        int a = stress::rand_int(0, n - 2), b = stress::rand_int(a + 1, n - 1);
        adj[perm[a]].push_back(perm[b]);
    }

    return adj;
}

struct TokenBrute{
    vector<vector<int>> adj;
    map<vector<int>, int> memo;

    int grundy(vector<int> tokens){
        sort(tokens.begin(), tokens.end());
        auto it = memo.find(tokens);
        if (it != memo.end()) return it->second;

        set<int> reach;
        for (int t = 0; t < (int)tokens.size(); t++){
            for (int v : adj[tokens[t]]){
                vector<int> next = tokens;
                next[t] = v;
                reach.insert(grundy(next));
            }
        }

        return memo[tokens] = brute_mex(reach);
    }
};

void test_dag(){
    for (int it = 0; it < stress::scaled(400); it++){
        int n = stress::rand_int(1, 7);
        TokenBrute brute;
        brute.adj = random_dag(n, stress::rand_int(0, 14));
        vector<int> g = grundy_dag(brute.adj);

        for (int q = 0; q < 20; q++){
            vector<int> tokens(stress::rand_int(1, 3));
            int x = 0;
            for (int& t : tokens){
                t = stress::rand_int(0, n - 1);
                x ^= g[t];
            }
            assert(x == brute.grundy(tokens));
        }
    }

    for (int it = 0; it < stress::scaled(20); it++){
        int n = stress::rand_int(500, 2000);
        TokenBrute brute;
        brute.adj = random_dag(n, stress::rand_int(0, 6 * n));
        vector<int> g = grundy_dag(brute.adj);
        for (int u = 0; u < n; u++) assert(g[u] == brute.grundy({u}));
    }

    /// Every edge u -> v for u < v: g[u] = n - 1 - u, values up to 199
    int n = 200;
    vector<vector<int>> complete(n);
    for (int u = 0; u < n; u++){
        for (int v = u + 1; v < n; v++) complete[u].push_back(v);
    }
    vector<int> g = grundy_dag(complete);
    for (int u = 0; u < n; u++) assert(g[u] == n - 1 - u);

    stress::expect_assert("misc/sprague_grundy.cpp", []{ grundy_dag({{1}, {2}, {0}}); });
    stress::expect_assert("misc/sprague_grundy.cpp", []{ grundy_dag({{}, {1}}); });
}

int main(){
    test_octal();
    test_random_table();
    test_treblecross();
    test_nim();
    test_mex();
    test_dag();
    return 0;
}
