#include "../common.h"

#define main library_main
#include "../../code_library/graphs/functional_graph.cpp"
#undef main

int walk(const vector<int>& next, int v, long long k){
    for (; k > 0; k--) v = next[v];

    return v;
}

/// Walking from v until a node repeats: the repeated node is v's entry and its first index is v's tail
pair<int, int> brute_entry(const vector<int>& next, int v){
    vector<int> seen(next.size(), -1);
    for (int i = 0; seen[v] == -1; i++) seen[v] = i, v = next[v];

    return {v, seen[v]};
}

void check_all(const vector<int>& next, bool all_pairs){
    int n = next.size();
    FunctionalGraph g(next);

    vector<int> len(n);
    for (int v = 0; v < n; v++){
        auto [e, t] = brute_entry(next, v);
        assert(g.entry[v] == e && g.tail[v] == t);
        assert(g.cycle[v] == g.cycle[e]);
        assert((g.pos[v] == -1) == (t != 0));
        len[v] = g.cycles[g.cycle[v]].size();
    }

    int total = 0;
    for (int c = 0; c < (int)g.cycles.size(); c++){
        const vector<int>& ring = g.cycles[c];
        total += ring.size();
        for (int i = 0; i < (int)ring.size(); i++){
            assert(g.cycle[ring[i]] == c && g.pos[ring[i]] == i);
            assert(next[ring[i]] == ring[(i + 1) % ring.size()]);
        }
    }
    int on_cycle = 0;
    for (int v = 0; v < n; v++) on_cycle += g.tail[v] == 0;
    assert(total == on_cycle);

    for (int v = 0; v < n; v++){
        for (int k = 0; k <= 2 * n + 2; k++) assert(g.jump(v, k) == walk(next, v, k));

        long long big = stress::rand_int(0, (long long)4e18);
        long long reduced = big < g.tail[v] ? big : g.tail[v] + (big - g.tail[v]) % len[v];
        assert(g.jump(v, big) == walk(next, v, reduced));
    }

    int pairs = all_pairs ? n * n : 4 * n;
    for (int it = 0; it < pairs; it++){
        int a = all_pairs ? it / n : stress::rand_int(0, n - 1);
        int b = all_pairs ? it % n : (it % 2 ? walk(next, a, stress::rand_int(0, 2 * n)) : stress::rand_int(0, n - 1));
        long long expected = -1;
        int v = a;
        for (int k = 0; k <= 2 * n; k++, v = next[v]){
            if (v == b){
                expected = k;
                break;
            }
        }
        assert(g.dist(a, b) == expected);
    }
}

vector<int> random_function(int n, int shape){
    vector<int> next(n);
    if (shape == 0){
        for (int& x : next) x = stress::rand_int(0, n - 1);
    }
    else if (shape == 1){
        iota(next.begin(), next.end(), 0);
        shuffle(next.begin(), next.end(), stress::rng());
    }
    else{
        /// A random tree of tails hanging off one small cycle, so tails get long
        int c = stress::rand_int(1, min(n, 3));
        vector<int> label(n);
        iota(label.begin(), label.end(), 0);
        shuffle(label.begin(), label.end(), stress::rng());
        for (int i = 0; i < c; i++) next[label[i]] = label[(i + 1) % c];
        for (int i = c; i < n; i++) next[label[i]] = label[stress::rand_int(max(0, i - 3), i - 1)];
    }

    return next;
}

int main(){
    FunctionalGraph empty(vector<int>{});
    assert(empty.cycles.empty());

    /// Every function on up to 5 nodes, every pair and every small k
    for (int n = 1; n <= 5; n++){
        int total = 1;
        for (int i = 0; i < n; i++) total *= n;
        for (int code = 0; code < total; code++){
            vector<int> next(n);
            for (int i = 0, c = code; i < n; i++, c /= n) next[i] = c % n;
            check_all(next, true);
        }
    }

    for (int it = 0; it < stress::scaled(600); it++){
        int n = stress::rand_int(1, 60);
        check_all(random_function(n, it % 3), n <= 30);
    }

    /// A path of 2e5 nodes into a self loop, jumps checked against the known closed form
    int n = 200000;
    vector<int> next(n);
    for (int v = 0; v < n; v++) next[v] = max(v - 1, 0);
    FunctionalGraph path(next);
    assert(path.tail[n - 1] == n - 1);
    for (int it = 0; it < stress::scaled(20000); it++){
        int v = stress::rand_int(0, n - 1);
        long long k = stress::rand_int(0, 2 * n);
        assert(path.jump(v, k) == max(0LL, v - k));
        assert(path.dist(v, k % n) == (k % n <= v ? v - k % n : -1));
    }

    stress::expect_assert("graphs/functional_graph.cpp", []{ FunctionalGraph bad({0, 2}); });
    stress::expect_assert("graphs/functional_graph.cpp", []{ FunctionalGraph bad({0, -1}); });

    return 0;
}
