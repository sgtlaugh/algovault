#include "../common.h"

#define main library_main
#include "../../code_library/dp/hamiltonian_dp.cpp"
#undef main

const long long NONE = LLONG_MAX;

long long order_cost(const vector<vector<long long>>& w, const vector<int>& order, bool cycle){
    int n = order.size();
    long long cost = 0;
    for (int k = 0; k + 1 < n; k++){
        if (w[order[k]][order[k + 1]] == NONE) return NONE;
        cost += w[order[k]][order[k + 1]];
    }
    if (!cycle) return cost;

    if (w[order[n - 1]][order[0]] == NONE) return NONE;
    return cost + w[order[n - 1]][order[0]];
}

/// Every permutation of the vertices, scored directly
long long brute(const vector<vector<long long>>& w, bool cycle){
    vector<int> order(w.size());
    iota(order.begin(), order.end(), 0);

    long long best = NONE;
    do best = min(best, order_cost(w, order, cycle));
    while (next_permutation(order.begin(), order.end()));

    return best;
}

void check(const vector<vector<long long>>& w, const pair<long long, vector<int>>& result, long long expected, bool cycle){
    auto [cost, order] = result;
    assert(cost == expected);
    if (expected == NONE){
        assert(order.empty());
        return;
    }

    vector<int> sorted = order, ids(w.size());
    sort(sorted.begin(), sorted.end());
    iota(ids.begin(), ids.end(), 0);
    assert(sorted == ids);
    assert(!cycle || order[0] == 0);
    assert(order_cost(w, order, cycle) == expected);
}

/// Random permutation gets weight light, every other ordered pair weight heavy, so the cheapest cycle is that permutation
void planted(int n, long long light, long long heavy){
    vector<int> perm(n);
    iota(perm.begin(), perm.end(), 0);
    shuffle(perm.begin() + 1, perm.end(), stress::rng());

    HamiltonianGraph g(n);
    for (int u = 0; u < n; u++){
        for (int v = 0; v < n; v++) g.add_directed_edge(u, v, heavy);
    }
    for (int k = 0; k < n; k++) g.add_directed_edge(perm[k], perm[(k + 1) % n], light);

    assert((g.shortest_cycle() == pair<long long, vector<int>>{light * n, perm}));
    assert(g.shortest_path().first == light * (n - 1));
}

int main(){
    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(1, 8), density = stress::rand_int(1, 100);
        long long range = it % 3 == 0 ? 4e17 : (it % 3 == 1 ? 3 : 100);
        long long lo = it % 2 ? -range : 0;

        HamiltonianGraph g(n);
        vector<vector<long long>> w(n, vector<long long>(n, NONE));
        int edges = stress::rand_int(0, 2 * n * n);
        for (int e = 0; e < edges; e++){
            if (stress::rand_int(1, 100) > density) continue;

            int u = stress::rand_int(0, n - 1), v = stress::rand_int(0, n - 1);
            long long weight = stress::rand_int(lo, range);
            if (stress::rand_int(0, 1)){
                g.add_directed_edge(u, v, weight);
                w[u][v] = min(w[u][v], weight);
            }
            else{
                g.add_edge(u, v, weight);
                w[u][v] = min(w[u][v], weight);
                w[v][u] = min(w[v][u], weight);
            }
        }

        check(w, g.shortest_path(), brute(w, false), false);
        check(w, g.shortest_cycle(), brute(w, true), true);
    }

    for (int it = 0; it < 20; it++) planted(stress::rand_int(9, 14), 1, 1000);
    planted(20, -400000000000000000LL, 400000000000000000LL);

    return 0;
}
