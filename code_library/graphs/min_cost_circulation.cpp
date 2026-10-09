/***
 *
 * Minimum Cost Circulation
 * Cheapest flow with zero net flow at every node, negative costs and negative cycles allowed
 *
 * Complexity: O(n^3 log(n C)), C = max |cost|, cost scaling with FIFO push-relabel (Goldberg and Tarjan)
 *
 * MinCostCirculation g(n);                  nodes 0 .. n-1
 * int id = g.add_edge(u, v, cap, cost);     cap >= 0, any cost, self loops and parallel edges allowed
 * long long cost = g.solve();               the minimum cost, 0 when no negative cycle exists
 * g.flow(id)                                flow on the edge in that optimal circulation
 *
 * Needs n <= 1e4 and |cost| <= 1e9: costs are scaled by n + 1 and potentials stay below about 7 n^2 C
 * The total cost sum(flow * cost) and the capacity sum into any node must fit in a long long as well
 *
 * Min cost max flow with negative cycles from s to t: add_edge(t, s, sum of caps out of s, -M), M = 1 + sum |cost|
 *     then max flow = flow on that edge and min cost = solve() + M * max flow, M must respect the cost limit above
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct MinCostCirculation{
    struct Arc{
        int to, rev;
        long long cap, cost;
    };

    static constexpr long long ALPHA = 4;

    int n;
    vector<vector<Arc>> adj;
    vector<pair<int, int>> pos;
    vector<long long> excess, potential;
    vector<int> cur;

    MinCostCirculation(int n) : n(n), adj(n), excess(n), potential(n), cur(n) {}

    int add_edge(int u, int v, long long cap, long long cost){
        assert(0 <= u && u < n && 0 <= v && v < n && cap >= 0);
        pos.push_back({u, (int)adj[u].size()});

        /// Scaled by n + 1: a cycle has at most n arcs, so 1-optimality in scaled units leaves no negative integer cycle
        cost *= n + 1;
        adj[u].push_back({v, (int)adj[v].size() + (u == v), cap, cost});
        adj[v].push_back({u, (int)adj[u].size() - 1, 0, -cost});

        return pos.size() - 1;
    }

    long long flow(int id) const{
        auto [u, i] = pos[id];
        const Arc& a = adj[u][i];
        return adj[a.to][a.rev].cap;
    }

    long long solve(){
        long long eps = 0;
        for (auto& arcs : adj){
            for (auto& a : arcs) eps = max(eps, abs(a.cost));
        }

        fill(potential.begin(), potential.end(), 0);
        while (eps > 1){
            eps = max(1LL, eps / ALPHA);
            refine(eps);
        }

        long long total = 0;
        for (int id = 0; id < (int)pos.size(); id++){
            auto [u, i] = pos[id];
            total += flow(id) * (adj[u][i].cost / (n + 1));
        }
        return total;
    }

    void push(int u, Arc& a, long long amount){
        a.cap -= amount, adj[a.to][a.rev].cap += amount;
        excess[u] -= amount, excess[a.to] += amount;
    }

    long long reduced_cost(int u, const Arc& a) const{
        return potential[u] + a.cost - potential[a.to];
    }

    /// Turns the circulation eps-optimal: saturate every negative reduced cost arc, then discharge the excesses FIFO
    void refine(long long eps){
        for (int u = 0; u < n; u++){
            for (auto& a : adj[u]){
                if (a.cap > 0 && reduced_cost(u, a) < 0) push(u, a, a.cap);
            }
        }

        queue<int> active;
        vector<bool> queued(n, false);
        for (int u = 0; u < n; u++){
            cur[u] = 0;
            if (excess[u] > 0) active.push(u), queued[u] = true;
        }

        while (!active.empty()){
            int u = active.front();
            active.pop(), queued[u] = false;
            while (excess[u] > 0){
                if (cur[u] == (int)adj[u].size()) relabel(u, eps);
                Arc& a = adj[u][cur[u]];
                if (a.cap == 0 || reduced_cost(u, a) >= 0){
                    cur[u]++;
                    continue;
                }

                push(u, a, min(excess[u], a.cap));
                if (excess[a.to] > 0 && !queued[a.to]) active.push(a.to), queued[a.to] = true;
            }
        }
    }

    /// An excess node always has a residual arc: its excess arrived along arcs whose reverses are residual
    void relabel(int u, long long eps){
        long long best = LLONG_MIN;
        for (auto& a : adj[u]){
            if (a.cap > 0) best = max(best, potential[a.to] - a.cost);
        }
        potential[u] = best - eps;
        cur[u] = 0;
    }
};

int main(){
    MinCostCirculation triangle(3);
    int e01 = triangle.add_edge(0, 1, 1, -1), e12 = triangle.add_edge(1, 2, 1, 1), e20 = triangle.add_edge(2, 0, 1, -1);
    assert(triangle.solve() == -1);
    assert(triangle.flow(e01) == 1 && triangle.flow(e12) == 1 && triangle.flow(e20) == 1);

    MinCostCirculation acyclic(3);
    int a01 = acyclic.add_edge(0, 1, 7, -5), a12 = acyclic.add_edge(1, 2, 7, -3), a02 = acyclic.add_edge(0, 2, 7, 4);
    assert(acyclic.solve() == 0);
    assert(acyclic.flow(a01) == 0 && acyclic.flow(a12) == 0 && acyclic.flow(a02) == 0);

    MinCostCirculation loops(1);
    int cheap = loops.add_edge(0, 0, 5, -2), dear = loops.add_edge(0, 0, 4, 3), zero = loops.add_edge(0, 0, 9, 0);
    assert(loops.solve() == -10);
    assert(loops.flow(cheap) == 5 && loops.flow(dear) == 0 && loops.flow(zero) == 0);

    MinCostCirculation empty(4);
    assert(empty.solve() == 0);

    /// Two negative cycles share arc 0 -> 1 of capacity 3: 2 units on 0-1-0 (-3 each) beat 0-1-2-0 (-2 each), 2-3-2 costs +1
    MinCostCirculation shared(4);
    int s01 = shared.add_edge(0, 1, 3, -4), s10 = shared.add_edge(1, 0, 2, 1), s12 = shared.add_edge(1, 2, 2, 1);
    int s20 = shared.add_edge(2, 0, 2, 1), s23 = shared.add_edge(2, 3, 5, 2), s32 = shared.add_edge(3, 2, 5, -1);
    assert(shared.solve() == -8);
    assert(shared.flow(s01) == 3 && shared.flow(s10) == 2 && shared.flow(s12) == 1 && shared.flow(s20) == 1);
    assert(shared.flow(s23) == 0 && shared.flow(s32) == 0);

    /// Min cost max flow recipe on the graph of mcmf_spfa.cpp (https://cp-algorithms.com/graph/edmonds_karp.html): cost 71, flow 10
    MinCostCirculation network(7);
    vector<array<int, 4>> arcs = {{1, 2, 7, 1}, {2, 3, 5, 2}, {3, 6, 8, 3}, {1, 4, 4, 4}, {4, 2, 3, 5}, {4, 5, 2, 4}, {2, 5, 3, 3}, {5, 3, 3, 2}, {5, 6, 5, 1}};
    long long big = 1;
    for (auto [u, v, cap, cost] : arcs) network.add_edge(u, v, cap, cost), big += abs(cost);
    int back = network.add_edge(6, 1, 11, -big);
    long long cost = network.solve();
    assert(network.flow(back) == 10 && cost + big * 10 == 71);

    /// Same recipe with an isolated negative cycle 1-2-1 worth 4 * (-3 + 1): max flow 1 at cost 5, total 5 - 8
    MinCostCirculation detour(4);
    detour.add_edge(1, 2, 4, -3), detour.add_edge(2, 1, 4, 1), detour.add_edge(0, 3, 1, 5);
    long long detour_big = 1 + 3 + 1 + 5;
    int detour_back = detour.add_edge(3, 0, 1, -detour_big);
    assert(detour.solve() + detour_big * detour.flow(detour_back) == -3 && detour.flow(detour_back) == 1);

    return 0;
}
