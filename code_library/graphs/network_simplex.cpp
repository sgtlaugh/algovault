/***
 *
 * Network Simplex
 * Cheapest flow meeting every node's supply or demand, on edges with lower and upper bounds and costs of either sign
 *
 * Complexity: O(n + m) per pivot, with no polynomial bound on the number of pivots
 *
 * NetworkSimplex g(n): empty graph on nodes 0 to n - 1, every supply 0
 * g.add_edge(u, v, cap, cost): adds u -> v carrying 0 to cap units, returns its id
 * g.add_edge(u, v, low, high, cost): adds u -> v carrying low to high units, low <= high of either sign, returns its id
 * g.add_supply(v, amount): node v must send amount more units out than in, a negative amount is a demand
 * g.solve(): the minimum total cost, nullopt when no flow meets every bound and supply
 * g.flow(id): flow on edge id after a feasible solve
 * g.potential[v]: optimal dual, u -> v has cost + potential[u] - potential[v] >= 0 if flow < high, <= 0 if flow > low
 *
 * Notes:
 *   - Requires __int128 for the total cost
 *   - Self loops, parallel edges and negative cycles are allowed
 *   - n max|cost| <= 1e18, and then |potential[v]| <= 2 n max|cost| + 1
 *   - sum |supply| + sum |low| + sum (high - low) must fit in long long
 *   - Several times faster than min_cost_circulation.cpp in practice, which instead guarantees a polynomial bound
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct NetworkSimplex{
    struct Edge{
        int from, to;
        long long low, cap, cost, flow;
    };

    static constexpr long long INF = numeric_limits<long long>::max();

    int n, m = 0;
    vector<Edge> edges;
    vector<long long> supply, potential;
    vector<int> parent, parent_edge, depth, first_child, next_sibling, prev_sibling, pending;
    int next_edge = 0;

    NetworkSimplex(int n) : n(n), supply(n) {}

    int add_edge(int u, int v, long long cap, long long cost){
        return add_edge(u, v, 0, cap, cost);
    }

    /// cap and flow are stored shifted down by low, so the simplex only sees 0 <= flow <= cap
    int add_edge(int u, int v, long long low, long long high, long long cost){
        assert(0 <= u && u < n && 0 <= v && v < n && low <= high);
        edges.resize(m);
        edges.push_back({u, v, low, high - low, cost, 0});
        return m++;
    }

    void add_supply(int v, long long amount){
        assert(0 <= v && v < n);
        supply[v] += amount;
    }

    optional<__int128> solve(){
        int root = n;
        edges.resize(m);
        vector<long long> b = supply;
        long long max_cost = 1;
        /// Negative cost edges start full, at their locally best bound: 2.4x faster on UVa 1659, 5x on a negative ring
        for (auto& e : edges){
            e.flow = e.cost < 0 ? e.cap : 0;
            b[e.from] -= e.low + e.flow, b[e.to] += e.low + e.flow;
            max_cost = max(max_cost, abs(e.cost));
        }

        /// A cycle through the root uses two root arcs, so big > n max_cost / 2 empties them when a flow is feasible
        long long big = n * max_cost + 1;
        potential.assign(n + 1, 0);
        parent.assign(n + 1, root), parent_edge.assign(n + 1, -1), depth.assign(n + 1, 1);
        first_child.assign(n + 1, -1), next_sibling.assign(n + 1, -1), prev_sibling.assign(n + 1, -1);
        parent[root] = -1, depth[root] = 0;

        /// Root arcs point so a zero flow one leads to the root: the tree starts strongly feasible, preventing cycling
        for (int v = 0; v < n; v++){
            if (b[v] >= 0) edges.push_back({v, root, 0, INF, big, b[v]}), potential[v] = -big;
            else edges.push_back({root, v, 0, INF, big, -b[v]}), potential[v] = big;
            parent_edge[v] = m + v;
            attach(v, root);
        }

        next_edge = 0;
        int block = max(10, (int)sqrt(edges.size()));
        while (true){
            int enter = entering_edge(block);
            if (enter == -1) break;
            pivot(enter);
        }
        potential.resize(n);

        /// Also catches a nonzero total supply: empty root arcs would mean conservation at every node
        for (int v = 0; v < n; v++){
            if (edges[m + v].flow > 0) return nullopt;
        }

        __int128 total = 0;
        for (int id = 0; id < m; id++) total += (__int128)edges[id].cost * flow(id);
        return total;
    }

    long long flow(int id) const{
        return edges[id].low + edges[id].flow;
    }

    long long reduced_cost(const Edge& e) const{
        return e.cost + potential[e.from] - potential[e.to];
    }

    /// Block search: scans block edges at a time from where the last search stopped, takes the worst violation found
    int entering_edge(int block){
        int total = edges.size(), enter = -1, in_block = 0;
        long long worst = 0;

        for (int scanned = 0; scanned < total; scanned++){
            const Edge& e = edges[next_edge];
            long long rc = reduced_cost(e), violation = max(e.flow < e.cap ? -rc : 0, e.flow > 0 ? rc : 0);
            if (violation > worst) worst = violation, enter = next_edge;

            next_edge = next_edge + 1 == total ? 0 : next_edge + 1;
            if (++in_block == block){
                if (enter != -1) break;
                in_block = 0;
            }
        }

        return enter;
    }

    /// Pushes flow around the cycle of enter and the tree, the last blocking arc after the apex leaves the tree
    void pivot(int enter){
        Edge& e = edges[enter];
        bool forward = reduced_cost(e) < 0;
        int p = forward ? e.from : e.to, q = forward ? e.to : e.from;

        int apex_p = p, apex_q = q;
        while (apex_p != apex_q){
            if (depth[apex_p] >= depth[apex_q]) apex_p = parent[apex_p];
            else apex_q = parent[apex_q];
        }
        int apex = apex_p;

        /// The cycle runs apex -> p, enter, q -> apex: strict on the p side and <= after keeps the last minimum
        long long delta = INF;
        int leave = -1;
        bool leave_on_p = false;
        for (int x = p; x != apex; x = parent[x]){
            long long r = residual(x, true);
            if (r < delta) delta = r, leave = x, leave_on_p = true;
        }
        long long r_enter = forward ? e.cap - e.flow : e.flow;
        if (r_enter <= delta) delta = r_enter, leave = -1;
        for (int x = q; x != apex; x = parent[x]){
            long long r = residual(x, false);
            if (r <= delta) delta = r, leave = x, leave_on_p = false;
        }

        if (delta > 0){
            e.flow += forward ? delta : -delta;
            for (int x = p; x != apex; x = parent[x]){
                Edge& t = edges[parent_edge[x]];
                t.flow += t.to == x ? delta : -delta;
            }
            for (int x = q; x != apex; x = parent[x]){
                Edge& t = edges[parent_edge[x]];
                t.flow += t.from == x ? delta : -delta;
            }
        }
        if (leave == -1) return;

        /// The path from the near end of enter up to leave flips direction, then hangs from the far end through enter
        int x = leave_on_p ? p : q, new_parent = leave_on_p ? q : p, new_edge = enter;
        while (true){
            int old_parent = parent[x], old_edge = parent_edge[x];
            detach(x);
            attach(x, new_parent);
            parent_edge[x] = new_edge;
            if (x == leave) break;
            new_parent = x, new_edge = old_edge, x = old_parent;
        }
        update_subtree(leave_on_p ? p : q);
    }

    /// Room to push along the tree arc between x and its parent, downward toward x or upward away from it
    long long residual(int x, bool down) const{
        const Edge& t = edges[parent_edge[x]];
        return (t.to == x) == down ? t.cap - t.flow : t.flow;
    }

    void detach(int x){
        if (prev_sibling[x] != -1) next_sibling[prev_sibling[x]] = next_sibling[x];
        else first_child[parent[x]] = next_sibling[x];
        if (next_sibling[x] != -1) prev_sibling[next_sibling[x]] = prev_sibling[x];
    }

    void attach(int x, int new_parent){
        parent[x] = new_parent, prev_sibling[x] = -1, next_sibling[x] = first_child[new_parent];
        if (first_child[new_parent] != -1) prev_sibling[first_child[new_parent]] = x;
        first_child[new_parent] = x;
    }

    /// Tree arcs have reduced cost 0, so each potential follows from its parent's
    void update_subtree(int top){
        pending.assign(1, top);
        while (!pending.empty()){
            int x = pending.back();
            pending.pop_back();
            const Edge& t = edges[parent_edge[x]];
            potential[x] = potential[parent[x]] + (t.to == x ? t.cost : -t.cost);
            depth[x] = depth[parent[x]] + 1;
            for (int c = first_child[x]; c != -1; c = next_sibling[c]) pending.push_back(c);
        }
    }
};

int main(){
    /// 4 units from 0 to 2: 3 along 0-1-2 at cost 2, 1 on 0-2 at cost 3, so potential[2] = potential[0] + 3
    NetworkSimplex g(3);
    g.add_supply(0, 4);
    g.add_supply(2, -4);
    int e01 = g.add_edge(0, 1, 3, 1), e12 = g.add_edge(1, 2, 3, 1), e02 = g.add_edge(0, 2, 5, 3);
    assert(g.solve() == 9);
    assert(g.flow(e01) == 3);
    assert(g.flow(e12) == 3);
    assert(g.flow(e02) == 1);
    assert(g.potential[2] - g.potential[0] == 3);

    /// No supplies, a cycle of cost 4 per unit forced to carry the lower bound of 2
    NetworkSimplex bounded(2);
    int forced = bounded.add_edge(0, 1, 2, 4, 5);
    bounded.add_edge(1, 0, 10, -1);
    assert(bounded.solve() == 8);
    assert(bounded.flow(forced) == 2);

    /// 3 units of demand but only 2 units of capacity
    NetworkSimplex starved(2);
    starved.add_supply(0, 3);
    starved.add_supply(1, -3);
    starved.add_edge(0, 1, 2, 1);
    assert(!starved.solve());
    return 0;
}
