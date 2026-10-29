/***
 *
 * Maximum flow with highest-label push-relabel (HLPP)
 * Max flow and min cut on nodes 0 to n - 1, a drop-in for maxflow.cpp's FlowGraph
 *
 * Complexity: O(n^2 sqrt(m)) for m edges
 *
 * HLPP<T = long long> g(n, src, sink): empty graph with edge capacities of type T, src != sink
 * g.add_directed_edge(u, v, cap): adds u -> v and returns its id into g.E
 * g.add_undirected_edge(u, v, cap): adds u -> v and v -> u, both with capacity cap
 * g.maxflow(): pushes as much flow as possible and returns the flow this call added
 * g.cut_edges(): ids into E of the edges leaving the source side, their capacities sum to the max flow
 * g.in_source_side(u): after maxflow(), true if src reaches u in the residual graph, the smallest min cut side
 * g.E[id].flow: flow on edge id after maxflow()
 *
 * Notes:
 *   - Capacities of edges touching src must sum to at most the max of T, as maxflow() pushes them at once, asserted
 *   - Unlike FlowGraph, maxflow() takes no limit
 *   - Each maxflow() call continues from the flow already pushed, so it can be called again after adding edges
 *   - Faster than FlowGraph where Dinic needs many phases: 3x on a 300 x 300 grid, 9x on a 100-layer graph
 *   - On par with FlowGraph on dense graphs, 1.1x to 1.3x slower on random sparse ones
 *   - Bipartite matching: 1.7x slower than FlowGraph on random graphs, 6x slower to 2.5x faster on Library Checker's
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T = long long>
struct HLPP{
    struct Edge{
        int u, v;
        T cap, flow;
    };

    int n, src, sink;
    int highest, max_level, work;
    vector<Edge> E;
    vector<int> start, to, rev, edge_id;
    vector<int> height, cur, active, next_active, prev_level, next_level, order;
    vector<T> res, excess;
    vector<bool> side;

    HLPP(int n, int src, int sink): n(n), src(src), sink(sink), side(n) {}

    int add_directed_edge(int u, int v, T cap){
        E.push_back({u, v, cap, 0});
        E.push_back({v, u, 0, 0});
        return E.size() - 2;
    }

    void add_undirected_edge(int u, int v, T cap){
        add_directed_edge(u, v, cap);
        add_directed_edge(v, u, cap);
    }

    T maxflow(){
        assert(src != sink);
        build_residual();
        excess.assign(n, 0);

        T emitted = 0;
        for (int i = start[src]; i < start[src + 1]; i++){
            bool overflow = __builtin_add_overflow(emitted, res[i], &emitted);
            assert(!overflow);
            excess[to[i]] += res[i], res[rev[i]] += res[i], res[i] = 0;
        }

        /// The first pass leaves a max preflow, whose leftover excess the second pass returns to src
        push_relabel(sink, src);
        push_relabel(src, sink);

        for (int i = 0; i < (int)res.size(); i++) E[edge_id[i]].flow = E[edge_id[i]].cap - res[i];
        find_source_side();
        return excess[sink];
    }

    vector<int> cut_edges() const{
        vector<int> ids;
        for (int id = 0; id < (int)E.size(); id++){
            if (E[id].cap > 0 && in_source_side(E[id].u) && !in_source_side(E[id].v)) ids.push_back(id);
        }

        return ids;
    }

    bool in_source_side(int u) const{
        return side[u];
    }

    /// Arcs of u sit at start[u] .. start[u + 1] - 1 in id order, counting sorted: adjacency lists cost 2x on sparse
    void build_residual(){
        int m = E.size();
        vector<int> pos(m);

        start.assign(n + 1, 0), to.resize(m), rev.resize(m), edge_id.resize(m), res.resize(m);

        for (auto& e : E) start[e.u + 1]++;
        for (int u = 0; u < n; u++) start[u + 1] += start[u];

        vector<int> fill_at(start.begin(), start.end() - 1);
        for (int id = 0; id < m; id++){
            int i = pos[id] = fill_at[E[id].u]++;
            edge_id[i] = id, to[i] = E[id].v, res[i] = E[id].cap - E[id].flow;
        }
        for (int i = 0; i < m; i++) rev[i] = pos[edge_id[i] ^ 1];
    }

    /// Moves every excess outside {target, blocked} into target; blocked keeps height n so nothing enters it
    void push_relabel(int target, int blocked){
        global_relabel(target, blocked);

        while (highest >= 0){
            int u = active[highest];
            if (u == -1){
                highest--;
                continue;
            }

            active[highest] = next_active[u];
            discharge(u, target);
            /// Global relabel once relabels scan 4n + m arcs: 1n to 8n plus m measure alike, 6n alone is 3x slower
            if (work > 4 * n + (int)res.size()) global_relabel(target, blocked);
        }
    }

    /// height[u] is the BFS distance from u to target in the residual graph, n when target is unreachable
    void global_relabel(int target, int blocked){
        height.assign(n, n), cur.resize(n), active.assign(n, -1), next_active.resize(n);
        prev_level.resize(2 * n), next_level.resize(2 * n);
        for (int k = n; k < 2 * n; k++) prev_level[k] = next_level[k] = k;

        height[target] = 0, order.assign(1, target);
        for (int k = 0; k < (int)order.size(); k++){
            int v = order[k];
            for (int i = start[v]; i < start[v + 1]; i++){
                int u = to[i];
                if (height[u] == n && u != blocked && res[rev[i]] > 0) height[u] = height[v] + 1, order.push_back(u);
            }
        }

        highest = -1, work = 0;
        for (int u : order){
            cur[u] = start[u], max_level = height[u];
            link_level(u);
            if (u != target && excess[u] > 0) activate(u), highest = height[u];
        }
    }

    /// Node n + k heads the circular list of the nodes at height k, which the gap heuristic empties in one sweep
    void link_level(int u){
        int head = n + height[u];
        prev_level[u] = head, next_level[u] = next_level[head];
        prev_level[next_level[head]] = u, next_level[head] = u;
    }

    void activate(int u){
        next_active[u] = active[height[u]], active[height[u]] = u;
    }

    void discharge(int u, int target){
        while (excess[u] > 0){
            if (cur[u] == start[u + 1]){
                relabel(u);
                if (height[u] == n) return;
                highest = height[u];
                continue;
            }

            int i = cur[u], v = to[i];
            if (res[i] > 0 && height[u] == height[v] + 1){
                T pushed = min(excess[u], res[i]);
                if (v != target && excess[v] == 0) activate(v);
                excess[u] -= pushed, excess[v] += pushed, res[i] -= pushed, res[rev[i]] += pushed;
            }
            else cur[u]++;
        }
    }

    /// Gap heuristic: when u leaves the last node at its height, every node above it can no longer reach target
    void relabel(int u){
        int old = height[u];
        prev_level[next_level[u]] = prev_level[u], next_level[prev_level[u]] = next_level[u];

        height[u] = n;
        for (int i = start[u]; i < start[u + 1]; i++){
            if (res[i] > 0 && height[to[i]] + 1 < height[u]) height[u] = height[to[i]] + 1, cur[u] = i;
        }
        work += start[u + 1] - start[u] + 12;

        if (next_level[n + old] == n + old){
            for (int k = old + 1; k <= max_level; k++){
                for (int v = next_level[n + k]; v != n + k; v = next_level[v]) height[v] = n;
                prev_level[n + k] = next_level[n + k] = n + k;
            }
            height[u] = n, max_level = old - 1;
        }
        else if (height[u] < n) link_level(u), max_level = max(max_level, height[u]);
    }

    void find_source_side(){
        fill(side.begin(), side.end(), false);
        order.assign(1, src);
        side[src] = true;

        for (int k = 0; k < (int)order.size(); k++){
            int u = order[k];
            for (int i = start[u]; i < start[u + 1]; i++){
                if (res[i] > 0 && !side[to[i]]) side[to[i]] = true, order.push_back(to[i]);
            }
        }
    }
};

int main(){
    /// Max flow 5 from 0 to 3, the two edges out of 0 are the min cut
    HLPP g(4, 0, 3);
    int first = g.add_directed_edge(0, 1, 3);
    int second = g.add_directed_edge(0, 2, 2);
    g.add_directed_edge(1, 2, 4);
    g.add_directed_edge(1, 3, 2);
    g.add_directed_edge(2, 3, 4);
    assert(g.maxflow() == 5);
    assert(g.E[first].flow == 3);
    assert((g.cut_edges() == vector<int>{first, second}));
    assert(g.in_source_side(0));
    assert(!g.in_source_side(1));

    int direct = g.add_directed_edge(0, 3, 4);
    assert(g.maxflow() == 4);
    assert(g.E[direct].flow == 4);

    HLPP<int> small(3, 0, 2);
    small.add_undirected_edge(0, 1, 7);
    small.add_directed_edge(1, 2, 5);
    assert(small.maxflow() == 5);

    return 0;
}
