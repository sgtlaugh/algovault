/***
 *
 * Flow with lower bounds (demands) on Dinic's algorithm
 * Every edge carries a flow f with low <= f <= high: feasible circulation, feasible s-t flow, max and min s-t flow
 *
 * Complexity: one Dinic run per call on n + 2 nodes and m + n + 2 edges, O(n^2 * m)
 *
 * LowerBoundFlow g(n); g.add_edge(u, v, low, high): directed edge u -> v, 0 <= low <= high, returns its index
 * g.circulation(): true if some flow satisfies every bound with conservation at every node
 * g.feasible_flow(s, t): value of some flow satisfying every bound, conservation everywhere but s and t
 * g.max_flow(s, t), g.min_flow(s, t): the largest / smallest such value
 * The three s-t calls return LowerBoundFlow::INFEASIBLE when no flow satisfies the bounds
 * g.flow(id): flow on edge id after a call that succeeded
 *
 * The value is the net flow out of s and can be negative when edges into s have lower bounds
 * Each call rebuilds the network from the stored edges, so one graph answers any number of queries
 * The sum of all upper bounds must fit in long long
 *
 * Reduction: edge u -> v [low, high] becomes capacity high - low, a super source S feeds each node's
 * incoming lower bounds minus outgoing ones and a super sink T drains the opposite, and a flow exists
 * iff max flow S -> T saturates S; for s-t flows, arcs t -> s and s -> t of unbounded capacity turn it into a circulation
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct LowerBoundFlow{
    static constexpr long long INFEASIBLE = numeric_limits<long long>::min();

    struct Bound{
        int u, v;
        long long low, high;
    };

    struct Arc{
        int to;
        long long cap;
    };

    int n;
    vector<Bound> bounds;
    vector<Arc> arcs;
    vector<vector<int>> adj;
    vector<int> dis, ptr;

    LowerBoundFlow(int n) : n(n), adj(n + 2), dis(n + 2), ptr(n + 2) {}

    int add_edge(int u, int v, long long low, long long high){
        assert(0 <= u && u < n && 0 <= v && v < n && 0 <= low && low <= high);
        bounds.push_back({u, v, low, high});
        return bounds.size() - 1;
    }

    bool circulation(){
        return build();
    }

    long long feasible_flow(int s, int t){
        assert(s != t);
        if (!build(s, t)) return INFEASIBLE;
        return detach_terminals();
    }

    /// Arc 2 * id is edge id's residual arc and arc 2 * id + 1 its reverse, whose capacity is the flow above low
    long long flow(int id) const{
        assert(0 <= id && 2 * id + 1 < (int)arcs.size());
        return bounds[id].low + arcs[2 * id + 1].cap;
    }

    long long max_flow(int s, int t){
        assert(s != t);
        if (!build(s, t)) return INFEASIBLE;

        long long value = detach_terminals();
        return value + dinic(s, t);
    }

    long long min_flow(int s, int t){
        assert(s != t);
        if (!build(s, t)) return INFEASIBLE;

        long long value = detach_terminals();
        return value - dinic(t, s);
    }

    void add_arc(int u, int v, long long cap){
        adj[u].push_back(arcs.size()), arcs.push_back({v, cap});
        adj[v].push_back(arcs.size()), arcs.push_back({u, 0});
    }

    bool bfs(int src, int sink){
        fill(dis.begin(), dis.end(), -1);
        vector<int> queue = {src};
        dis[src] = 0;

        for (int i = 0; i < (int)queue.size(); i++){
            int u = queue[i];
            for (int id : adj[u]){
                if (dis[arcs[id].to] == -1 && arcs[id].cap > 0) dis[arcs[id].to] = dis[u] + 1, queue.push_back(arcs[id].to);
            }
        }
        return dis[sink] != -1;
    }

    /// Iterative: level graphs as deep as n would overflow an 8 MB stack with a recursive DFS
    long long blocking_flow(int src, int sink){
        long long total = 0;
        vector<int> path;
        int u = src;

        while (true){
            if (u == sink){
                long long f = LLONG_MAX;
                for (int id : path) f = min(f, arcs[id].cap);
                for (int id : path) arcs[id].cap -= f, arcs[id ^ 1].cap += f;
                total += f;

                int keep = 0;
                while (arcs[path[keep]].cap > 0) keep++;
                path.resize(keep);
                u = keep ? arcs[path.back()].to : src;
                continue;
            }

            int& i = ptr[u];
            while (i < (int)adj[u].size() && (arcs[adj[u][i]].cap <= 0 || dis[arcs[adj[u][i]].to] != dis[u] + 1)) i++;
            if (i < (int)adj[u].size()){
                path.push_back(adj[u][i]);
                u = arcs[adj[u][i]].to;
                continue;
            }

            if (path.empty()) break;
            u = arcs[path.back() ^ 1].to;
            path.pop_back();
            ptr[u]++;
        }
        return total;
    }

    /// Feeds lower bounds from S = n to T = n + 1, with the s-t arcs last so detach_terminals finds them at the back
    bool build(int s = -1, int t = -1){
        int super_src = n, super_sink = n + 1;
        arcs.clear();
        for (auto& list : adj) list.clear();

        vector<long long> excess(n, 0);
        long long total_high = 0;
        for (auto& b : bounds){
            add_arc(b.u, b.v, b.high - b.low);
            excess[b.v] += b.low, excess[b.u] -= b.low;
            total_high += b.high;
        }

        long long demand = 0;
        for (int x = 0; x < n; x++){
            if (excess[x] > 0) add_arc(super_src, x, excess[x]), demand += excess[x];
            if (excess[x] < 0) add_arc(x, super_sink, -excess[x]);
        }
        if (s != -1) add_arc(t, s, total_high), add_arc(s, t, total_high);

        return dinic(super_src, super_sink) == demand;
    }

    /// Removes the t -> s and s -> t arcs, returning the value they carried; S and T stay saturated so no path crosses them
    long long detach_terminals(){
        int m = arcs.size();
        long long value = arcs[m - 3].cap - arcs[m - 1].cap;
        for (int id = m - 4; id < m; id++) arcs[id].cap = 0;
        return value;
    }

    long long dinic(int src, int sink){
        long long flow = 0;

        while (bfs(src, sink)){
            fill(ptr.begin(), ptr.end(), 0);
            flow += blocking_flow(src, sink);
        }
        return flow;
    }
};

int main(){
    LowerBoundFlow cycle(3);
    int a = cycle.add_edge(0, 1, 1, 3), b = cycle.add_edge(1, 2, 2, 4), c = cycle.add_edge(2, 0, 0, 5);
    assert(cycle.circulation());
    assert(cycle.flow(a) == cycle.flow(b) && cycle.flow(b) == cycle.flow(c) && 2 <= cycle.flow(a) && cycle.flow(a) <= 3);

    LowerBoundFlow stuck(2);
    stuck.add_edge(0, 1, 2, 3), stuck.add_edge(1, 0, 0, 1);
    assert(!stuck.circulation());

    /***
     * s = 0, t = 3: 0 -> 1 [0, 4], 0 -> 2 [0, 2], 1 -> 3 [3, 3], 2 -> 3 [1, 5], 1 -> 2 [0, 1]
     * Max 6: 0 -> 1 = 4 splits into 1 -> 3 = 3 and 1 -> 2 = 1, which joins 0 -> 2 = 2 so 2 -> 3 = 3
     * Min 4: 1 -> 3 is pinned at 3 and 2 -> 3 needs at least 1
    ***/
    LowerBoundFlow g(4);
    int e01 = g.add_edge(0, 1, 0, 4), e02 = g.add_edge(0, 2, 0, 2), e13 = g.add_edge(1, 3, 3, 3), e23 = g.add_edge(2, 3, 1, 5), e12 = g.add_edge(1, 2, 0, 1);
    assert(g.max_flow(0, 3) == 6);
    assert(g.flow(e01) == 4 && g.flow(e02) == 2 && g.flow(e13) == 3 && g.flow(e23) == 3 && g.flow(e12) == 1);
    assert(g.min_flow(0, 3) == 4);
    assert(g.flow(e13) == 3 && g.flow(e23) == 1 && g.flow(e12) == 0 && g.flow(e01) == 3 && g.flow(e02) == 1);
    long long any = g.feasible_flow(0, 3);
    assert(4 <= any && any <= 6);
    assert(!g.circulation());

    LowerBoundFlow backwards(2);
    backwards.add_edge(1, 0, 2, 2), backwards.add_edge(0, 1, 0, 1);
    assert(backwards.max_flow(0, 1) == -1 && backwards.min_flow(0, 1) == -2);
    assert(backwards.max_flow(1, 0) == 2 && backwards.min_flow(1, 0) == 1);

    LowerBoundFlow bottleneck(3);
    bottleneck.add_edge(0, 1, 0, 1), bottleneck.add_edge(1, 2, 2, 3);
    assert(bottleneck.feasible_flow(0, 2) == LowerBoundFlow::INFEASIBLE);
    assert(bottleneck.max_flow(0, 2) == LowerBoundFlow::INFEASIBLE && bottleneck.min_flow(0, 2) == LowerBoundFlow::INFEASIBLE);

    LowerBoundFlow empty(2);
    assert(empty.circulation() && empty.max_flow(0, 1) == 0 && empty.min_flow(1, 0) == 0);

    LowerBoundFlow loop(2);
    int self = loop.add_edge(1, 1, 4, 7);
    assert(loop.circulation() && loop.flow(self) >= 4 && loop.max_flow(0, 1) == 0);

    return 0;
}
