/***
 *
 * Reachability
 * Which vertices of a directed graph reach which: the full closure, per-vertex counts or offline queries
 *
 * Complexity: O((n + m) n / 64) per call, closure() uses n^2 / 8 bytes of memory, counts() and reaches() O(n + m + q)
 *
 * Reachability g(n): an empty directed graph on vertices 0 to n - 1
 * g.add_edge(u, v): adds the edge u -> v
 * g.closure(): rows of n bits, v is reachable from u iff res[u][v >> 6] >> (v & 63) & 1
 * g.counts(): res[u] is the number of vertices reachable from u
 * g.reaches(queries): res[i] is 1 iff queries[i].second is reachable from queries[i].first
 *
 * Notes:
 *   - Every vertex reaches itself, so counts()[u] >= 1 and reaches({{u, u}}) is {1}
 *   - closure() fits n up to about 25000 in 128 MB; counts() and reaches() take any n
 *   - counts() and reaches() make n / 64 passes over the m edges: at n = 200000, m = 400000 that takes 2.5-5 s
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// BEGIN COPY scc from code_library/graphs/scc.cpp
struct SCC{
    int n, count = 0;
    vector<vector<int>> adj, rev;
    vector<int> comp;

    SCC(int n) : n(n), adj(n), rev(n), comp(n, -1) {}

    void add_edge(int u, int v){
        adj[u].push_back(v);
        rev[v].push_back(u);
    }

    void run(){
        vector<int> order, it(n, 0), stack;
        vector<char> visited(n, 0);

        for (int s = 0; s < n; s++){
            if (visited[s]) continue;
            stack.push_back(s);
            visited[s] = 1;
            while (!stack.empty()){
                int u = stack.back();
                if (it[u] < (int)adj[u].size()){
                    int v = adj[u][it[u]++];
                    if (!visited[v]) visited[v] = 1, stack.push_back(v);
                }
                else{
                    order.push_back(u);
                    stack.pop_back();
                }
            }
        }

        count = 0;
        fill(comp.begin(), comp.end(), -1);
        for (int i = n - 1; i >= 0; i--){
            int s = order[i];
            if (comp[s] != -1) continue;
            stack.push_back(s);
            comp[s] = count;
            while (!stack.empty()){
                int u = stack.back();
                stack.pop_back();
                for (int v : rev[u]){
                    if (comp[v] == -1) comp[v] = count, stack.push_back(v);
                }
            }
            count++;
        }
    }

    vector<vector<int>> condensation() const{
        vector<vector<int>> dag(count);

        for (int u = 0; u < n; u++){
            for (int v : adj[u]){
                if (comp[u] != comp[v]) dag[comp[u]].push_back(comp[v]);
            }
        }
        for (auto& out : dag){
            sort(out.begin(), out.end());
            out.erase(unique(out.begin(), out.end()), out.end());
        }

        return dag;
    }
};
/// END COPY scc

struct Reachability{
    int n;
    SCC scc;

    Reachability(int n) : n(n), scc(n) {}

    void add_edge(int u, int v){
        scc.add_edge(u, v);
    }

    vector<vector<uint64_t>> closure(){
        scc.run();
        vector<vector<int>> dag = scc.condensation();
        vector<int> first(scc.count, -1);
        for (int u = n - 1; u >= 0; u--) first[scc.comp[u]] = u;

        /// Row first[c] gathers component c's reach, successors first, then is copied to c's other vertices
        vector<vector<uint64_t>> res(n, vector<uint64_t>((n + 63) / 64, 0));
        for (int u = 0; u < n; u++) res[first[scc.comp[u]]][u >> 6] |= 1ULL << (u & 63);
        for (int c = scc.count - 1; c >= 0; c--){
            vector<uint64_t>& row = res[first[c]];
            for (int d : dag[c]){
                for (int w = 0; w < (int)row.size(); w++) row[w] |= res[first[d]][w];
            }
        }
        for (int u = 0; u < n; u++) res[u] = res[first[scc.comp[u]]];

        return res;
    }

    vector<int> counts(){
        scc.run();
        vector<int> reached(scc.count, 0);
        sweep(by_component(), [&](int, const vector<uint64_t>& mask, int hi){
            /// Skipping zero masks made counts() 1.4x faster: popcount is a library call without -mpopcnt
            for (int c = 0; c <= hi; c++){
                if (mask[c]) reached[c] += __builtin_popcountll(mask[c]);
            }
        });

        vector<int> res(n);
        for (int u = 0; u < n; u++) res[u] = reached[scc.comp[u]];
        return res;
    }

    vector<char> reaches(const vector<pair<int, int>>& queries){
        scc.run();
        vector<int> targets = by_component(), pos(n);
        for (int i = 0; i < n; i++) pos[targets[i]] = i;

        vector<vector<int>> by_chunk((n + 63) / 64);
        for (int i = 0; i < (int)queries.size(); i++) by_chunk[pos[queries[i].second] >> 6].push_back(i);

        vector<char> res(queries.size(), 0);
        sweep(targets, [&](int chunk, const vector<uint64_t>& mask, int){
            for (int i : by_chunk[chunk]){
                auto [u, v] = queries[i];
                res[i] = mask[scc.comp[u]] >> (pos[v] & 63) & 1;
            }
        });

        return res;
    }

    /// Chunks of vertices sorted by component span few components, so sweep starts low
    vector<int> by_component() const{
        vector<int> res(n);
        iota(res.begin(), res.end(), 0);
        sort(res.begin(), res.end(), [&](int a, int b){ return scc.comp[a] < scc.comp[b]; });
        return res;
    }

    /// After scc.run(), f(chunk, mask, hi): bit i of mask[c] says c reaches targets[64 * chunk + i], 0 above hi
    template <typename F>
    void sweep(const vector<int>& targets, F f) const{
        vector<vector<int>> dag = scc.condensation();
        /// One flat edge array sorted by source descending: reaches() ran 3x faster than walking dag
        vector<pair<int, int>> edges;
        for (int c = scc.count - 1; c >= 0; c--){
            for (int d : dag[c]) edges.push_back({c, d});
        }
        vector<uint64_t> mask(scc.count);

        for (int chunk = 0; 64 * chunk < n; chunk++){
            fill(mask.begin(), mask.end(), 0);
            int hi = 0;
            for (int i = 64 * chunk; i < min(n, 64 * chunk + 64); i++){
                int c = scc.comp[targets[i]];
                mask[c] |= 1ULL << (i & 63), hi = max(hi, c);
            }

            /// Edges go from lower to higher components, so none from above hi leads into the chunk
            int e = partition_point(edges.begin(), edges.end(), [&](pair<int, int> x){ return x.first > hi; }) - edges.begin();
            for (; e < (int)edges.size(); e++) mask[edges[e].first] |= mask[edges[e].second];
            f(chunk, mask, hi);
        }
    }
};

int main(){
    /***
     * 0 -> 1 -> 2 -> 0 is a cycle, 2 -> 3 -> 4, and 5 -> 4
     * {0, 1, 2} reach everything but 5, 3 reaches {3, 4}, 4 only itself, 5 reaches {4, 5}
    ***/
    Reachability g(6);
    for (auto [u, v] : vector<pair<int, int>>{{0, 1}, {1, 2}, {2, 0}, {2, 3}, {3, 4}, {5, 4}}) g.add_edge(u, v);

    assert((g.counts() == vector<int>{5, 5, 5, 2, 1, 2}));
    assert((g.reaches({{1, 0}, {0, 4}, {4, 3}, {5, 4}, {3, 5}, {4, 4}}) == vector<char>{1, 1, 0, 1, 0, 1}));

    auto reach = g.closure();
    assert(reach[2][0] == 0b011111);
    assert(reach[5][0] == 0b110000);
    return 0;
}
