/***
 *
 * Functional Graph
 * Cycles, tails and long jumps in a graph where every node has exactly one outgoing edge, v -> next[v]
 *
 * Complexity: O(n) to build and memory, O(log n) per dist and jump
 *
 * FunctionalGraph g(next): nodes 0..n-1 with next[v] in [0, n), self loops allowed
 * g.dist(a, b): the fewest steps from a to b, 0 when a = b, -1 if b is never reached
 * g.jump(v, k): the node reached from v after k >= 0 steps, k a long long
 * g.cycles[c]: the nodes of cycle c in walking order, so next[cycles[c][i]] = cycles[c][i + 1] cyclically
 * g.cycle[v]: the cycle that v ends on, cycles are numbered from 0
 * g.entry[v]: the first cycle node reached from v, v itself on a cycle
 * g.pos[v]: v's index in cycles[cycle[v]] for a cycle node, -1 for a tail node
 * g.tail[v]: the number of steps from v to entry[v], 0 on a cycle
 *
 * Notes:
 *   - The cycle length of v is cycles[cycle[v]].size(), and the cycle position it first reaches is pos[entry[v]]
 *   - Iterative, so tails of length 1e6 are safe
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct FunctionalGraph{
    int n;
    vector<int> next;
    vector<vector<int>> cycles;
    vector<int> cycle, entry, pos, tail;
    vector<int> tin, depth_start, by_depth;  /// tail nodes in preorder, grouped by tail: by_depth[depth_start[d]..]

    FunctionalGraph(const vector<int>& next) : n(next.size()), next(next), cycle(n, -1), entry(n), pos(n, -1), tail(n){
        for (int v : next) assert(0 <= v && v < n);

        find_cycles();
        order_tails();
    }

    long long dist(int a, int b) const{
        if (cycle[a] != cycle[b]) return -1;

        if (pos[b] != -1){
            long long len = cycles[cycle[b]].size();
            return tail[a] + (pos[b] - pos[entry[a]] + len) % len;
        }
        if (tail[a] < tail[b]) return -1;

        int k = tail[a] - tail[b];
        return jump(a, k) == b ? k : -1;
    }

    int jump(int v, long long k) const{
        assert(k >= 0);
        if (k < tail[v]){
            /// The ancestor at depth d is the last node of depth d whose preorder index is at most v's
            int d = tail[v] - k;
            auto first = by_depth.begin() + depth_start[d], last = by_depth.begin() + depth_start[d + 1];
            auto it = upper_bound(first, last, tin[v], [&](int t, int u){ return t < tin[u]; });
            return *prev(it);
        }

        const vector<int>& ring = cycles[cycle[v]];
        return ring[(pos[entry[v]] + (k - tail[v]) % ring.size()) % ring.size()];  /// reduced first: k may be LLONG_MAX
    }

    void find_cycles(){
        /// 0 unvisited, 1 on the current walk, 2 finished
        vector<char> state(n, 0);
        vector<int> walk;
        for (int s = 0; s < n; s++){
            walk.clear();
            int v = s;
            while (state[v] == 0) state[v] = 1, walk.push_back(v), v = next[v];

            if (state[v] == 1){
                vector<int> ring(find(walk.begin(), walk.end(), v), walk.end());
                for (int i = 0; i < (int)ring.size(); i++) cycle[ring[i]] = cycles.size(), pos[ring[i]] = i;
                cycles.push_back(ring);
            }
            for (int u : walk) state[u] = 2;
        }
    }

    void order_tails(){
        /// Reverse tree edges u -> next[u] for every tail node u, as flat adjacency lists
        vector<int> start(n + 1, 0), child(n);
        for (int u = 0; u < n; u++) if (pos[u] == -1) start[next[u] + 1]++;
        for (int v = 0; v < n; v++) start[v + 1] += start[v];
        vector<int> slot(start.begin(), start.end() - 1);
        for (int u = 0; u < n; u++) if (pos[u] == -1) child[slot[next[u]]++] = u;

        vector<int> preorder, stack;
        tin.assign(n, 0);
        for (int r = 0; r < n; r++){
            if (pos[r] == -1) continue;
            entry[r] = r, stack.push_back(r);
            while (!stack.empty()){
                int v = stack.back();
                stack.pop_back();
                tin[v] = preorder.size(), preorder.push_back(v);
                for (int i = start[v]; i < start[v + 1]; i++){
                    int u = child[i];
                    cycle[u] = cycle[v], entry[u] = entry[v], tail[u] = tail[v] + 1;
                    stack.push_back(u);
                }
            }
        }

        int max_tail = 0;
        for (int t : tail) max_tail = max(max_tail, t);
        depth_start.assign(max_tail + 2, 0);
        for (int v : preorder) depth_start[tail[v] + 1]++;
        for (int d = 0; d <= max_tail; d++) depth_start[d + 1] += depth_start[d];
        by_depth.resize(n);
        vector<int> fill_at(depth_start.begin(), depth_start.end() - 1);
        for (int v : preorder) by_depth[fill_at[tail[v]]++] = v;
    }
};

int main(){
    /// 0 -> 1 -> 2 -> 3 -> 1 is a tail 0 into the cycle 1, 2, 3, and 4 -> 5 -> 4 is a second cycle
    FunctionalGraph g({1, 2, 3, 1, 5, 4});
    assert(g.cycles.size() == 2);
    assert((g.cycles[g.cycle[0]] == vector<int>{1, 2, 3}));
    assert(g.entry[0] == 1);
    assert(g.tail[0] == 1);
    assert(g.pos[0] == -1);
    assert(g.pos[3] == 2);

    /// From 0: 1 step to 1, then k - 1 more around the 3-cycle, so k = 10 lands on 1 and k = 10^18 + 1 on 2
    assert(g.jump(0, 10) == 1);
    assert(g.jump(0, 1000000000000000001LL) == 2);
    assert(g.dist(0, 3) == 3);
    assert(g.dist(2, 1) == 2);
    assert(g.dist(1, 0) == -1);
    assert(g.dist(0, 4) == -1);

    return 0;
}
