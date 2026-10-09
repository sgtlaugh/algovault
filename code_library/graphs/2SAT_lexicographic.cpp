/***
 *
 * Lexicographically first 2SAT assignment (1 based index for variables)
 * Each variable can have two possible values: true or false
 * Variables must satisfy a system of constraints on pairs of variables
 * Positive variables are identified from 1 to n, while negative variables from -1 to -n
 *
 * Among all satisfying assignments, returns the one that sets x1 true if possible, then x2, and so on
 * Use 2SAT_kosaraju.cpp when any satisfying assignment will do, it is always linear
 *
 * Fixes variables greedily in index order, a literal can be set iff it does not imply its own negation
 * Easy inputs run in O(n + m), unsatisfiable inputs in O(n + m) once the greedy wastes O(n + m) work,
 * and the worst case is O(n * m / 64) with bitset reachability, instead of O(n * m) for the plain greedy
 *
***/

#include <bits/stdc++.h>

using namespace std;

struct Graph{
    int n;
    vector<vector<int>> adj;  /// node 2x is x true, 2x + 1 is x false
    vector<char> value;

    Graph(int n) : n(n), adj(2 * n + 2) {}

    inline int node(int x){
        assert(x != 0 && abs(x) <= n);
        return x > 0 ? 2 * x : -2 * x + 1;
    }

    /// Add implication, if a then b
    inline void add_implication(int a, int b){
        int u = node(a), v = node(b);
        adj[u].push_back(v), adj[v ^ 1].push_back(u ^ 1);
    }

    inline void add_or(int a, int b){
        add_implication(-a, b);
    }

    inline void add_xor(int a, int b){
        add_or(a, b);
        add_or(-a, -b);
    }

    inline void add_and(int a, int b){
        add_or(a, b);
        add_or(a, -b);
        add_or(-a, b);
    }

    /// Force variable x to be true (if x is negative, force !x to be true)
    inline void force_true(int x){
        adj[node(x) ^ 1].push_back(node(x));  /// !x -> x is its own contrapositive, add it once
    }

    /// Force variable x to be false (if x is negative, force !x to be false)
    inline void force_false(int x){
        force_true(-x);
    }

    /// Iterative Tarjan, components are numbered in reverse topological order (sinks first)
    int scc(vector<int>& comp){
        int N = adj.size(), timer = 0, C = 0;
        vector<int> low(N), num(N, -1), it(N, 0), st, call;
        vector<char> on(N, 0);
        comp.assign(N, -1);

        for (int s = 2; s < N; s++){
            if (num[s] != -1) continue;
            call.assign(1, s), num[s] = low[s] = timer++, st.push_back(s), on[s] = 1;
            while (!call.empty()){
                int u = call.back();
                if (it[u] < (int)adj[u].size()){
                    int v = adj[u][it[u]++];
                    if (num[v] == -1) num[v] = low[v] = timer++, st.push_back(v), on[v] = 1, call.push_back(v);
                    else if (on[v]) low[u] = min(low[u], num[v]);
                    continue;
                }
                call.pop_back();
                if (!call.empty()) low[call.back()] = min(low[call.back()], low[u]);
                if (low[u] == num[u]){
                    int w;
                    do { w = st.back(), st.pop_back(), on[w] = 0, comp[w] = C; } while (w != u);
                    C++;
                }
            }
        }

        return C;
    }

    /// Returns true if the system is 2-satisfiable, the assignment is then available with get_assignment()
    bool is_satisfiable(){
        int N = adj.size(), C = 0;
        long long E = 0;
        for (auto& v : adj) E += v.size();

        /// Wasted (failed or aborted) work allowed before SCC pays for itself, then a per query budget that keeps
        /// a bitset batch of 64 queries, which costs O(n + m), amortized to the same O((n + m) / 64) per query
        long long allowance = N + E, budget = max(64LL, (N + E) / 64);

        vector<int> comp, stk, trail;
        vector<vector<int>> dag;
        vector<unsigned long long> mask;
        vector<signed char> reach(n + 1, -1);  /// reach[x]: whether node 2x reaches 2x + 1, -1 while unknown
        value.assign(N, 0);

        auto mark_closure = [&](int l){
            stk.assign(1, l);
            while (!stk.empty()){
                int u = stk.back();
                stk.pop_back();
                if (value[u]) continue;
                value[u] = 1;
                for (int w : adj[u]) if (!value[w]) stk.push_back(w);
            }
        };

        /// 1 if l was consistent and its closure is kept, 0 on a contradiction, -1 when over the limit
        auto attempt = [&](int l, long long limit){
            long long work = 0;
            int status = 1;
            trail.clear(), stk.assign(1, l);
            while (!stk.empty()){
                int u = stk.back();
                stk.pop_back();
                if (value[u]) continue;
                if (value[u ^ 1]){
                    status = 0;
                    break;
                }
                if (++work > limit){
                    status = -1;
                    break;
                }
                value[u] = 1, trail.push_back(u);
                for (int w : adj[u]) if (!value[w]) stk.push_back(w);  /// fixed nodes are closed under implication, stop there
            }
            if (status != 1){
                for (int u : trail) value[u] = 0;
                allowance -= work;
            }
            return status;
        };

        /// One bitset pass over the condensation answers up to 64 pending queries at once
        auto answer_batch = [&](int from){
            if (dag.empty()){
                dag.resize(C), mask.resize(C);
                for (int u = 2; u < N; u++){
                    for (int v : adj[u]) if (comp[u] != comp[v]) dag[comp[u]].push_back(comp[v]);
                }
            }

            vector<int> qs;
            for (int x = from; x <= n && qs.size() < 64; x++){
                if (reach[x] == -1 && comp[2 * x] > comp[2 * x + 1]) qs.push_back(x);
            }
            fill(mask.begin(), mask.end(), 0);
            for (int j = 0; j < (int)qs.size(); j++) mask[comp[2 * qs[j]]] |= 1ULL << j;
            for (int c = C - 1; c >= 0; c--){
                if (mask[c]) for (int d : dag[c]) mask[d] |= mask[c];
            }
            for (int j = 0; j < (int)qs.size(); j++) reach[qs[j]] = mask[comp[2 * qs[j] + 1]] >> j & 1;
        };

        for (int x = 1; x <= n; x++){
            int t = 2 * x;
            if (value[t] || value[t ^ 1]) continue;

            /// A literal whose component comes after its negation's cannot reach the negation
            if (!comp.empty() && comp[t] < comp[t ^ 1]){
                mark_closure(t);
                continue;
            }

            if (reach[x] == -1){
                int status = attempt(t, comp.empty() ? max(allowance, 0LL) : budget);
                if (status == 1) continue;

                if (status == 0 && comp.empty()){  /// t implies !t, if !t also implies t the system is unsatisfiable
                    int other = attempt(t ^ 1, max(allowance, 0LL));
                    if (other == 1) continue;
                    if (other == 0) return false;
                }

                if (comp.empty()){  /// the waste has paid for SCC, redo x with the shortcut and batches available
                    C = scc(comp);
                    for (int y = 1; y <= n; y++) if (comp[2 * y] == comp[2 * y + 1]) return false;
                    x--;
                    continue;
                }

                if (status == 0) reach[x] = 1;
                else answer_batch(x);
            }
            mark_closure(reach[x] ? t ^ 1 : t);
        }

        return true;
    }

    /***
     * returns variables set to true in a vector, in increasing order
     * only call when is_satisfiable()
    ***/

    vector <int> get_assignment(){
        vector <int> set_nodes;
        for (int x = 1; x <= n; x++){
            if (value[2 * x]) set_nodes.push_back(x);
        }

        return set_nodes;
    }
};

int main(){
    auto g = Graph(5);
    assert(g.is_satisfiable());
    assert((g.get_assignment() == vector<int>{1, 2, 3, 4, 5}));

    g.add_implication(-1, -2);  /// if 1 is false then 2 is false
    g.add_implication(-2, -3);  /// if 2 is false then 3 is false
    g.add_implication(-3, 4);   /// if 3 is false then 4 is true
    g.add_implication(-3, -5);  /// if 3 is false then 5 is false
    g.add_xor(1, 4);            /// exactly one of 1 or 4 must be true
    g.add_or(4, 5);             /// either of 4 or 5 must be true

    assert(g.is_satisfiable());
    assert((g.get_assignment() == vector<int>{1, 2, 3, 5}));

    g.force_false(1);           /// 1 false forces 2 and 3 false, so 4 true and 5 false
    assert(g.is_satisfiable());
    assert((g.get_assignment() == vector<int>{4}));

    g.add_xor(1, 5);
    assert(!g.is_satisfiable());
    return 0;
}
