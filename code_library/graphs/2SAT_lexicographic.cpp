/***
 *
 * Lexicographically first 2SAT assignment
 * Among all satisfying assignments, finds the one that sets x1 true if possible, then x2, and so on
 *
 * Complexity:
 *   - O(n + m) when no greedy attempt fails
 *   - O(n (n + m)) worst case, n variables implying one literal of out-degree m cost O(n m), 7 s at n = m = 1.6e5
 *   - O(n + m) memory, plus 32 bytes per component once a batch runs
 *
 * TwoSat g(n): no constraints on variables 1 to n, literal x means x is true and -x means x is false
 * g.add_and(a, b): requires both literals a and b
 * g.add_implication(a, b): requires b whenever a holds
 * g.add_or(a, b): requires at least one of a and b
 * g.add_var(): adds variable n + 1 and returns its index
 * g.add_xor(a, b): requires exactly one of a and b
 * g.at_most_one(literals): allows at most one of the literals to hold, a literal listed twice counts twice
 * g.force_false(x): requires literal x to be false
 * g.force_true(x): requires literal x to be true
 * g.get_assignment(): the true variables, ascending, valid after is_satisfiable() returns true
 * g.is_satisfiable(): whether the constraints can all hold, finding the first assignment when they can
 * g.value(x): whether literal x holds in that assignment, valid after is_satisfiable() returns true
 *
 * Notes:
 *   - Use 2SAT_kosaraju.cpp when any satisfying assignment will do, it is always linear
 *   - at_most_one adds k - 2 variables through add_var() and 3k - 5 clauses for k >= 2 literals
 *   - The other variables still get the first assignment, as if at_most_one added no variables
 *   - n literals that each imply their negation along one chain take one batch, 0.2 s at n = 5e5
 *   - n / 2 such literals with no ancestor relation, sharing one path of length n / 2, take 9 to 14 s at n = 5e5
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct TwoSat{
    /// Queries per bitset batch, the worst case runs 2.2x faster than with 64, for B / 8 bytes per component
    static constexpr int B = 256;

    int n;
    vector<vector<int>> adj;  /// node 2x is x true, 2x + 1 is x false
    vector<char> truth;

    TwoSat(int n) : n(n), adj(2 * n + 2) {}

    inline void add_and(int a, int b){
        force_true(a);
        force_true(b);
    }

    inline void add_implication(int a, int b){
        int u = node(a), v = node(b);
        adj[u].push_back(v), adj[v ^ 1].push_back(u ^ 1);
    }

    inline void add_or(int a, int b){
        add_implication(-a, b);
    }

    int add_var(){
        n++;
        adj.resize(2 * n + 2);
        return n;
    }

    inline void add_xor(int a, int b){
        add_or(a, b);
        add_or(-a, -b);
    }

    /// prev is a literal meaning "one of lits[0..i-1] is true", it forbids lits[i] and carries over to the next prefix
    void at_most_one(const vector<int>& lits){
        int k = lits.size();
        if (k <= 1) return;

        int prev = lits[0];
        for (int i = 1; i + 1 < k; i++){
            int cur = add_var();
            add_implication(prev, -lits[i]);
            add_implication(prev, cur);
            add_implication(lits[i], cur);
            prev = cur;
        }

        add_implication(prev, -lits[k - 1]);
    }

    inline void force_false(int x){
        force_true(-x);
    }

    inline void force_true(int x){
        adj[node(x) ^ 1].push_back(node(x));  /// !x -> x is its own contrapositive, add it once
    }

    vector<int> get_assignment(){
        vector<int> set_nodes;
        for (int x = 1; x <= n; x++){
            if (truth[2 * x]) set_nodes.push_back(x);
        }

        return set_nodes;
    }

    bool is_satisfiable(){
        int nodes = adj.size(), components = 0;
        long long edges = 0;
        for (auto& v : adj) edges += v.size();

        /// Failed work allowed before SCC pays for itself, then a per query budget of one O(n + m) batch over B queries
        long long allowance = nodes + edges, budget = max((long long)B, (nodes + edges) / B);

        vector<int> comp, stk, trail, rep, order;
        size_t next = 0;
        int ahead = 1;
        vector<vector<int>> dag;
        vector<bitset<B>> mask;
        vector<signed char> reach(n + 1, -1);  /// reach[x]: whether node 2x reaches 2x + 1, -1 while unknown
        truth.assign(nodes, 0);

        auto mark_closure = [&](int l){
            stk.assign(1, l);
            while (!stk.empty()){
                int u = stk.back();
                stk.pop_back();
                if (truth[u]) continue;
                truth[u] = 1;
                for (int w : adj[u]) if (!truth[w]) stk.push_back(w);
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
                if (truth[u]) continue;
                if (truth[u ^ 1]){
                    status = 0;
                    break;
                }
                if (++work > limit){
                    status = -1;
                    break;
                }
                truth[u] = 1, trail.push_back(u);
                /// fixed nodes are closed under implication, stop there
                for (int w : adj[u]) if (!truth[w]) stk.push_back(w);
            }
            if (status != 1){
                for (int u : trail) truth[u] = 0;
                allowance -= work;
            }
            return status;
        };

        /// One bitset pass answers the first B / 2 undecided queries, the current x always among them, then B sinks
        auto answer_batch = [&](){
            if (dag.empty()){
                dag.resize(components), mask.resize(components), rep.resize(components);
                for (int u = 2; u < nodes; u++){
                    rep[comp[u]] = u;
                    for (int v : adj[u]) if (comp[u] != comp[v]) dag[comp[u]].push_back(comp[v]);
                }
                for (int y = 1; y <= n; y++) if (comp[2 * y] > comp[2 * y + 1]) order.push_back(y);
                sort(order.begin(), order.end(), [&](int a, int b){ return comp[2 * a] < comp[2 * b]; });
            }

            /// A skipped variable stays skipped, so both scans resume where the last batch stopped
            vector<int> qs;
            int lo = components, hi = -1;
            auto take = [&](int y){
                if (reach[y] != -1 || truth[2 * y] || truth[2 * y + 1] || comp[2 * y] < comp[2 * y + 1]) return;
                reach[y] = 0;  /// queued, the pass below overwrites it
                qs.push_back(y), lo = min(lo, comp[2 * y + 1]), hi = max(hi, comp[2 * y]);
            };
            for (; ahead <= n && qs.size() < B / 2; ahead++) take(ahead);
            for (; next < order.size() && qs.size() < B; next++) take(order[next]);

            fill(mask.begin() + lo, mask.begin() + hi + 1, bitset<B>());
            for (int j = 0; j < (int)qs.size(); j++) mask[comp[2 * qs[j]]][j] = 1;
            /// A path from an undecided literal to its negation never enters a decided component, so those are skipped
            for (int c = hi; c >= lo; c--){
                if (mask[c].none() || truth[rep[c]] || truth[rep[c] ^ 1]) continue;
                for (int d : dag[c]) if (d >= lo) mask[d] |= mask[c];
            }

            /// y implying !y makes !y hold in every solution, so it is fixed right away and decides the ancestors of y
            for (int j = 0; j < (int)qs.size(); j++){
                reach[qs[j]] = mask[comp[2 * qs[j] + 1]][j];
                if (reach[qs[j]]) mark_closure(2 * qs[j] + 1);
            }
        };

        for (int x = 1; x <= n; x++){
            int t = 2 * x;
            if (truth[t] || truth[t ^ 1]) continue;

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
                    components = scc(comp);
                    for (int y = 1; y <= n; y++) if (comp[2 * y] == comp[2 * y + 1]) return false;
                    x--;
                    continue;
                }

                if (status == 0) reach[x] = 1;
                else answer_batch();
            }
            mark_closure(reach[x] ? t ^ 1 : t);
        }

        return true;
    }

    bool value(int x){
        return truth[node(x)];
    }

    inline int node(int x){
        assert(x != 0 && abs(x) <= n);
        return x > 0 ? 2 * x : -2 * x + 1;
    }

    /// Iterative Tarjan, components are numbered in reverse topological order (sinks first)
    int scc(vector<int>& comp){
        int nodes = adj.size(), timer = 0, components = 0;
        vector<int> low(nodes), num(nodes, -1), it(nodes, 0), st, call;
        vector<char> on(nodes, 0);
        comp.assign(nodes, -1);

        for (int s = 2; s < nodes; s++){
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
                    do { w = st.back(), st.pop_back(), on[w] = 0, comp[w] = components; } while (w != u);
                    components++;
                }
            }
        }

        return components;
    }
};

int main(){
    /// x1 true makes x4 false through the xor, and x4 false forces x3, x5 and then x2
    TwoSat g(5);
    g.add_implication(-1, -2);
    g.add_implication(-2, -3);
    g.add_implication(-3, 4);
    g.add_implication(-3, -5);
    g.add_xor(1, 4);
    g.add_or(4, 5);
    assert(g.is_satisfiable());
    assert((g.get_assignment() == vector<int>{1, 2, 3, 5}));

    /// x1 false forces x2 and x3 false, so x4 true and x5 false
    g.force_false(1);
    assert(g.is_satisfiable());
    assert((g.get_assignment() == vector<int>{4}));
    assert(!g.value(1));
    assert(g.value(-1));
    assert(!g.value(5));

    g.add_xor(1, 5);
    assert(!g.is_satisfiable());

    /// at_most_one takes variable 4 for itself, and the new variable 5 needs x2, which loses to x1
    TwoSat h(3);
    h.at_most_one({1, 2, 3});
    int y = h.add_var();
    assert(y == 5);
    h.add_implication(y, 2);
    assert(h.is_satisfiable());
    assert((h.get_assignment() == vector<int>{1, 4}));

    return 0;
}
