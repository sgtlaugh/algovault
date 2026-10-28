/***
 *
 * 2SAT with Kosaraju's algorithm
 * Finds a true or false value for every variable that satisfies a system of constraints on pairs of variables
 *
 * Complexity: O(n + m), n = number of variables, m = number of clauses
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
 * g.is_satisfiable(): whether the constraints can all hold, finding an assignment when they can
 * g.value(x): whether literal x holds in that assignment, valid after is_satisfiable() returns true
 *
 * Notes:
 *   - at_most_one adds k - 2 variables via add_var() and 3k - 5 clauses for k >= 2, not k^2 / 2 pairwise clauses
 *   - Variables from add_var() are numbered after the existing ones, so earlier indices never change
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct TwoSat{
    int n, t;
    vector<pair<int, int>> edges;  /// signed literals, mapped to nodes 1..2n only in build() since add_var() changes n
    vector<int> adj_start, adj_to, rev_start, rev_to;
    vector<int> visited, comp, order, dfs_t;

    TwoSat(int n): n(n) {}

    inline void add_and(int a, int b){
        force_true(a);
        force_true(b);
    }

    inline void add_implication(int a, int b){
        assert(a != 0 && abs(a) <= n && b != 0 && abs(b) <= n);
        edges.push_back({a, b}), edges.push_back({-b, -a});
    }

    inline void add_or(int a, int b){
        add_implication(-a, b);
    }

    int add_var(){
        return ++n;
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
        assert(x != 0 && abs(x) <= n);
        edges.push_back({x, -x});  /// x -> !x is its own contrapositive, add it once
    }

    inline void force_true(int x){
        assert(x != 0 && abs(x) <= n);
        edges.push_back({-x, x});  /// !x -> x is its own contrapositive, add it once
    }

    vector<int> get_assignment(){
        vector<int> res;
        for (int x = 1; x <= n; x++){
            if (value(x)) res.push_back(x);
        }

        return res;
    }

    bool is_satisfiable(){
        build();
        for (int i = 1; i <= n; i++){
            if (comp[i] == comp[i + n]) return false;
        }

        return true;
    }

    /// If x implies !x then comp[!x] is discovered first, so x is false
    inline bool value(int x){
        if (x < 0) return !value(-x);
        return comp[x] < comp[x + n];
    }

    /// Components are numbered in discovery order, which is a reverse topological order of the implication graph
    void build(){
        int m = 2 * n + 2, c = 0;

        flatten(adj_start, adj_to, false);
        flatten(rev_start, rev_to, true);
        dfs_t.assign(m, 0), order.assign(m, 0), comp.assign(m, 0);

        vector<pair<int, int>> topsort_stack;
        visited.assign(m, 0);
        t = 0;
        for (int i = 2 * n; i >= 1; i--){
            if (!visited[i]) topsort(i, topsort_stack);
            order[dfs_t[i]] = i;
        }

        vector<int> dfs_stack;
        visited.assign(m, 0);
        for (int i = 2 * n; i >= 1; i--){
            int x = order[i];
            if (!visited[x]) dfs(x, c++, dfs_stack);
        }
    }

    void flatten(vector<int>& start, vector<int>& to, bool reverse){
        start.assign(2 * n + 2, 0), to.resize(edges.size());
        for (auto [a, b]: edges) start[node(reverse ? b : a)]++;
        for (int i = 1; i <= 2 * n + 1; i++) start[i] += start[i - 1];
        for (auto [a, b]: edges){
            int u = node(a), v = node(b);
            if (reverse) swap(u, v);
            to[--start[u]] = v;
        }
    }

    /// Literal x is node x, literal -x is node n + x
    inline int node(int x){
        return x > 0 ? x : n - x;
    }

    void topsort(int s, vector<pair<int, int>>& stack){
        visited[s] = true, stack.push_back({s, rev_start[s]});
        while (!stack.empty()){
            auto& [u, e] = stack.back();
            if (e < rev_start[u + 1]){
                int v = rev_to[e++];
                if (!visited[v]) visited[v] = true, stack.push_back({v, rev_start[v]});
            }
            else{
                dfs_t[u] = ++t;
                stack.pop_back();
            }
        }
    }

    /// Labels everything reachable from s, the visiting order does not matter here
    void dfs(int s, int c, vector<int>& stack){
        comp[s] = c, visited[s] = true, stack.push_back(s);
        while (!stack.empty()){
            int u = stack.back();
            stack.pop_back();
            for (int e = adj_start[u]; e < adj_start[u + 1]; e++){
                int v = adj_to[e];
                if (!visited[v]) comp[v] = c, visited[v] = true, stack.push_back(v);
            }
        }
    }
};

int main(){
    /// 2 is forced false, so the xor makes 4 true, the implication makes 1 false and its contrapositive 3 false
    TwoSat g(4);
    g.add_implication(1, 2);
    g.add_implication(-2, -3);
    g.force_false(2);
    g.add_xor(2, 4);
    g.add_or(1, 4);
    assert(g.is_satisfiable());
    assert((g.get_assignment() == vector<int>{4}));
    assert(g.value(4));
    assert(!g.value(-4));

    /// One of 1 and 3 holds, so -2 and 4 cannot, and the 4 literals took variables 5 and 6
    TwoSat h(4);
    h.at_most_one({1, -2, 3, 4});
    assert(h.n == 6);
    h.add_or(1, 3);
    assert(h.is_satisfiable());
    assert(h.value(2));
    assert(!h.value(4));
    assert(h.value(1) != h.value(3));

    h.force_true(4);
    assert(!h.is_satisfiable());

    return 0;
}
