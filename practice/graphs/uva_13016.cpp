/***
 *
 * 2SAT with Kosaraju's algorithm (1 based index for variables)
 * Each variable can have two possible values: true or false
 * Variables must satisfy a system of constraints on pairs of variables
 * After is_satisfiable() returns true, value(x) gives a satisfying assignment
 *
 * Complexity: O(n + m), n = number of variables, m = number of clauses
 *
 * Literals are signed variable indices: x means x is true, -x means x is false
 * at_most_one(literals) allows at most one of k literals to be true (a literal listed twice counts twice)
 * It uses the prefix encoding: for k >= 2, k - 2 auxiliary variables and 3k - 5 clauses instead of k^2 / 2 pairwise clauses
 * Auxiliary variables come from add_var(), numbered after the existing ones, so earlier indices never change
 *
 * Both DFS passes are iterative, so long implication paths cannot overflow the stack
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

    /// Literal x is node x, literal -x is node n + x
    inline int node(int x){
        return x > 0 ? x : n - x;
    }

    /// Adds a fresh variable and returns its index
    int add_var(){
        return ++n;
    }

    /// Add implication, if a then b
    inline void add_implication(int a, int b){
        assert(a != 0 && abs(a) <= n && b != 0 && abs(b) <= n);
        edges.push_back({a, b}), edges.push_back({-b, -a});
    }

    inline void add_or(int a, int b){
        add_implication(-a, b);
    }

    inline void add_xor(int a, int b){
        add_or(a, b);
        add_or(-a, -b);
    }

    inline void add_and(int a, int b){
        force_true(a);
        force_true(b);
    }

    /// Force variable x to be true (if x is negative, force !x to be true)
    inline void force_true(int x){
        assert(x != 0 && abs(x) <= n);
        edges.push_back({-x, x});  /// !x -> x is its own contrapositive, add it once
    }

    /// Force variable x to be false (if x is negative, force !x to be false)
    inline void force_false(int x){
        assert(x != 0 && abs(x) <= n);
        edges.push_back({x, -x});  /// x -> !x is its own contrapositive, add it once
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

    /// Stable counting sort of the edges by source, so the DFS visits neighbours in insertion order
    void flatten(vector<int>& start, vector<int>& to, bool reverse){
        start.assign(2 * n + 2, 0), to.resize(edges.size());
        for (auto [a, b]: edges) start[node(reverse ? b : a)]++;
        for (int i = 1; i <= 2 * n + 1; i++) start[i] += start[i - 1];
        for (int i = (int)edges.size() - 1; i >= 0; i--){
            int u = node(edges[i].first), v = node(edges[i].second);
            if (reverse) swap(u, v);
            to[--start[u]] = v;
        }
    }

    /// Post-order DFS on the reverse graph, each stack entry keeps the next edge to scan so neighbours go in insertion order
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

    /// Components are numbered in discovery order, which is a reverse topological order of the implication graph
    void build(){
        int i, x, c = 0, m = 2 * n + 2;
        flatten(adj_start, adj_to, false), flatten(rev_start, rev_to, true);
        dfs_t.assign(m, 0), order.assign(m, 0), comp.assign(m, 0);

        vector<pair<int, int>> topsort_stack;
        visited.assign(m, 0);
        for (i = 2 * n, t = 0; i >= 1; i--){
            if (!visited[i]) topsort(i, topsort_stack);
            order[dfs_t[i]] = i;
        }

        vector<int> dfs_stack;
        visited.assign(m, 0);
        for (i = 2 * n; i >= 1; i--){
            x = order[i];
            if (!visited[x]) dfs(x, c++, dfs_stack);
        }
    }

    bool is_satisfiable(){
        build();
        for (int i = 1; i <= n; i++){
            if (comp[i] == comp[i + n]) return false;
        }
        return true;
    }

    /// Value of literal x in a satisfying assignment, valid only after is_satisfiable() returned true
    /// If x implies !x then comp[!x] is discovered first, so x is false
    inline bool value(int x){
        if (x < 0) return !value(-x);
        return comp[x] < comp[x + n];
    }
};

struct Reply{
    int a, b, player, count;
};

/// Variable x is true when suspect x is in the hand of player 2, the three circle suspects are in neither hand
bool admissible(const vector<Reply>& replies, const array<bool, 27>& circle){
    TwoSat sat(26);
    for (auto [a, b, player, count]: replies){
        vector<int> held;
        for (int x: {a, b}){
            if (!circle[x]) held.push_back(player == 2 ? x : -x);
        }

        if (count > (int)held.size()) return false;
        if (held.size() == 2 && count == 1) sat.add_xor(held[0], held[1]);
        else{
            for (int lit: held){
                if (count) sat.force_true(lit);
                else sat.force_false(lit);
            }
        }
    }
    return sat.is_satisfiable();
}

int main(){
    int n;
    while (scanf("%d", &n) == 1){
        vector<Reply> replies(n);
        for (auto& r: replies){
            char s[8];
            scanf("%s %d %d", s, &r.player, &r.count);
            r.a = s[0] - 'A' + 1, r.b = s[1] - 'A' + 1;
        }

        int res = 0;
        for (int i = 1; i <= 26; i++){
            for (int j = i + 1; j <= 26; j++){
                for (int k = j + 1; k <= 26; k++){
                    array<bool, 27> circle{};
                    circle[i] = circle[j] = circle[k] = true;
                    res += admissible(replies, circle);
                }
            }
        }
        printf("%d\n", res);
    }
    return 0;
}
