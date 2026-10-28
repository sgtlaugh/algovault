#include "../common.h"

#define main library_main
#include "../../code_library/graphs/2SAT_lexicographic.cpp"
#undef main

typedef pair<int, int> Clause;  /// a or b, a == b forces a single literal

void add_all(TwoSat& g, const vector<Clause>& clauses){
    for (auto [a, b] : clauses){
        if (a == b) g.force_true(a);
        else g.add_or(a, b);
    }
}

/// The plain greedy, an independent O(n * m) reference that is fast enough at a few thousand variables
pair<bool, vector<int>> greedy(int n, const vector<Clause>& clauses){
    vector<vector<int>> adj(2 * n + 2);
    auto node = [](int x){ return x > 0 ? 2 * x : -2 * x + 1; };
    for (auto [a, b] : clauses) adj[node(-a)].push_back(node(b)), adj[node(-b)].push_back(node(a));

    vector<char> on(2 * n + 2, 0);
    vector<int> trail, stk;
    auto attempt = [&](int s){
        trail.clear(), stk.assign(1, s);
        while (!stk.empty()){
            int u = stk.back();
            stk.pop_back();
            if (on[u ^ 1]) return false;
            if (on[u]) continue;
            on[u] = 1, trail.push_back(u);
            for (int v : adj[u]) stk.push_back(v);
        }

        return true;
    };

    for (int x = 1; x <= n; x++){
        if (on[2 * x] || on[2 * x + 1]) continue;
        if (!attempt(2 * x)){
            for (int u : trail) on[u] = 0;
            if (!attempt(2 * x + 1)) return {false, {}};
        }
    }

    vector<int> res;
    for (int x = 1; x <= n; x++) if (on[2 * x]) res.push_back(x);
    return {true, res};
}

int random_literal(int n){
    int x = stress::rand_int(1, n);
    return stress::rand_int(0, 1) ? x : -x;
}

/// Exhaustive over the caller's variables, at_most_one helpers excluded: the first assignment preferring true
void check_brute(TwoSat& g, const vector<int>& vars, const vector<Clause>& clauses, const vector<vector<int>>& groups){
    int k = vars.size(), best = -1;
    vector<int> bit(g.n + 1, -1);
    for (int i = 0; i < k; i++) bit[vars[i]] = i;
    auto holds = [&](int mask, int x){ return (mask >> bit[abs(x)] & 1) == (x > 0); };
    auto key = [&](int mask){ int r = 0; for (int i = 0; i < k; i++) r = r * 2 + (mask >> i & 1); return r; };
    for (int mask = 0; mask < (1 << k); mask++){
        bool ok = true;
        for (auto [a, b] : clauses) ok &= holds(mask, a) || holds(mask, b);
        for (auto& lits : groups){
            int count = 0;
            for (int x : lits) count += holds(mask, x);
            ok &= count <= 1;
        }
        if (ok && (best == -1 || key(mask) > key(best))) best = mask;
    }

    assert(g.is_satisfiable() == (best != -1));
    if (best == -1) return;
    vector<int> expected, got;
    for (int i = 0; i < k; i++) if (best >> i & 1) expected.push_back(vars[i]);
    for (int x : g.get_assignment()) if (bit[x] != -1) got.push_back(x);
    assert(got == expected);
    for (int x : vars) assert(g.value(x) == holds(best, x));
    for (int x = 1; x <= g.n; x++) assert(g.value(-x) != g.value(x));
}

void check_greedy(int n, const vector<Clause>& clauses){
    TwoSat g(n);
    add_all(g, clauses);
    auto [sat, expected] = greedy(n, clauses);
    assert(g.is_satisfiable() == sat);
    if (sat) assert(g.get_assignment() == expected);
}

/// x_i = true walks the whole chain p, so attempts blow the budget and SCC, the shortcut and bitset batches all run
vector<Clause> adversarial(int n, int k, bool contradiction){
    vector<Clause> clauses;
    /// SCC starts once failed attempts (~k + n each) waste ~2 (vars + clauses), a reaching prefix also hits true x
    long long size = 2LL * (n + k) + 2LL * (k + 2 * n);
    int reaching = max(1LL, min<long long>(n / 2, size / (k + n) + stress::rand_int(-3, 2)));
    if (stress::rand_int(0, 2) == 0) reaching = stress::rand_int(1, n / 2);

    for (int j = 1; j < k; j++) clauses.push_back({-(n + j), n + j + 1});
    for (int i = 1; i <= n; i++){
        /// past the reaching prefix x_i may skip x_i -> p_1, and then must not reach !x_i
        if (i <= reaching || stress::rand_int(0, 1)) clauses.push_back({-i, n + 1});
        /// p_k -> !x_i keeps them off the topological shortcut
        clauses.push_back({-(n + k), -i});
        /// a few forced false x fail within the budget after SCC
        if (i > reaching && stress::rand_int(0, 9) == 0) clauses.push_back({-i, -i});
    }
    /// noise only adds !x_a or !x_b, which never forces an x true
    for (int extra = stress::rand_int(0, n / 4); extra; extra--){
        clauses.push_back({-stress::rand_int(1, n), -stress::rand_int(1, n)});
    }
    if (contradiction) clauses.push_back({n, n});  /// the last x, reached only after the waste has paid for SCC
    return clauses;
}

/// Segments of shuffled signed literals, false in a hidden assignment except the hub paths p
vector<Clause> chains(int n){
    int tail = stress::rand_int(0, n / 4), vars = tail + n;
    vector<int> lit(n);
    for (int i = 0; i < n; i++) lit[i] = stress::rand_int(0, 1) ? tail + i + 1 : -(tail + i + 1);
    shuffle(lit.begin(), lit.end(), stress::rng());

    vector<Clause> clauses;
    vector<char> hidden(vars + n + 1, 0), plain(vars + 1, 1);  /// hidden[x]: x is true, room for one decoy per segment
    /// 1..tail is a true chain fixed first that chain literals point into, so batches run into decided components
    for (int x = 1; x <= tail; x++){
        hidden[x] = 1;
        if (x < tail) clauses.push_back({-x, x + 1});
    }

    for (int start = 0, decoy = vars; start < n;){
        int end = min(n, start + (int)stress::rand_int(1, 1200));
        if (stress::rand_int(0, 1)){  /// positive and ascending, so each x fails at full cost and fixes little
            for (int i = start; i < end; i++) lit[i] = abs(lit[i]), plain[lit[i]] = 0;
            sort(lit.begin() + start, lit.begin() + end);
        }
        if (stress::rand_int(0, 2) || end - start < 2){  /// a chain l_1 -> ... -> l_k -> !l_k
            for (int i = start; i < end; i++){
                hidden[abs(lit[i])] = lit[i] < 0;
                if (i + 1 < end) clauses.push_back({-lit[i], lit[i + 1]});
                if (tail && stress::rand_int(0, 7) == 0) clauses.push_back({-lit[i], (int)stress::rand_int(1, tail)});
            }
            clauses.push_back({-lit[end - 1], -lit[end - 1]});
            if (stress::rand_int(0, 1)){  /// a forced false decoy d -> !l_k, so Tarjan meets the negated chain first
                decoy++;
                clauses.push_back({-decoy, -decoy}), clauses.push_back({-decoy, -lit[end - 1]});
            }
        }
        else {  /// a hub: s -> p_1 -> ... -> p_k -> !s, no s an ancestor of another, so each s is its own query
            int mid = stress::rand_int(start + 1, end - 1);
            for (int i = mid; i < end; i++){
                hidden[abs(lit[i])] = lit[i] > 0;
                if (i + 1 < end) clauses.push_back({-lit[i], lit[i + 1]});
            }
            for (int i = start; i < mid; i++){
                hidden[abs(lit[i])] = lit[i] < 0;
                clauses.push_back({-lit[i], lit[mid]}), clauses.push_back({-lit[end - 1], -lit[i]});
            }
        }
        start = end;
    }

    /// Noise only between plain variables: noise into a negated chain would fix its whole prefix
    vector<int> plain_vars;
    for (int x = 1; x <= vars; x++) if (plain[x]) plain_vars.push_back(x);
    /// drawing from plain_vars instead of rejecting other picks ends the loop when no variable is plain
    for (int extra = plain_vars.empty() ? 0 : stress::rand_int(0, n / 4); extra;){
        int a = plain_vars[stress::rand_int(0, plain_vars.size() - 1)];
        int b = plain_vars[stress::rand_int(0, plain_vars.size() - 1)];
        if (stress::rand_int(0, 1)) a = -a;
        if (stress::rand_int(0, 1)) b = -b;
        if ((a > 0) == hidden[abs(a)] || (b > 0) == hidden[abs(b)]) clauses.push_back({a, b}), extra--;
    }

    return clauses;
}

/// The batch must answer a source g_i ahead of the decoys y_j: taking y_j until g_i comes up costs k / B passes over R
vector<Clause> decoys(int k, int len){
    int lq = stress::rand_int(len / 12, len / 6), a = (6 * k + 4 * len) / lq + stress::rand_int(4, 12);
    int g0 = 1, q0 = g0 + a, y0 = q0 + lq, r0 = y0 + k;

    vector<Clause> clauses = {{-1, -(r0 + len)}};
    /// variable 1 fails first and its Tarjan DFS numbers the y_j and R low, so the y_j lead the sinks first order
    for (int j = 1; j <= k; j++) clauses.push_back({-1, y0 + j}), clauses.push_back({-(y0 + j), r0 + 1});
    /// sources g_i over a shared path Q burn the waste allowance, so SCC starts on some g_i
    for (int i = 1; i <= a; i++) clauses.push_back({-(g0 + i), q0 + 1}), clauses.push_back({-(q0 + lq), -(g0 + i)});
    for (int i = 1; i < lq; i++) clauses.push_back({-(q0 + i), q0 + i + 1});
    for (int i = 1; i < len; i++) clauses.push_back({-(r0 + i), r0 + i + 1});
    return clauses;
}

int main(){
    for (long long it = 0; it < stress::scaled(4000); it++){
        int n = stress::rand_int(1, 10), m = stress::rand_int(0, 3 * n);
        TwoSat g(n);
        vector<int> vars(n);
        iota(vars.begin(), vars.end(), 1);
        vector<Clause> clauses;
        vector<vector<int>> groups;
        auto literal = [&](){
            int x = vars[stress::rand_int(0, vars.size() - 1)];
            return stress::rand_int(0, 1) ? x : -x;
        };
        for (int c = 0; c < m; c++){
            int a = literal(), b = literal(), kind = stress::rand_int(0, 7);
            if (kind == 0) g.add_or(a, b), clauses.push_back({a, b});
            else if (kind == 1) g.add_implication(a, b), clauses.push_back({-a, b});
            else if (kind == 2) g.add_xor(a, b), clauses.push_back({a, b}), clauses.push_back({-a, -b});
            else if (kind == 3) g.add_and(a, b), clauses.push_back({a, a}), clauses.push_back({b, b});
            else if (kind == 4) g.force_true(a), clauses.push_back({a, a});
            else if (kind == 5) g.force_false(a), clauses.push_back({-a, -a});
            else if (kind == 6){
                vector<int> lits(stress::rand_int(0, 5));
                for (auto& x : lits) x = literal();
                g.at_most_one(lits), groups.push_back(lits);
            }
            /// a variable added after at_most_one helpers must still lose only to the variables before it
            else if (vars.size() < 10) vars.push_back(g.add_var());
        }
        check_brute(g, vars, clauses, groups);
    }

    for (long long it = 0; it < stress::scaled(60); it++){
        int mode = it % 6, n = stress::rand_int(200, 1500);
        vector<Clause> clauses;
        if (mode < 2) clauses = adversarial(n, stress::rand_int(n / 2, 2 * n), mode == 1);
        else if (mode >= 4) clauses = chains(2 * n);
        else {
            int vars = 4 * n, m = mode == 2 ? 9 * vars / 10 : 3 * vars;  /// under the 2SAT threshold, and dense
            vector<int> hidden(vars + 1);
            for (auto& h : hidden) h = stress::rand_int(0, 1);
            while ((int)clauses.size() < m){
                int a = random_literal(vars), b = random_literal(vars);
                if (mode == 2 || (a > 0) == hidden[abs(a)] || (b > 0) == hidden[abs(b)]) clauses.push_back({a, b});
            }
        }

        int vars = 0;
        for (auto [a, b] : clauses) vars = max({vars, abs(a), abs(b)});
        check_greedy(vars, clauses);
    }

    for (int it = 0; it < 2; it++){
        auto clauses = decoys(stress::rand_int(50000, 60000), stress::rand_int(50000, 60000));
        int vars = 0;
        for (auto [a, b] : clauses) vars = max({vars, abs(a), abs(b)});
        check_greedy(vars, clauses);
    }

    /// x_i -> p_1 -> ... -> p_k -> !x_i: every x_i fails at full cost, so only the waste allowance keeps this linear
    int n = 30000;
    TwoSat chain(2 * n);
    for (int j = 1; j < n; j++) chain.add_implication(n + j, n + j + 1);
    for (int i = 1; i <= n; i++) chain.add_implication(i, n + 1), chain.add_implication(2 * n, -i);
    clock_t start = clock();  /// CPU time, parallel sanitizer jobs inflate wall time
    assert(chain.is_satisfiable());
    assert(clock() - start < 2 * CLOCKS_PER_SEC);
    vector<int> path(n);
    iota(path.begin(), path.end(), n + 1);
    assert(chain.get_assignment() == path);

    /// Solving again after more constraints must start from scratch
    TwoSat g(3);
    g.add_or(-1, -2);
    assert(g.is_satisfiable() && (g.get_assignment() == vector<int>{1, 3}));
    g.force_false(1);
    assert(g.is_satisfiable() && (g.get_assignment() == vector<int>{2, 3}));
    g.force_true(1);
    assert(!g.is_satisfiable());
    return 0;
}
