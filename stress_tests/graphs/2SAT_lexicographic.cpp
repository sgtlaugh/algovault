#include "../common.h"

#define main library_main
#include "../../code_library/graphs/2SAT_lexicographic.cpp"
#undef main

typedef pair<int, int> Clause;  /// a or b, a == b forces a single literal

void add_all(Graph& g, const vector<Clause>& clauses){
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

/// Exhaustive: satisfiable iff some assignment works, and the answer is the lexicographically first one preferring true
void check_brute(Graph& g, int n, const vector<Clause>& clauses){
    int best = -1;
    auto key = [&](int mask){ int r = 0; for (int i = 0; i < n; i++) r = r * 2 + (mask >> i & 1); return r; };
    for (int mask = 0; mask < (1 << n); mask++){
        bool ok = true;
        for (auto [a, b] : clauses){
            bool va = (mask >> (abs(a) - 1) & 1) == (a > 0), vb = (mask >> (abs(b) - 1) & 1) == (b > 0);
            ok &= va || vb;
        }
        if (ok && (best == -1 || key(mask) > key(best))) best = mask;
    }

    assert(g.is_satisfiable() == (best != -1));
    if (best == -1) return;
    vector<int> expected;
    for (int x = 1; x <= n; x++) if (best >> (x - 1) & 1) expected.push_back(x);
    assert(g.get_assignment() == expected);
}

void check_greedy(int n, const vector<Clause>& clauses){
    Graph g(n);
    add_all(g, clauses);
    auto [sat, expected] = greedy(n, clauses);
    assert(g.is_satisfiable() == sat);
    if (sat) assert(g.get_assignment() == expected);
}

/// x_i = true walks the whole chain p, so attempts blow the budget and SCC, the topological shortcut and bitset batches
/// all run; a random length prefix reaches !x_i and is wasted work, so SCC starts at varying points, later x may skip
/// x_i -> p_1 and then must not reach !x_i, while p_k -> !x_i keeps them off the shortcut
/// A few forced false x fail within the budget after SCC, noise only adds !x_a or !x_b, which never forces an x true
vector<Clause> adversarial(int n, int k, bool contradiction){
    vector<Clause> clauses;
    /// The SCC switch happens once the waste of failed attempts (~k + n each) passes ~2 (vars + clauses), aim the
    /// reaching prefix near that point so the switch also lands on non-reaching x, whose answer is true
    long long size = 2LL * (n + k) + 2LL * (k + 2 * n);
    int reaching = max(1LL, min<long long>(n / 2, size / (k + n) + stress::rand_int(-3, 2)));
    if (stress::rand_int(0, 2) == 0) reaching = stress::rand_int(1, n / 2);

    for (int j = 1; j < k; j++) clauses.push_back({-(n + j), n + j + 1});
    for (int i = 1; i <= n; i++){
        if (i <= reaching || stress::rand_int(0, 1)) clauses.push_back({-i, n + 1});
        clauses.push_back({-(n + k), -i});
        if (i > reaching && stress::rand_int(0, 9) == 0) clauses.push_back({-i, -i});
    }
    for (int extra = stress::rand_int(0, n / 4); extra; extra--) clauses.push_back({-stress::rand_int(1, n), -stress::rand_int(1, n)});
    if (contradiction) clauses.push_back({n, n});  /// the last x, reached only after the waste has paid for SCC
    return clauses;
}

/// Shuffled signed literals form segments, half of them made positive ascending, false in a hidden assignment except p:
/// a chain l_1 -> ... -> l_k -> !l_k, half the time with a forced false decoy d -> !l_k so Tarjan meets the negated
/// chain first, or a hub, sources s -> p_1 -> ... -> p_k -> !s forcing each s false with no s an ancestor of another,
/// so each takes its own batch query. Variables 1..tail are a true chain fixed first, which chain literals point into,
/// so batches run into decided components
vector<Clause> chains(int n){
    int tail = stress::rand_int(0, n / 4), vars = tail + n;
    vector<int> lit(n);
    for (int i = 0; i < n; i++) lit[i] = stress::rand_int(0, 1) ? tail + i + 1 : -(tail + i + 1);
    shuffle(lit.begin(), lit.end(), stress::rng());

    vector<Clause> clauses;
    vector<char> hidden(vars + n + 1, 0), plain(vars + 1, 1);  /// hidden[x]: x is true, room for one decoy per segment
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
        if (stress::rand_int(0, 2) || end - start < 2){
            for (int i = start; i < end; i++){
                hidden[abs(lit[i])] = lit[i] < 0;
                if (i + 1 < end) clauses.push_back({-lit[i], lit[i + 1]});
                if (tail && stress::rand_int(0, 7) == 0) clauses.push_back({-lit[i], (int)stress::rand_int(1, tail)});
            }
            clauses.push_back({-lit[end - 1], -lit[end - 1]});
            if (stress::rand_int(0, 1)){
                decoy++;
                clauses.push_back({-decoy, -decoy}), clauses.push_back({-decoy, -lit[end - 1]});
            }
        }
        else {
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

    for (int extra = stress::rand_int(0, n / 4); extra;){
        int a = random_literal(vars), b = random_literal(vars);
        if (!plain[abs(a)] || !plain[abs(b)]) continue;  /// noise into a negated chain would fix its whole prefix
        if ((a > 0) == hidden[abs(a)] || (b > 0) == hidden[abs(b)]) clauses.push_back({a, b}), extra--;
    }
    return clauses;
}

/// Variable 1 fails first, its Tarjan DFS numbers the y_j and the path R they all imply low, so the y_j lead the sinks
/// first order with ranges spanning R. Sources g_i over a shared path Q burn the waste allowance, SCC starts on some
/// g_i, whose batch must answer it ahead of the y_j: taking y_j until x comes up costs k / B passes over R
vector<Clause> decoys(int k, int len){
    int lq = stress::rand_int(len / 12, len / 6), a = (6 * k + 4 * len) / lq + stress::rand_int(4, 12);
    int g0 = 1, q0 = g0 + a, y0 = q0 + lq, r0 = y0 + k;

    vector<Clause> clauses = {{-1, -(r0 + len)}};
    for (int j = 1; j <= k; j++) clauses.push_back({-1, y0 + j}), clauses.push_back({-(y0 + j), r0 + 1});
    for (int i = 1; i <= a; i++) clauses.push_back({-(g0 + i), q0 + 1}), clauses.push_back({-(q0 + lq), -(g0 + i)});
    for (int i = 1; i < lq; i++) clauses.push_back({-(q0 + i), q0 + i + 1});
    for (int i = 1; i < len; i++) clauses.push_back({-(r0 + i), r0 + i + 1});
    return clauses;
}

int main(){
    for (long long it = 0; it < stress::scaled(4000); it++){
        int n = stress::rand_int(1, 10), m = stress::rand_int(0, 3 * n);
        vector<Clause> clauses;
        Graph g(n);
        for (int c = 0; c < m; c++){
            int a = random_literal(n), b = random_literal(n), kind = stress::rand_int(0, 5);
            if (kind == 0) g.add_or(a, b), clauses.push_back({a, b});
            else if (kind == 1) g.add_implication(a, b), clauses.push_back({-a, b});
            else if (kind == 2) g.add_xor(a, b), clauses.push_back({a, b}), clauses.push_back({-a, -b});
            else if (kind == 3) g.add_and(a, b), clauses.push_back({a, b}), clauses.push_back({a, -b}), clauses.push_back({-a, b});
            else if (kind == 4) g.force_true(a), clauses.push_back({a, a});
            else g.force_false(a), clauses.push_back({-a, -a});
        }
        check_brute(g, n, clauses);
    }

    for (long long it = 0; it < stress::scaled(60); it++){
        int mode = it % 6, n = stress::rand_int(200, 1500);
        vector<Clause> clauses;
        if (mode < 2) clauses = adversarial(n, stress::rand_int(n / 2, 2 * n), mode == 1);
        else if (mode >= 4) clauses = chains(2 * n);
        else {
            int vars = 4 * n, m = mode == 2 ? 9 * vars / 10 : 3 * vars;  /// below the satisfiability threshold, and dense
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

    /// Solving again after more constraints must start from scratch
    Graph g(3);
    g.add_or(-1, -2);
    assert(g.is_satisfiable() && (g.get_assignment() == vector<int>{1, 3}));
    g.force_false(1);
    assert(g.is_satisfiable() && (g.get_assignment() == vector<int>{2, 3}));
    g.force_true(1);
    assert(!g.is_satisfiable());
    return 0;
}
