#include "../common.h"

#define main library_main
#include "../../code_library/graphs/2SAT_kosaraju.cpp"
#undef main

/// Random mixes of every constraint type, at_most_one included, against exhaustive search over all 2^n assignments
/// of the original variables, and the returned assignment must satisfy every constraint
void random_systems(){
    for (long long it = 0; it < stress::scaled(60000); it++){
        int n = stress::rand_int(1, 7), m = stress::rand_int(0, 12);
        Graph g(n);
        vector<array<int, 3>> constraints;  /// kind, a, b
        vector<vector<int>> groups;
        int vars = n;
        size_t edges = 0;

        auto literal = [&](){ int v = stress::rand_int(1, n); return stress::rand_int(0, 1) ? v : -v; };
        int midway = stress::rand_int(0, m);  /// an early query must not stop later constraints from counting
        for (int i = 0; i < m; i++){
            if (i == midway) g.is_satisfiable();
            int kind = stress::rand_int(0, 6), a = literal(), b = literal();
            if (kind == 0) g.add_implication(a, b), edges += 2;
            if (kind == 1) g.add_or(a, b), edges += 2;
            if (kind == 2) g.add_xor(a, b), edges += 4;
            if (kind == 3) g.add_and(a, b), edges += 2;
            if (kind == 4) g.force_true(a), edges += 1;
            if (kind == 5) g.force_false(a), edges += 1;
            if (kind == 6){
                vector<int> lits(stress::rand_int(0, 7));
                for (auto& x : lits) x = literal();
                g.at_most_one(lits);
                groups.push_back(lits);
                int k = lits.size();
                if (k >= 2) vars += k - 2, edges += 2 * (3 * k - 5);
            }
            if (kind != 6) constraints.push_back({kind, a, b});
        }
        assert(g.n == vars && g.edges.size() == edges);

        auto satisfies = [&](auto value){
            for (auto& c : constraints){
                bool a = value(c[1]), b = value(c[2]);
                if (c[0] == 0 && a && !b) return false;
                if (c[0] == 1 && !a && !b) return false;
                if (c[0] == 2 && a == b) return false;
                if (c[0] == 3 && !(a && b)) return false;
                if (c[0] == 4 && !a) return false;
                if (c[0] == 5 && a) return false;
            }
            for (auto& lits : groups){
                int count = 0;
                for (int x : lits) count += value(x);
                if (count > 1) return false;
            }
            return true;
        };

        bool satisfiable = false;
        for (int mask = 0; mask < (1 << n) && !satisfiable; mask++){
            satisfiable = satisfies([&](int x){ bool v = mask >> (abs(x) - 1) & 1; return x > 0 ? v : !v; });
        }

        assert(g.is_satisfiable() == satisfiable);
        if (satisfiable){
            for (int x = 1; x <= g.n; x++) assert(g.value(-x) == !g.value(x));
            assert(satisfies([&](int x){ return g.value(x); }));
        }
    }
}

/// One long at_most_one over shuffled distinct literals with a planted true literal: the rest must come out false,
/// and planting a second true literal must make it unsatisfiable
void long_groups(){
    for (long long it = 0; it < stress::scaled(20); it++){
        int n = stress::rand_int(2, 5000);
        vector<int> lits(n);
        for (int i = 0; i < n; i++) lits[i] = stress::rand_int(0, 1) ? i + 1 : -(i + 1);
        shuffle(lits.begin(), lits.end(), stress::rng());

        Graph g(n);
        g.at_most_one(lits);
        int a = stress::rand_int(0, n - 1), b = stress::rand_int(0, n - 2);
        if (b >= a) b++;
        g.force_true(lits[a]);
        assert(g.n == n + max(0, n - 2) && g.is_satisfiable());
        for (int i = 0; i < n; i++) assert(g.value(lits[i]) == (i == a));

        g.force_true(lits[b]);
        assert(!g.is_satisfiable());
    }
}

int main(){
    random_systems();
    long_groups();
    return 0;
}
