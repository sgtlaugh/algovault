#include "../common.h"

#define main library_main
#include "../../code_library/graphs/2SAT_kosaraju.cpp"
#undef main

/// Random mixes of every constraint type against exhaustive search over all 2^n assignments,
/// and the returned assignment must satisfy every constraint
int main(){
    for (long long it = 0; it < stress::scaled(60000); it++){
        int n = stress::rand_int(1, 7), m = stress::rand_int(0, 12);
        Graph g(n);
        vector<array<int, 3>> constraints;  /// kind, a, b

        auto literal = [&](){ int v = stress::rand_int(1, n); return stress::rand_int(0, 1) ? v : -v; };
        int midway = stress::rand_int(0, m);  /// an early query must not stop later constraints from counting
        for (int i = 0; i < m; i++){
            if (i == midway) g.is_satisfiable();
            int kind = stress::rand_int(0, 5), a = literal(), b = literal();
            if (kind == 0) g.add_implication(a, b);
            if (kind == 1) g.add_or(a, b);
            if (kind == 2) g.add_xor(a, b);
            if (kind == 3) g.add_and(a, b);
            if (kind == 4) g.force_true(a);
            if (kind == 5) g.force_false(a);
            constraints.push_back({kind, a, b});
        }

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
            return true;
        };

        bool satisfiable = false;
        for (int mask = 0; mask < (1 << n) && !satisfiable; mask++){
            satisfiable = satisfies([&](int x){ bool v = mask >> (abs(x) - 1) & 1; return x > 0 ? v : !v; });
        }

        assert(g.is_satisfiable() == satisfiable);
        if (satisfiable){
            for (int x = 1; x <= n; x++) assert(g.value(-x) == !g.value(x));
            assert(satisfies([&](int x){ return g.value(x); }));
        }
    }
    return 0;
}
