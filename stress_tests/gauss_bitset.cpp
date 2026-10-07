#include "common.h"

#define main library_main
#include "../code_library/gauss_bitset.cpp"
#undef main

bool satisfies(const vector<bitset<MAX>>& equations, int m, const bitset<MAX>& x){
    for (const auto& e : equations){
        int v = 0;
        for (int j = 0; j < m; j++) v ^= e[j] & x[j];
        if (v != e[m]) return false;
    }
    return true;
}

int main(){
    for (long long it = 0; it < stress::scaled(3000); it++){
        int m = stress::rand_int(0, 12), n = stress::rand_int(0, 14), density = stress::rand_int(1, 3);
        vector<bitset<MAX>> equations(n);
        for (auto& e : equations){
            for (int j = 0; j <= m; j++) e[j] = stress::rand_int(0, density) == 0;
        }
        /// Duplicated and combined rows make dependent systems common
        if (n >= 2 && stress::rand_int(0, 1)) equations[stress::rand_int(0, n - 1)] = equations[0] ^ equations[n - 1];

        long long solutions = 0;
        for (int mask = 0; mask < (1 << m); mask++) solutions += satisfies(equations, m, bitset<MAX>(mask));

        bitset<MAX> res;
        res.set();  /// stale bits must be cleared by gauss
        int f_var = gauss(m, equations, res);
        if (!solutions) assert(f_var == -1);
        else{
            assert(f_var >= 0 && solutions == 1LL << f_var);
            assert(satisfies(equations, m, res));
            assert((res >> m).none());
        }
    }

    /// Large systems: a planted solution must be recovered as a valid one with the right rank
    for (long long it = 0; it < stress::scaled(30); it++){
        int m = stress::rand_int(1, MAX - 1), n = stress::rand_int(1, 300), rank = stress::rand_int(0, min(n, m));
        bitset<MAX> planted;
        for (int j = 0; j < m; j++) planted[j] = stress::rand_int(0, 1);

        vector<bitset<MAX>> basis(rank), equations(n);
        for (int r = 0; r < rank; r++){
            basis[r][stress::rand_int(0, m - 1)] = 1;  /// may still be dependent, the rank is only an upper bound
            for (int j = 0; j < m; j++) if (stress::rand_int(0, 1)) basis[r][j] = 1;
        }
        for (auto& e : equations){
            for (int r = 0; r < rank; r++) if (stress::rand_int(0, 1)) e ^= basis[r];
            e[m] = (e & planted).count() & 1;
        }

        bitset<MAX> res;
        int f_var = gauss(m, equations, res);
        assert(f_var >= m - rank && satisfies(equations, m, res));
    }
    return 0;
}
