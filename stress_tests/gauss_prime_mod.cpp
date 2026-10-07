#include "common.h"

#define main library_main
#include "../code_library/gauss_prime_mod.cpp"
#undef main

bool satisfies(const vector<vector<int>>& equations, const vector<int>& x, int mod){
    int m = x.size();
    for (const auto& e : equations){
        long long v = 0;
        for (int j = 0; j < m; j++) v = (v + (long long)e[j] * x[j]) % mod;
        if (v != e[m]) return false;
    }
    return true;
}

int main(){
    const int small_primes[] = {2, 3, 5, 7};
    for (long long it = 0; it < stress::scaled(3000); it++){
        int mod = small_primes[stress::rand_int(0, 3)];
        int m = stress::rand_int(1, mod == 2 ? 10 : mod == 3 ? 6 : 4), n = stress::rand_int(1, 7);
        vector<vector<int>> equations(n, vector<int>(m + 1));
        int density = stress::rand_int(0, 2);
        for (auto& e : equations) for (auto& x : e) x = stress::rand_int(0, density) ? stress::rand_int(0, mod - 1) : 0;
        if (n >= 2 && stress::rand_int(0, 1)){
            int r = stress::rand_int(1, n - 1), c = stress::rand_int(1, mod - 1);
            for (int j = 0; j <= m; j++) equations[r][j] = (equations[r][j] + c * equations[0][j]) % mod;  /// a dependent row
        }

        long long solutions = 0, total = 1;
        for (int j = 0; j < m; j++) total *= mod;
        vector<int> x(m);
        for (long long code = 0; code < total; code++){
            long long c = code;
            for (int j = 0; j < m; j++) x[j] = c % mod, c /= mod;
            solutions += satisfies(equations, x, mod);
        }

        vector<int> res = {42, 42};  /// stale contents must be replaced
        int f_var = gauss(equations, res, mod);
        if (!solutions) assert(f_var == -1);
        else{
            long long expected = 1;
            for (int j = 0; j < f_var; j++) expected *= mod;
            assert(f_var >= 0 && solutions == expected);
            assert((int)res.size() == m && satisfies(equations, res, mod));
            for (int v : res) assert(0 <= v && v < mod);
        }
    }

    /// A large prime: rank r systems built as products of n x r and r x m random factors, full rank r almost surely
    const int mod = 1000000007;
    for (long long it = 0; it < stress::scaled(150); it++){
        int n = stress::rand_int(1, it % 10 ? 10 : 60), m = stress::rand_int(1, it % 10 ? 10 : 60), r = stress::rand_int(0, min(n, m));
        vector<vector<long long>> left(n, vector<long long>(r)), right(r, vector<long long>(m));
        for (auto& row : left) for (auto& v : row) v = stress::rand_int(0, mod - 1);
        for (auto& row : right) for (auto& v : row) v = stress::rand_int(0, mod - 1);
        vector<int> planted(m);
        for (auto& v : planted) v = stress::rand_int(0, mod - 1);

        vector<vector<int>> equations(n, vector<int>(m + 1));
        for (int i = 0; i < n; i++){
            for (int j = 0; j < m; j++){
                long long v = 0;
                for (int k = 0; k < r; k++) v = (v + left[i][k] * right[k][j]) % mod;
                equations[i][j] = v;
            }
            long long v = 0;
            for (int j = 0; j < m; j++) v = (v + (long long)equations[i][j] * planted[j]) % mod;
            equations[i][m] = v;
        }

        if (r < n){
            /// The other n - 1 >= r rows span the row space almost surely, so shifting any RHS makes the system inconsistent
            vector<vector<int>> probe = equations;
            int row = stress::rand_int(0, n - 1);
            probe[row][m] = (probe[row][m] + stress::rand_int(1, mod - 1)) % mod;
            vector<int> res;
            assert(gauss(probe, res, mod) == -1);
        }

        vector<int> res;
        int f_var = gauss(equations, res, mod);
        assert(f_var == m - r && satisfies(equations, res, mod));
    }
    return 0;
}
