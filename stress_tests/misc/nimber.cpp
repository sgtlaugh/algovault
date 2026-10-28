#include "../common.h"

#define main library_main
#include "../../code_library/misc/nimber.cpp"
#undef main

unsigned long long brute_product(unsigned long long a, unsigned long long b);

/// 2^m * 2^n = 2^(m ^ n) times (F ^ F / 2) for each Fermat 2-power F = 2^(2^t) with bit t in both m and n
unsigned long long brute_power_product(int m, int n){
    static unsigned long long memo[64][64] = {};
    if (memo[m][n]) return memo[m][n];

    unsigned long long res = 1ULL << (m ^ n);
    for (int t = 0; t < 6; t++){
        if ((m & n) >> t & 1) res = brute_product(res, (1ULL << (1 << t)) | (1ULL << ((1 << t) - 1)));
    }

    return memo[m][n] = res;
}

unsigned long long brute_product(unsigned long long a, unsigned long long b){
    unsigned long long res = 0;
    for (int i = 0; i < 64; i++){
        if (!(a >> i & 1)) continue;
        for (int j = 0; j < 64; j++){
            if (b >> j & 1) res ^= brute_power_product(i, j);
        }
    }

    return res;
}

unsigned long long random_nimber(){
    int bits = stress::rand_int(0, 4) ? 64 : (int)stress::rand_int(1, 64);
    unsigned long long x = stress::rng()();
    return bits == 64 ? x : x & ((1ULL << bits) - 1);
}

int main(){
    /// The mex definition itself, a * b = mex{a' * b ^ a * b' ^ a' * b' : a' < a, b' < b}, on every pair below 64
    const int M = 64;
    vector<vector<int>> mex_table(M, vector<int>(M));
    for (int a = 0; a < M; a++){
        for (int b = 0; b < M; b++){
            vector<char> seen(256, 0);  /// 64 is not a Fermat 2-power, products of values below 64 reach 255
            for (int x = 0; x < a; x++){
                for (int y = 0; y < b; y++) seen[mex_table[x][b] ^ mex_table[a][y] ^ mex_table[x][y]] = 1;
            }
            while (seen[mex_table[a][b]]) mex_table[a][b]++;
            assert(nim_product(a, b) == (unsigned long long)mex_table[a][b]);
            assert(brute_product(a, b) == (unsigned long long)mex_table[a][b]);
        }
    }

    /// Every pair of bytes, the whole base table
    for (unsigned long long a = 0; a < 256; a++){
        for (unsigned long long b = 0; b < 256; b++) assert(nim_product(a, b) == brute_product(a, b));
    }

    /// Every pair of powers of two, which hits each table row and every split level
    for (int i = 0; i < 64; i++){
        for (int j = 0; j < 64; j++) assert(nim_product(1ULL << i, 1ULL << j) == brute_power_product(i, j));
    }

    /// Fermat 2-powers square to 3/2 of themselves
    for (int t = 0; t < 6; t++){
        unsigned long long f = 1ULL << (1 << t);
        assert(nim_product(f, f) == (f | f >> 1));
    }

    /// Random values of every bit length against the brute force
    for (long long it = 0; it < stress::scaled(10000); it++){
        unsigned long long a = random_nimber(), b = random_nimber();
        if (it % 50 == 0) a = ~0ULL - stress::rand_int(0, 3);
        assert(nim_product(a, b) == brute_product(a, b));
    }

    /// Field axioms on random 64-bit values
    for (long long it = 0; it < stress::scaled(20000); it++){
        unsigned long long a = random_nimber(), b = random_nimber(), c = random_nimber();
        assert(nim_product(a, b) == nim_product(b, a));
        assert(nim_product(nim_product(a, b), c) == nim_product(a, nim_product(b, c)));
        assert(nim_product(a, b ^ c) == (nim_product(a, b) ^ nim_product(a, c)));
        assert(nim_product(a, 1) == a);
    }

    /// Inverses, the largest values included
    for (long long it = 0; it < stress::scaled(300); it++){
        unsigned long long a = it < 4 ? ~0ULL - it : random_nimber();
        if (a == 0) continue;
        assert(nim_product(a, nim_inverse(a)) == 1);
    }
    stress::expect_assert("misc/nimber.cpp", []{ nim_inverse(0); });

    return 0;
}
