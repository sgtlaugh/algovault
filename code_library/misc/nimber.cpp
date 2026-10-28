/***
 *
 * Nimber Arithmetic
 * Nim product and inverse of 64-bit nimbers, the field GF(2^64) of Conway's nim arithmetic
 *
 * Complexity: O(1) per nim_product (46 table lookups), 126 products per nim_inverse, 78 KB of tables
 *
 * nim_product(a, b): a times b in nim arithmetic, commutative, associative and distributive over the nim sum a ^ b
 * nim_inverse(a): the x with nim_product(a, x) == 1, a != 0
 *
 * Notes:
 *   - Nimbers below 2^(2^k) are closed under both operations, a subfield: GF(4) below 4, GF(16) below 16
 *   - Products of 2^(2^k) for distinct k are ordinary products, and 2^(2^k) squared is 3 2^(2^k - 1): 4 * 4 = 6
 *   - The tables are built on the first call, about 0.3 ms
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// With x = 2^k, x * x is the nim sum of x and 2^(k - 1), so a product of 2k-bit values needs 3 k-bit products, not 4
struct NimberTables{
    unsigned char small[256][256] = {};
    unsigned long long half[6][256]; /// rows 0..1 for k = 16 and 2..5 for k = 32: byte i of c, times 2^(k - 1)

    NimberTables(){
        small[1][1] = 1;
        for (int h = 2; h <= 8; h *= 2){
            for (int a = 0; a < (1 << h); a++){
                for (int b = 0; b < (1 << h); b++){
                    if (max(a, b) >= (1 << (h / 2))) small[a][b] = split(a, b, h);
                }
            }
        }

        for (int k = 16; k <= 32; k *= 2){
            for (int i = 0; i < k / 8; i++){
                for (int v = 0; v < 256; v++){
                    half[k / 8 - 2 + i][v] = split((unsigned long long)v << (8 * i), 1ULL << (k - 1), k);
                }
            }
        }
    }

    unsigned long long product(unsigned long long a, unsigned long long b, int h) const{
        return h <= 8 ? small[a][b] : split(a, b, h);
    }

    /// a, b < 2^h, and every table entry it reads for width h / 2 is already built
    unsigned long long split(unsigned long long a, unsigned long long b, int h) const{
        int k = h / 2;
        unsigned long long mask = (1ULL << k) - 1;
        unsigned long long ah = a >> k, al = a & mask, bh = b >> k, bl = b & mask;

        unsigned long long low = product(al, bl, k), high = product(ah, bh, k), mid = product(ah ^ al, bh ^ bl, k);
        return ((mid ^ low) << k) ^ low ^ times_half(high, k);
    }

    /// c * 2^(k - 1) for c < 2^k, linear in c, so one lookup per byte of c
    unsigned long long times_half(unsigned long long c, int k) const{
        if (k <= 8) return small[c][1 << (k - 1)];

        unsigned long long res = 0;
        for (int i = 0; i < k / 8; i++) res ^= half[k / 8 - 2 + i][c >> (8 * i) & 255];
        return res;
    }
};

unsigned long long nim_product(unsigned long long a, unsigned long long b){
    static const NimberTables tables;
    return tables.product(a, b, 64);
}

/// The multiplicative group has order 2^64 - 1, so a^(2^64 - 2) is the inverse
unsigned long long nim_inverse(unsigned long long a){
    assert(a != 0);

    unsigned long long res = 1;
    for (int i = 1; i < 64; i++){
        a = nim_product(a, a);
        res = nim_product(res, a);
    }

    return res;
}

int main(){
    /// The field GF(4) on {0, 1, 2, 3}: 2 * 2 = 3, so 2 * 3 = 3 ^ 2 and 3 * 3 = 3 ^ 2 ^ 2 ^ 1
    assert(nim_product(2, 2) == 3);
    assert(nim_product(2, 3) == 1);
    assert(nim_product(3, 3) == 2);

    /// Distinct 2^(2^k) multiply as integers, equal ones square to 3/2 of themselves, the rest follows by distributing
    assert(nim_product(4, 2) == 8);
    assert(nim_product(1ULL << 32, 1ULL << 31) == 1ULL << 63);  /// 2^31 = 2 * 4 * 16 * 256 * 65536
    assert(nim_product(16, 16) == 24);
    assert(nim_product(2, 5) == 10);  /// 2 * (4 ^ 1) = 8 ^ 2

    /// Division by y is a product with nim_inverse(y)
    assert(nim_inverse(2) == 3);

    unsigned long long x = 1000000007, y = 998244353;
    assert(nim_product(nim_product(x, y), nim_inverse(y)) == x);

    return 0;
}
