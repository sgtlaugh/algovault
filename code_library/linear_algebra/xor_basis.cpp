/***
 *
 * XOR Basis
 * Linear basis over GF(2) of 64-bit values: the span is the set of XORs of all subsets of inserted values
 *
 * Complexity: O(64) per insert and query, O(64^2) for intersect, 520 bytes per instance
 *
 * Values are any unsigned 64-bit integers, including 2^64 - 1
 * The basis is kept in reduced row echelon form: basis[i] is 0 or has leading bit i,
 * and no other basis vector has bit i set
 *
 *   insert(x)        adds x, returns false if x was already in the span
 *   contains(x)      whether some subset XORs to x
 *   max_xor(x)       max of x ^ s over the span (max subset XOR for x = 0)
 *   min_xor(x)       min of x ^ s over the span
 *   kth(k)           k-th smallest distinct span value, 0-indexed, requires k < 2^rank (kth(0) = 0)
 *   count_less(x)    number of distinct span values < x
 *   intersect(other) basis of the intersection of both spans (Zassenhaus)
 *
 * The span has 2^rank distinct values, and each one is the XOR of exactly 2^(n - rank)
 * of the 2^n subsets of the n inserted values
 *
 * Example:
 *   XorBasis b;
 *   b.insert(6), b.insert(10), b.insert(12);  // span {0, 6, 10, 12}, the last insert returns false
 *   b.kth(2) == 10, b.count_less(11) == 3, b.max_xor(1) == 13
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

using ull = unsigned long long;

struct XorBasis{
    ull basis[64] = {};
    int rank = 0;

    bool insert(ull x){
        for (int i = 63; i >= 0; i--){
            if (x >> i & 1) x ^= basis[i];
        }
        if (!x) return false;

        /// x now has no pivot bits, so clearing its leading bit from higher vectors keeps the form reduced
        int top = 63 - __builtin_clzll(x);
        for (int i = top + 1; i < 64; i++){
            if (basis[i] >> top & 1) basis[i] ^= x;
        }

        basis[top] = x, rank++;
        return true;
    }

    bool contains(ull x) const{
        return min_xor(x) == 0;
    }

    ull max_xor(ull x = 0) const{
        for (int i = 63; i >= 0; i--) x = max(x, x ^ basis[i]);
        return x;
    }

    ull min_xor(ull x) const{
        for (int i = 63; i >= 0; i--) x = min(x, x ^ basis[i]);
        return x;
    }

    ull kth(ull k) const{
        assert(rank == 64 || (k >> rank) == 0);
        ull res = 0;
        for (int i = 0, j = 0; i < 64; i++){
            if (!basis[i]) continue;
            if (k >> j & 1) res ^= basis[i];
            j++;
        }
        return res;
    }

    ull count_less(ull x) const{
        ull res = 0, cur = 0;
        int below = rank;

        /// cur is the span value matching x on the bits above i, below counts the pivots under i
        for (int i = 63; i >= 0; i--){
            bool want = x >> i & 1, differ = (cur ^ x) >> i & 1;
            if (basis[i]){
                below--;
                if (want) res += 1ULL << below;
                if (differ) cur ^= basis[i];
            }
            else if (differ) return want ? res + (1ULL << below) : res;
        }
        return res;
    }

    XorBasis intersect(const XorBasis& other) const{
        XorBasis res;
        ull head[64] = {}, tail[64] = {};

        /// Rows are (u, u) for this basis and (v, 0) for the other; a row reduced to (0, w) has w in both spans
        auto add_row = [&](ull h, ull t){
            for (int i = 63; i >= 0; i--){
                if (!(h >> i & 1)) continue;
                if (!head[i]){
                    head[i] = h, tail[i] = t;
                    return;
                }
                h ^= head[i], t ^= tail[i];
            }
            res.insert(t);
        };

        for (ull u : basis){
            if (u) add_row(u, u);
        }
        for (ull v : other.basis){
            if (v) add_row(v, 0);
        }
        return res;
    }
};

int main(){
    XorBasis b;
    assert(b.insert(6));
    assert(b.insert(10));
    assert(!b.insert(12));          /// 6 ^ 10, the span is {0, 6, 10, 12}
    assert(b.rank == 2);
    assert(b.contains(12));
    assert(!b.contains(4));
    assert(b.max_xor() == 12);      /// max subset XOR
    assert(b.max_xor(1) == 13);     /// 1 ^ 12
    assert(b.min_xor(5) == 3);      /// 5 ^ 6
    assert(b.kth(2) == 10);
    assert(b.count_less(11) == 3);  /// 0, 6 and 10

    XorBasis other;
    other.insert(12), other.insert(1);
    XorBasis common = b.intersect(other);  /// {0, 6, 10, 12} and {0, 1, 12, 13} share {0, 12}
    assert(common.rank == 1);
    assert(common.contains(12));
    return 0;
}
