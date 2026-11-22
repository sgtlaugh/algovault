/***
 *
 * XOR Basis
 * Linear basis over GF(2) of 64-bit values: the span is the set of XORs of all subsets of inserted values
 *
 * Complexity:
 *   - XorBasis: O(64) per insert and query, O(64^2) for intersect, 520 bytes per instance
 *   - RangeXorBasis: O(BITS) per push and query, 12 BITS + 8 bytes per instance, 4 less for odd BITS
 *
 * XorBasis b: an empty basis, its span is {0}
 * b.insert(x): adds x, false if x was already in the span
 * b.contains(x): whether some subset XORs to x
 * b.count_less(x): the number of distinct span values below x
 * b.intersect(other): the basis of the intersection of both spans
 * b.kth(k): the k-th smallest distinct span value, 0-indexed, k < 2^rank, kth(0) = 0
 * b.max_xor(x = 0): the max of x ^ s over the span, the max subset XOR for x = 0
 * b.min_xor(x): the min of x ^ s over the span
 * b.rank: the number of nonzero basis vectors
 *
 * RangeXorBasis<BITS = 64> rb: an empty array, values below 2^BITS
 * rb.push(x): appends x as a[size()]
 * rb.contains(l, x): whether some subset of a[l..size() - 1] XORs to x, 0 <= l
 * rb.max_xor(l, x = 0): the max of x ^ s over the subset XORs s of a[l..size() - 1]
 * rb.rank(l): the rank of a[l..size() - 1], its span has 2^rank(l) values
 * rb.size(): the number of values pushed
 *
 * Notes:
 *   - Values are any unsigned 64-bit integers, including 2^64 - 1
 *   - XorBasis keeps basis[i] 0 or with leading bit i, and no other basis vector has bit i set
 *   - The span has 2^rank values, each the XOR of exactly 2^(n - rank) of the 2^n subsets of the n inserted values
 *   - RangeXorBasis answers any [l, r] offline with queries sorted by r, asking each once size() = r + 1
 *   - Online [l, r] queries keep a copy of rb after each push: 776 bytes each at BITS = 64, so narrow values pass BITS
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

    ull max_xor(ull x = 0) const{
        for (int i = 63; i >= 0; i--) x = max(x, x ^ basis[i]);
        return x;
    }

    ull min_xor(ull x) const{
        for (int i = 63; i >= 0; i--) x = min(x, x ^ basis[i]);
        return x;
    }
};

template<int BITS = 64>
struct RangeXorBasis{
    ull basis[BITS] = {};
    int pos[BITS] = {}, count = 0;

    /// basis[i] is 0 or has leading bit i and lies in the span of a[pos[i]..], so suffix l uses those with pos[i] >= l
    void push(ull x){
        assert(x >> (BITS - 1) <= 1);  /// x < 2^BITS without a shift by 64
        int p = count++;
        for (int i = BITS - 1; i >= 0; i--){
            if (!(x >> i & 1)) continue;
            if (!basis[i]){
                basis[i] = x, pos[i] = p;
                return;
            }

            /// The newer vector stays, the older one moves down, so every suffix keeps a basis of its own span
            if (pos[i] < p) swap(basis[i], x), swap(pos[i], p);
            x ^= basis[i];
        }
    }

    bool contains(int l, ull x) const{
        for (int i = BITS - 1; i >= 0; i--){
            if (!(x >> i & 1)) continue;
            if (!basis[i] || pos[i] < l) return false;
            x ^= basis[i];
        }

        return x == 0;
    }

    ull max_xor(int l, ull x = 0) const{
        for (int i = BITS - 1; i >= 0; i--){
            if (basis[i] && pos[i] >= l) x = max(x, x ^ basis[i]);
        }

        return x;
    }

    int rank(int l) const{
        int res = 0;
        for (int i = 0; i < BITS; i++) res += basis[i] && pos[i] >= l;
        return res;
    }

    int size() const{
        return count;
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

    /// a = 5, 3, 8, 1: a[0..3] spans all 16 values, a[1..3] = {3, 8, 1} spans 8, a[2..3] = {8, 1} spans {0, 1, 8, 9}
    RangeXorBasis<4> rb;
    rb.push(5), rb.push(3), rb.push(8), rb.push(1);
    assert(rb.max_xor(0) == 15);
    assert(rb.max_xor(1) == 11);     /// 8 ^ 3
    assert(rb.max_xor(2) == 9);      /// 8 ^ 1
    assert(rb.max_xor(2, 6) == 15);  /// 6 ^ 9
    assert(rb.contains(1, 10));      /// 3 ^ 8 ^ 1
    assert(!rb.contains(2, 3));
    assert(rb.rank(0) == 4);
    assert(rb.rank(2) == 2);
    assert(rb.size() == 4);
    return 0;
}
