/***
 *
 * DSU with Potentials
 * Union-find that also knows how much each value differs from the others in its set
 *
 * Complexity: amortized O(α(n)) per operation, O(n) memory
 *
 * PotentialDSU dsu(n, mod = 0): nodes 0..n, all values unknown
 * dsu.relate(a, b, d): record x[b] = x[a] + d, false if that contradicts what is known (nothing changes then)
 * dsu.diff(a, b): x[b] - x[a], a and b must be in the same set
 * dsu.same(a, b): true if a and b are in the same set, so diff(a, b) is known
 * dsu.size(a): number of nodes in a's set
 *
 * Example, the heights of A, B, C and D:
 *     dsu.relate(A, B, 3);    // true, B is 3 taller than A
 *     dsu.relate(B, C, 2);    // true, C is 2 taller than B
 *     dsu.diff(A, C);         // 5, so C is 5 taller than A
 *     dsu.diff(C, A);         // -5
 *     dsu.relate(A, C, 5);    // true, already implied
 *     dsu.relate(C, A, 1);    // false, A would be 1 taller than C, but it is 5 shorter
 *     dsu.same(A, D);         // false, nothing relates D yet
 *
 * mod = 0: plain integers, every implied difference must fit in long long
 * mod > 0: differences modulo mod, mod <= 1e18; mod = 2 is parity, relate(a, b, 1) means "a and b differ"
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct PotentialDSU{
    long long mod;
    vector<int> parent, sz;
    vector<long long> up;  /// x[i] = x[parent[i]] + up[i], so x[i] = x[root] + up[i] right after find(i)

    PotentialDSU(int n, long long mod = 0) : mod(mod), parent(n + 1), sz(n + 1, 1), up(n + 1, 0){
        assert(0 <= mod && mod <= (long long)1e18);
        iota(parent.begin(), parent.end(), 0);
    }

    int find(int i){
        if (parent[i] == i) return i;
        int root = find(parent[i]);
        up[i] = norm(up[i] + up[parent[i]]);
        return parent[i] = root;
    }

    bool relate(int a, int b, long long d){
        int ra = find(a), rb = find(b);
        d = norm(d);
        if (ra == rb) return norm(up[b] - up[a]) == d;

        /// x[rb] = x[b] - up[b] = x[ra] + up[a] + d - up[b]; d + up[a] first, it is an implied difference
        long long w = norm(d + up[a] - up[b]);
        if (sz[ra] < sz[rb]) swap(ra, rb), w = norm(-w);
        parent[rb] = ra, up[rb] = w, sz[ra] += sz[rb];
        return true;
    }

    long long diff(int a, int b){
        assert(find(a) == find(b));
        return norm(up[b] - up[a]);
    }

    bool same(int a, int b){
        return find(a) == find(b);
    }

    int size(int a){
        return sz[find(a)];
    }

    long long norm(long long x) const{
        if (!mod) return x;
        x %= mod;
        return x < 0 ? x + mod : x;
    }
};

int main(){
    /// The header example, A = 0, B = 1, C = 2, D = 3
    auto heights = PotentialDSU(3);
    assert(heights.relate(0, 1, 3) && heights.relate(1, 2, 2));
    assert(heights.diff(0, 2) == 5 && heights.diff(2, 0) == -5);
    assert(heights.relate(0, 2, 5) && !heights.relate(2, 0, 1) && !heights.same(0, 3));
    assert(heights.diff(0, 2) == 5 && heights.size(0) == 3 && heights.size(3) == 1);

    /// AtCoder ABC087 D samples: x[R] = x[L] + D
    auto line = PotentialDSU(3);
    assert(line.relate(1, 2, 1) && line.relate(2, 3, 1) && line.relate(1, 3, 2));
    assert(line.diff(1, 3) == 2 && line.diff(3, 1) == -2 && line.diff(2, 2) == 0);

    auto liar = PotentialDSU(3);
    assert(liar.relate(1, 2, 1) && liar.relate(2, 3, 1) && !liar.relate(1, 3, 5));
    assert(liar.diff(1, 3) == 2 && liar.size(2) == 3);

    auto huge = PotentialDSU(3);
    assert(huge.relate(0, 1, (long long)1e18) && huge.relate(2, 1, -(long long)1e18));
    assert(huge.diff(0, 2) == (long long)2e18 && huge.diff(2, 0) == -(long long)2e18);

    /// Parity: an odd cycle contradicts, an even cycle does not
    auto parity = PotentialDSU(7, 2);
    assert(parity.relate(1, 2, 1) && parity.relate(2, 3, 1) && !parity.relate(3, 1, 1));
    assert(parity.diff(1, 3) == 0 && parity.diff(3, 2) == 1);
    assert(parity.relate(4, 5, 1) && parity.relate(5, 6, 1) && parity.relate(6, 7, 1) && parity.relate(7, 4, 1));
    assert(parity.relate(4, 6, -2) && !parity.relate(5, 7, 3) && !parity.same(1, 4));

    /// Rock-paper-scissors: x[b] = x[a] + 1 (mod 3) when a beats b
    auto rps = PotentialDSU(3, 3);
    assert(rps.relate(1, 2, 1) && rps.relate(2, 3, 1) && rps.relate(3, 1, 1) && !rps.relate(1, 3, 1));
    assert(rps.diff(1, 3) == 2 && rps.diff(3, 1) == 1);

    auto wide = PotentialDSU(2, (long long)1e18);
    assert(wide.relate(0, 1, -1) && wide.diff(0, 1) == (long long)1e18 - 1 && wide.diff(1, 0) == 1);
    assert(wide.relate(1, 0, 1 - 3 * (long long)1e18) && !wide.relate(1, 0, 2));

    return 0;
}
