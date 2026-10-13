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
    /// The header example
    int A = 0, B = 1, C = 2, D = 3;
    PotentialDSU heights(3);
    assert(heights.relate(A, B, 3));
    assert(heights.relate(B, C, 2));
    assert(heights.diff(A, C) == 5);
    assert(heights.diff(C, A) == -5);
    assert(heights.relate(A, C, 5));
    assert(!heights.relate(C, A, 1));
    assert(!heights.same(A, D));
    assert(heights.size(A) == 3);

    /// Parity: an odd cycle contradicts
    PotentialDSU parity(3, 2);
    assert(parity.relate(1, 2, 1));
    assert(parity.relate(2, 3, 1));
    assert(!parity.relate(3, 1, 1));
    assert(parity.diff(1, 3) == 0);

    /// An implied difference can exceed every single relation
    PotentialDSU huge(2);
    assert(huge.relate(0, 1, (long long)1e18));
    assert(huge.relate(2, 1, -(long long)1e18));
    assert(huge.diff(0, 2) == (long long)2e18);

    return 0;
}
