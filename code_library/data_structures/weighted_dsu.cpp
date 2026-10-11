/***
 *
 * Weighted Disjoint Set Union
 * Unknowns x[0..n] with relations x[b] - x[a] = d, over the integers (mod = 0) or modulo mod
 * Nodes can be 0 or 1 based, allocates nodes 0..n
 *
 * Complexity: amortized O(α(n)) per operation, where α(n) is the inverse Ackermann function, O(n) memory
 *
 * connect(a, b, d): add the relation, false if it contradicts the earlier ones, which leaves the state unchanged
 * diff(a, b): x[b] - x[a] for connected a and b, reduced to [0, mod) when mod > 0
 * Parity (bipartiteness, same / different constraints) is mod = 2 with d = 1 for different
 * Requires 0 <= mod <= 1e18, and with mod = 0 every difference implied by the relations must fit in long long
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct WeightedDSU{
    long long mod;
    vector<int> counter, parent;
    vector<long long> pot;  /// x[i] - x[parent[i]]

    WeightedDSU(int n, long long mod = 0) : mod(mod), counter(n + 1, 1), parent(n + 1), pot(n + 1, 0){
        assert(0 <= mod && mod <= (long long)1e18);
        iota(parent.begin(), parent.end(), 0);
    }

    /// Recursion depth is O(log n) thanks to union by size
    int find_root(int i){
        if (i == parent[i]) return i;

        int root = find_root(parent[i]);
        pot[i] = reduce(pot[i] + pot[parent[i]]);
        return parent[i] = root;
    }

    bool connect(int a, int b, long long d){
        int ra = find_root(a), rb = find_root(b);
        d = reduce(d);
        if (ra == rb) return reduce(pot[b] - pot[a]) == d;

        /// d + pot[a] = x[b] - x[ra] is an implied difference, pot[a] - pot[b] alone is not and can overflow
        long long shift = reduce(d + pot[a] - pot[b]);
        if (counter[rb] > counter[ra]) swap(ra, rb), shift = reduce(-shift);
        parent[rb] = ra, pot[rb] = shift, counter[ra] += counter[rb];
        return true;
    }

    long long diff(int a, int b){
        assert(find_root(a) == find_root(b));
        return reduce(pot[b] - pot[a]);
    }

    bool is_connected(int a, int b){
        return find_root(a) == find_root(b);
    }

    int component_size(int i){
        return counter[find_root(i)];
    }

    long long reduce(long long x) const{
        if (!mod) return x;
        x %= mod;
        return x < 0 ? x + mod : x;
    }
};

int main(){
    /// AtCoder ABC087 D samples: x[R] - x[L] = D
    auto line = WeightedDSU(3);
    assert(line.connect(1, 2, 1) && line.connect(2, 3, 1) && line.connect(1, 3, 2));
    assert(line.diff(1, 3) == 2 && line.diff(3, 1) == -2 && line.diff(2, 2) == 0);

    auto liar = WeightedDSU(3);
    assert(liar.connect(1, 2, 1) && liar.connect(2, 3, 1) && !liar.connect(1, 3, 5));
    assert(liar.diff(1, 3) == 2 && liar.component_size(2) == 3);

    auto huge = WeightedDSU(3);
    assert(huge.connect(0, 1, (long long)1e18) && huge.connect(2, 1, -(long long)1e18));
    assert(huge.diff(0, 2) == (long long)2e18 && huge.diff(2, 0) == -(long long)2e18);

    /// Parity: an odd cycle contradicts, an even cycle does not
    auto parity = WeightedDSU(7, 2);
    assert(parity.connect(1, 2, 1) && parity.connect(2, 3, 1) && !parity.connect(3, 1, 1));
    assert(parity.diff(1, 3) == 0 && parity.diff(3, 2) == 1);
    assert(parity.connect(4, 5, 1) && parity.connect(5, 6, 1) && parity.connect(6, 7, 1) && parity.connect(7, 4, 1));
    assert(parity.connect(4, 6, -2) && !parity.connect(5, 7, 3) && !parity.is_connected(1, 4));

    /// Rock-paper-scissors: x[b] - x[a] = 1 (mod 3) when a beats b
    auto rps = WeightedDSU(3, 3);
    assert(rps.connect(1, 2, 1) && rps.connect(2, 3, 1) && rps.connect(3, 1, 1) && !rps.connect(1, 3, 1));
    assert(rps.diff(1, 3) == 2 && rps.diff(3, 1) == 1);

    auto wide = WeightedDSU(2, (long long)1e18);
    assert(wide.connect(0, 1, -1) && wide.diff(0, 1) == (long long)1e18 - 1 && wide.diff(1, 0) == 1);
    assert(wide.connect(1, 0, 1 - 3 * (long long)1e18) && !wide.connect(1, 0, 2));

    return 0;
}
