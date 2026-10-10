/***
 *
 * Disjoint Set Union and its variants
 * Plain, with rollback, with weighted relations (potentials, parity) and partially persistent
 * Nodes can be 0 or 1 based, every struct allocates nodes 0..n
 * connect returns true if a and b were in different components
 *
 * Complexity:
 *   - DSU, WeightedDSU: amortized O(α(n)) per operation, where α(n) is the inverse Ackermann function
 *   - RollbackDSU, PartiallyPersistentDSU: O(log n) per operation, union by size without path compression
 *   - O(n) memory for all of them
 *
***/

#include <stdio.h>
#include <bits/stdtr1c++.h>

using namespace std;

/***
 * Path compression and union by size
***/

struct DSU{
    vector<int> counter, parent;

    DSU(int n){
        parent.resize(n + 1);
        counter.resize(n + 1, 1);
        for (int i = 0; i <= n; i++) parent[i] = i;
    }

    int find_root(int i){
        while (i != parent[i]){
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        return parent[i];
    }

    bool connect(int a, int b){
        int c = find_root(a), d = find_root(b);
        if (c == d) return false;
        if (counter[c] > counter[d]) swap(c, d);
        parent[c] = d;
        counter[d] += counter[c], counter[c] = 0;
        return true;
    }

    bool is_connected(int a, int b){
        return find_root(a) == find_root(b);
    }

    int component_size(int i){
        int p = find_root(i);
        return counter[p];
    }
};

/***
 * Undoable unions, for offline dynamic connectivity and divide and conquer over time
 * time(): number of successful connects not yet rolled back
 * rollback(t): undo every union made after time() returned t, O(1) per undone union
 * offline_dynamic_connectivity.cpp carries a checked copy of this struct
***/

/// BEGIN SHARED rollback_dsu
struct RollbackDSU{
    vector<int> counter, parent, history;

    RollbackDSU(int n) : counter(n + 1, 1), parent(n + 1){
        iota(parent.begin(), parent.end(), 0);
    }

    int find_root(int i){
        while (i != parent[i]) i = parent[i];
        return i;
    }

    bool connect(int a, int b){
        a = find_root(a), b = find_root(b);
        if (a == b) return false;

        if (counter[a] > counter[b]) swap(a, b);
        parent[a] = b, counter[b] += counter[a];
        history.push_back(a);
        return true;
    }

    bool is_connected(int a, int b){
        return find_root(a) == find_root(b);
    }

    int component_size(int i){
        return counter[find_root(i)];
    }

    int time(){
        return history.size();
    }

    void rollback(int t){
        assert(0 <= t && t <= time());
        while ((int)history.size() > t){
            int a = history.back();
            history.pop_back();
            counter[parent[a]] -= counter[a], parent[a] = a;
        }
    }
};
/// END SHARED rollback_dsu

/***
 * Unknowns x[0..n] with relations x[b] - x[a] = d, over the integers (mod = 0) or modulo mod
 * connect(a, b, d): add the relation, false if it contradicts the earlier ones, which leaves the state unchanged
 * diff(a, b): x[b] - x[a] for connected a and b, reduced to [0, mod) when mod > 0
 * Parity (bipartiteness, same / different constraints) is mod = 2 with d = 1 for different
 * Requires 0 <= mod <= 1e18, and with mod = 0 every difference implied by the relations must fit in long long
***/

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

/***
 * The k-th call to connect happens at time k, now is the number of calls so far
 * find_root(i, t), is_connected(a, b, t), component_size(i, t): the state after the first t calls, t >= 0
 * Any t >= now gives the current state
 * The first time a and b were connected is a binary search over t on is_connected
***/

struct PartiallyPersistentDSU{
    int now = 0;
    vector<int> parent, joined;  /// joined[i]: time i stopped being a root, INT_MAX while it still is one
    vector<vector<pair<int, int>>> sizes;  /// (time, size) of root i, times increasing

    PartiallyPersistentDSU(int n) : parent(n + 1), joined(n + 1, INT_MAX), sizes(n + 1, {{0, 1}}){
        iota(parent.begin(), parent.end(), 0);
    }

    int find_root(int i, int t){
        assert(t >= 0);
        while (i != parent[i] && joined[i] <= t) i = parent[i];
        return i;
    }

    bool connect(int a, int b){
        now++;
        a = find_root(a, now), b = find_root(b, now);
        if (a == b) return false;

        if (sizes[a].back().second > sizes[b].back().second) swap(a, b);
        parent[a] = b, joined[a] = now;
        sizes[b].push_back({now, sizes[a].back().second + sizes[b].back().second});
        return true;
    }

    bool is_connected(int a, int b, int t){
        return find_root(a, t) == find_root(b, t);
    }

    int component_size(int i, int t){
        auto& history = sizes[find_root(i, t)];
        return prev(upper_bound(history.begin(), history.end(), make_pair(t, INT_MAX)))->second;
    }
};

int main(){
    auto dsu = DSU(8);
    dsu.connect(1, 2);
    dsu.connect(1, 3);
    dsu.connect(4, 5);
    dsu.connect(5, 6);
    dsu.connect(6, 7);

    assert(dsu.component_size(3) == 3);
    assert(dsu.component_size(5) == 4);
    assert(!dsu.is_connected(3, 5));

    assert(dsu.connect(2, 5));
    assert(!dsu.connect(7, 1));

    assert(dsu.component_size(2) == 7);
    assert(dsu.component_size(6) == 7);
    assert(dsu.is_connected(3, 5));

    /// Rollback
    auto rb = RollbackDSU(6);
    assert(rb.connect(1, 2));
    int checkpoint = rb.time();
    assert(rb.connect(3, 4) && rb.connect(2, 3) && !rb.connect(1, 4));

    assert(rb.time() == 3 && checkpoint == 1);
    assert(rb.component_size(4) == 4 && rb.is_connected(1, 4));

    rb.rollback(checkpoint);
    assert(rb.time() == 1);
    assert(rb.component_size(1) == 2 && rb.component_size(3) == 1 && rb.component_size(4) == 1);
    assert(!rb.is_connected(1, 3) && !rb.is_connected(3, 4));

    assert(rb.connect(4, 5) && rb.component_size(5) == 2);
    rb.rollback(0);
    for (int i = 0; i <= 6; i++) assert(rb.component_size(i) == 1);

    /// Weighted, AtCoder ABC087 D samples: x[R] - x[L] = D
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

    /// Partially persistent
    auto pp = PartiallyPersistentDSU(5);
    assert(pp.connect(1, 2) && pp.connect(3, 4) && !pp.connect(2, 1) && pp.connect(2, 4));
    assert(pp.now == 4);

    assert(!pp.is_connected(1, 2, 0) && pp.is_connected(1, 2, 1));
    assert(!pp.is_connected(1, 4, 3) && pp.is_connected(1, 4, 4) && pp.is_connected(3, 2, 100));
    assert(pp.find_root(1, 0) == 1 && pp.find_root(5, 4) == 5);

    assert(pp.component_size(1, 0) == 1 && pp.component_size(1, 1) == 2 && pp.component_size(4, 1) == 1);
    assert(pp.component_size(4, 3) == 2 && pp.component_size(3, 4) == 4 && pp.component_size(1, 100) == 4);
    assert(pp.component_size(5, 4) == 1);
    assert(pp.is_connected(3, 2, INT_MAX) && pp.component_size(1, INT_MAX) == 4 && !pp.is_connected(1, 5, INT_MAX));

    return 0;
}
