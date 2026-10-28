/***
 *
 * DSU with Potentials
 * Union-find that also knows how much each value differs from the others in its set
 *
 * Complexity: amortized O(α(n)) per operation
 *
 * PotentialDSU dsu(n, mod = 0): nodes 0..n, all values unknown
 * dsu.relate(a, b, d): record x[b] = x[a] + d, false if that contradicts what is known, which leaves the state unchanged
 * dsu.diff(a, b): x[b] - x[a], a and b must be in the same set
 * dsu.same(a, b): true if a and b are in the same set, so diff(a, b) is known
 * dsu.size(a): number of nodes in a's set
 *
 * mod = 0 works on plain integers, every implied difference must fit in long long
 * mod > 0 works modulo mod, mod <= 1e18
 * mod = 2 is parity, relate(a, b, 1) means a and b differ
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
    int n;

    scanf("%d", &n);

    map<long long, int> rows, cols;
    auto id = [&](map<long long, int>& ids, long long key){
        auto it = ids.find(key);
        if (it != ids.end()) return it->second;
        int next = rows.size() + cols.size();
        return ids[key] = next;
    };

    PotentialDSU dsu(2 * n, 2);
    bool ok = true;
    while (n--){
        long long x, y;
        int a, b;
        scanf("%lld %lld %d %d", &x, &y, &a, &b);
        int r = id(rows, x), c = id(cols, y);
        if (ok) ok = dsu.relate(r, c, a != b);
        puts(ok ? "Yes" : "No");
    }

    return 0;
}
