#include "../common.h"

#define main library_main
#include "../../code_library/linear_algebra/determinant.cpp"
#undef main

/// Exact determinant by Laplace expansion along the first row, in __int128
__int128 laplace(const vector<vector<long long>>& a){
    int n = a.size();
    if (n == 0) return 1;
    __int128 res = 0;

    for (int c = 0; c < n; c++){
        vector<vector<long long>> minor;
        for (int r = 1; r < n; r++){
            vector<long long> row;
            for (int k = 0; k < n; k++){
                if (k != c) row.push_back(a[r][k]);
            }
            minor.push_back(row);
        }
        __int128 sub = laplace(minor) * a[0][c];
        res += c % 2 ? -sub : sub;
    }

    return res;
}

/// Exact determinant by fraction free Bareiss elimination, fine while the minors stay small
__int128 bareiss(vector<vector<__int128>> a){
    int n = a.size(), sign = 1;
    __int128 prev = 1;

    for (int i = 0; i < n; i++){
        if (a[i][i] == 0){
            int p = i + 1;
            while (p < n && a[p][i] == 0) p++;
            if (p == n) return 0;
            swap(a[i], a[p]), sign = -sign;
        }
        for (int j = i + 1; j < n; j++){
            for (int k = i + 1; k < n; k++) a[j][k] = (a[j][k] * a[i][i] - a[j][i] * a[i][k]) / prev;
        }
        prev = a[i][i];
    }

    return sign * (n ? a[n - 1][n - 1] : 1);
}

long long reduce(__int128 x, long long m){
    x %= m;
    return (long long)(x < 0 ? x + m : x);
}

/// Spanning trees by trying every (n - 1)-edge subset and checking it is acyclic with a union find
long long brute_spanning_trees(int n, const vector<pair<int, int>>& edges){
    int e = edges.size();
    long long total = 0;

    for (int mask = 0; mask < (1 << e); mask++){
        if (__builtin_popcount(mask) != n - 1) continue;
        vector<int> par(n);
        iota(par.begin(), par.end(), 0);
        function<int(int)> find = [&](int x){ return par[x] == x ? x : par[x] = find(par[x]); };
        bool ok = true;
        for (int i = 0; i < e && ok; i++){
            if (!(mask >> i & 1)) continue;
            int a = find(edges[i].first), b = find(edges[i].second);
            if (a == b) ok = false;
            par[a] = b;
        }
        total += ok;
    }

    return total;
}

/// Arborescences by trying every (n - 1)-edge subset: one parent per non-root node, and every parent chain reaches root
long long brute_arborescences(int n, const vector<pair<int, int>>& edges, int root){
    int e = edges.size();
    long long total = 0;

    for (int mask = 0; mask < (1 << e); mask++){
        if (__builtin_popcount(mask) != n - 1) continue;
        vector<int> par(n, -1), indeg(n, 0);
        for (int i = 0; i < e; i++){
            if (mask >> i & 1) indeg[edges[i].second]++, par[edges[i].second] = edges[i].first;
        }
        bool ok = indeg[root] == 0;
        for (int v = 0; v < n && ok; v++){
            if (v != root && indeg[v] != 1) ok = false;
        }
        for (int v = 0; v < n && ok; v++){
            int x = v, steps = 0;
            while (x != root && steps <= n) x = par[x], steps++;
            if (x != root) ok = false;
        }
        total += ok;
    }

    return total;
}

long long power_mod(long long b, long long e, long long m){
    long long r = 1 % m;
    for (b %= m; e; e >>= 1, b = (__int128)b * b % m){
        if (e & 1) r = (__int128)r * b % m;
    }
    return r;
}

long long pick_modulus(long long it){
    if (it % 3 == 0) return stress::rand_int(1, 30);
    if (it % 3 == 1) return (1LL << 62) - 1;
    return stress::rand_int(1, (1LL << 62) - 1);
}

/// Every multiset of at most max_edges ordered pairs on n <= 3 nodes, every root, so the tiny corner cases are covered deterministically
void exhaustive_matrix_tree(int max_edges){
    for (int n = 1; n <= 3; n++){
        vector<pair<int, int>> edges;
        function<void(int)> extend = [&](int from){
            SpanningTrees g(n);
            Arborescences d(n);
            for (auto [u, v] : edges) g.add_edge(u, v), d.add_edge(u, v);
            for (long long m : {(1LL << 62) - 1, 2LL}){
                assert(g.count(m) == brute_spanning_trees(n, edges) % m);
                for (int root = 0; root < n; root++) assert(d.count(root, m) == brute_arborescences(n, edges, root) % m);
            }

            if ((int)edges.size() == max_edges) return;
            for (int p = from; p < n * n; p++){
                edges.push_back({p / n, p % n});
                extend(p);
                edges.pop_back();
            }
        };
        extend(0);
    }
}

void stress_matrix_tree(){
    const long long BIG = (1LL << 62) - 1;
    exhaustive_matrix_tree(4);

    for (long long it = 0; it < stress::scaled(6000); it++){
        int n = stress::rand_int(1, 6), e = stress::rand_int(0, 12);
        vector<pair<int, int>> edges;
        for (int i = 0; i < e; i++){
            int u = stress::rand_int(0, n - 1), v = stress::rand_int(0, n - 1);
            if (it % 5 == 0 && stress::rand_int(0, 3) == 0) v = u;
            edges.push_back({u, v});
        }

        SpanningTrees g(n);
        Arborescences d(n);
        for (auto [u, v] : edges) g.add_edge(u, v), d.add_edge(u, v);
        long long m = pick_modulus(it);
        assert(g.count(m) == brute_spanning_trees(n, edges) % m);
        int root = stress::rand_int(0, n - 1);
        assert(d.count(root, m) == brute_arborescences(n, edges, root) % m);
    }

    for (long long it = 0; it < stress::scaled(300); it++){
        int n = stress::rand_int(1, 40), e = stress::rand_int(0, 3 * n);
        SpanningTrees g(n);
        Arborescences d(n);
        for (int i = 0; i < e; i++){
            int u = stress::rand_int(0, n - 1), v = stress::rand_int(0, n - 1);
            g.add_edge(u, v), d.add_edge(u, v), d.add_edge(v, u);
        }
        long long m = pick_modulus(it);
        assert(g.count(m) == d.count(stress::rand_int(0, n - 1), m));
    }

    for (int n = 1; n <= 50; n++){
        SpanningTrees g(n);
        Arborescences d(n);
        for (int u = 0; u < n; u++){
            for (int v = 0; v < n; v++){
                if (u < v) g.add_edge(u, v);
                if (u != v) d.add_edge(u, v);
            }
        }
        long long cayley = n == 1 ? 1 : power_mod(n, n - 2, BIG);
        assert(g.count(1000000007) == (n == 1 ? 1 : power_mod(n, n - 2, 1000000007)));
        assert(d.count(n - 1, 998244353) == (n == 1 ? 1 : power_mod(n, n - 2, 998244353)));
        if (n <= 17) assert(g.count(BIG) == cayley && d.count(0, BIG) == cayley);
    }

    for (int a = 1; a <= 12; a++){
        for (int b = 1; b <= 12; b++){
            SpanningTrees g(a + b);
            for (int u = 0; u < a; u++){
                for (int v = 0; v < b; v++) g.add_edge(u, a + v);
            }
            long long m = 1000000007;
            assert(g.count(m) == power_mod(a, b - 1, m) * power_mod(b, a - 1, m) % m);
        }
    }
}

int main(){
    for (long long it = 0; it < stress::scaled(20000); it++){
        int n = stress::rand_int(0, 5);
        long long range = it % 3 == 0 ? 3 : (it % 3 == 1 ? 1000 : 1000000);
        vector<vector<long long>> a(n, vector<long long>(n));
        for (auto& row : a){
            for (auto& x : row) x = stress::rand_int(-range, range);
        }
        if (n >= 2 && it % 7 == 0) a[1] = a[0];
        __int128 exact = laplace(a);

        long long m = it % 4 == 0 ? stress::rand_int(1, 30) : (it % 4 == 1 ? (1LL << 62) - 1 : stress::rand_int(1, (1LL << 62) - 1));
        assert(determinant_mod(a, m) == reduce(exact, m));
        if (exact >= -2258704744122758558LL && exact <= 2258704744122758558LL) assert(determinant(a) == (long long)exact);
    }

    for (long long it = 0; it < stress::scaled(300); it++){
        int n = stress::rand_int(6, 18);
        vector<vector<long long>> a(n, vector<long long>(n));
        vector<vector<__int128>> b(n, vector<__int128>(n));
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++) a[i][j] = b[i][j] = stress::rand_int(-2, 2);
        }

        __int128 exact = bareiss(b);
        assert(determinant(a) == (long long)exact);
        long long m = stress::rand_int(1, 1000000007);
        assert(determinant_mod(a, m) == reduce(exact, m));
    }

    const long long HALF = 2258704744122758558LL;
    for (long long it = 0; it < stress::scaled(2000); it++){
        long long v = it < 6 ? vector<long long>{HALF, -HALF, HALF - 1, 1 - HALF, 0, 1}[it] : stress::rand_int(-HALF, HALF);
        assert(determinant({{v}}) == v);
        assert(determinant({{v, 0}, {0, 1}}) == v && determinant({{0, v}, {1, 0}}) == -v);
    }

    stress_matrix_tree();

    return 0;
}
