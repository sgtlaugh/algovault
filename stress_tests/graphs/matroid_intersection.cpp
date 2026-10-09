#include "../common.h"

#define main library_main
#include "../../code_library/graphs/matroid_intersection.cpp"
#undef main

using ull = unsigned long long;

enum Kind {GRAPHIC, PARTITION, XOR};

struct Spec{
    Kind kind;
    int v = 0;
    vector<array<int, 2>> edges;
    vector<int> color, cap;
    vector<ull> vec;
};

/// Plain Gaussian elimination: independent iff no vector reduces to 0 against the earlier ones
bool xor_independent(const vector<ull>& values){
    ull basis[64] = {};
    for (ull a : values){
        for (int b = 63; b >= 0 && a; b--){
            if (!(a >> b & 1)) continue;
            if (!basis[b]){
                basis[b] = a;
                break;
            }
            a ^= basis[b];
        }
        if (a == 0) return false;
    }
    return true;
}

/// Independence from the definitions: no cycle, no color over its cap, no vector in the span of the others
bool independent(const Spec& s, const vector<int>& subset){
    if (s.kind == GRAPHIC){
        vector<int> root(s.v);
        iota(root.begin(), root.end(), 0);
        function<int(int)> find = [&](int a){ return root[a] == a ? a : root[a] = find(root[a]); };
        for (int e : subset){
            int a = find(s.edges[e][0]), b = find(s.edges[e][1]);
            if (a == b) return false;
            root[a] = b;
        }
        return true;
    }
    if (s.kind == PARTITION){
        vector<int> used(s.cap.size(), 0);
        for (int x : subset){
            if (++used[s.color[x]] > s.cap[s.color[x]]) return false;
        }
        return true;
    }

    vector<ull> values;
    for (int x : subset) values.push_back(s.vec[x]);
    return xor_independent(values);
}

ull random_vector(int bits){
    if (stress::rand_int(0, 3) == 0) return (ull)stress::rand_int(0, 7) << 61;
    return stress::rand_int(0, (1LL << bits) - 1);
}

Spec random_spec(int n){
    Spec s;
    s.kind = (Kind)stress::rand_int(0, 2);
    if (s.kind == GRAPHIC){
        s.v = stress::rand_int(1, 7);
        for (int i = 0; i < n; i++) s.edges.push_back({(int)stress::rand_int(0, s.v - 1), (int)stress::rand_int(0, s.v - 1)});
    }
    else if (s.kind == PARTITION){
        int colors = stress::rand_int(1, 5);
        for (int i = 0; i < n; i++) s.color.push_back(stress::rand_int(0, colors - 1));
        for (int c = 0; c < colors; c++) s.cap.push_back(stress::rand_int(0, 3));
    }
    else{
        int bits = stress::rand_int(1, 5);
        for (int i = 0; i < n; i++) s.vec.push_back(random_vector(bits));
    }
    return s;
}

template<typename F>
auto with_oracle(const Spec& s, F solve){
    if (s.kind == GRAPHIC){
        GraphicMatroid g(s.v);
        for (auto [a, b] : s.edges) g.add_edge(a, b);
        return solve(g);
    }
    if (s.kind == PARTITION){
        PartitionMatroid p(s.color, s.cap);
        return solve(p);
    }
    XorMatroid x(s.vec);
    return solve(x);
}

vector<int> run(int n, const Spec& a, const Spec& b){
    vector<int> res = with_oracle(a, [&](auto& m1){
        return with_oracle(b, [&](auto& m2){ return matroid_intersection(n, m1, m2); });
    });

    for (int i = 0; i < (int)res.size(); i++){
        assert(0 <= res[i] && res[i] < n);
        assert(i == 0 || res[i - 1] < res[i]);
    }
    assert(independent(a, res) && independent(b, res));
    return res;
}

/// The driver never asks can_exchange(y, x) when can_add(x) holds, so the oracle contract is checked on its own
void check_oracle(int n, const Spec& s){
    vector<int> order(n), chosen;
    iota(order.begin(), order.end(), 0);
    shuffle(order.begin(), order.end(), stress::rng());
    int limit = stress::rand_int(0, n);
    for (int x : order){
        chosen.push_back(x);
        if ((int)chosen.size() > limit || !independent(s, chosen)) chosen.pop_back();
    }
    sort(chosen.begin(), chosen.end());

    with_oracle(s, [&](auto& m){
        m.build(chosen);
        for (int x = 0; x < n; x++){
            if (binary_search(chosen.begin(), chosen.end(), x)) continue;
            vector<int> added = chosen;
            added.push_back(x);
            assert(m.can_add(x) == independent(s, added));
            for (int i = 0; i < (int)chosen.size(); i++){
                vector<int> swapped = added;
                swapped.erase(swapped.begin() + i);
                assert(m.can_exchange(chosen[i], x) == independent(s, swapped));
            }
        }
    });
}

int brute(int n, const Spec& a, const Spec& b){
    int best = 0;
    for (int mask = 0; mask < (1 << n); mask++){
        if (__builtin_popcount(mask) <= best) continue;
        vector<int> subset;
        for (int i = 0; i < n; i++){
            if (mask >> i & 1) subset.push_back(i);
        }
        if (independent(a, subset) && independent(b, subset)) best = subset.size();
    }
    return best;
}

/// Encodes "at most one element per label" as each matroid kind, so bipartite matching has a Kuhn reference at scale
Spec label_spec(Kind kind, const vector<int>& label, int labels){
    Spec s;
    s.kind = kind;
    if (kind == PARTITION){
        s.color = label, s.cap.assign(labels, 1);
        return s;
    }
    if (kind == GRAPHIC){
        /// Label l is edge l of a random spanning tree on labels + 1 vertices, a path half of the time for deep trees
        bool path = stress::rand_int(0, 1);
        vector<array<int, 2>> tree;
        for (int l = 0; l < labels; l++) tree.push_back({l + 1, path ? l : (int)stress::rand_int(0, l)});
        s.v = labels + 1;
        for (int l : label) s.edges.push_back(tree[l]);
        return s;
    }

    vector<ull> basis;
    while ((int)basis.size() < labels){
        basis.push_back(stress::rng()());
        if (!xor_independent(basis)) basis.pop_back();
    }
    for (int l : label) s.vec.push_back(basis[l]);
    return s;
}

int kuhn(int n_left, int n_right, const vector<array<int, 2>>& edges){
    vector<vector<int>> adj(n_left);
    for (auto [u, w] : edges) adj[u].push_back(w);
    vector<int> match(n_right, -1), seen;
    function<bool(int)> try_kuhn = [&](int u){
        for (int w : adj[u]){
            if (seen[w]) continue;
            seen[w] = 1;
            if (match[w] == -1 || try_kuhn(match[w])){
                match[w] = u;
                return true;
            }
        }
        return false;
    };

    int res = 0;
    for (int u = 0; u < n_left; u++){
        seen.assign(n_right, 0);
        res += try_kuhn(u);
    }
    return res;
}

void check_matching(int n_left, int n_right, const vector<array<int, 2>>& edges){
    vector<int> left, right;
    for (auto [u, w] : edges) left.push_back(u), right.push_back(w);
    Kind kl = (Kind)stress::rand_int(0, n_left <= 64 ? 2 : 1), kr = (Kind)stress::rand_int(0, n_right <= 64 ? 2 : 1);
    Spec a = label_spec(kl, left, n_left), b = label_spec(kr, right, n_right);
    assert((int)run(edges.size(), a, b).size() == kuhn(n_left, n_right, edges));
}

int main(){
    for (long long it = 0; it < stress::scaled(6000); it++){
        int n = it < 13 ? it : stress::rand_int(0, 12);
        Spec a = random_spec(n), b = random_spec(n);
        assert((int)run(n, a, b).size() == brute(n, a, b));
    }

    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(0, 14);
        check_oracle(n, random_spec(n));
    }

    for (long long it = 0; it < stress::scaled(300); it++){
        int n_left = stress::rand_int(1, 70), n_right = stress::rand_int(1, 70), m = stress::rand_int(0, 250);
        vector<array<int, 2>> edges;
        for (int i = 0; i < m; i++) edges.push_back({(int)stress::rand_int(0, n_left - 1), (int)stress::rand_int(0, n_right - 1)});
        check_matching(n_left, n_right, edges);
    }

    /// Staircase: left i -> right i and i + 1 with the greedy-friendly edge listed second, forcing long augmenting paths
    for (int n : {64, 200}){
        vector<array<int, 2>> stairs;
        for (int i = 0; i < n; i++){
            if (i + 1 < n) stairs.push_back({i, i + 1});
            stairs.push_back({i, i});
        }
        check_matching(n, n, stairs);
        check_matching(n, n, stairs);
    }

    /// Rainbow forests at scale: no fast reference, so the two oracle orders must agree on the size
    for (long long it = 0; it < stress::scaled(20); it++){
        int v = stress::rand_int(2, 150), m = stress::rand_int(1, 500), colors = stress::rand_int(1, 200);
        Spec g{GRAPHIC, v, {}, {}, {}, {}}, p{PARTITION, 0, {}, {}, {}, {}};
        for (int i = 0; i < m; i++){
            g.edges.push_back({(int)stress::rand_int(0, v - 1), (int)stress::rand_int(0, v - 1)});
            p.color.push_back(stress::rand_int(0, colors - 1));
        }
        for (int c = 0; c < colors; c++) p.cap.push_back(stress::rand_int(0, 2));
        assert(run(m, g, p).size() == run(m, p, g).size());
    }

    return 0;
}
