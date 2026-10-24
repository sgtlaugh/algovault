#include "../common.h"

#define main library_main
#include "../../code_library/graphs/weighted_matroid_intersection.cpp"
#undef main

#define main hungarian_main
#include "../../code_library/graphs/hungarian_algorithm.cpp"
#undef main

using ull = unsigned long long;

enum Kind{GRAPHIC, PARTITION, XOR};

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

long long total(const vector<long long>& weight, const vector<int>& subset){
    long long sum = 0;
    for (int x : subset) sum += weight[x];
    return sum;
}

/// Checks the answer is a sorted common independent set and returns its weight
long long run(int n, const vector<long long>& weight, const Spec& a, const Spec& b){
    vector<int> res = with_oracle(a, [&](auto& m1){
        return with_oracle(b, [&](auto& m2){ return weighted_matroid_intersection(n, weight, m1, m2); });
    });

    for (int i = 0; i < (int)res.size(); i++){
        assert(0 <= res[i] && res[i] < n);
        assert(i == 0 || res[i - 1] < res[i]);
    }
    assert(independent(a, res) && independent(b, res));
    return total(weight, res);
}

/// The max weight over every common independent set, and the max weight over the largest ones
pair<long long, long long> brute(int n, const vector<long long>& weight, const Spec& a, const Spec& b){
    long long best = 0, best_largest = 0;
    int largest = 0;
    for (int mask = 0; mask < (1 << n); mask++){
        vector<int> subset;
        for (int i = 0; i < n; i++){
            if (mask >> i & 1) subset.push_back(i);
        }
        if (!independent(a, subset) || !independent(b, subset)) continue;

        long long w = total(weight, subset);
        best = max(best, w);
        if ((int)subset.size() > largest) largest = subset.size(), best_largest = w;
        else if ((int)subset.size() == largest) best_largest = max(best_largest, w);
    }
    return {best, best_largest};
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

/// Max weight bipartite matching on a complete bipartite graph with weights >= 0 is the max cost assignment
void check_assignment(int n_left, int n_right, long long max_weight){
    vector<vector<long long>> mat(n_left, vector<long long>(n_right));
    vector<int> left, right;
    vector<long long> weight;
    for (int i = 0; i < n_left; i++){
        for (int j = 0; j < n_right; j++){
            mat[i][j] = stress::rand_int(0, max_weight);
            left.push_back(i), right.push_back(j), weight.push_back(mat[i][j]);
        }
    }

    Kind kl = (Kind)stress::rand_int(0, 2), kr = (Kind)stress::rand_int(0, 2);
    Spec a = label_spec(kl, left, n_left), b = label_spec(kr, right, n_right);
    assert(run(weight.size(), weight, a, b) == hungarian(mat, false).first);
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

/// Weight 1 on every edge of a sparse bipartite graph: the answer is a maximum matching, found by Kuhn
void check_matching(int n_left, int n_right, const vector<array<int, 2>>& edges){
    vector<int> left, right;
    for (auto [u, w] : edges) left.push_back(u), right.push_back(w);
    Kind kl = (Kind)stress::rand_int(0, n_left <= 64 ? 2 : 1), kr = (Kind)stress::rand_int(0, n_right <= 64 ? 2 : 1);
    Spec a = label_spec(kl, left, n_left), b = label_spec(kr, right, n_right);
    assert(run(edges.size(), vector<long long>(edges.size(), 1), a, b) == kuhn(n_left, n_right, edges));
}

int main(){
    for (long long it = 0; it < stress::scaled(6000); it++){
        int n = it < 11 ? it : stress::rand_int(0, 10);
        long long max_weight = it % 3 == 0 ? 1000000000000LL : 6;
        vector<long long> weight(n);
        for (auto& w : weight) w = stress::rand_int(-max_weight, max_weight);
        /// Equal weights make every path equally cheap, so only the fewest-elements tie-break keeps chosen independent
        if (it % 3 == 1) fill(weight.begin(), weight.end(), stress::rand_int(1, 3));
        Spec a = random_spec(n), b = random_spec(n);
        auto [best, best_largest] = brute(n, weight, a, b);
        assert(run(n, weight, a, b) == best);

        /// The header's recipe for the best among the largest sets
        long long big = 1;
        for (long long w : weight) big += 2 * abs(w);
        vector<long long> shifted = weight;
        for (auto& w : shifted) w += big;
        vector<int> res = with_oracle(a, [&](auto& m1){
            return with_oracle(b, [&](auto& m2){ return weighted_matroid_intersection(n, shifted, m1, m2); });
        });
        assert(total(weight, res) == best_largest);
    }

    /// Small weights tie often, where a path that is cheapest but not shortest breaks independence
    for (long long it = 0; it < stress::scaled(150); it++){
        check_assignment(stress::rand_int(1, 12), stress::rand_int(1, 12), it % 2 ? 3 : 1000000000);
    }

    for (long long it = 0; it < stress::scaled(100); it++){
        int n_left = stress::rand_int(1, 40), n_right = stress::rand_int(1, 40), m = stress::rand_int(0, 120);
        vector<array<int, 2>> edges;
        for (int i = 0; i < m; i++) edges.push_back({(int)stress::rand_int(0, n_left - 1), (int)stress::rand_int(0, n_right - 1)});
        check_matching(n_left, n_right, edges);
    }

    /// Staircase: left i -> right i and i + 1, the greedy-friendly edge listed second, forcing long augmenting paths
    for (int n : {16, 64}){
        vector<array<int, 2>> stairs;
        for (int i = 0; i < n; i++){
            if (i + 1 < n) stairs.push_back({i, i + 1});
            stairs.push_back({i, i});
        }
        check_matching(n, n, stairs);
    }

    /// Weighted rainbow forests at scale: no fast reference, so the two oracle orders must agree on the weight
    for (long long it = 0; it < stress::scaled(10); it++){
        int v = stress::rand_int(2, 60), m = stress::rand_int(1, 150), colors = stress::rand_int(1, 60);
        Spec g{GRAPHIC, v, {}, {}, {}, {}}, p{PARTITION, 0, {}, {}, {}, {}};
        vector<long long> weight(m);
        for (int i = 0; i < m; i++){
            g.edges.push_back({(int)stress::rand_int(0, v - 1), (int)stress::rand_int(0, v - 1)});
            p.color.push_back(stress::rand_int(0, colors - 1));
            weight[i] = stress::rand_int(-5, 20);
        }
        for (int c = 0; c < colors; c++) p.cap.push_back(stress::rand_int(0, 2));
        assert(run(m, weight, g, p) == run(m, weight, p, g));
    }

    return 0;
}
