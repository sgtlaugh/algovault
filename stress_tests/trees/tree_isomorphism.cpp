#include "../common.h"
#include <sys/wait.h>
#include <unistd.h>

#define main library_main
#include "../../code_library/trees/tree_isomorphism.cpp"
#undef main

using Edges = vector<pair<int, int>>;

IsoTree build(int n, const Edges& edges){
    IsoTree t(n);
    for (auto [u, v] : edges) t.add_edge(u, v);
    return t;
}

Edges prufer_decode(int n, const vector<int>& seq){
    if (n == 1) return {};
    vector<int> deg(n, 1);
    for (int x : seq) deg[x]++;

    Edges edges;
    for (int x : seq){
        int leaf = 0;
        while (deg[leaf] != 1) leaf++;
        edges.push_back({leaf, x});
        deg[leaf]--, deg[x]--;
    }

    vector<int> last;
    for (int v = 0; v < n; v++){
        if (deg[v] == 1) last.push_back(v);
    }
    edges.push_back({last[0], last[1]});
    return edges;
}

/// Random tree of a given shape, labels not shuffled
Edges random_edges(int n, int shape){
    Edges edges;
    for (int i = 1; i < n; i++){
        int p;
        if (shape == 0) p = stress::rand_int(0, i - 1);
        else if (shape == 1) p = i - 1;
        else if (shape == 2) p = 0;
        else if (shape == 3) p = i % 2 ? i - 1 : max(0, i - 2);
        else p = stress::rand_int(max(0, i - 3), i - 1);
        edges.push_back({i, p});
    }
    return edges;
}

/// Applies the permutation, flips edge directions and shuffles the edge order
Edges relabel(const Edges& edges, const vector<int>& perm){
    Edges res;
    for (auto [u, v] : edges){
        if (stress::rand_int(0, 1)) swap(u, v);
        res.push_back({perm[u], perm[v]});
    }
    shuffle(res.begin(), res.end(), stress::rng());
    return res;
}

vector<int> random_perm(int n){
    vector<int> perm(n);
    iota(perm.begin(), perm.end(), 0);
    shuffle(perm.begin(), perm.end(), stress::rng());
    return perm;
}

/// Brute force: some relabeling maps a's edge set onto b's (and ra onto rb when rooted)
bool brute_isomorphic(int n, const Edges& a, const Edges& b, int ra = -1, int rb = -1){
    vector<int> mask(n, 0), perm(n);
    for (auto [u, v] : b) mask[u] |= 1 << v, mask[v] |= 1 << u;
    iota(perm.begin(), perm.end(), 0);

    do{
        if (ra != -1 && perm[ra] != rb) continue;
        bool ok = true;
        for (auto [u, v] : a){
            if (!(mask[perm[u]] >> perm[v] & 1)){
                ok = false;
                break;
            }
        }
        if (ok) return true;
    } while (next_permutation(perm.begin(), perm.end()));
    return false;
}

/// Every labeled tree on n nodes via Prufer sequences: distinct unrooted ids must number OEIS A000055(n), rooted ids A000081(n)
void check_exhaustive(int n){
    const int unrooted_count[] = {0, 1, 1, 1, 2, 3, 6, 11};
    const int rooted_count[] = {0, 1, 1, 2, 4, 9, 20, 48};
    TreeCanonizer canon;
    set<int> unrooted, rooted;
    vector<int> seq(max(0, n - 2), 0);

    while (true){
        IsoTree t = build(n, prufer_decode(n, seq));
        unrooted.insert(canon.unrooted_id(t));
        for (int r = 0; r < n; r++) rooted.insert(canon.rooted_id(t, r));

        int i = 0;
        while (i < (int)seq.size() && ++seq[i] == n) seq[i++] = 0;
        if (i == (int)seq.size()) break;
    }

    assert((int)unrooted.size() == unrooted_count[n]);
    assert((int)rooted.size() == rooted_count[n]);
}

/// Pairs on n <= 8 nodes against the all relabelings brute force, half of them relabeled copies so both answers occur
void check_small_pair(int n){
    Edges a = random_edges(n, stress::rand_int(0, 4));
    a = relabel(a, random_perm(n));
    Edges b = stress::rand_int(0, 1) ? relabel(a, random_perm(n)) : relabel(random_edges(n, stress::rand_int(0, 4)), random_perm(n));
    IsoTree ta = build(n, a), tb = build(n, b);
    TreeCanonizer canon;

    assert((canon.unrooted_id(ta) == canon.unrooted_id(tb)) == brute_isomorphic(n, a, b));
    int ra = stress::rand_int(0, n - 1), rb = stress::rand_int(0, n - 1);
    assert((canon.rooted_id(ta, ra) == canon.rooted_id(tb, rb)) == brute_isomorphic(n, a, b, ra, rb));
}

/// Reference canonical form: "(" + sorted child strings + ")", and centers as the nodes of minimum eccentricity
struct NaiveTree{
    int n;
    vector<vector<int>> adj;

    NaiveTree(int n, const Edges& edges) : n(n), adj(n){
        for (auto [u, v] : edges) adj[u].push_back(v), adj[v].push_back(u);
    }

    vector<int> dist_from(int s) const{
        vector<int> dist(n, -1), queue = {s};
        dist[s] = 0;
        for (int i = 0; i < (int)queue.size(); i++){
            for (int v : adj[queue[i]]){
                if (dist[v] == -1) dist[v] = dist[queue[i]] + 1, queue.push_back(v);
            }
        }
        return dist;
    }

    vector<string> subtree_strings(int root) const{
        vector<int> dist = dist_from(root), order(n);
        iota(order.begin(), order.end(), 0);
        sort(order.begin(), order.end(), [&](int x, int y){ return dist[x] > dist[y]; });

        vector<string> str(n);
        for (int v : order){
            vector<string> kids;
            for (int c : adj[v]){
                if (dist[c] == dist[v] + 1) kids.push_back(str[c]);
            }
            sort(kids.begin(), kids.end());
            str[v] = "(";
            for (auto& k : kids) str[v] += k;
            str[v] += ")";
        }
        return str;
    }

    vector<int> centers() const{
        vector<int> ecc(n), res;
        for (int v = 0; v < n; v++){
            vector<int> d = dist_from(v);
            ecc[v] = *max_element(d.begin(), d.end());
        }
        int best = *min_element(ecc.begin(), ecc.end());
        for (int v = 0; v < n; v++){
            if (ecc[v] == best) res.push_back(v);
        }
        return res;
    }

    string unrooted_string() const{
        string best;
        for (int c : centers()){
            string s = subtree_strings(c)[c];
            if (best.empty() || s < best) best = s;
        }
        return best;
    }
};

/// Ids and canonical strings must be in bijection across every subtree of a batch of trees sharing one canonizer
void check_batch(int count, int max_n){
    TreeCanonizer canon;
    map<string, int> rooted_by_string, unrooted_by_string;
    map<int, string> string_by_rooted, string_by_unrooted;
    auto bind = [](map<string, int>& to_id, map<int, string>& to_string, const string& s, int id){
        assert(to_id.try_emplace(s, id).first->second == id);
        assert(to_string.try_emplace(id, s).first->second == s);
    };

    Edges base = random_edges(max_n, stress::rand_int(0, 4));
    for (int it = 0; it < count; it++){
        int n = stress::rand_int(1, max_n);
        Edges edges = stress::rand_int(0, 2) ? random_edges(n, stress::rand_int(0, 4)) : Edges(base.begin(), base.begin() + n - 1);
        if (n > 2 && stress::rand_int(0, 1)){
            auto& [child, parent] = edges[stress::rand_int(0, n - 2)];
            parent = stress::rand_int(0, child - 1);
        }
        edges = relabel(edges, random_perm(n));

        IsoTree t = build(n, edges);
        NaiveTree naive(n, edges);
        assert(t.centers() == naive.centers());
        bind(unrooted_by_string, string_by_unrooted, naive.unrooted_string(), canon.unrooted_id(t));

        int root = stress::rand_int(0, n - 1);
        vector<int> ids = canon.subtree_ids(t, root);
        vector<string> strs = naive.subtree_strings(root);
        for (int v = 0; v < n; v++) bind(rooted_by_string, string_by_rooted, strs[v], ids[v]);
    }
}

/// A relabeled copy keeps every id, the same tree with one leaf moved changes it
void check_large(int n, int shape){
    Edges a = random_edges(n, shape);
    vector<int> perm = random_perm(n);
    IsoTree ta = build(n, a), tb = build(n, relabel(a, perm));
    TreeCanonizer canon;

    assert(canon.unrooted_id(ta) == canon.unrooted_id(tb));
    int r = stress::rand_int(0, n - 1);
    assert(canon.rooted_id(ta, r) == canon.rooted_id(tb, perm[r]));
    vector<int> ca = ta.centers(), cb = tb.centers();
    for (int& c : ca) c = perm[c];
    sort(ca.begin(), ca.end());
    assert(ca == cb);
}

/// Path against the path with its last node hung from node 1, centers of the path are its middle
void check_path(int n){
    Edges path = random_edges(n, 1), moved = path;
    moved.back().second = 1;
    vector<int> perm = random_perm(n);
    IsoTree tp = build(n, relabel(path, perm)), tm = build(n, relabel(moved, random_perm(n)));
    TreeCanonizer canon;

    assert(canon.unrooted_id(tp) != canon.unrooted_id(tm));
    assert(canon.rooted_id(tp, perm[0]) == canon.rooted_id(tp, perm[n - 1]));
    assert(canon.rooted_id(tp, perm[0]) != canon.rooted_id(tp, perm[n / 2]));
    vector<int> expected = {perm[(n - 1) / 2], perm[n / 2]};
    sort(expected.begin(), expected.end());
    expected.erase(unique(expected.begin(), expected.end()), expected.end());
    assert(tp.centers() == expected);
}

/// Every call must abort on a non-tree, the alarm turns an endless traversal into a failure instead of a hang
void check_rejects(int n, const Edges& edges, int call){
    pid_t pid = fork();
    assert(pid >= 0);
    if (pid == 0){
        alarm(5);
        assert(freopen("/dev/null", "w", stderr));
        IsoTree t = build(n, edges);
        TreeCanonizer canon;
        if (call == 0) t.centers();
        else if (call == 1) canon.rooted_id(t, 0);
        else canon.unrooted_id(t);
        _exit(0);
    }

    int status;
    assert(waitpid(pid, &status, 0) == pid);
    assert(WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT);
}

int main(){
    for (int call = 0; call < 3; call++){
        check_rejects(4, {{0, 1}, {1, 2}, {2, 0}}, call);
        check_rejects(5, {{0, 1}, {1, 2}, {2, 3}, {3, 1}}, call);
        check_rejects(3, {{0, 1}, {1, 2}, {2, 0}}, call);
        check_rejects(3, {{0, 1}, {0, 1}, {1, 2}}, call);
        check_rejects(3, {{0, 1}}, call);
        check_rejects(2, {{0, 0}, {0, 1}}, call);
        check_rejects(6, {{0, 1}, {1, 2}, {3, 4}, {4, 5}, {5, 3}}, call);
    }

    for (int n = 1; n <= 7; n++) check_exhaustive(n);

    for (long long it = 0; it < stress::scaled(1500); it++) check_small_pair(it % 8 + 1);

    for (long long it = 0; it < stress::scaled(150); it++) check_batch(40, it % 3 ? 12 : 60);

    for (long long it = 0; it < stress::scaled(50); it++) check_large(stress::rand_int(1, 5000), it % 5);

    for (int n : {4, 5, 6, 7, 10, 11}) check_path(n);

    check_large(100000, 0);
    check_large(200000, 2);
    check_path(200000);

    return 0;
}
