#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/mo.cpp"
#undef main

/// Distinct values and sum of every range against direct scans
void check_array(int n, int q, int max_v, int verify = 1 << 30){
    vector<int> a(n);
    for (auto& x : a) x = stress::rand_int(0, max_v);
    vector<pair<int, int>> queries(q);
    for (auto& [l, r] : queries) l = stress::rand_int(0, n - 1), r = stress::rand_int(l, n - 1);

    vector<int> freq(max_v + 1, 0);
    int distinct = 0, window = 0;
    long long sum = 0;
    vector<int> got_distinct(q);
    vector<long long> got_sum(q);

    for (int round = 0; round < 2; round++){
        mo(n, queries,
            [&](int i){ distinct += freq[a[i]]++ == 0, sum += a[i], window++; },
            [&](int i){ distinct -= --freq[a[i]] == 0, sum -= a[i], window--; },
            [&](int qi){ got_distinct[qi] = distinct, got_sum[qi] = sum; });
        assert(window == 0 && sum == 0 && distinct == 0);
        for (int i = 0; i < q && i < verify; i++){
            auto [l, r] = queries[i];
            assert(got_sum[i] == accumulate(a.begin() + l, a.begin() + r + 1, 0LL));
            assert(got_distinct[i] == (int)set<int>(a.begin() + l, a.begin() + r + 1).size());
        }
    }
}

/// Hashes (position, value) so a window sum of hashes pins down exactly which cells, at which version, are inside
unsigned long long cell_hash(int i, int v){
    unsigned long long x = (unsigned long long)i * 1000003 + v + 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

/// Distinct values and cell hash of every timed range against replaying the updates on a fresh copy
/// check_calls holds the callback count to 2 * the header bound, random inputs reach about 1.16 at large n
void check_updates(int n, int q, int u, int max_v, int verify = 1 << 30, bool check_calls = false){
    vector<int> a(n), update_pos(u), value(u);
    for (auto& x : a) x = stress::rand_int(0, max_v);
    for (int j = 0; j < u; j++) update_pos[j] = stress::rand_int(0, n - 1), value[j] = stress::rand_int(0, max_v);
    vector<array<int, 3>> queries(q);
    for (auto& [l, r, t] : queries) l = stress::rand_int(0, n - 1), r = stress::rand_int(l, n - 1), t = stress::rand_int(0, u);
    const vector<int> initial = a, initial_value = value;

    vector<int> freq(max_v + 1, 0), got_distinct(q);
    vector<unsigned long long> got_hash(q);
    int distinct = 0;
    unsigned long long hash = 0;
    double block = cbrt((double)n * n);

    for (int round = 0; round < 2; round++){
        long long calls = 0;
        mo_with_updates(n, queries, update_pos,
            [&](int i){ distinct += freq[a[i]]++ == 0, hash += cell_hash(i, a[i]), calls++; },
            [&](int i){ distinct -= --freq[a[i]] == 0, hash -= cell_hash(i, a[i]), calls++; },
            [&](int j){ swap(a[update_pos[j]], value[j]), calls++; },
            [&](int qi){ got_distinct[qi] = distinct, got_hash[qi] = hash; });
        assert(distinct == 0 && hash == 0 && a == initial && value == initial_value);
        assert(!check_calls || calls <= 2 * ((q + u) * block + block * block));

        for (int i = 0; i < q && i < verify; i++){
            auto [l, r, t] = queries[i];
            vector<int> b = initial;
            for (int j = 0; j < t; j++) b[update_pos[j]] = initial_value[j];

            unsigned long long h = 0;
            for (int k = l; k <= r; k++) h += cell_hash(k, b[k]);
            assert(got_hash[i] == h && got_distinct[i] == (int)set<int>(b.begin() + l, b.begin() + r + 1).size());
        }
    }
}

/// max(value * count) and cell hash of every range against direct scans, with add-only state and a contiguity check
/// check_calls holds the add count to 2 * the header bound, random inputs reach about 0.99 at large n
void check_rollback(int n, int q, int max_v, int verify = 1 << 30, bool check_calls = false){
    vector<int> a(n);
    for (auto& x : a) x = stress::rand_int(1, max_v);
    vector<pair<int, int>> queries(q);
    for (auto& [l, r] : queries) l = stress::rand_int(0, n - 1), r = stress::rand_int(l, n - 1);

    struct Saved{
        int pos, wl, wr;
        long long best;
        unsigned long long hash;
    };
    vector<int> freq(max_v + 1, 0);
    vector<Saved> history;
    vector<size_t> marks;
    int wl = 0, wr = -1;
    long long best = 0;
    unsigned long long hash = 0;
    vector<long long> got_best(q);
    vector<unsigned long long> got_hash(q);

    for (int round = 0; round < 2; round++){
        long long calls = 0;
        mo_with_rollback(n, queries,
            [&](int i){
                assert(wl > wr || i == wl - 1 || i == wr + 1);
                calls++;
                history.push_back({i, wl, wr, best, hash});
                if (wl > wr) wl = wr = i;
                else wl = min(wl, i), wr = max(wr, i);
                best = max(best, (long long)a[i] * ++freq[a[i]]), hash += cell_hash(i, a[i]);
            },
            [&](){
                assert(marks.size() < 2);
                marks.push_back(history.size());
            },
            [&](){
                assert(!marks.empty());
                for (; history.size() > marks.back(); history.pop_back()){
                    Saved s = history.back();
                    freq[a[s.pos]]--, wl = s.wl, wr = s.wr, best = s.best, hash = s.hash;
                }
                marks.pop_back();
            },
            [&](int qi){ got_best[qi] = best, got_hash[qi] = hash; });
        assert(marks.empty() && history.empty() && wl > wr && best == 0 && hash == 0);
        assert(!check_calls || calls <= 2 * (n * sqrt((double)q) + q));

        for (int i = 0; i < q && i < verify; i++){
            auto [l, r] = queries[i];
            vector<int> count(max_v + 1, 0);
            long long expected = 0;
            unsigned long long h = 0;
            for (int k = l; k <= r; k++) expected = max(expected, (long long)a[k] * ++count[a[k]]), h += cell_hash(k, a[k]);
            assert(got_best[i] == expected && got_hash[i] == h);
        }
    }
}

vector<pair<int, int>> make_tree(int n, int shape){
    vector<int> label(n);
    iota(label.begin(), label.end(), 0);
    shuffle(label.begin(), label.end(), stress::rng());

    vector<pair<int, int>> edges;
    for (int i = 1; i < n; i++){
        int p = shape == 0 ? stress::rand_int(0, i - 1) : shape == 1 ? i - 1 : shape == 2 ? 0 : (i - 1) / 2;
        edges.push_back({label[i], label[p]});
    }
    return edges;
}

/// Multiset of values on each tree path, nodes or edges, against climbing parents
void check_tree(int n, int q, int shape, int verify = 1 << 30){
    auto edges = make_tree(n, shape);
    int root = stress::rand_int(0, n - 1);
    TreePathMo tree(n);
    for (auto [u, v] : edges) tree.add_edge(u, v);
    tree.build(root);

    vector<int> value(n);
    for (auto& x : value) x = stress::rand_int(0, 20);
    vector<pair<int, int>> paths(q);
    for (auto& [u, v] : paths) u = stress::rand_int(0, n - 1), v = stress::rand_int(0, n - 1);

    for (bool edge_mode : {false, true}){
        vector<int> freq(21, 0);
        int distinct = 0;
        long long sum = 0;
        vector<int> got_distinct(q);
        vector<long long> got_sum(q);

        tree.run(paths,
            [&](int v){ distinct += freq[value[v]]++ == 0, sum += value[v]; },
            [&](int v){ distinct -= --freq[value[v]] == 0, sum -= value[v]; },
            [&](int qi){ got_distinct[qi] = distinct, got_sum[qi] = sum; }, edge_mode);
        assert(sum == 0 && distinct == 0);

        for (int i = 0; i < q && i < verify; i++){
            int u = paths[i].first, v = paths[i].second, l;
            for (int x = u, y = v; ; ){
                if (x == y){
                    l = x;
                    break;
                }
                if (tree.depth[x] >= tree.depth[y]) x = tree.up[0][x];
                else y = tree.up[0][y];
            }

            vector<int> nodes;
            for (int x = u; x != l; x = tree.up[0][x]) nodes.push_back(x);
            for (int x = v; x != l; x = tree.up[0][x]) nodes.push_back(x);
            if (!edge_mode) nodes.push_back(l);

            long long s = 0;
            set<int> vals;
            for (int x : nodes) s += value[x], vals.insert(value[x]);
            assert(got_sum[i] == s && got_distinct[i] == (int)vals.size());
        }
    }
}

int main(){
    for (long long it = 0; it < stress::scaled(2000); it++){
        check_array(stress::rand_int(1, 60), stress::rand_int(1, 60), it % 2 ? 3 : 50);
        check_tree(stress::rand_int(1, 40), stress::rand_int(1, 40), it % 4);
        check_updates(stress::rand_int(1, 60), stress::rand_int(1, 60), stress::rand_int(0, 30), it % 2 ? 3 : 50);
        check_rollback(stress::rand_int(1, 60), stress::rand_int(1, 60), it % 2 ? 3 : 50);
    }

    check_array(100000, 100000, 1000, 300);
    check_updates(20000, 20000, 20000, 1000, 200, true);
    check_rollback(50000, 50000, 1000, 300, true);
    check_tree(100000, 2000, 1, 200);  /// a path, deep enough to break a recursive Euler tour on an 8 MB stack
    check_tree(30000, 30000, 0, 2000);

    return 0;
}
