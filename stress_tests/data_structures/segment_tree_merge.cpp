#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/segment_tree_merge.cpp"
#undef main

/// Every node reachable from a live root holds cnt > 0 equal to its children's sum, and no node leaks from the pool
void check_pool(const SegmentTreeMerge& st, const vector<int>& roots){
    long long reachable = 0;
    vector<array<int, 3>> stack;
    for (int root : roots) if (root) stack.push_back({root, 0, st.n - 1});

    while (!stack.empty()){
        auto [cur, lo, hi] = stack.back();
        stack.pop_back();
        reachable++;
        const auto& node = st.nodes[cur];
        assert(node.cnt > 0);
        if (lo == hi){
            assert(!node.l && !node.r);
            continue;
        }

        int mid = lo + (hi - lo) / 2;
        assert(node.cnt == st.nodes[node.l].cnt + st.nodes[node.r].cnt);
        if (node.l) stack.push_back({node.l, lo, mid});
        if (node.r) stack.push_back({node.r, mid + 1, hi});
    }
    assert(reachable == (long long)st.nodes.size() - 1 - (long long)st.free_nodes.size());
}

long long live_nodes(const SegmentTreeMerge& st){
    return (long long)st.nodes.size() - (long long)st.free_nodes.size();
}

/// A root-to-leaf path over [0, n) with halving ranges has ceil(log2 n) + 1 nodes
int path_nodes(int n){
    int depth = 0;
    while ((1LL << depth) < n) depth++;
    return depth + 1;
}

int random_key(int n){
    int pick = stress::rand_int(0, 9);
    if (pick == 0) return INT_MIN;
    if (pick == 1) return INT_MAX;
    return stress::rand_int(-1, n + 1);
}

long long kth_of(const multiset<int>& s, long long k){
    return *next(s.begin(), k);
}

/// Random inserts, merges, splits and queries over several sets against one std::multiset per set
void check(int n, int sets, int ops, int check_every){
    SegmentTreeMerge st(n);
    vector<int> roots(sets, 0);
    vector<multiset<int>> ref(sets);
    int bound = path_nodes(n);

    for (int op = 0; op < ops; op++){
        int i = stress::rand_int(0, sets - 1), j = stress::rand_int(0, sets - 1);
        int type = stress::rand_int(0, 9);
        if (type <= 3){
            int x = stress::rand_int(0, n - 1), c = stress::rand_int(1, 3);
            if (stress::rand_int(0, 4) == 0) x = stress::rand_int(0, 1) ? 0 : n - 1;
            long long before = live_nodes(st);
            size_t pool = st.nodes.size(), freed = st.free_nodes.size();
            st.insert(roots[i], x, c);
            assert(live_nodes(st) - before <= bound);
            if ((int)freed >= bound) assert(st.nodes.size() == pool);  /// freed nodes are reused before the pool grows
            for (int k = 0; k < c; k++) ref[i].insert(x);
        }
        else if (type == 4 && i != j){
            roots[i] = st.merge(roots[i], roots[j]);
            roots[j] = 0;
            ref[i].insert(ref[j].begin(), ref[j].end());
            ref[j].clear();
        }
        else if (type == 5 && i != j){
            long long k = stress::rand_int(0, ref[i].size() + 1);
            long long before = live_nodes(st);
            int rest = st.split_k(roots[i], k);
            assert(live_nodes(st) - before <= bound);
            roots[j] = st.merge(roots[j], rest);
            while ((long long)ref[i].size() > k){
                ref[j].insert(*prev(ref[i].end()));
                ref[i].erase(prev(ref[i].end()));
            }
        }
        else if (type == 6 && i != j){
            int x = random_key(n);
            long long before = live_nodes(st);
            int rest = st.split_key(roots[i], x);
            assert(live_nodes(st) - before <= bound);
            roots[j] = st.merge(roots[j], rest);
            auto from = ref[i].lower_bound(x);
            ref[j].insert(from, ref[i].end());
            ref[i].erase(from, ref[i].end());
        }
        else if (type == 7 && !ref[i].empty()){
            long long k = stress::rand_int(0, ref[i].size() - 1);
            assert(st.kth(roots[i], k) == kth_of(ref[i], k));
        }
        else{
            int lo = random_key(n), hi = random_key(n);
            if (lo > hi) swap(lo, hi);
            long long expected = distance(ref[i].lower_bound(lo), ref[i].upper_bound(hi));
            assert(st.count(roots[i], lo, hi) == expected);
        }

        assert(st.size(roots[i]) == (long long)ref[i].size() && st.size(roots[j]) == (long long)ref[j].size());
        if (op % check_every == 0) check_pool(st, roots);
    }

    check_pool(st, roots);
    for (int i = 0; i < sets; i++){
        assert(st.size(roots[i]) == (long long)ref[i].size());
        long long k = 0;
        for (int x : ref[i]) assert(st.kth(roots[i], k++) == x);
    }
}

int main(){
    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = it < 64 ? it + 1 : stress::rand_int(1, it % 2 ? 40 : 1000000000);
        check(n, stress::rand_int(1, 6), stress::rand_int(1, 80), 1);
    }

    for (int n : {1, 2, 3, 127, 128, 129}) check(n, 5, 3000, 7);
    check(1000000000, 40, 30000, 5000);
    check(200000, 20, 30000, 5000);

    /// Counts and ranks past the int range
    SegmentTreeMerge wide(1000000000);
    int w = 0;
    wide.insert(w, 999999999);
    wide.insert(w, 0, 1000000000000LL);
    assert(wide.size(w) == 1000000000001LL && wide.count(w, 0, 0) == 1000000000000LL);
    assert(wide.kth(w, 999999999999LL) == 0 && wide.kth(w, 1000000000000LL) == 999999999);

    return 0;
}
