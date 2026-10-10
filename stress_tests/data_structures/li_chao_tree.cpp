#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/li_chao_tree.cpp"
#undef main

struct Segment{
    long long k, b, l, r;
};

/// Min or max over the segments covering x, evaluated in __int128 so the brute force cannot overflow
template <bool MAXIMIZE>
long long brute_query(const vector<Segment>& segments, long long x){
    bool found = false;
    __int128 res = 0;
    for (const auto& s : segments){
        if (x < s.l || x > s.r) continue;
        __int128 y = (__int128)s.k * x + s.b;
        if (!found || (MAXIMIZE ? y > res : y < res)) res = y, found = true;
    }

    if (!found) return LiChaoTree<long long, MAXIMIZE>::NONE;
    assert(res >= LLONG_MIN && res <= LLONG_MAX);
    return (long long)res;
}

/// Random lines and segments, segment ends may fall outside [lo, hi] or come reversed
template <bool MAXIMIZE>
void run(long long lo, long long hi, long long max_k, long long max_b, int ops, bool exhaustive){
    LiChaoTree<long long, MAXIMIZE> tree(lo, hi);
    vector<Segment> segments;
    vector<long long> probes = {lo, hi, lo + (hi - lo) / 2};
    int lines_only = 0;

    auto random_x = [&](){
        long long span = max(1LL, (hi - lo) / 8);
        return stress::rand_int(lo - span / 2, hi + span / 2);
    };

    for (int op = 0; op < ops; op++){
        long long k = stress::rand_int(0, 3) ? stress::rand_int(-max_k, max_k) : stress::rand_int(-2, 2);
        long long b = stress::rand_int(-max_b, max_b);

        if (stress::rand_int(0, 2)){
            tree.add_line(k, b);
            segments.push_back({k, b, lo, hi});
            lines_only++;
        }
        else{
            long long l = random_x(), r = random_x();
            if (stress::rand_int(0, 3)) tie(l, r) = make_pair(min(l, r), max(l, r));
            tree.add_segment(k, b, l, r);
            segments.push_back({k, b, max(l, lo), min(r, hi)});
            for (long long x : {l - 1, l, r, r + 1}){
                if (lo <= x && x <= hi) probes.push_back(x);
            }
        }

        if (exhaustive){
            for (long long x = lo; x <= hi; x++) assert(tree.query(x) == brute_query<MAXIMIZE>(segments, x));
        }
        else if (op % 16 == 0){
            for (int q = 0; q < 8; q++) probes.push_back(stress::rand_int(lo, hi));
            for (long long x : probes) assert(tree.query(x) == brute_query<MAXIMIZE>(segments, x));
            probes = {lo, hi, lo + (hi - lo) / 2};
        }
    }

    /// A segment splits into at most 2 nodes per level, each insert adds at most 1 node, plus 2 split paths
    int depth = 1;
    while (depth < 64 && ((hi - lo) >> depth) > 0) depth++;
    long long segment_ops = ops - lines_only;
    assert((long long)tree.nodes.size() <= lines_only + segment_ops * 4 * (depth + 1));
    if (lines_only == ops) assert((int)tree.nodes.size() <= ops);
}

/// Path length of a root to leaf walk over [lo, hi], every add_line copies at most this many nodes
int max_path(long long lo, long long hi){
    int len = 1;
    for (unsigned long long size = (unsigned long long)(hi - lo) + 1; size > 1; size = (size + 1) / 2) len++;
    return len;
}

/// Each version picks a random parent version, the brute force walks the parent chain of the queried version
template <bool MAXIMIZE>
void run_persistent(long long lo, long long hi, long long max_k, long long max_b, int ops, bool exhaustive){
    PersistentLiChaoTree<long long, MAXIMIZE> tree(lo, hi);
    vector<int> parent = {-1};
    vector<Segment> lines = {{0, 0, lo, hi}};
    int path = max_path(lo, hi);

    auto brute = [&](int version, long long x){
        vector<Segment> chain;
        for (int v = version; v > 0; v = parent[v]) chain.push_back(lines[v]);
        return brute_query<MAXIMIZE>(chain, x);
    };

    auto check = [&](int version, long long x){
        assert(tree.query(version, x) == brute(version, x));
    };

    for (int op = 0; op < ops; op++){
        int last = (int)parent.size() - 1;
        int choice = stress::rand_int(0, 3);
        int from = choice == 0 ? 0 : choice == 1 ? last : (int)stress::rand_int(0, last);
        long long k = stress::rand_int(0, 3) ? stress::rand_int(-max_k, max_k) : stress::rand_int(-2, 2);
        long long b = stress::rand_int(-max_b, max_b);

        size_t before = tree.nodes.size();
        int version = tree.add_line(from, k, b);
        assert(version == (int)parent.size());
        assert(tree.nodes.size() - before <= (size_t)path);
        parent.push_back(from);
        lines.push_back({k, b, lo, hi});

        if (exhaustive){
            int old = stress::rand_int(0, version);
            for (long long x = lo; x <= hi; x++) check(version, x), check(old, x), check(from, x);
        }
        else if (op % 16 == 0){
            for (int q = 0; q < 12; q++){
                int v = stress::rand_int(0, version);
                for (long long x : {lo, hi, lo + (hi - lo) / 2, stress::rand_int(lo, hi)}) check(v, x);
            }
        }
    }
}

int main(){
    const long long E9 = 1000000000LL, E18 = 1000000000000000000LL;

    for (long long it = 0; it < stress::scaled(1500); it++){
        long long lo = stress::rand_int(-20, 20), hi = lo + stress::rand_int(0, 40);
        int ops = stress::rand_int(1, 25);
        if (it % 2) run<false>(lo, hi, 6, 30, ops, true);
        else run<true>(lo, hi, 6, 30, ops, true);
    }

    for (long long it = 0; it < stress::scaled(40); it++){
        run<false>(-E9, E9, E9, E18, 2000, false);
        run<true>(-E9, E9, E9, E18, 2000, false);
        run<false>(-E18, E18, 4, E18, 2000, false);
        run<true>(-E18, E18, 4, E18, 2000, false);
    }

    for (long long it = 0; it < stress::scaled(1500); it++){
        long long lo = stress::rand_int(-20, 20), hi = lo + stress::rand_int(0, 40);
        int ops = stress::rand_int(1, 30);
        if (it % 2) run_persistent<false>(lo, hi, 6, 30, ops, true);
        else run_persistent<true>(lo, hi, 6, 30, ops, true);
    }

    for (long long it = 0; it < stress::scaled(20); it++){
        run_persistent<false>(-E9, E9, E9, E18, 2000, false);
        run_persistent<true>(-E9, E9, E9, E18, 2000, false);
        run_persistent<false>(-E18, E18, 4, E18, 2000, false);
        run_persistent<true>(-E18, E18, 4, E18, 2000, false);
    }

    for (long long it = 0; it < stress::scaled(200); it++){
        LiChaoTree<long long> tree(-E9, E9);
        vector<Segment> segments;
        int n = stress::rand_int(1, 400);
        for (int i = 0; i < n; i++){
            long long k = stress::rand_int(-E9, E9), b = stress::rand_int(-E18, E18);
            tree.add_line(k, b);
            segments.push_back({k, b, -E9, E9});
        }

        assert((int)tree.nodes.size() <= n);
        for (long long x : {-E9, -E9 + 1, -1LL, 0LL, 1LL, E9 - 1, E9}) assert(tree.query(x) == brute_query<false>(segments, x));
    }

    /// A line that wins or loses at both ends of a node's range settles there, so neither kind descends or adds a node
    for (long long it = 0; it < stress::scaled(200); it++){
        LiChaoTree<long long> settled(-E9, E9);
        LiChaoTree<long long, true> settled_max(-E9, E9);
        long long k = stress::rand_int(-E9, E9), base = stress::rand_int(-E18 / 2, E18 / 2);
        long long low = LLONG_MAX, high = LLONG_MIN;
        for (int i = 0; i < 50; i++){
            long long b = base + stress::rand_int(-E18 / 4, E18 / 4);
            settled.add_line(k, b), settled_max.add_line(k, b);
            low = min(low, b), high = max(high, b);
        }

        assert(settled.nodes.size() == 1 && settled_max.nodes.size() == 1);
        assert(settled.query(1) == k + low && settled_max.query(-1) == -k + high);
    }

    /// Fractional slopes, every other test uses long long
    LiChaoTree<double> real(0, 10);
    real.add_line(0.5, 0.25), real.add_line(-0.5, 4.0);
    assert(abs(real.query(3) - 1.75) < 1e-9 && abs(real.query(10) + 1.0) < 1e-9);

    /// Splitting a negative range with a truncating (l + r) / 2 breaks this case, random tests hit it about once in 1e5
    LiChaoTree<long long> negative(-26, 13);
    for (auto [k, b] : vector<pair<long long, long long>>{{-4, -16}, {1, -16}, {-6, -12}, {6, 4}, {-1, 10}}) negative.add_line(k, b);
    assert(negative.query(-1) == -17);

    return 0;
}
