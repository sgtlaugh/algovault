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

    return 0;
}
