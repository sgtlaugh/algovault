#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/historic_segment_tree.cpp"
#undef main

/// The reference keeps every cell's current value, the largest value it ever held, and its sum over ticks
struct History{
    vector<long long> a, peak, total;

    History(const vector<long long>& init) : a(init), peak(init), total(init.size(), 0) {}

    void add(int l, int r, long long x){
        for (int i = l; i <= r; i++) a[i] += x, peak[i] = max(peak[i], a[i]);
    }

    void tick(){
        for (size_t i = 0; i < a.size(); i++) total[i] += a[i];
    }
};

void check(HistoricSegmentTree& st, const History& h, int l, int r){
    long long sum = 0, mx = LLONG_MIN, peak = LLONG_MIN, total = 0;
    for (int i = l; i <= r; i++) sum += h.a[i], mx = max(mx, h.a[i]), peak = max(peak, h.peak[i]), total += h.total[i];
    assert(st.query_sum(l, r) == sum);
    assert(st.query_max(l, r) == mx);
    assert(st.query_historic_max(l, r) == peak);
    assert(st.query_historic_sum(l, r) == total);
}

/// burst > 1 runs that many updates between queries, so tags stack several levels deep before any push
void run(int n, int ops, long long lim, int burst){
    vector<long long> init(n, 0);
    bool from_array = stress::rand_int(0, 1);
    if (from_array) for (auto& x : init) x = stress::rand_int(-lim, lim);
    HistoricSegmentTree st = from_array ? HistoricSegmentTree(init) : HistoricSegmentTree(n);
    History h(init);

    for (int it = 0; it < ops; it++){
        for (int b = 0; b < burst; b++){
            if (stress::rand_int(0, 3) == 0){
                st.tick();
                h.tick();
                continue;
            }

            int l = stress::rand_int(0, n - 1), r = stress::rand_int(l, n - 1);
            if (stress::rand_int(0, 2) == 0) l = 0, r = n - 1;
            long long x = stress::rand_int(-lim, lim);
            st.add(l, r, x);
            h.add(l, r, x);
        }

        int l = stress::rand_int(0, n - 1), r = stress::rand_int(l, n - 1);
        check(st, h, l, r);
    }

    for (int l = 0; l < n; l++){
        for (int r = l; r < n; r++) check(st, h, l, r);
    }
}

int main(){
    /// Every small size, mixing single-step and bursty update patterns
    for (int it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(1, 12);
        run(n, stress::rand_int(1, 40), stress::rand_int(0, 1) ? 5 : 1000000, stress::rand_int(0, 1) ? 1 : 8);
    }

    for (int it = 0; it < stress::scaled(60); it++){
        int n = stress::rand_int(50, 300);
        run(n, 400, 1000000, stress::rand_int(1, 20));
    }

    /// 2e5 cells ticked 1000 times: the root's historic sum reaches 1e17, every node's stays in range
    HistoricSegmentTree big(200000);
    for (int it = 0; it < 1000; it++){
        big.add(0, 199999, it % 2 ? -1000000000 : 1000000000);
        big.tick();
    }
    assert(big.query_historic_max(0, 199999) == 1000000000);
    assert(big.query_historic_sum(0, 199999) == 100000000000000000LL);
    assert(big.query_historic_sum(7, 7) == 500000000000LL);
    assert(big.query_sum(0, 199999) == 0);

    /// Every value, range sum and historic sum stays in range while the pending tick terms cancel past 1e19
    const long long X = 1000000000000000LL;
    HistoricSegmentTree lifted(vector<long long>(4, -X));
    lifted.add(0, 3, X);
    for (int it = 0; it < 20000; it++) lifted.tick();
    lifted.add(1, 2, -X);
    lifted.tick();
    assert(lifted.query_historic_sum(0, 3) == -2 * X);
    assert(lifted.query_historic_sum(1, 1) == -X);
    assert(lifted.query_historic_max(0, 3) == 0);
    assert(lifted.query_sum(0, 3) == -2 * X);

    /// 20000 ticks pushed at once onto a child holding a pending add of 1e15
    HistoricSegmentTree split(vector<long long>{-X, -X, 0, 0});
    split.add(0, 1, X);
    for (int it = 0; it < 20000; it++) split.tick();
    assert(split.query_historic_sum(0, 0) == 0);
    assert(split.query_historic_sum(0, 3) == 0);

    /// A sum change of 1.6e19 between two in-range sums: add times len overflows on its own
    const long long Y = 4000000000000000000LL;
    HistoricSegmentTree swing(vector<long long>{-Y, -Y});
    swing.add(0, 1, 2 * Y);
    assert(swing.query_sum(0, 0) == Y);
    assert(swing.query_max(0, 1) == Y);

    stress::expect_assert("data_structures/historic_segment_tree.cpp", []{ HistoricSegmentTree bad(0); });
    stress::expect_assert("data_structures/historic_segment_tree.cpp", []{ HistoricSegmentTree st(5); st.add(3, 2, 1); });
    stress::expect_assert("data_structures/historic_segment_tree.cpp", []{ HistoricSegmentTree st(5); st.query_historic_sum(0, 5); });

    return 0;
}
