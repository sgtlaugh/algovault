#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/interval_set.cpp"
#undef main

/// Cell i stands for [coord(i), coord(i + 1)); a boolean array of cells is the reference
void simulate(int n, long long base, long long step, int ops){
    auto coord = [&](int i){ return base + i * step; };
    IntervalSet<long long> s;
    vector<bool> on(n, false);

    for (int op = 0; op < ops; op++){
        int l = stress::rand_int(0, n), r = stress::rand_int(0, n);
        if (stress::rand_int(0, 5) == 0) r = l;
        bool adding = stress::rand_int(0, 1);
        if (adding) s.add(coord(l), coord(r));
        else s.remove(coord(l), coord(r));
        for (int i = l; i < r; i++) on[i] = adding;

        map<long long, long long> runs;
        long long total = 0;
        for (int i = 0; i < n; i++){
            if (!on[i] || (i > 0 && on[i - 1])) continue;
            int j = i;
            while (j < n && on[j]) j++;
            runs[coord(i)] = coord(j);
            total += coord(j) - coord(i);
        }
        assert(s.seg == runs);
        assert(s.covered_length() == total);

        for (int i = 0; i < n; i++){
            assert(s.contains(coord(i)) == on[i]);
            assert(s.contains(coord(i) + step - 1) == on[i]);
        }
        assert(!s.contains(coord(n)) && !s.contains(coord(0) - 1));

        for (int q = 0; q < 10; q++){
            int a = stress::rand_int(0, n), b = stress::rand_int(a, n);
            bool all = true;
            for (int i = a; i < b; i++) all = all && on[i];
            assert(s.covers(coord(a), coord(b)) == all);

            if (step > 1 && a < b){
                long long ql = coord(a) + stress::rand_int(0, step - 1), qr = coord(b) - stress::rand_int(0, step - 1);
                assert(s.covers(ql, qr) == (qr <= ql || all));
            }
        }
    }
}

/// Large random workload in int, total stays <= 1e9: checks the canonical form and the running total
void large_run(int ops){
    IntervalSet<int> s;
    for (int op = 0; op < ops; op++){
        int l = stress::rand_int(0, 1000000000);
        int len = stress::rand_int(0, op % 3 ? 1000 : 100000000);
        int r = (int)min(1000000000LL, (long long)l + len);
        if (stress::rand_int(0, 2)) s.add(l, r);
        else s.remove(l, r);
    }

    long long total = 0;
    int last = INT_MIN;
    bool first = true;
    for (auto [l, r] : s.seg){
        assert(l < r && (first || last < l));
        total += r - l, last = r, first = false;
    }
    assert(total == s.covered_length());
}

/// IntervalMap against an array of cell values, cell i standing for [coord(i), coord(i + 1)) as in simulate
void simulate_map(int n, long long base, long long step, int ops){
    auto coord = [&](int i){ return base + i * step; };
    long long start = stress::rand_int(-5, 5);
    IntervalMap<long long, long long> m(coord(0), coord(n), start);
    vector<long long> val(n, start);

    for (int op = 0; op < ops; op++){
        int l = stress::rand_int(0, n), r = stress::rand_int(0, n);
        if (l > r) swap(l, r);
        if (stress::rand_int(0, 5) == 0) r = l;
        int kind = stress::rand_int(0, 3);

        if (kind == 0){
            long long v = stress::rand_int(-5, 5);
            m.assign(coord(l), coord(r), v);
            for (int i = l; i < r; i++) val[i] = v;
        }
        else if (kind == 1){
            long long d = stress::rand_int(-3, 3);
            m.for_each(coord(l), coord(r), [&](long long, long long, long long& v){ v += d; });
            for (int i = l; i < r; i++) val[i] += d;
        }
        else{
            long long sum = 0, prev_end = coord(l);
            int visits = 0;
            m.for_each(coord(l), coord(r), [&](long long a, long long b, long long v){
                assert(a == prev_end && a < b && (b - a) % step == 0);
                for (long long c = a; c < b; c += step) assert(val[(c - base) / step] == v);
                sum += (b - a) / step * v, prev_end = b, visits++;
            });
            assert(prev_end == coord(r) && (visits > 0) == (l < r));

            long long want = 0;
            for (int i = l; i < r; i++) want += val[i];
            assert(sum == want);
        }

        long long last = LLONG_MIN;
        for (auto [k, v] : m.seg){
            assert(last < k && coord(0) <= k && k <= coord(n) && (k - base) % step == 0);
            last = k;
        }
        assert(m.seg.begin()->first == coord(0) && last == coord(n));

        for (int i = 0; i < n; i++){
            assert(m.get(coord(i)) == val[i]);
            assert(m.get(coord(i) + step - 1) == val[i]);
        }
    }
}

/// Pieces split by for_each stay put, so a fragmented map keeps its exact values; assign then collapses them
void fragment_then_collapse(){
    int n = 2000;
    IntervalMap<int, long long> m(0, n, 0);
    for (int i = 0; i < n; i++) m.for_each(i, i + 1, [&](int, int, long long& v){ v = i; });
    assert((int)m.seg.size() == n + 1);
    for (int i = 0; i < n; i++) assert(m.get(i) == i);

    m.assign(0, n, 7);
    assert((m.seg == map<int, long long>{{0, 7}, {n, 0}}));
}

int main(){
    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = it < 300 ? it % 12 + 1 : stress::rand_int(1, 60);
        if (it % 3 == 0) simulate(n, stress::rand_int(-1000, 1000), 1, 40);
        else if (it % 3 == 1) simulate(n, stress::rand_int(-1000000, 1000000), stress::rand_int(2, 1000), 40);
        else simulate(n, -4000000000000000000LL, 8000000000000000000LL / n, 40);
    }

    large_run(stress::scaled(300000));

    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = it < 300 ? it % 12 + 1 : stress::rand_int(1, 60);
        if (it % 3 == 0) simulate_map(n, stress::rand_int(-1000, 1000), 1, 40);
        else if (it % 3 == 1) simulate_map(n, stress::rand_int(-1000000, 1000000), stress::rand_int(2, 1000), 40);
        else simulate_map(n, -4000000000000000000LL, 8000000000000000000LL / n, 40);
    }
    fragment_then_collapse();

    IntervalMap<int, int> bounded(0, 10, 1);
    stress::expect_assert("data_structures/interval_set.cpp", [&]{ bounded.assign(-1, 3, 2); });
    stress::expect_assert("data_structures/interval_set.cpp", [&]{ bounded.for_each(5, 11, [](int, int, int&) {}); });
    stress::expect_assert("data_structures/interval_set.cpp", []{ IntervalMap<int, int> empty(4, 4, 0); });

    /// Floating T, which simulate never uses: cuts at fractional dyadic points keep their exact endpoints
    IntervalSet<double> reals;
    reals.add(0.5, 1.75), reals.add(1.25, 2.0), reals.add(3.0, 3.5);
    assert((reals.seg == map<double, double>{{0.5, 2.0}, {3.0, 3.5}}) && reals.covered_length() == 2.0);
    reals.remove(1.5, 3.25);
    assert((reals.seg == map<double, double>{{0.5, 1.5}, {3.25, 3.5}}) && reals.covered_length() == 1.25);
    assert(reals.contains(1.25) && !reals.contains(1.5) && reals.covers(3.25, 3.5) && !reals.covers(1.25, 3.3));

    return 0;
}
