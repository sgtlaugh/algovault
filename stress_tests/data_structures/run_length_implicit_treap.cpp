#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/run_length_implicit_treap.cpp"
#undef main

/// Runs that are arithmetic progressions a, a + d, a + 2d, ... with range add and range sum
struct Progression{
    long long a, d;
};

struct ProgressionRuns{
    using Value = Progression;
    using Agg = long long;
    using Lazy = long long;

    static Agg run_agg(Value x, long long count){ return x.a * count + x.d * (count * (count - 1) / 2); }
    static Agg combine(Agg a, Agg b){ return a + b; }
    static pair<Value, Value> cut(Value x, long long k){ return {x, {x.a + k * x.d, x.d}}; }
    static Value apply(Value x, Lazy lz){ return {x.a + lz, x.d}; }
    static Agg apply(Agg a, Lazy lz, long long count){ return a + lz * count; }
    static Lazy compose(Lazy old_lz, Lazy new_lz){ return old_lz + new_lz; }
    static Value reversed(Value x, long long count){ return {x.a + (count - 1) * x.d, -x.d}; }
    static Agg reversed(Agg a){ return a; }
};

/// An Agg that reads differently backwards, so a wrong reversed(agg) or a lazy applied in the wrong order shows
struct Ends{
    long long first, last, sum;
};

struct EndsRuns{
    using Value = long long;
    using Agg = Ends;
    using Lazy = long long;

    static Agg run_agg(Value x, long long count){ return {x, x, x * count}; }
    static Agg combine(Agg a, Agg b){ return {a.first, b.last, a.sum + b.sum}; }
    static Value apply(Value x, Lazy lz){ return x + lz; }
    static Agg apply(Agg a, Lazy lz, long long count){ return {a.first + lz, a.last + lz, a.sum + lz * count}; }
    static Lazy compose(Lazy old_lz, Lazy new_lz){ return old_lz + new_lz; }
    static Value reversed(Value x, long long){ return x; }
    static Agg reversed(Agg a){ return {a.last, a.first, a.sum}; }
};

/// Polynomial hashes read forwards and backwards: a combine that does not commute, for the descents, and an Agg
/// whose reverse differs, which only a descent can see (queries recompute every Agg on their split paths)
/// No Lazy, so the treap must work without one
struct Hash{
    long long h, hr, pw;
};

struct HashRuns{
    static constexpr long long MOD = 1000000007, B = 131;
    using Value = long long;
    using Agg = Hash;

    static Agg run_agg(Value x, long long count){
        Agg a{0, 0, 1};
        for (long long i = 0; i < count; i++) a.h = (a.h * B + x) % MOD, a.pw = a.pw * B % MOD;
        a.hr = a.h;   /// a run of equal values reads the same backwards
        return a;
    }
    static Agg combine(Agg a, Agg b){ return {(a.h * b.pw + b.h) % MOD, (b.hr * a.pw + a.hr) % MOD, a.pw * b.pw % MOD}; }
    static Value reversed(Value x, long long){ return x; }
    static Agg reversed(Agg a){ return {a.hr, a.h, a.pw}; }
};

/// Only the required parts, no Lazy, reversed or cut: insert, erase, assign, query and the descents still work
struct MinRuns{
    using Value = long long;
    using Agg = long long;

    static Agg run_agg(Value x, long long){ return x; }
    static Agg combine(Agg a, Agg b){ return min(a, b); }
};

/// Sums of 1e18 long runs do not fit in long long
struct WideSumRuns{
    using Value = long long;
    using Agg = __int128;
    using Lazy = long long;

    static Agg run_agg(Value x, long long count){ return (__int128)x * count; }
    static Agg combine(Agg a, Agg b){ return a + b; }
    static Value apply(Value x, Lazy lz){ return x + lz; }
    static Agg apply(Agg a, Lazy lz, long long count){ return a + (__int128)lz * count; }
    static Lazy compose(Lazy old_lz, Lazy new_lz){ return old_lz + new_lz; }
    static Value reversed(Value x, long long){ return x; }
    static Agg reversed(Agg a){ return a; }
};

vector<pair<long long, long long>> random_runs(int count, long long max_len, long long lo, long long hi){
    vector<pair<long long, long long>> runs;
    for (int i = 0; i < count; i++) runs.push_back({stress::rand_int(lo, hi), stress::rand_int(1, max_len)});
    return runs;
}

long long brute_sum(const vector<long long>& a, long long l, long long r){
    return accumulate(a.begin() + l, a.begin() + r + 1, 0LL);
}

/// Every operation of the default SumRuns against a plain array, values kept non-negative so sums are monotone
void check_sum_runs(int ops, long long max_len){
    vector<pair<long long, long long>> runs = random_runs(stress::rand_int(0, 4), max_len, 0, 5);
    vector<long long> a;
    for (auto [v, len] : runs) a.insert(a.end(), len, v);
    bool from_runs = stress::rand_int(0, 1);   /// both constructors
    if (!from_runs) fill(a.begin(), a.end(), 7);
    RunLengthImplicitTreap<SumRuns> t = from_runs ? RunLengthImplicitTreap<SumRuns>(runs) : RunLengthImplicitTreap<SumRuns>(a.size(), 7);

    for (int op = 0; op < ops; op++){
        int n = a.size(), kind = stress::rand_int(0, 9);
        if (n == 0 || kind == 0){
            long long pos = stress::rand_int(0, n), v = stress::rand_int(0, 5), len = stress::rand_int(1, max_len);
            t.insert(pos, v, len);
            a.insert(a.begin() + pos, len, v);
        }
        else{
            long long l = stress::rand_int(0, n - 1), r = stress::rand_int(l, min<long long>(n - 1, l + 2 * max_len));
            if (kind == 1 && n > 1){
                t.erase(l, r);
                a.erase(a.begin() + l, a.begin() + r + 1);
            }
            else if (kind == 2){
                long long v = stress::rand_int(0, 5);
                t.assign(l, r, v);
                fill(a.begin() + l, a.begin() + r + 1, v);
            }
            else if (kind == 3){
                long long lz = stress::rand_int(0, 3);
                t.update(l, r, lz);
                for (long long i = l; i <= r; i++) a[i] += lz;
            }
            else if (kind == 4){
                t.reverse(l, r);
                std::reverse(a.begin() + l, a.begin() + r + 1);
            }
            else if (kind == 5){
                assert(t.query(l, r) == brute_sum(a, l, r));
                assert(t.get(l) == a[l]);
            }
            else if (kind == 6){
                long long k = stress::rand_int(0, 40), expected = l - 1;
                while (expected + 1 < n && brute_sum(a, l, expected + 1) <= k) expected++;
                assert(t.max_right(l, [k](long long s){ return s <= k; }) == expected);
                expected = r + 1;
                while (expected - 1 >= 0 && brute_sum(a, expected - 1, r) <= k) expected--;
                assert(t.min_left(r, [k](long long s){ return s <= k; }) == expected);
            }
            else if (kind == 7){
                vector<long long> seen;
                t.for_each_run(l, r, [&](long long v, long long len){
                    assert(len >= 1);
                    seen.insert(seen.end(), len, v);
                });
                assert(seen == vector<long long>(a.begin() + l, a.begin() + r + 1));
            }
            else{
                assert(t.max_right(n, [](long long){ return false; }) == n - 1);   /// empty range from the end
                assert(t.min_left(-1, [](long long){ return false; }) == 0);
            }
        }
        assert(t.size() == (long long)a.size());
    }
}

/// Arithmetic progression runs: cut must restart the right part, reverse must flip the step
void check_progression_runs(int ops){
    vector<long long> a;
    RunLengthImplicitTreap<ProgressionRuns> t(vector<pair<Progression, long long>>{});
    auto push_run = [&](long long pos, Progression p, long long len){
        vector<long long> run;
        for (long long i = 0; i < len; i++) run.push_back(p.a + i * p.d);
        a.insert(a.begin() + pos, run.begin(), run.end());
    };

    for (int op = 0; op < ops; op++){
        int n = a.size(), kind = stress::rand_int(0, 7);
        if (n == 0 || kind == 0){
            long long pos = stress::rand_int(0, n), len = stress::rand_int(1, 8);
            Progression p{stress::rand_int(-20, 20), stress::rand_int(-5, 5)};
            t.insert(pos, p, len);
            push_run(pos, p, len);
            continue;
        }

        long long l = stress::rand_int(0, n - 1), r = stress::rand_int(l, n - 1);
        if (kind == 1){
            t.erase(l, r);
            a.erase(a.begin() + l, a.begin() + r + 1);
        }
        else if (kind == 2){
            Progression p{stress::rand_int(-20, 20), stress::rand_int(-5, 5)};
            t.assign(l, r, p);
            a.erase(a.begin() + l, a.begin() + r + 1);
            push_run(l, p, r - l + 1);
        }
        else if (kind == 3){
            long long lz = stress::rand_int(-4, 4);
            t.update(l, r, lz);
            for (long long i = l; i <= r; i++) a[i] += lz;
        }
        else if (kind == 4){
            t.reverse(l, r);
            std::reverse(a.begin() + l, a.begin() + r + 1);
        }
        else if (kind == 5){
            assert(t.query(l, r) == brute_sum(a, l, r));
            assert(t.get(r).a == a[r]);
        }
        else{
            vector<long long> seen;
            t.for_each_run(l, r, [&](Progression p, long long len){
                for (long long i = 0; i < len; i++) seen.push_back(p.a + i * p.d);
            });
            assert(seen == vector<long long>(a.begin() + l, a.begin() + r + 1));
        }
        assert(t.size() == (long long)a.size());
    }
}

/// First and last values change places under reverse, and negative lazies interleave with reversals and cuts
void check_ends_runs(int ops){
    vector<long long> a(stress::rand_int(1, 12), 0);
    RunLengthImplicitTreap<EndsRuns> t(a.size(), 0);
    for (int op = 0; op < ops; op++){
        int n = a.size(), kind = stress::rand_int(0, 6);
        long long l = stress::rand_int(0, n - 1), r = stress::rand_int(l, n - 1);
        if (kind == 0){
            long long v = stress::rand_int(-9, 9), len = stress::rand_int(1, 4);
            t.insert(l, v, len);
            a.insert(a.begin() + l, len, v);
        }
        else if (kind == 1){
            long long lz = stress::rand_int(-5, 5);
            t.update(l, r, lz);
            for (long long i = l; i <= r; i++) a[i] += lz;
        }
        else if (kind == 2){
            t.reverse(l, r);
            std::reverse(a.begin() + l, a.begin() + r + 1);
        }
        else if (kind == 3){
            long long v = stress::rand_int(-9, 9);
            t.assign(l, r, v);
            fill(a.begin() + l, a.begin() + r + 1, v);
        }
        else if (kind == 4){
            Ends e = t.query(l, r);
            assert(e.first == a[l] && e.last == a[r] && e.sum == brute_sum(a, l, r));
        }
        else{
            /// A query around a freshly reversed range reads its root's aggregate without splitting through it
            t.reverse(l, r);
            std::reverse(a.begin() + l, a.begin() + r + 1);
            long long x = stress::rand_int(0, l), y = stress::rand_int(r, n - 1);
            Ends e = t.query(x, y);
            assert(e.first == a[x] && e.last == a[y] && e.sum == brute_sum(a, x, y));
        }
    }
}

/// Descents under a merge that does not commute: the longest stretch from l (or ending at r) that matches a target
/// equal to the sequence except at one position, pred is monotone because a mismatch never goes away
void check_hash_descents(int ops){
    vector<long long> a;
    RunLengthImplicitTreap<HashRuns> t(vector<pair<long long, long long>>{});
    for (int op = 0; op < ops; op++){
        int n = a.size();
        if (n == 0 || stress::rand_int(0, 2) == 0){
            long long pos = stress::rand_int(0, n), v = stress::rand_int(0, 2), len = stress::rand_int(1, 5);
            t.insert(pos, v, len);
            a.insert(a.begin() + pos, len, v);
            continue;
        }

        long long l = stress::rand_int(0, n - 1), r = stress::rand_int(0, n - 1);
        int kind = stress::rand_int(0, 5);
        if (kind <= 1){
            long long x = min(l, r), y = max(l, r), v = stress::rand_int(0, 2);
            if (kind == 0) t.assign(x, y, v), fill(a.begin() + x, a.begin() + y + 1, v);
            else t.reverse(x, y), std::reverse(a.begin() + x, a.begin() + y + 1);
            continue;
        }

        vector<long long> target = a;
        long long j = stress::rand_int(0, n - 1);
        if (stress::rand_int(0, 3)) target[j] = (target[j] + 1) % 3;
        vector<long long> pre(n + 1, 0), pw(n + 1, 1);
        for (int i = 0; i < n; i++){
            pre[i + 1] = (pre[i] * HashRuns::B + target[i]) % HashRuns::MOD;
            pw[i + 1] = pw[i] * HashRuns::B % HashRuns::MOD;
        }
        auto segment = [&](long long x, long long len){
            return ((pre[x + len] - pre[x] * pw[len]) % HashRuns::MOD + HashRuns::MOD) % HashRuns::MOD;
        };
        auto length = [&](const Hash& a){ return (long long)(find(pw.begin(), pw.end(), a.pw) - pw.begin()); };

        long long expected = l;
        while (expected < n && a[expected] == target[expected]) expected++;
        assert(t.max_right(l, [&](const Hash& s){
            long long len = length(s);
            return l + len <= n && s.h == segment(l, len);
        }) == expected - 1);

        expected = r;
        while (expected >= 0 && a[expected] == target[expected]) expected--;
        assert(t.min_left(r, [&](const Hash& s){
            long long len = length(s);
            return len <= r + 1 && s.h == segment(r - len + 1, len);
        }) == expected + 1);
    }
}

/// Descents through arithmetic progression runs, where the part of a run tried must be the right one at either end
void check_progression_descents(int ops){
    vector<long long> a;
    RunLengthImplicitTreap<ProgressionRuns> t(vector<pair<Progression, long long>>{});
    for (int op = 0; op < ops; op++){
        int n = a.size(), kind = stress::rand_int(0, 3);
        if (n == 0 || kind == 0){
            long long pos = stress::rand_int(0, n), len = stress::rand_int(1, 8);
            Progression p{stress::rand_int(0, 9), stress::rand_int(0, 4)};   /// non-negative values keep sums monotone
            t.insert(pos, p, len);
            vector<long long> run;
            for (long long i = 0; i < len; i++) run.push_back(p.a + i * p.d);
            a.insert(a.begin() + pos, run.begin(), run.end());
            continue;
        }

        long long l = stress::rand_int(0, n - 1), r = stress::rand_int(l, n - 1);
        if (kind == 1){
            t.reverse(l, r);
            std::reverse(a.begin() + l, a.begin() + r + 1);
            continue;
        }

        long long k = stress::rand_int(0, 150), expected = l - 1;
        while (expected + 1 < n && brute_sum(a, l, expected + 1) <= k) expected++;
        assert(t.max_right(l, [k](long long s){ return s <= k; }) == expected);
        expected = r + 1;
        while (expected - 1 >= 0 && brute_sum(a, expected - 1, r) <= k) expected--;
        assert(t.min_left(r, [k](long long s){ return s <= k; }) == expected);
    }
}

/// A Runs type with only Value, Agg, run_agg and combine (a minimum) against a plain array
void check_min_runs(int ops){
    vector<long long> a;
    RunLengthImplicitTreap<MinRuns> t(vector<pair<long long, long long>>{});
    for (int op = 0; op < ops; op++){
        int n = a.size(), kind = stress::rand_int(0, 4);
        if (n == 0 || kind == 0){
            long long pos = stress::rand_int(0, n), x = stress::rand_int(0, 20), count = stress::rand_int(1, 4);
            t.insert(pos, x, count);
            a.insert(a.begin() + pos, count, x);
            continue;
        }

        long long l = stress::rand_int(0, n - 1), r = stress::rand_int(l, n - 1);
        if (kind == 1 && n > 1){
            t.erase(l, r);
            a.erase(a.begin() + l, a.begin() + r + 1);
        }
        else if (kind == 2){
            long long x = stress::rand_int(0, 20);
            t.assign(l, r, x);
            fill(a.begin() + l, a.begin() + r + 1, x);
        }
        else{
            assert(t.query(l, r) == *min_element(a.begin() + l, a.begin() + r + 1));
            long long k = stress::rand_int(0, 20), expected = l - 1;
            while (expected + 1 < n && a[expected + 1] >= k) expected++;
            assert(t.max_right(l, [k](long long m){ return m >= k; }) == expected);
            expected = r + 1;
            while (expected - 1 >= 0 && a[expected - 1] >= k) expected--;
            assert(t.min_left(r, [k](long long m){ return m >= k; }) == expected);
        }
        assert(t.size() == (long long)a.size());
    }
}

/// Reference for huge lengths: the sequence as a plain list of runs, cut by a linear scan
struct RunList{
    vector<pair<long long, long long>> runs;   /// {value, len}

    long long size() const{
        long long s = 0;
        for (auto& [v, len] : runs) s += len;
        return s;
    }

    /// Index of the first run starting at pos, after cutting the run that straddles it
    int cut(long long pos){
        long long start = 0;
        for (int i = 0; i < (int)runs.size(); i++){
            if (start == pos) return i;
            if (pos < start + runs[i].second){
                long long v = runs[i].first, len = runs[i].second;
                runs[i].second = pos - start;
                runs.insert(runs.begin() + i + 1, {v, len - (pos - start)});
                return i + 1;
            }
            start += runs[i].second;
        }
        return runs.size();
    }

    __int128 sum(long long l, long long r){
        int i = cut(l), j = cut(r + 1);
        __int128 s = 0;
        for (int k = i; k < j; k++) s += (__int128)runs[k].first * runs[k].second;
        return s;
    }

    /// Largest r with sum(l, r) <= k for non-negative values
    long long max_right(long long l, __int128 k){
        int i = cut(l);
        long long r = l - 1;
        for (int k2 = i; k2 < (int)runs.size(); k2++){
            auto [v, len] = runs[k2];
            if ((__int128)v * len <= k){
                k -= (__int128)v * len, r += len;
                continue;
            }
            return r + (long long)(k / v);   /// v > 0 here, as v * len > k >= 0
        }
        return r;
    }
};

/// Runs up to 1e17 long and positions up to 1e18 against the run list
void check_huge(int ops){
    RunList ref;
    RunLengthImplicitTreap<WideSumRuns> t(vector<pair<long long, long long>>{});
    const long long BIG = 100000000000000000LL;
    for (int op = 0; op < ops; op++){
        long long n = ref.size();
        int kind = stress::rand_int(0, 6);
        if (n == 0 || (kind == 0 && n < 10 * BIG)){
            long long pos = stress::rand_int(0, n), v = stress::rand_int(0, 3), len = stress::rand_int(1, BIG);
            t.insert(pos, v, len);
            int i = ref.cut(pos);
            ref.runs.insert(ref.runs.begin() + i, {v, len});
            continue;
        }

        long long l = stress::rand_int(0, n - 1), r = stress::rand_int(l, n - 1);
        if (stress::rand_int(0, 1)) r = min(n - 1, l + stress::rand_int(0, 3));   /// short ranges deep inside long runs
        if (kind == 1 && n > BIG){
            t.erase(l, r);
            int i = ref.cut(l), j = ref.cut(r + 1);
            ref.runs.erase(ref.runs.begin() + i, ref.runs.begin() + j);
        }
        else if (kind == 2){
            long long v = stress::rand_int(0, 3);
            t.assign(l, r, v);
            int i = ref.cut(l), j = ref.cut(r + 1);
            ref.runs.erase(ref.runs.begin() + i, ref.runs.begin() + j);
            ref.runs.insert(ref.runs.begin() + i, {v, r - l + 1});
        }
        else if (kind == 3){
            t.update(l, r, 1);
            int i = ref.cut(l), j = ref.cut(r + 1);
            for (int k = i; k < j; k++) ref.runs[k].first++;
        }
        else if (kind == 4){
            t.reverse(l, r);
            int i = ref.cut(l), j = ref.cut(r + 1);
            std::reverse(ref.runs.begin() + i, ref.runs.begin() + j);
        }
        else if (kind == 5){
            assert(t.query(l, r) == ref.sum(l, r));
        }
        else{
            __int128 k = (__int128)stress::rand_int(0, 3) * stress::rand_int(0, 3 * BIG);
            assert(t.max_right(l, [k](__int128 s){ return s <= k; }) == ref.max_right(l, k));
        }
        assert(t.size() == ref.size());
    }
}

/// Every operation allocates at most 3 nodes, the bound the header promises, and a long run stays fast
void check_node_bound_and_speed(){
    const int q = 200000;
    const long long N = 1000000000000000000LL;
    RunLengthImplicitTreap<WideSumRuns> t(N, 1);
    long long ops = 0;
    auto start = chrono::steady_clock::now();
    for (int i = 0; i < q; i++){
        long long n = t.size(), l = stress::rand_int(0, n - 1), r = stress::rand_int(l, n - 1);
        int kind = stress::rand_int(0, 4);
        if (kind == 0) t.assign(l, r, stress::rand_int(0, 3));
        else if (kind == 1) t.update(l, r, 1);
        else if (kind == 2) t.reverse(l, r);
        else if (kind == 3) t.insert(l, stress::rand_int(0, 3), stress::rand_int(1, 1000));
        else t.query(l, r);
        ops++;
        assert((long long)t.nodes.size() <= 1 + 3 * ops);
    }
    double seconds = chrono::duration<double>(chrono::steady_clock::now() - start).count();
    fprintf(stderr, "2e5 operations on a 1e18 sequence: %.2fs, %zu nodes\n", seconds, t.nodes.size());

    /// Queries alone only cut: every piece comes from the one initial run, which must not leave them chained
    RunLengthImplicitTreap<WideSumRuns> cuts(N, 1);
    for (int i = 0; i < q; i++){
        long long l = stress::rand_int(0, N - 1), r = min(N - 1, l + stress::rand_int(0, 4));
        assert(cuts.query(l, r) == r - l + 1);
    }
    vector<pair<int, int>> stack = {{cuts.root, 1}};
    int depth = 0;
    while (!stack.empty()){
        auto [x, d] = stack.back();
        stack.pop_back();
        if (x == -1) continue;
        depth = max(depth, d);
        stack.push_back({cuts.nodes[x].l, d + 1}), stack.push_back({cuts.nodes[x].r, d + 1});
    }
    fprintf(stderr, "2e5 queries cutting one run: %zu nodes, depth %d\n", cuts.nodes.size(), depth);
    assert(depth <= 120);
}

int main(){
    for (long long it = 0; it < stress::scaled(3000); it++) check_sum_runs(60, it % 3 ? 4 : 1);
    for (long long it = 0; it < stress::scaled(3000); it++) check_progression_runs(60);
    for (long long it = 0; it < stress::scaled(3000); it++) check_ends_runs(60);
    for (long long it = 0; it < stress::scaled(2000); it++) check_hash_descents(60);
    for (long long it = 0; it < stress::scaled(2000); it++) check_progression_descents(60);
    for (long long it = 0; it < stress::scaled(2000); it++) check_min_runs(60);
    for (long long it = 0; it < stress::scaled(1500); it++) check_huge(60);
    check_node_bound_and_speed();
    return 0;
}
