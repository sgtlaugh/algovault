#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/segment_tree_beats.cpp"
#undef main

void update(SegmentTreeBeats& st, vector<long long>& a, int type, int l, int r, long long x){
    if (type == 0){
        st.chmin(l, r, x);
        for (int i = l - 1; i < r; i++) a[i] = min(a[i], x);
    }
    else if (type == 1){
        st.chmax(l, r, x);
        for (int i = l - 1; i < r; i++) a[i] = max(a[i], x);
    }
    else if (type == 2){
        st.add(l, r, x);
        for (int i = l - 1; i < r; i++) a[i] += x;
    }
    else{
        st.assign(l, r, x);
        for (int i = l - 1; i < r; i++) a[i] = x;
    }
}

/// Sums are compared mod 2^64 because node sums may wrap while the documented queried sums still fit
void check(SegmentTreeBeats& st, const vector<long long>& a, int l, int r){
    __int128 sum = 0;
    for (int i = l - 1; i < r; i++) sum += a[i];
    assert((unsigned long long)st.query_sum(l, r) == (unsigned long long)sum);
    assert(st.query_min(l, r) == *min_element(a.begin() + l - 1, a.begin() + r));
    assert(st.query_max(l, r) == *max_element(a.begin() + l - 1, a.begin() + r));
}

void run(int n, int ops, long long lo, long long hi){
    vector<long long> a(n, 0);
    bool from_array = stress::rand_int(0, 1);
    if (from_array) for (auto& x : a) x = stress::rand_int(lo, hi);
    SegmentTreeBeats st = from_array ? SegmentTreeBeats(a) : SegmentTreeBeats(n);

    for (int op = 0; op < ops; op++){
        int type = stress::rand_int(0, 6);
        int l = stress::rand_int(1, n), r = stress::rand_int(l, n);
        long long x = stress::rand_int(lo, hi);

        if (type < 4) update(st, a, type, l, r, x);
        else check(st, a, l, r);
    }
}

/// Two updates with |x| up to 4.5e18 from zeros spend the whole 9e18 budget, so a chmin / chmax delta
/// times a count of 2+ leaves the signed range and only the unsigned sum arithmetic survives it
void run_extreme(int n){
    const long long C = 4500000000000000000LL;
    vector<long long> a(n, 0);
    SegmentTreeBeats st(n);

    for (int op = 0; op < 2; op++){
        int l = stress::rand_int(1, n), r = stress::rand_int(l, n);
        update(st, a, stress::rand_int(0, 3), l, r, stress::rand_int(-C, C));
        for (int ql = 1; ql <= n; ql++){
            for (int qr = ql; qr <= n; qr++) check(st, a, ql, qr);
        }
    }
}

/// n = 1e6 is the classic beats size (HDU 5306 Gorgeous Sequence, 256 MB), so the tree must leave room for the input
/// Most ranges are short to keep the brute force cheap, a few span nearly everything to reach the deep nodes
void run_million(){
    const int n = 1000000;
    vector<long long> a(n);
    for (auto& x : a) x = stress::rand_int(-1000000000, 1000000000);
    SegmentTreeBeats st(a);
    assert(st.tree.size() * sizeof(SegmentTreeBeats::Node) <= (200u << 20));

    for (int op = 0; op < 100000; op++){
        int l = stress::rand_int(1, n), r = min(n, l + (int)stress::rand_int(0, 64));
        if (op % 2000 == 0) l = stress::rand_int(1, 10), r = stress::rand_int(n - 10, n);
        int type = stress::rand_int(0, 6);
        long long x = stress::rand_int(-1000000000, 1000000000);

        if (type < 4) update(st, a, type, l, r, x);
        else check(st, a, l, r);
    }
    check(st, a, 1, n);
}

int main(){
    for (long long it = 0; it < stress::scaled(20000); it++) run(stress::rand_int(1, 10), 60, -4, 4);
    for (long long it = 0; it < stress::scaled(3000); it++) run(stress::rand_int(1, 200), 300, -1000, 1000);
    for (long long it = 0; it < stress::scaled(20); it++) run(stress::rand_int(1000, 3000), 3000, -1000000000, 1000000000);

    /// max |a[i]| + sum of |x| reaches up to 31 * 2.9e17 ~= 9e18, the documented bound, and node sums wrap
    const long long B = 290000000000000000LL;
    for (long long it = 0; it < stress::scaled(3000); it++) run(stress::rand_int(1, 40), 30, -B, B);
    for (long long it = 0; it < stress::scaled(20000); it++) run_extreme(stress::rand_int(1, 8));
    /// The build visits every node, so a tree sized one level short reads out of range for some n here
    for (int n = 1; n <= 4100; n++){
        vector<long long> a(n);
        for (auto& x : a) x = stress::rand_int(-4, 4);
        SegmentTreeBeats st(a);
        check(st, a, 1, n);
    }
    run_million();

    return 0;
}
