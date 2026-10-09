/***
 *
 * Interval Set
 * Set of disjoint half-open intervals [l, r) with merging add, splitting remove and coverage queries
 *
 * Complexity: O(log n) amortized per add and remove, O(log n) per query, n = number of stored intervals
 *
 * IntervalSet<T> s:
 *   add(l, r): covers [l, r), merging with every stored interval it overlaps or touches
 *   remove(l, r): uncovers [l, r), cutting stored intervals that stick out of it
 *   contains(x): whether the point x is covered
 *   covers(l, r): whether all of [l, r) is covered (true for an empty range)
 *   covered_length(): total length covered, kept as a running sum
 *
 * s.seg is the map l -> r of the stored intervals in increasing order; it never holds two intervals
 * that touch, so it is the canonical form of the covered set
 * Empty ranges (l >= r) are ignored by add and remove
 * T is an integer or floating point type; the total covered length must fit in T, so long long
 * coordinates in [-4e18, 4e18] are safe
 * With floating T the running sum drifts over long add/remove workloads, so covered_length() can
 * end slightly off (e.g. a tiny nonzero value once everything is removed)
 * For closed integer intervals [l, r] call add(l, r + 1)
 *
 * Example:
 *   IntervalSet<int> s;
 *   s.add(1, 4), s.add(6, 9), s.add(4, 6);   // seg = {[1, 9)}
 *   s.remove(3, 5);                          // seg = {[1, 3), [5, 9)}, covered_length() = 6
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T = long long>
struct IntervalSet{
    map<T, T> seg;
    T total = 0;

    void add(T l, T r){
        if (!(l < r)) return;

        auto it = seg.upper_bound(l);
        if (it != seg.begin() && !(prev(it)->second < l)) --it;

        while (it != seg.end() && !(r < it->first)){
            l = min(l, it->first);
            r = max(r, it->second);
            total -= it->second - it->first;
            it = seg.erase(it);
        }

        seg.emplace_hint(it, l, r);
        total += r - l;
    }

    bool contains(T x) const{
        auto it = seg.upper_bound(x);
        return it != seg.begin() && x < prev(it)->second;
    }

    bool covers(T l, T r) const{
        if (!(l < r)) return true;
        auto it = seg.upper_bound(l);
        return it != seg.begin() && !(prev(it)->second < r);
    }

    T covered_length() const{
        return total;
    }

    void remove(T l, T r){
        if (!(l < r)) return;

        auto it = seg.upper_bound(l);
        if (it != seg.begin() && l < prev(it)->second) --it;

        while (it != seg.end() && it->first < r){
            T a = it->first, b = it->second;
            total -= b - a;
            it = seg.erase(it);

            if (a < l) seg.emplace(a, l), total += l - a;
            if (r < b) seg.emplace(r, b), total += b - r;
        }
    }
};

int main(){
    IntervalSet<int> s;
    assert(s.covered_length() == 0 && !s.contains(0) && s.covers(5, 5) && !s.covers(0, 1));

    s.add(1, 4), s.add(6, 9);
    assert((s.seg == map<int, int>{{1, 4}, {6, 9}}));
    assert(s.covered_length() == 6);
    assert(s.contains(1) && s.contains(3) && !s.contains(4) && !s.contains(0) && s.contains(8) && !s.contains(9));
    assert(s.covers(1, 4) && s.covers(2, 3) && !s.covers(3, 7) && !s.covers(0, 2));

    s.add(4, 6);
    assert((s.seg == map<int, int>{{1, 9}}));
    assert(s.covered_length() == 8 && s.covers(1, 9));

    s.add(5, 5), s.add(7, 3), s.remove(6, 6);
    assert((s.seg == map<int, int>{{1, 9}}));

    s.remove(3, 5);
    assert((s.seg == map<int, int>{{1, 3}, {5, 9}}));
    assert(s.covered_length() == 6 && !s.contains(3) && !s.contains(4) && s.contains(5));

    s.add(-10, -7), s.add(20, 30);
    s.remove(-8, 25);
    assert((s.seg == map<int, int>{{-10, -8}, {25, 30}}));
    assert(s.covered_length() == 7);

    s.add(-9, 26);
    assert((s.seg == map<int, int>{{-10, 30}}));
    s.remove(-10, 30);
    assert(s.seg.empty() && s.covered_length() == 0);

    IntervalSet<long long> big;
    big.add(-4000000000000000000LL, 0), big.add(0, 4000000000000000000LL);
    assert(big.seg.size() == 1 && big.covered_length() == 8000000000000000000LL);
    big.remove(-1, 1);
    assert(big.covered_length() == 7999999999999999998LL && !big.contains(0) && big.contains(1));

    IntervalSet<double> reals;
    reals.add(0.5, 1.75), reals.add(1.25, 2.0), reals.add(3.0, 3.5);
    assert(reals.seg.size() == 2 && abs(reals.covered_length() - 2.0) < 1e-9);
    assert(reals.contains(1.9) && !reals.contains(2.0) && reals.covers(0.5, 2.0) && !reals.covers(1.5, 3.25));

    return 0;
}
