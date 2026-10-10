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
    s.add(1, 4), s.add(6, 9);
    assert((s.seg == map<int, int>{{1, 4}, {6, 9}}));
    assert(s.covered_length() == 6);
    assert(s.contains(3));
    assert(!s.contains(4));   /// half-open: 4 is not in [1, 4)
    assert(!s.covers(3, 7));  /// the gap [4, 6) is uncovered

    s.add(4, 6);  /// touches both sides, so all three merge
    assert((s.seg == map<int, int>{{1, 9}}));

    s.remove(3, 5);
    assert((s.seg == map<int, int>{{1, 3}, {5, 9}}));
    assert(s.covered_length() == 6);
    return 0;
}
