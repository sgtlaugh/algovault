/***
 *
 * Interval Set
 * Set of disjoint half-open intervals [l, r) that merges on add, splits on remove and answers coverage queries
 *
 * Complexity: amortized O(log n) per add and remove, O(log n) per query, n = number of stored intervals
 *
 * IntervalSet<T> s: empty set, T is long long by default
 * s.add(l, r): cover [l, r), merging with every stored interval it overlaps or touches
 * s.remove(l, r): uncover [l, r), cutting stored intervals that stick out of it
 * s.contains(x): true if the point x is covered
 * s.covered_length(): total covered length, kept as a running sum
 * s.covers(l, r): true if all of [l, r) is covered, always true for an empty range
 * s.seg: map l -> r of the stored intervals in increasing order, never two that touch
 *
 * Empty ranges, l >= r, are ignored by add and remove
 * For closed integer intervals [l, r] call add(l, r + 1)
 * The total covered length must fit in T, long long coordinates in [-4e18, 4e18] are safe
 * With floating point T the running sum drifts, so covered_length() can end slightly off
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

    bool contains(T x) const{
        auto it = seg.upper_bound(x);
        return it != seg.begin() && x < prev(it)->second;
    }

    T covered_length() const{
        return total;
    }

    bool covers(T l, T r) const{
        if (!(l < r)) return true;
        auto it = seg.upper_bound(l);
        return it != seg.begin() && !(prev(it)->second < r);
    }
};

int main(){
    /// Add, merge and split
    IntervalSet<int> s;
    s.add(1, 4);
    s.add(6, 9);
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
