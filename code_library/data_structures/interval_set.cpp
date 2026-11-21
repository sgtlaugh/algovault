/***
 *
 * Interval Set and Map
 * Disjoint half-open intervals [l, r): a set that merges and splits on cover changes, and a map that assigns values
 *
 * Complexity, n = number of stored intervals or pieces:
 *   - IntervalSet: amortized O(log n) per add and remove, O(log n) per query
 *   - IntervalMap assign: amortized O(log n), it removes every piece it covers
 *   - IntervalMap for_each: O(log n + k) for the k pieces it visits, get: O(log n)
 *
 * IntervalSet<T> s: empty set, T is long long by default
 * s.add(l, r): cover [l, r), merging with every stored interval it overlaps or touches
 * s.remove(l, r): uncover [l, r), cutting stored intervals that stick out of it
 * s.contains(x): true if the point x is covered
 * s.covered_length(): total covered length, kept as a running sum
 * s.covers(l, r): true if all of [l, r) is covered, always true for an empty range
 * s.seg: map l -> r of the stored intervals in increasing order, never two that touch
 *
 * IntervalMap<K, V> m(lo, hi, v): every point of [lo, hi) holds v, lo < hi, K and V are long long by default
 * m.assign(l, r, v): every point of [l, r) now holds v, lo <= l, r <= hi
 * m.for_each(l, r, f): calls f(a, b, value) for each piece [a, b) of [l, r) in order, value is a reference f may change
 * m.get(x): the value at the point x, lo <= x < hi
 * m.seg: map l -> value of the piece [l, next key), its last key is hi and starts no piece
 *
 * Notes:
 *   - Empty ranges, l >= r, are ignored by add, remove, assign and for_each
 *   - For closed integer intervals [l, r] call add(l, r + 1)
 *   - The total covered length must fit in T, long long coordinates in [-4e18, 4e18] are safe
 *   - With floating point T the running sum drifts, so covered_length() can end slightly off
 *   - IntervalMap is fast only while assignments keep the piece count low, as on random data (Chtholly tree)
 *   - Without assignments the pieces only multiply, so for_each over many of them costs O(n) each time
 *   - IntervalMap never merges equal neighbouring pieces, and V must be default constructible
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

template <typename K = long long, typename V = long long>
struct IntervalMap{
    map<K, V> seg;

    IntervalMap(K lo, K hi, V v){
        assert(lo < hi);
        seg[lo] = v, seg[hi] = V();
    }

    void assign(K l, K r, V v){
        if (!(l < r)) return;
        auto end = split(r), it = split(l);
        seg.erase(it, end);
        seg.emplace_hint(end, l, v);
    }

    template <typename F>
    void for_each(K l, K r, F f){
        if (!(l < r)) return;
        auto end = split(r);
        for (auto it = split(l); it != end; it++) f(it->first, next(it)->first, it->second);
    }

    V get(K x) const{
        return prev(seg.upper_bound(x))->second;
    }

    /// Makes x a key without erasing anything, so iterators from earlier splits stay valid
    typename map<K, V>::iterator split(K x){
        assert(!(x < seg.begin()->first) && !(prev(seg.end())->first < x));
        auto it = prev(seg.upper_bound(x));
        if (!(it->first < x)) return it;
        return seg.emplace_hint(next(it), x, it->second);
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

    /// Assign then add 1 on [4, 7): pieces 0 on [0, 2), 7 on [2, 4), 8 on [4, 5), 1 on [5, 7), 0 on [7, 10)
    IntervalMap<int, int> m(0, 10, 0);
    m.assign(2, 5, 7);
    m.for_each(4, 7, [](int, int, int& v){ v += 1; });
    assert(m.get(4) == 8);
    assert(m.get(6) == 1);
    assert(m.get(9) == 0);

    long long total = 0;
    m.for_each(0, 10, [&](int a, int b, int v){ total += (long long)(b - a) * v; });
    assert(total == 24);

    return 0;
}
