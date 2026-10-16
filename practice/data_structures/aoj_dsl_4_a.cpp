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
 * s.covers(l, r): true if all of [l, r) is covered, always true for an empty range
 * s.covered_length(): total covered length, kept as a running sum
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
    int n;
    scanf("%d", &n);
    vector<array<long long, 4>> rects(n);
    vector<long long> xs;
    for (auto& r: rects){
        scanf("%lld %lld %lld %lld", &r[0], &r[1], &r[2], &r[3]);
        xs.push_back(r[0]);
        xs.push_back(r[2]);
    }
    sort(xs.begin(), xs.end());
    xs.erase(unique(xs.begin(), xs.end()), xs.end());

    /// Inside one x strip the covered set is a union of y intervals, rebuilt per strip: O(n^2 log n) for n = 2000
    long long area = 0;
    for (int i = 0; i + 1 < (int)xs.size(); i++){
        IntervalSet<long long> strip;
        for (auto& r: rects){
            if (r[0] <= xs[i] && xs[i + 1] <= r[2]) strip.add(r[1], r[3]);
        }
        area += strip.covered_length() * (xs[i + 1] - xs[i]);
    }

    printf("%lld\n", area);
    return 0;
}
