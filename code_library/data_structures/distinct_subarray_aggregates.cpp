/***
 *
 * Distinct Subarray Aggregates
 * Push values one by one; after each push, lists every distinct gcd/or/and of a subarray ending at the newest value
 *
 * Complexity: O(log A) calls of op per push for gcd of values in [0, A], O(B) for or/and on B-bit values
 *
 * DistinctSubarrayAggregates<T, Op> agg(op): no values yet, folding with op
 * DistinctSubarrayAggregates<T, Op> agg: the same with a default constructed Op, like bit_or<T>
 * agg.push(x): append x as a[i], i = number of earlier pushes
 * agg.segments(): {value, left, right} sorted by left, op(a[l..i]) == value exactly for l in [left, right]
 *
 * Notes:
 *   - op must be associative, commutative and idempotent, and every strict change must lose a bit or a factor
 *   - gcd, bitwise or and bitwise and qualify, sum does not, min and max work but give up to n segments per push
 *   - The segments cover l = 0..i, at most floor(log2 A) + 2 of them for gcd and B + 1 for or/and
 *   - Subarrays per value, as in CF 475D: after each push, cnt[s.value] += s.right - s.left + 1 for every segment
 *   - A lambda op is passed with its type:
 *         auto gcd_op = [](long long x, long long y){ return gcd(x, y); };
 *         DistinctSubarrayAggregates<long long, decltype(gcd_op)> agg(gcd_op);
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T, typename Op>
struct DistinctSubarrayAggregates{
    struct Segment{
        T value;
        int left, right;
    };

    Op combine;
    int n = 0;
    vector<Segment> segs;

    explicit DistinctSubarrayAggregates(Op combine = Op()) : combine(combine) {}

    void push(T x){
        for (auto& s : segs) s.value = combine(s.value, x);
        segs.push_back({x, n, n});
        n++;

        int k = 0;
        for (int j = 1; j < (int)segs.size(); j++){
            if (segs[j].value == segs[k].value) segs[k].right = segs[j].right;
            else segs[++k] = segs[j];
        }
        segs.resize(k + 1);
    }

    const vector<Segment>& segments() const{
        return segs;
    }
};

int main(){
    /// gcd of the subarrays of {6, 4, 2, 3}, counted per value over every push
    auto gcd_op = [](long long x, long long y){ return gcd(x, y); };
    DistinctSubarrayAggregates<long long, decltype(gcd_op)> g(gcd_op);
    map<long long, long long> cnt;
    for (long long x : {6, 4, 2, 3}){
        g.push(x);
        for (auto& s : g.segments()) cnt[s.value] += s.right - s.left + 1;
    }

    auto& segs = g.segments();  /// gcd(a[l..3]) is 1 for l <= 2, 3 for l = 3
    assert(segs.size() == 2);
    assert(segs[0].value == 1);
    assert(segs[0].left == 0);
    assert(segs[0].right == 2);
    assert(segs[1].value == 3);
    assert(segs[1].left == 3);
    assert(segs[1].right == 3);
    assert((cnt == map<long long, long long>{{1, 3}, {2, 4}, {3, 1}, {4, 1}, {6, 1}}));

    /// Bitwise or with a default constructed functor
    DistinctSubarrayAggregates<long long, bit_or<long long>> ors;
    ors.push(1);
    ors.push(2);
    ors.push(4);
    ors.push(1);
    auto& or_segs = ors.segments();  /// a[l..3] or to 7 for l <= 1, 5 for l = 2, 1 for l = 3
    assert(or_segs.size() == 3);
    assert(or_segs[0].value == 7);
    assert(or_segs[0].right == 1);
    assert(or_segs[1].value == 5);
    assert(or_segs[2].value == 1);

    return 0;
}
