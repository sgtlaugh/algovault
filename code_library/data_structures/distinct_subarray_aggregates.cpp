/***
 *
 * Distinct Subarray Aggregates
 * For every right end i, the distinct values of op(a[l..i]) over all l <= i, each with the range of starts l giving it
 *
 * Complexity: O(n log A) calls of op over n pushes, O(log A) memory
 *
 * op must be associative, commutative and idempotent, and every strict change must lose a bit or a factor:
 * gcd of values in [0, A], bitwise or, bitwise and qualify; sum does not, a value can repeat across segments
 * min and max stay correct but give up to n segments per push
 * Under that condition equal values of op(a[l..i]) form one contiguous range of l, so each value appears once
 *
 * DistinctSubarrayAggregates<T, Op> agg(op);
 * agg.push(x) appends x as a[i], i = number of earlier pushes, starting from 0
 * agg.segments() then lists {value, left, right}: op(a[l..i]) == value exactly for l in [left, right]
 * Segments are sorted by left ascending, cover l = 0..i, adjacent values differ, the last one ends at right = i
 * At most floor(log2 A) + 2 segments for gcd of values in [0, A] with A >= 1, a single segment when every value is 0
 * At most B + 1 segments for or/and on B-bit values
 *
 * Count of subarrays per value (CF 475D): after each push, cnt[s.value] += s.right - s.left + 1 for every segment
 * Longest subarray with a given value: right end i, start s.left of the segment holding that value
 *
 * auto gcd_op = [](long long x, long long y){ return gcd(x, y); };
 * DistinctSubarrayAggregates<long long, decltype(gcd_op)> agg(gcd_op);
 * DistinctSubarrayAggregates<unsigned long long, bit_or<unsigned long long>> ors;
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

    Op op;
    int n = 0;
    vector<Segment> segs;

    explicit DistinctSubarrayAggregates(Op op = Op()) : op(op) {}

    void push(T x){
        for (auto& s : segs) s.value = op(s.value, x);
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
    using Tuples = vector<tuple<long long, int, int>>;
    auto as_tuples = [](const auto& agg){
        Tuples res;
        for (auto& s : agg.segments()) res.emplace_back(s.value, s.left, s.right);
        return res;
    };

    /// Count of subarrays per gcd value of {6, 4, 2, 3}
    auto gcd_op = [](long long x, long long y){ return gcd(x, y); };
    DistinctSubarrayAggregates<long long, decltype(gcd_op)> g(gcd_op);
    map<long long, long long> cnt;
    for (long long x : {6, 4, 2, 3}){
        g.push(x);
        for (auto& s : g.segments()) cnt[s.value] += s.right - s.left + 1;
    }
    assert((as_tuples(g) == Tuples{{1, 0, 2}, {3, 3, 3}}));    /// gcd(a[l..3]) is 1 for l <= 2, 3 for l = 3
    assert((cnt == map<long long, long long>{{1, 3}, {2, 4}, {3, 1}, {4, 1}, {6, 1}}));

    DistinctSubarrayAggregates<long long, bit_or<long long>> ors;
    for (long long x : {1, 2, 4, 1}) ors.push(x);
    assert((as_tuples(ors) == Tuples{{7, 0, 1}, {5, 2, 2}, {1, 3, 3}}));
    return 0;
}
