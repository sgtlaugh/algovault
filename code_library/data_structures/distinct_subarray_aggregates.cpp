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
    auto gcd_op = [](long long x, long long y){ return gcd(x, y); };

    DistinctSubarrayAggregates<long long, decltype(gcd_op)> g(gcd_op);
    assert(g.segments().empty());
    map<long long, long long> cnt;
    vector<Tuples> expected_gcd = {
        {{6, 0, 0}},
        {{2, 0, 0}, {4, 1, 1}},
        {{2, 0, 2}},
        {{1, 0, 2}, {3, 3, 3}},
    };
    vector<long long> a = {6, 4, 2, 3};
    for (int i = 0; i < 4; i++){
        g.push(a[i]);
        assert(as_tuples(g) == expected_gcd[i]);
        for (auto& s : g.segments()) cnt[s.value] += s.right - s.left + 1;
    }
    assert((cnt == map<long long, long long>{{1, 3}, {2, 4}, {3, 1}, {4, 1}, {6, 1}}));

    DistinctSubarrayAggregates<long long, decltype(gcd_op)> zeros(gcd_op);
    zeros.push(0);
    zeros.push(0);
    assert((as_tuples(zeros) == Tuples{{0, 0, 1}}));
    zeros.push(5);
    assert((as_tuples(zeros) == Tuples{{5, 0, 2}}));

    DistinctSubarrayAggregates<long long, decltype(gcd_op)> big(gcd_op);
    big.push(1LL << 60);
    big.push(3LL << 59);
    assert((as_tuples(big) == Tuples{{1LL << 59, 0, 0}, {3LL << 59, 1, 1}}));

    DistinctSubarrayAggregates<long long, bit_or<long long>> ors;
    for (long long x : {1, 2, 4}) ors.push(x);
    assert((as_tuples(ors) == Tuples{{7, 0, 0}, {6, 1, 1}, {4, 2, 2}}));
    ors.push(1);
    assert((as_tuples(ors) == Tuples{{7, 0, 1}, {5, 2, 2}, {1, 3, 3}}));

    DistinctSubarrayAggregates<long long, bit_and<long long>> ands;
    for (long long x : {12, 10, 14}) ands.push(x);
    assert((as_tuples(ands) == Tuples{{8, 0, 0}, {10, 1, 1}, {14, 2, 2}}));
    ands.push(1);
    assert((as_tuples(ands) == Tuples{{0, 0, 2}, {1, 3, 3}}));

    DistinctSubarrayAggregates<unsigned long long, bit_or<unsigned long long>> widest;
    for (int b = 0; b < 64; b++) widest.push(1ULL << b);
    widest.push(0);
    assert(widest.segments().size() == 65);
    assert(widest.segments()[0].value == ~0ULL && widest.segments()[64].value == 0);

    DistinctSubarrayAggregates<long long, decltype(gcd_op)> deepest(gcd_op);
    for (int b = 0; b <= 62; b++) deepest.push(1LL << b);
    deepest.push(0);
    assert(deepest.segments().size() == 64);

    return 0;
}
