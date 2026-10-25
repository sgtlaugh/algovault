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
 * op must be associative, commutative and idempotent, and every strict change must lose a bit or a factor
 * gcd, bitwise or and bitwise and qualify, sum does not
 * min and max stay correct but give up to n segments per push
 * The segments cover l = 0..i, at most floor(log2 A) + 2 of them for gcd and B + 1 for or/and
 * Subarrays per value, as in CF 475D: after each push, cnt[s.value] += s.right - s.left + 1 for every segment
 *
 * auto gcd_op = [](long long x, long long y){ return gcd(x, y); };
 * DistinctSubarrayAggregates<long long, decltype(gcd_op)> agg(gcd_op);
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
    int t;

    scanf("%d", &t);
    auto gcd_op = [](long long x, long long y){ return __gcd(x, y); };
    while (t--){
        int n;
        scanf("%d", &n);
        DistinctSubarrayAggregates<long long, decltype(gcd_op)> agg(gcd_op);
        long long res = 0;
        for (int i = 0; i < n; i++){
            long long x;
            scanf("%lld", &x);
            agg.push(x);
            for (auto& s : agg.segments()) res = max(res, s.value * (i - s.left + 1));
        }

        printf("%lld\n", res);
    }
    return 0;
}
