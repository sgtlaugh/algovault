/***
 *
 * Continued Fractions
 * Continued fraction expansion, convergents, best rational approximation and fraction binary search
 *
 * Complexity: O(log q) for the expansion, convergents and best approximation, O(log N) predicate calls for the search
 *
 * continued_fraction(p, q): [a0; a1, ..., ak] with p / q = a0 + 1 / (a1 + 1 / (... + 1 / ak))
 *     a0 = floor(p / q) can be negative, a1..ak >= 1 and ak >= 2 when k >= 1, so the expansion is unique
 *     p / q need not be reduced, requires q >= 1, any p works
 *
 * convergents(a): the fractions [a0; a1, ..., ai] for every i, reduced, the last one is p / q reduced
 *     Requires a non-empty expansion, never overflows on continued_fraction(p, q) for p > LLONG_MIN
 *
 * best_approximation(p, q, N): the fraction closest to p / q among all fractions with denominator <= N
 *     Ties go to the smaller denominator, then to the smaller value: (1, 2, 1) gives 0 / 1, (5, 12, 3) gives 1 / 2
 *     Requires |p| <= 1e18, 1 <= q <= 1e18, 1 <= N <= 1e18 and __int128
 *
 * fraction_search(f, N): the smallest fraction p / q with 0 <= p <= N and 1 <= q <= N for which f(p, q) is true
 *     f must be monotone over the non-negative rationals: false up to some point, true from there on
 *     Returns {1, 0} if f is false for every such fraction
 *     f is only called with reduced fractions inside the bounds, at most 4 log2(N) + 10 times, N <= 1e18
 *     f(p, q) = p * p >= 2 * q * q with N = 10 gives 10 / 7, the smallest fraction >= sqrt(2) with p, q <= 10
 *
 * Fractions are {numerator, denominator} pairs, the denominator is always positive
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

vector<long long> continued_fraction(long long p, long long q){
    assert(q > 0);
    vector<long long> a;

    while (q){
        long long d = p / q, r = p % q;
        if (r < 0) d--, r += q;
        a.push_back(d);
        p = q, q = r;
    }

    return a;
}

vector<pair<long long, long long>> convergents(const vector<long long>& a){
    assert(!a.empty());
    vector<pair<long long, long long>> res;
    long long pp = 1, pq = 0, cp = a[0], cq = 1;
    res.push_back({cp, cq});

    for (size_t i = 1; i < a.size(); i++){
        tie(pp, pq, cp, cq) = make_tuple(cp, cq, a[i] * cp + pp, a[i] * cq + pq);
        res.push_back({cp, cq});
    }

    return res;
}

pair<long long, long long> best_approximation(long long p, long long q, long long n){
    assert(q > 0 && n > 0);
    auto a = continued_fraction(p, q);
    long long pp = 1, pq = 0, cp = a[0], cq = 1;

    for (size_t i = 1; i < a.size(); i++){
        if (a[i] <= (n - pq) / cq){
            tie(pp, pq, cp, cq) = make_tuple(cp, cq, a[i] * cp + pp, a[i] * cq + pq);
            continue;
        }

        /// The convergent and the largest semiconvergent are the neighbours of p / q among denominators <= n
        long long t = (n - pq) / cq, sp = pp + t * cp, sq = pq + t * cq;
        __int128 dc = (__int128)p * cq - (__int128)cp * q, ds = (__int128)p * sq - (__int128)sp * q;
        if (dc < 0) dc = -dc;
        if (ds < 0) ds = -ds;

        __int128 lhs = dc * sq, rhs = ds * cq;
        if (lhs != rhs) return lhs < rhs ? make_pair(cp, cq) : make_pair(sp, sq);
        if (cq != sq) return cq < sq ? make_pair(cp, cq) : make_pair(sp, sq);
        return (__int128)cp * sq < (__int128)sp * cq ? make_pair(cp, cq) : make_pair(sp, sq);
    }

    return {cp, cq};
}

template<typename F>
pair<long long, long long> fraction_search(F f, long long n){
    assert(n > 0);
    if (f(0LL, 1LL)) return {0, 1};

    /// Invariant: f(lp / lq) is false, f(rp / rq) is true, with 1 / 0 standing for infinity
    long long lp = 0, lq = 1, rp = 1, rq = 0;
    while (lp + rp <= n && lq + rq <= n){
        bool side = f(lp + rp, lq + rq);
        long long& bp = side ? rp : lp;
        long long& bq = side ? rq : lq;
        long long dp = side ? lp : rp, dq = side ? lq : rq;

        auto valid = [&](long long k){
            if (dp && k > (n - bp) / dp) return false;
            if (dq && k > (n - bq) / dq) return false;
            return f(bp + k * dp, bq + k * dq) == side;
        };

        long long k = 1, step = 1;
        while (valid(k + step)) k += step, step *= 2;
        for (step /= 2; step; step /= 2){
            if (valid(k + step)) k += step;
        }

        bp += k * dp, bq += k * dq;
    }

    return {rp, rq};
}

int main(){
    long long p, q;

    while (scanf("%lld %lld", &p, &q) == 2){
        auto a = continued_fraction(p, q);
        printf("[%lld", a[0]);
        for (size_t i = 1; i < a.size(); i++) printf("%c%lld", i == 1 ? ';' : ',', a[i]);
        puts("]");
    }
    return 0;
}
