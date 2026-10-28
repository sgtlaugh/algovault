/***
 *
 * Ternary Search and Golden-Section Search
 * Optimum of a unimodal function: integer argmax, real argmin
 *
 * Complexity: ternary_search_max O(log(hi - lo)) calls of f (2 per step)
 *             golden_section_min O(log((hi - lo) / tol)) calls of f (1 per step, plus 1 per 70+ steps for drift)
 *
 * ternary_search_max(lo, hi, f): smallest i in [lo, hi] maximizing f(i), for any comparable return type
 *   Requires f(lo) < ... < f(p) >= f(p + 1) >= ... >= f(hi): strictly increasing up to the peak p, then
 *   non-increasing, so plateaus are allowed at the peak and on the falling side but not on the rising side
 *   (0, 0, 5, 0 is not searchable in O(log) by any method). -4e18 <= lo <= hi <= 4e18
 *   Plateaus on the rising side instead, f(lo) <= ... <= f(p) > ... > f(hi): search g(i) = f(lo + hi - i)
 *   and map back with i -> lo + hi - i, which gives the largest argmax
 *   Minimum: search -f (strictly decreasing, then non-decreasing)
 *
 * golden_section_min(lo, hi, f, eps = 1e-9): x in [lo, hi] minimizing a real f, |lo|, |hi| <= 1e300 so hi - lo stays finite
 *   Requires f strictly decreasing, then strictly increasing, with an optional flat bottom
 *   Stops once hi - lo <= eps * max(1, |lo| + |hi|), so eps is relative for large coordinates; keep eps >= 1e-12
 *   Near a smooth minimum f changes by about (x - x*)^2, so double rounding limits x to about 1e-8 relative
 *   accuracy whatever eps is; f(x) itself is accurate to rounding. A kink at c is located to eps * max(1, |c|)
 *   Maximum: minimize -f
 *
 * Example:
 *   long long i = ternary_search_max(0, n - 1, [&](long long i){ return a[i]; });
 *   double x = golden_section_min(-1000, 1000, [](double x){ return 4 + x + 0.3 * x * x; });
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template<typename F>
long long ternary_search_max(long long lo, long long hi, F f){
    assert(lo <= hi);

    /// f(i) < f(i + 1) holds exactly for i < p, so this is a binary search for the first i where it fails
    while (lo < hi){
        long long mid = lo + (hi - lo) / 2;
        if (f(mid) < f(mid + 1)) lo = mid + 1;
        else hi = mid;
    }

    return lo;
}

template<typename F>
double golden_section_min(double lo, double hi, F f, double eps = 1e-9){
    assert(lo <= hi);

    const double r = (sqrt(5.0) - 1) / 2;
    double x1 = hi - r * (hi - lo), x2 = lo + r * (hi - lo);
    double f1 = f(x1), f2 = f(x2);

    while (hi - lo > eps * max(1.0, abs(lo) + abs(hi))){
        /// The reused point's offset from its golden position grows by 1/r per step from rounding, and left alone it
        /// first stalls the shrink and then breaks lo < x1 < x2 < hi, so a comparison discards the minimum
        if (f1 < f2){
            hi = x2, x2 = x1, f2 = f1;
            x1 = hi - r * (hi - lo), f1 = f(x1);
            double ideal = lo + r * (hi - lo);
            if (abs(x2 - ideal) > (hi - lo) / 64) x2 = ideal, f2 = f(x2);
        }
        else{
            lo = x1, x1 = x2, f1 = f2;
            x2 = lo + r * (hi - lo), f2 = f(x2);
            double ideal = hi - r * (hi - lo);
            if (abs(x1 - ideal) > (hi - lo) / 64) x1 = ideal, f1 = f(x1);
        }
    }

    return (lo + hi) / 2;
}

int main(){
    double t;
    int n;

    while (scanf("%lf %d", &t, &n) == 2){
        vector<double> run(n), cycle(n);
        for (int i = 0; i < n; i++) scanf("%lf %lf", &run[i], &cycle[i]);

        auto margin = [&](double r){
            double best = HUGE_VAL;
            for (int i = 0; i + 1 < n; i++) best = min(best, 3600 * (r / run[i] + (t - r) / cycle[i]));
            return 3600 * (r / run[n - 1] + (t - r) / cycle[n - 1]) - best;
        };

        double r = golden_section_min(0, t, margin, 1e-12);
        if (margin(r) <= 1e-6) printf("The cheater can win by %.0f seconds with r = %.2fkm and k = %.2fkm.\n", fabs(margin(r)), r, t - r);
        else puts("The cheater cannot win.");
    }

    return 0;
}
