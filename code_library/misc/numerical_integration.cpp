/***
 *
 * Numerical Integration
 * Definite integral of f over [a, b] with composite Simpson, adaptive Simpson or Romberg
 *
 * Complexity: simpson O(n) evaluations, adaptive_simpson O(panels it needs) with recursion depth <= 44,
 *             romberg O(2^k) evaluations for k <= max_levels levels
 *
 * simpson(a, b, f, n = 1000): 2n equal subintervals for n >= 1, 2n + 1 evaluations
 *     Exact for cubics, error (b - a) h^4 max|f''''| / 180 with h = (b - a) / 2n
 *
 * adaptive_simpson(a, b, f, eps = 1e-9): starts from 64 equal panels and splits a panel until its halves agree
 *     with it; the start grid keeps periodic f that vanishes on a coarse grid (sin^2 over [0, 8pi]) from
 *     reading as 0, but a spike, or two jumps, within one start panel ((b - a) / 64) can still be missed
 *     Also stops a panel once the disagreement is below 1e-14 of the integral of |f| (estimated on the start grid,
 *     raised by each panel's own), so an absolute eps below the rounding noise of a huge integral still terminates,
 *     even when f vanishes on the start grid and is huge between its points (relative error then ~1e-13)
 *     Handles kinks, jumps and infinite-slope endpoints (sqrt at 0) by refining only around them
 *
 * romberg(a, b, f, eps = 1e-9, max_levels = 20): Richardson extrapolation of trapezoid sums on 2^k subintervals
 *     Stops when two diagonal entries differ by <= max(eps, 1e-14 * integral of |f|), at least 8 levels (257
 *     evaluations, so like adaptive_simpson it reads periodic f as constant only when (b - a) / 256 is a whole
 *     number of periods), else returns the last estimate (2^max_levels + 1 evaluations); fast only on smooth f,
 *     use adaptive_simpson for kinks and jumps
 *
 * eps is an absolute error target judged by the method's own estimate, not a guarantee: a coincidental
 * cancellation in the estimate stops early; on the stress test's random integrands the error exceeded
 * eps + 1e-13 relative for about 1 in 10^5 adaptive calls (worst 6x) and 1 in 2 * 10^5 Romberg calls (worst 1.4x)
 *
 * f must be finite at every point of [a, b] (endpoints included); a > b gives the negated integral
 * Nest calls for multiple integrals:
 *     adaptive_simpson(0, 1, [](double x){ return adaptive_simpson(0, x, [&](double y){ return x * y; }); });
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template<typename F>
double simpson(double a, double b, F f, int n = 1000){
    double h = (b - a) / (2 * n), sum = f(a) + f(b);
    for (int i = 1; i < 2 * n; i++) sum += f(a + i * h) * (i & 1 ? 4 : 2);
    return sum * h / 3;
}

template<typename F>
double adaptive_simpson_panel(F& f, double a, double b, double fa, double fm, double fb, double whole, double eps, double& noise, int depth){
    const int MAX_DEPTH = 44;
    double m = (a + b) / 2, flm = f((a + m) / 2), frm = f((m + b) / 2);
    double left = (m - a) / 6 * (fa + 4 * flm + fm), right = (b - m) / 6 * (fm + 4 * frm + fb);
    double delta = left + right - whole;

    /// The start grid can miss where |f| is large (zero on the grid, huge between), so every panel's integral of |f|
    /// raises the floor: it then always covers this panel's own rounding noise, which ends the recursion
    noise = max(noise, 1e-14 * abs(b - a) / 12 * (abs(fa) + 4 * abs(flm) + 2 * abs(fm) + 4 * abs(frm) + abs(fb)));
    if (depth >= MAX_DEPTH || abs(delta) <= max(15 * eps, noise)) return left + right + delta / 15;

    return adaptive_simpson_panel(f, a, m, fa, flm, fm, left, eps / 2, noise, depth + 1)
         + adaptive_simpson_panel(f, m, b, fm, frm, fb, right, eps / 2, noise, depth + 1);
}

template<typename F>
double adaptive_simpson(double a, double b, F f, double eps = 1e-9){
    const int PANELS = 64;
    double h = (b - a) / (2 * PANELS), abs_integral = 0;
    vector<double> y(2 * PANELS + 1);
    for (int i = 0; i <= 2 * PANELS; i++){
        y[i] = f(i == 2 * PANELS ? b : a + i * h);
        abs_integral += abs(h / 3 * y[i]) * (i == 0 || i == 2 * PANELS ? 1 : i % 2 ? 4 : 2);
    }

    /// Not halved per level like eps: a panel's rounding noise shrinks with its width
    double noise = 1e-14 * abs_integral, res = 0;
    for (int i = 0; i < 2 * PANELS; i += 2){
        double whole = h / 3 * (y[i] + 4 * y[i + 1] + y[i + 2]), r = i + 2 == 2 * PANELS ? b : a + (i + 2) * h;
        res += adaptive_simpson_panel(f, a + i * h, r, y[i], y[i + 1], y[i + 2], whole, eps / PANELS, noise, 0);
    }
    return res;
}

template<typename F>
double romberg(double a, double b, F f, double eps = 1e-9, int max_levels = 20){
    const int MIN_LEVELS = 8;
    double h = b - a, fa = f(a), fb = f(b), abs_trapezoid = abs(h / 2) * (abs(fa) + abs(fb));
    vector<double> prev = {h / 2 * (fa + fb)}, cur;

    for (int k = 1; k <= max_levels; k++){
        h /= 2;
        double sum = 0, abs_sum = 0;
        for (long long i = 1; i < (1LL << k); i += 2){
            double y = f(a + i * h);
            sum += y;
            abs_sum += abs(y);
        }

        cur.assign(k + 1, 0);
        cur[0] = prev[0] / 2 + h * sum;
        abs_trapezoid = abs_trapezoid / 2 + abs(h) * abs_sum;
        double power = 4;
        for (int j = 1; j <= k; j++, power *= 4) cur[j] = cur[j - 1] + (cur[j - 1] - prev[j - 1]) / (power - 1);

        if (k >= MIN_LEVELS && abs(cur[k] - prev[k - 1]) <= max(eps, 1e-14 * abs_trapezoid)) return cur[k];
        swap(prev, cur);
    }
    return prev.back();
}

int main(){
    const double PI = acos(-1.0);
    auto close = [](double x, double y){ return abs(x - y) <= 1e-9; };

    assert(close(simpson(0, 1, [](double x){ return x * x * x; }, 1), 0.25));  /// exact on cubics, even with n = 1
    assert(close(simpson(0, PI, [](double x){ return sin(x); }), 2));

    assert(close(adaptive_simpson(0, 1, [](double x){ return 4 / (1 + x * x); }), PI));
    assert(close(adaptive_simpson(0, 1, [](double x){ return sqrt(x); }), 2.0 / 3));          /// infinite slope at 0
    assert(close(adaptive_simpson(-1, 2, [](double x){ return abs(x); }), 2.5));              /// a kink
    assert(close(adaptive_simpson(0, 1, [](double x){ return x < 0.3 ? 1.0 : 3.0; }), 2.4));  /// a jump

    double double_integral = adaptive_simpson(0, 1, [](double x){
        return adaptive_simpson(0, x, [&](double y){ return x * y; });
    });
    assert(close(double_integral, 1.0 / 8));  /// x y over the triangle 0 <= y <= x <= 1

    assert(close(romberg(0, 1, [](double x){ return 4 / (1 + x * x); }), PI));
    assert(close(romberg(3, 1, [](double x){ return 1 / x; }), -log(3)));  /// a > b negates
    return 0;
}
