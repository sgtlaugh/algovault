/***
 *
 * Numerical Integration
 * Definite integral of f over [a, b] with composite Simpson, adaptive Simpson or Romberg
 *
 * Complexity:
 *   - adaptive_simpson: O(panels it needs) evaluations, recursion depth <= 44
 *   - romberg: 2^k + 1 evaluations for k levels, 8 <= k <= max_levels
 *   - simpson: 2n + 1 evaluations
 *
 * adaptive_simpson(a, b, f, eps = 1e-9): splits each of 64 equal start panels until its halves agree with it
 * romberg(a, b, f, eps = 1e-9, max_levels = 20): Richardson extrapolation of trapezoid sums on 2^k subintervals
 * simpson(a, b, f, n = 1000): composite Simpson on 2n equal subintervals for n >= 1, exact for cubics
 *
 * Notes:
 *   - f must be finite at every point of [a, b], endpoints included, and a > b gives the negated integral
 *   - eps is absolute and judged by the method's own error estimate, so a coincidental cancellation stops early
 *   - On random integrands the error exceeded eps + 1e-13 relative in about 1 in 10^5 calls, by at most 6x
 *   - adaptive_simpson and romberg also stop at 1e-14 of the integral of |f|, so an eps below rounding still terminates
 *   - adaptive_simpson refines only around kinks, jumps and infinite-slope endpoints like sqrt at 0
 *   - adaptive_simpson can miss a spike, or two jumps, within one start panel of width (b - a) / 64
 *   - romberg reads periodic f as constant when (b - a) / 256 is a whole number of periods
 *   - romberg returns its last estimate after max_levels, it is fast only on smooth f: use adaptive_simpson for kinks
 *   - simpson errs by (b - a) h^4 max|f''''| / 180 with h = (b - a) / 2n
 *   - Nest calls for multiple integrals:
 *         adaptive_simpson(0, 1, [](double x){ return adaptive_simpson(0, x, [&](double y){ return x * y; }); });
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template<typename F>
double adaptive_simpson_panel(F& f, double a, double b, double fa, double fm, double fb, double whole, double eps, double& noise, int depth){
    const int MAX_DEPTH = 44;
    double m = (a + b) / 2, flm = f((a + m) / 2), frm = f((m + b) / 2);
    double left = (m - a) / 6 * (fa + 4 * flm + fm), right = (b - m) / 6 * (fm + 4 * frm + fb);
    double delta = left + right - whole;

    /// The start grid can miss where |f| is large, so each panel's integral of |f| lifts the floor over its rounding
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

template<typename F>
double simpson(double a, double b, F f, int n = 1000){
    double h = (b - a) / (2 * n), sum = f(a) + f(b);
    for (int i = 1; i < 2 * n; i++) sum += f(a + i * h) * (i & 1 ? 4 : 2);
    return sum * h / 3;
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
