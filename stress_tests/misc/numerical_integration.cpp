#include "../common.h"

#define main library_main
#include "../../code_library/misc/numerical_integration.cpp"
#undef main

struct Family {
    function<double(double)> f, antiderivative;
    double lo, hi, abs_bound;
    bool smooth;
};

double rand_real(double lo, double hi){
    return lo + (hi - lo) * (stress::rand_int(0, 1LL << 52) / double(1LL << 52));
}

double poly_eval(const vector<double>& c, double x){
    double res = 0;
    for (int i = (int)c.size() - 1; i >= 0; i--) res = res * x + c[i];
    return res;
}

double poly_antiderivative(const vector<double>& c, double x){
    double res = 0;
    for (int i = (int)c.size() - 1; i >= 0; i--) res = (res + c[i] / (i + 1)) * x;
    return res;
}

vector<double> random_poly(int degree, double bound){
    vector<double> c(degree + 1);
    for (auto& v: c) v = rand_real(-bound, bound);
    return c;
}

/// Random integrands with closed-form antiderivatives: abs_bound bounds |f| on [lo, hi], scaling the tolerances
Family random_family(){
    int kind = stress::rand_int(0, 7);
    double scale = pow(10.0, stress::rand_int(-6, 12));

    if (kind == 0){
        auto c = random_poly(stress::rand_int(0, 8), 5);
        for (auto& v: c) v *= scale;
        double bound = 0;
        for (int i = 0; i < (int)c.size(); i++) bound += abs(c[i]) * pow(3, i);
        return {[=](double x){ return poly_eval(c, x); }, [=](double x){ return poly_antiderivative(c, x); }, -3, 3, bound, true};
    }
    if (kind == 1){
        double k = rand_real(-3, 3);
        if (abs(k) < 1e-3) k = 1;
        return {[=](double x){ return scale * exp(k * x); }, [=](double x){ return scale * exp(k * x) / k; }, -3, 3, scale * exp(3 * abs(k)), true};
    }
    if (kind == 2){
        double w = rand_real(0.1, 20), phase = rand_real(0, 7);
        return {[=](double x){ return scale * sin(w * x + phase); }, [=](double x){ return -scale * cos(w * x + phase) / w; }, -3, 3, scale, true};
    }
    if (kind == 3){
        double c = rand_real(-2, 2), s = rand_real(0.3, 2);
        return {[=](double x){ return scale / (1 + (x - c) * (x - c) / (s * s)); },
                [=](double x){ return scale * s * atan((x - c) / s); }, -3, 3, scale, true};
    }
    if (kind == 4){
        double base = rand_real(-3, 2);
        return {[=](double x){ return scale * sqrt(max(0.0, x - base)); },
                [=](double x){ return scale * 2.0 / 3 * pow(max(0.0, x - base), 1.5); }, base, 3, scale * sqrt(3 - base), false};
    }
    if (kind == 5){
        double c = rand_real(-3, 3);
        return {[=](double x){ return scale * abs(x - c); },
                [=](double x){ return scale * (x - c) * abs(x - c) / 2; }, -3, 3, scale * 6, false};
    }

    if (kind == 6){
        double c = rand_real(-3, 3), left = rand_real(-5, 5), right = rand_real(-5, 5);
        return {[=](double x){ return scale * (x < c ? left : right); },
                [=](double x){ return scale * (x < c ? left * (x - c) : right * (x - c)); }, -3, 3, scale * 5, false};
    }

    /// floor(k x) has up to 6k jumps, each needing its own deep refinement: int_0^t floor = n t - n (n + 1) / 2, n = floor(t)
    /// k <= 10 keeps the jumps further apart than a start panel ((b - a) / 64), two jumps in one panel can cancel
    double k = rand_real(1, 10);
    return {[=](double x){ return scale * floor(k * x); },
            [=](double x){ double n = floor(k * x); return scale * (n * k * x - n * (n + 1) / 2) / k; }, -3, 3, scale * 3 * k + scale, false};
}

/// Every integrand is checked against its closed-form antiderivative: Simpson exactly on cubics and with
/// the exact error term (b - a) h^4 / 180 * 24 c4 on quartics, adaptive Simpson and Romberg within eps and
/// stopping at the first check where they are exact (adaptive on cubics, Romberg up to degree 9)
/// (their estimates can be fooled, so a rare miss up to 1000 eps is allowed, a frequent one is not)
int main(){
    for (int c0 = -1; c0 <= 1; c0++){
        for (int c1 = -1; c1 <= 1; c1++){
            for (int c2 = -1; c2 <= 1; c2++){
                for (int c3 = -1; c3 <= 1; c3++){
                    for (int a = -2; a <= 2; a++){
                        for (int b = -2; b <= 2; b++){
                            vector<double> c = {double(c0), double(c1), double(c2), double(c3)};
                            auto f = [&](double x){ return poly_eval(c, x); };
                            double exact = poly_antiderivative(c, b) - poly_antiderivative(c, a);
                            for (int n = 1; n <= 3; n++) assert(abs(simpson(a, b, f, n) - exact) < 1e-12);
                            long long calls = 0;
                            assert(abs(adaptive_simpson(a, b, [&](double x){ calls++; return f(x); }) - exact) < 1e-12);
                            assert(calls == 257);
                            assert(abs(romberg(a, b, f) - exact) < 1e-12);
                        }
                    }
                }
            }
        }
    }

    /// R(k, k) is exact up to degree 2k + 1, so on degree <= 9 Romberg stops at its first check (k = 8, 257 evaluations)
    for (long long it = 0; it < stress::scaled(3000); it++){
        auto c = random_poly(stress::rand_int(0, 9), 5);
        double a = rand_real(-2, 2), b = rand_real(-2, 2);
        long long calls = 0;
        double res = romberg(a, b, [&](double x){ calls++; return poly_eval(c, x); });
        assert(calls == 257 && abs(res - (poly_antiderivative(c, b) - poly_antiderivative(c, a))) < 1e-9);
    }

    for (long long it = 0; it < stress::scaled(3000); it++){
        auto c = random_poly(4, 10);
        double a = rand_real(-10, 10), b = rand_real(-10, 10);
        int n = stress::rand_int(1, 60);
        if (it % 2) c[4] = 0;

        double h = (b - a) / (2 * n), m = max({1.0, abs(a), abs(b)});
        double exact = poly_antiderivative(c, b) - poly_antiderivative(c, a) + c[4] * (b - a) * pow(h, 4) * 24 / 180;
        double bound = abs(c[0]) + abs(c[1]) * m + abs(c[2]) * m * m + abs(c[3]) * pow(m, 3) + abs(c[4]) * pow(m, 4);
        double tol = 1e-13 * (abs(b - a) * n + m) * bound;
        assert(abs(simpson(a, b, [&](double x){ return poly_eval(c, x); }, n) - exact) <= tol);
    }

    long long iterations = stress::scaled(4000), misses = 0;
    for (long long it = 0; it < iterations; it++){
        Family fam = random_family();
        double a = rand_real(fam.lo, fam.hi), b = rand_real(fam.lo, fam.hi);
        if (it % 5 == 0) a = fam.lo;
        /// the default absolute eps against integrands up to ~1e13 only terminates through the rounding floor
        double eps = it % 4 ? pow(10.0, -stress::rand_int(4, 11)) * fam.abs_bound : 1e-9;
        double fa = fam.antiderivative(a), fb = fam.antiderivative(b), exact = fb - fa;

        /// the last term is the reference's own cancellation error when F(a) and F(b) are close and large
        double tol = eps + 1e-13 * fam.abs_bound * abs(b - a) + 1e-15 * (abs(fa) + abs(fb));

        long long calls = 0;
        auto counted = [&](double x){ assert(++calls <= 100000); return fam.f(x); };
        double err = abs(adaptive_simpson(a, b, counted, eps) - exact);
        assert(err <= 1000 * tol);
        misses += err > tol;

        if (!fam.smooth) continue;
        calls = 0;
        err = abs(romberg(a, b, counted, eps) - exact);
        assert(err <= 1000 * tol);
        misses += err > tol;
    }
    assert(misses <= 2 + iterations / 1000);

    /// The delta / 15 Richardson term keeps the adaptive error far below eps on smooth f (worst 0.0034 eps over
    /// 1.2e5 samples, plain Simpson on the halves reaches 0.1 eps)
    for (long long it = 0; it < stress::scaled(4000); it++){
        Family fam = random_family();
        if (!fam.smooth) continue;
        double a = rand_real(fam.lo, fam.hi), b = rand_real(fam.lo, fam.hi), eps = 1e-6 * fam.abs_bound;
        double exact = fam.antiderivative(b) - fam.antiderivative(a);
        assert(abs(adaptive_simpson(a, b, fam.f, eps) - exact) <= eps / 100);
    }

    /// c + cos(2 pi p t) over p whole periods reads as constant to samples spaced a whole number of periods apart:
    /// Romberg's first check (2^8 + 1 samples) and adaptive Simpson's first split (257 samples) both see it for p < 256
    for (int p = 1; p < 256; p++){
        double a = rand_real(-3, 3), b = a + rand_real(0.5, 6), c = rand_real(-2, 2), phase = rand_real(0, 7);
        auto f = [&](double x){ return c + cos(2 * acos(-1.0) * p * (x - a) / (b - a) + phase); };
        assert(abs(romberg(a, b, f) - c * (b - a)) <= 1e-8);
        assert(abs(adaptive_simpson(a, b, f) - c * (b - a)) <= 1e-8);
    }

    /// Integrands that vanish on the whole start grid but are huge between its points: a rounding floor taken from the
    /// start grid alone is 0 and the panels recurse toward depth 44 (5.7e6 calls at scale 1e6, over 2e8 at 1e8); with
    /// the floor raised by the panels it costs no more than f shifted off the grid (2.7e5 against 4.9e5 calls, m = 384)
    for (long long it = 0; it < stress::scaled(300); it++){
        double scale = pow(10.0, stress::rand_int(0, 14)), exact;
        long long calls = 0;
        function<double(double)> f;
        if (it % 2){
            double c = rand_real(0.01, 0.99), w = rand_real(0.002, 0.005);
            f = [=](double x){ double t = (x - c) / w; return abs(t) < 1 ? scale * (1 - t * t) * (1 - t * t) : 0.0; };
            exact = scale * w * 16 / 15;
        }
        else{
            int m = 128 * (2 * stress::rand_int(0, 1) + 1);
            f = [=](double x){ double v = sin(m * acos(-1.0) * x); return scale * v * v; };
            exact = scale / 2;
        }
        double res = adaptive_simpson(0, 1, [&](double x){ assert(++calls <= 1000000); return f(x); });
        assert(abs(res - exact) <= 1e-9 + 1e-12 * exact);
    }

    for (long long it = 0; it < stress::scaled(100); it++){
        int i = stress::rand_int(0, 4), j = stress::rand_int(0, 4);
        double x0 = rand_real(-2, 2), x1 = rand_real(-2, 2), y0 = rand_real(-2, 2), y1 = rand_real(-2, 2);
        double res = adaptive_simpson(x0, x1, [&](double x){
            return adaptive_simpson(y0, y1, [&](double y){ return pow(x, i) * pow(y, j); });
        });
        double exact = (pow(x1, i + 1) - pow(x0, i + 1)) / (i + 1) * (pow(y1, j + 1) - pow(y0, j + 1)) / (j + 1);
        assert(abs(res - exact) <= 1e-7);
    }

    return 0;
}
