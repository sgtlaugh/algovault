// LINK: -lquadmath
#include "common.h"
#include <quadmath.h>

#define main library_main
#include "../code_library/circle_circle_intersection_area.cpp"
#undef main

/// Lens area from the segment formula with the kite area taken by Heron's formula, a different route than the library's sines
__float128 lens_area(__float128 x0, __float128 y0, __float128 r0, __float128 x1, __float128 y1, __float128 r1){
    __float128 d = sqrtq((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1)), pi = acosq(-1);
    if (d >= r0 + r1) return 0;
    if (d <= fabsq(r0 - r1)) return pi * fminq(r0, r1) * fminq(r0, r1);

    __float128 a0 = acosq((d * d + r0 * r0 - r1 * r1) / (2 * d * r0)), a1 = acosq((d * d + r1 * r1 - r0 * r0) / (2 * d * r1));
    __float128 kite = sqrtq((-d + r0 + r1) * (d + r0 - r1) * (d - r0 + r1) * (d + r0 + r1)) / 2;
    return r0 * r0 * a0 + r1 * r1 * a1 - kite;
}

/// Midpoint rule over x of the overlap of the two vertical chords, crude but sharing nothing with either formula
long double integrate(const Circle& a, const Circle& b){
    long double lo = max(a.centre.x - a.radius, b.centre.x - b.radius), hi = min(a.centre.x + a.radius, b.centre.x + b.radius), res = 0;
    if (lo >= hi) return 0;
    const int steps = 20000;
    long double h = (hi - lo) / steps;
    for (int i = 0; i < steps; i++){
        long double x = lo + (i + 0.5L) * h;
        long double ha = sqrtl(max(0.0L, a.radius * a.radius - (x - a.centre.x) * (x - a.centre.x)));
        long double hb = sqrtl(max(0.0L, b.radius * b.radius - (x - b.centre.x) * (x - b.centre.x)));
        long double top = min(a.centre.y + ha, b.centre.y + hb), bottom = max(a.centre.y - ha, b.centre.y - hb);
        if (top > bottom) res += (top - bottom) * h;
    }
    return res;
}

int main(){
    /// 2 * acos(0) is the double overload, about 1e-4 off for a contained circle of radius 1e6
    assert(fabsl(PI - (long double)acosq(-1)) <= 1e-18L);
    Circle huge(Point(0, 0), 2e6), inner(Point(1, 1), 1e6);
    assert(fabsl(intersection_area(huge, inner) - (long double)lens_area(0, 0, 2e6, 1, 1, 1e6)) <= 1e-6L);

    for (long long it = 0; it < stress::scaled(200000); it++){
        int range = it % 3 ? 10 : 1000;
        long double x0 = stress::rand_int(-range, range), y0 = stress::rand_int(-range, range), r0 = stress::rand_int(0, range);
        long double x1 = stress::rand_int(-range, range), y1 = stress::rand_int(-range, range), r1 = stress::rand_int(0, range);
        if (it % 4 == 0) x0 += stress::rand_int(0, 999) / 1000.0L, r1 += stress::rand_int(0, 999) / 1000.0L;

        Circle a(Point(x0, y0), r0), b(Point(x1, y1), r1);
        long double res = intersection_area(a, b), swapped = intersection_area(b, a);
        long double expected = lens_area(x0, y0, r0, x1, y1, r1);
        long double scale = max(1.0L, max(r0, r1) * max(r0, r1));

        assert(!std::isnan(res) && res == res);
        assert(fabsl(res - expected) <= 1e-9L * scale);
        assert(fabsl(res - swapped) <= 1e-9L * scale);
        assert(res >= -1e-9L * scale && res <= PI * min(r0, r1) * min(r0, r1) + 1e-9L * scale);

        if (it % 500 == 0) assert(fabsl(res - integrate(a, b)) <= 1e-4L * scale);
    }

    /// Tangent and nearly tangent pairs, where the acos arguments sit right at +-1
    for (long long it = 0; it < stress::scaled(20000); it++){
        long double r0 = stress::rand_int(1, 1000), r1 = stress::rand_int(1, 1000), t = stress::rand_int(-1000, 1000) * 1e-12L;
        long double d = stress::rand_int(0, 1) ? r0 + r1 + t : fabsl(r0 - r1) + t;
        long double angle = stress::rand_int(0, 359) * PI / 180;
        Circle a(Point(0, 0), r0), b(Point(d * cosl(angle), d * sinl(angle)), r1);

        long double res = intersection_area(a, b), scale = max(r0, r1) * max(r0, r1);
        assert(!std::isnan(res));
        assert(fabsl(res - (long double)lens_area(0, 0, r0, b.centre.x, b.centre.y, r1)) <= 1e-6L * scale);
    }
    return 0;
}
