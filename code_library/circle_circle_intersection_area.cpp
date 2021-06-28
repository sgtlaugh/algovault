/***
 *
 * Calculates the intersection area of two circles
 *
***/

#include <bits/stdc++.h>

using namespace std;

typedef long double fType;

const fType EPS = 1e-9;
const fType PI  = 2 * acos(0);

struct Point{
    fType x, y;

    Point() {}
    Point(fType x, fType y) : x(x), y(y) {}

    fType get_dis(const Point& other){
        fType d1 = x - other.x, d2 = y - other.y;
        return sqrt(d1 * d1 + d2 * d2);
    }
};

struct Circle{
    Point centre;
    fType radius;

    Circle(){}
    Circle(Point centre, fType radius) : centre(centre), radius(radius) {}

    fType get_area(){
        return PI * radius * radius;
    }
};

/// Returns true if circle B is inside of circle A
bool is_inside(Circle A, Circle B){
    fType x = A.centre.get_dis(B.centre) + B.radius;
    return (abs(x - A.radius) < EPS) || (x < A.radius);
}

fType intersection_area(Circle A, Circle B){
    if (is_inside(A, B)) return B.get_area();
    if (is_inside(B, A)) return A.get_area();

    fType x = A.radius + B.radius;
    fType y = A.centre.get_dis(B.centre);
    if ((abs(x - y) < EPS) || (y > x)) return 0;

    fType c = y, r0 = A.radius, r1 = B.radius;
    x = (r1 * r1 - r0 * r0 + c * c) / (2 * r1 * c);
    fType CBD = acos(x) * 2;
    y = (r0 * r0 - r1 * r1 + c * c) / (2 * r0 * c);
    fType CAD = acos(y) * 2;

    fType res = CBD * r1 * r1 - r1 * r1 * sin(CBD) + CAD * r0 * r0 - r0 * r0 * sin(CAD);
    return 0.5 * res;
}

int main(){
    Circle A = Circle(Point(0, 0), 5);
    Circle B = Circle(Point(4, 0), 7);
    assert(abs(intersection_area(A, B) - 62.692413539573) < EPS);

    B = Circle(Point(10, 0), 5);
    assert(abs(intersection_area(A, B) - 0) < EPS);

    B = Circle(Point(0, 0), 10);
    assert(abs(intersection_area(A, B) - 0) < 78.539816339745);
    assert(abs(intersection_area(B, A) - 0) < 78.539816339745);

    return 0;
}
