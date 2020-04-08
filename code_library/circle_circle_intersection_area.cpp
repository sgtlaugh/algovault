#include <stdio.h>
#include <bits/stdtr1c++.h>

using namespace std;


const long double EPS = 1e-9;
const long double PI  = 2.0 * acos(0.0);

struct Point{
    long double x, y;

    Point() {}
    Point(long double x, long double y) : x(x), y(y) {}

    long double get_dis(const Point& other){
        long double d1 = x - other.x, d2 = y - other.y;
        return sqrt(d1 * d1 + d2 * d2);
    }
};

struct Circle{
    Point centre;
    long double radius;

    Circle(){}
    Circle(Point centre, long double radius) : centre(centre), radius(radius) {}

    long double get_area(){
        return PI * radius * radius;
    }
};

/// return true if circle B is inside of circle A
bool is_inside(Circle A, Circle B){
    long double x = A.centre.get_dis(B.centre) + B.radius;
    return (abs(x - A.radius) < EPS) || (x < A.radius);
}

/// returns the intersection area of circle A and circle B
long double intersection_area(Circle A, Circle B){
    if (is_inside(A, B)) return B.get_area();
    if (is_inside(B, A)) return A.get_area();

    long double x = A.radius + B.radius;
    long double y = A.centre.get_dis(B.centre);
    if ((abs(x - y) < EPS) || (y > x)) return 0.0;

    long double c = y, r0 = A.radius, r1 = B.radius;
    x = (r1 * r1 - r0 * r0 + c * c) / (2.0 * r1 * c);
    long double CBD = acos(x) * 2.0;
    y = (r0 * r0 - r1 * r1 + c * c) / (2.0 * r0 * c);
    long double CAD = acos(y) * 2.0;

    long double res = (0.5 * CBD * r1 * r1) - (0.5 * r1 * r1 * sin(CBD)) + (0.5 * CAD * r0 * r0) - (0.5 * r0 * r0 * sin(CAD));
    return res;
}

int main(){
    Circle A = Circle(Point(0, 0), 5);
    Circle B = Circle(Point(4, 0), 7);
    printf("%0.5f\n", (double)intersection_area(A, B));   /// 62.69241

    B = Circle(Point(10, 0), 5);
    printf("%0.5f\n", (double)intersection_area(A, B));   /// 0.00000

    B = Circle(Point(0, 0), 10);
    printf("%0.5f\n", (double)intersection_area(A, B));   /// 78.53982
    printf("%0.5f\n", (double)intersection_area(B, A));   /// 78.53982

    return 0;
}
