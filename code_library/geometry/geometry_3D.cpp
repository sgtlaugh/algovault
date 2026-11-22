/***
 *
 * 3D Geometry
 * Distances between points, lines, segments, planes and triangles in space, line-plane intersection and volume
 *
 * Complexity: O(1) per call, O(number of face vertices) for polyhedron_volume
 *
 * Point3 of doubles with +, -, * and / by a scalar; a triangle is any 3 points a, b, c, a plane 3 non-collinear ones
 * cross(u, v): the cross product u x v
 * dist(p, q): Euclidean distance between two points
 * dist_point_line(p, a, b): distance from p to the line through a != b
 * dist_point_plane(p, a, b, c): signed distance to the plane, positive on the side of cross(b - a, c - a)
 * dist_point_segment(p, a, b): distance from p to segment ab, a == b allowed
 * dist_point_triangle(p, a, b, c): distance from p to the filled triangle abc, degenerate triangles allowed
 * dist_segment_segment(a, b, c, d): distance between segments ab and cd, zero-length ones allowed
 * dist_segment_triangle(p, q, a, b, c): distance between segment pq and the filled triangle abc, about 0 when they meet
 * dot(u, v): the dot product
 * line_plane_intersection(p, q, a, b, c): {1, point} one point, {0, p} parallel, {-1, p} in the plane, p != q
 * norm(u): length of u
 * polyhedron_volume(points, faces): volume of a closed polyhedron, each face counter-clockwise seen from outside
 *
 * Notes:
 *   - Every result is a double; distances hold to 1e-7 at coordinates up to 1e6, nearly parallel segments included
 *   - Near a triangle whose smallest angle has sine s the error grows to about 1e-10 / s, so 5e-7 at s = 1e-4
 *   - A triangle with s below EPS counts as its three edges, so collinear points with rounding noise work
 *   - A line is parallel to the plane when the sine of their angle is below EPS = 1e-9
 *   - Two triangles share a point iff some edge of one is within the tolerance of the other: dist_segment_triangle
 *   - For a convex hull of integer points use convex_hull_3D.cpp, whose volume() is exact
 *   - Defines EPS = 1e-9 as geometry.cpp does: pasted next to it, keep one EPS
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

const double EPS = 1e-9;

struct Point3{
    double x, y, z;

    Point3 operator+(const Point3& p) const{
        return {x + p.x, y + p.y, z + p.z};
    }

    Point3 operator-(const Point3& p) const{
        return {x - p.x, y - p.y, z - p.z};
    }

    Point3 operator*(double k) const{
        return {x * k, y * k, z * k};
    }

    Point3 operator/(double k) const{
        return {x / k, y / k, z / k};
    }
};

Point3 cross(const Point3& u, const Point3& v){
    return {u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x};
}

double dot(const Point3& u, const Point3& v){
    return u.x * v.x + u.y * v.y + u.z * v.z;
}

double norm(const Point3& u){
    return sqrt(dot(u, u));
}

double dist(const Point3& p, const Point3& q){
    return norm(p - q);
}

double dist_point_line(const Point3& p, const Point3& a, const Point3& b){
    return norm(cross(b - a, p - a)) / norm(b - a);
}

double dist_point_plane(const Point3& p, const Point3& a, const Point3& b, const Point3& c){
    Point3 n = cross(b - a, c - a);
    return dot(n, p - a) / norm(n);
}

double dist_point_segment(const Point3& p, const Point3& a, const Point3& b){
    Point3 ab = b - a;
    double len2 = dot(ab, ab), t = len2 > 0 ? dot(p - a, ab) / len2 : 0;
    return dist(p, a + ab * max(0.0, min(1.0, t)));
}

/// The sign tests need no tolerance: a point misjudged near an edge gets that edge's distance from the caller anyway
bool projects_inside(const Point3& p, const Point3& a, const Point3& b, const Point3& c){
    Point3 n = cross(b - a, c - a);
    if (dot(n, n) <= EPS * EPS * dot(b - a, b - a) * dot(c - a, c - a)) return false;  /// sine below EPS: n is noise
    return dot(cross(b - a, p - a), n) >= 0 && dot(cross(c - b, p - b), n) >= 0 && dot(cross(a - c, p - c), n) >= 0;
}

double dist_point_triangle(const Point3& p, const Point3& a, const Point3& b, const Point3& c){
    if (projects_inside(p, a, b, c)) return fabs(dist_point_plane(p, a, b, c));
    return min({dist_point_segment(p, a, b), dist_point_segment(p, b, c), dist_point_segment(p, c, a)});
}

double dist_segment_segment(const Point3& a, const Point3& b, const Point3& c, const Point3& d){
    Point3 u = b - a, v = d - c, r = a - c;
    double uu = dot(u, u), vv = dot(v, v), uv = dot(u, v), ur = dot(u, r), vr = dot(v, r);
    if (uu == 0) return dist_point_segment(a, c, d);
    if (vv == 0) return dist_point_segment(c, a, b);

    /// (u x v).(v x r) equals uv vr - vv ur without its cancellation, which cost 1e-6 on long near-parallel segments
    Point3 w = cross(u, v);
    double den = dot(w, w), s = den > 0 ? max(0.0, min(1.0, dot(w, cross(v, r)) / den)) : 0;
    double t = (uv * s + vr) / vv;
    if (t < 0) t = 0, s = max(0.0, min(1.0, -ur / uu));
    else if (t > 1) t = 1, s = max(0.0, min(1.0, (uv - ur) / uu));

    return dist(a + u * s, c + v * t);
}

double dist_segment_triangle(const Point3& p, const Point3& q, const Point3& a, const Point3& b, const Point3& c){
    Point3 n = cross(b - a, c - a);
    double sp = dot(n, p - a), sq = dot(n, q - a);

    if ((sp < 0 && sq > 0) || (sp > 0 && sq < 0)){
        if (projects_inside(p + (q - p) * (sp / (sp - sq)), a, b, c)) return 0;
    }

    double res = min(dist_point_triangle(p, a, b, c), dist_point_triangle(q, a, b, c));
    return min({res, dist_segment_segment(p, q, a, b), dist_segment_segment(p, q, b, c), dist_segment_segment(p, q, c, a)});
}

pair<int, Point3> line_plane_intersection(const Point3& p, const Point3& q, const Point3& a, const Point3& b, const Point3& c){
    Point3 n = cross(b - a, c - a), d = q - p;
    double den = dot(n, d);
    if (fabs(den) > EPS * norm(n) * norm(d)) return {1, p + d * (dot(n, a - p) / den)};

    bool in_plane = fabs(dist_point_plane(p, a, b, c)) <= EPS * max(norm(p - a), norm(d));
    return {in_plane ? -1 : 0, p};
}

double polyhedron_volume(const vector<Point3>& points, const vector<vector<int>>& faces){
    double res = 0;
    for (const vector<int>& f : faces){
        /// Fanned from the face's first vertex: each term is 6 times a signed tetrahedron against the origin
        for (int i = 1; i + 1 < (int)f.size(); i++) res += dot(points[f[0]], cross(points[f[i]], points[f[i + 1]]));
    }

    return res / 6;
}

int main(){
    /// Unit square at z = 0 split into two triangles; the point (0.5, 0.5, 2) hangs 2 above it
    Point3 a{0, 0, 0}, b{1, 0, 0}, c{1, 1, 0}, p{0.5, 0.5, 2};
    assert(fabs(dist_point_plane(p, a, b, c) - 2) < 1e-12);
    assert(fabs(dist_point_triangle(p, a, b, c) - 2) < 1e-12);
    assert(fabs(dist_point_line(p, a, b) - sqrt(4.25)) < 1e-12);
    assert(fabs(dist_point_segment({3, 0, 0}, a, b) - 2) < 1e-12);

    /// Skew segments: the x axis piece and a vertical piece at (0.5, 1), 1 apart at (0.5, 0, 0) and (0.5, 1, 0)
    assert(fabs(dist_segment_segment(a, b, {0.5, 1, -1}, {0.5, 1, 1}) - 1) < 1e-12);

    /// A segment piercing the triangle touches it, one floating 1 above it does not
    assert(dist_segment_triangle({0.8, 0.2, -1}, {0.8, 0.2, 1}, a, b, c) == 0);
    assert(fabs(dist_segment_triangle({0, 0, 1}, {1, 1, 1}, a, b, c) - 1) < 1e-12);

    /// The vertical line through (0.25, 0.75) meets the plane z = 0 there, a horizontal one at z = 1 never does
    pair<int, Point3> hit = line_plane_intersection({0.25, 0.75, 5}, {0.25, 0.75, 3}, a, b, c);
    assert(hit.first == 1);
    assert(dist(hit.second, {0.25, 0.75, 0}) < 1e-12);
    assert(line_plane_intersection({0, 0, 1}, {1, 0, 1}, a, b, c).first == 0);
    assert(line_plane_intersection({0, 0, 0}, {1, 1, 0}, a, b, c).first == -1);

    /// Corner tetrahedron of the unit cube: volume 1/6
    vector<Point3> tet = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    vector<vector<int>> faces = {{0, 2, 1}, {0, 1, 3}, {0, 3, 2}, {1, 2, 3}};
    assert(fabs(polyhedron_volume(tet, faces) - 1.0 / 6) < 1e-12);

    return 0;
}
