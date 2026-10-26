#include "../common.h"

#define main library_main
#include "../../code_library/misc/local_minimization.cpp"
#undef main

using Point = vector<double>;

double rand_real(double lo, double hi){
    return lo + (hi - lo) * stress::rand_int(0, 1 << 30) / (1 << 30);
}

Point rand_point(int d, double lo, double hi){
    Point p(d);
    for (double& v : p) v = rand_real(lo, hi);
    return p;
}

double dist(const Point& a, const Point& b){
    double s = 0;
    for (size_t i = 0; i < a.size(); i++) s += (a[i] - b[i]) * (a[i] - b[i]);
    return sqrt(s);
}

/// Every run: no worse than the start or than tol relative above the reference minimum, in at most 2500 d^2 calls
template<typename F>
void check_min(Point x0, F f, double step, double fstar, double tol = 1e-9, double eps = 1e-12){
    long long calls = 0;
    auto counted = [&](const Point& y){ calls++; return f(y); };
    Point x = local_minimum(x0, counted, step, eps);
    assert(calls <= 2500 * (long long)(x0.size() * x0.size()));
    assert(x.size() == x0.size());
    assert(f(x) <= f(x0));
    assert(f(x) <= fstar + tol * max(1.0, abs(fstar)));
}

/// (x - c)^T Q (x - c) + offset with Q's eigenvalues spread over 1..cond along random axes
void check_quadratic(int d, double cond){
    vector<Point> axes(d);
    for (int i = 0; i < d; i++){
        axes[i] = rand_point(d, -1, 1);
        for (int j = 0; j < i; j++){
            double dot = 0;
            for (int k = 0; k < d; k++) dot += axes[i][k] * axes[j][k];
            for (int k = 0; k < d; k++) axes[i][k] -= dot * axes[j][k];
        }
        double len = dist(axes[i], Point(d, 0));
        for (double& v : axes[i]) v /= len;
    }

    Point c = rand_point(d, -100, 100), lambda(d);
    for (int i = 0; i < d; i++) lambda[i] = d == 1 ? cond : pow(cond, (double)i / (d - 1));
    double offset = rand_real(-1000, 1000);
    auto f = [&](const Point& x){
        double s = offset;
        for (int i = 0; i < d; i++){
            double y = 0;
            for (int k = 0; k < d; k++) y += axes[i][k] * (x[k] - c[k]);
            s += lambda[i] * y * y;
        }
        return s;
    };
    check_min(rand_point(d, -100, 100), f, pow(10, stress::rand_int(-1, 2)), offset);
}

/// Curved valley, its only minimum is 0 at (1, ..., 1) for 2 and 3 variables
void check_rosenbrock(int d){
    auto f = [&](const Point& x){
        double s = 0;
        for (int i = 0; i + 1 < d; i++) s += 100 * pow(x[i + 1] - x[i] * x[i], 2) + pow(1 - x[i], 2);
        return s;
    };
    check_min(rand_point(d, -2, 2), f, 0.5, 0);
}

/// Sum of distances against Weiszfeld's iteration, or a data point when that is the minimum
void check_fermat(int n, int dim){
    vector<Point> p(n);
    bool collinear = dim == 2 && stress::rand_int(0, 5) == 0;
    for (int i = 0; i < n; i++){
        p[i] = rand_point(dim, 0, 100);
        if (i > 0 && stress::rand_int(0, 4) == 0) p[i] = p[stress::rand_int(0, i - 1)];
        if (collinear) p[i][1] = p[i][0];  /// the minimum is then a median point, where f has a kink
    }
    auto f = [&](const Point& x){
        double s = 0;
        for (auto& q : p) s += dist(x, q);
        return s;
    };

    Point y(dim, 0);
    for (auto& q : p) for (int k = 0; k < dim; k++) y[k] += q[k] / n;
    for (int it = 0; it < 100000; it++){
        Point num(dim, 0);
        double den = 0;
        bool on_point = false;
        for (auto& q : p){
            double r = dist(y, q);
            if (r < 1e-12) on_point = true;
            for (int k = 0; k < dim; k++) num[k] += q[k] / max(r, 1e-12);
            den += 1 / max(r, 1e-12);
        }
        for (int k = 0; k < dim; k++) num[k] /= den;
        if (on_point || dist(num, y) < 1e-13) break;
        y = num;
    }
    double fstar = f(y);
    for (auto& q : p) fstar = min(fstar, f(q));
    check_min(rand_point(dim, 0, 100), f, 10, fstar);
}

/// Center and radius of the ball with every point of s on its sphere, centered in their affine hull
bool ball_through(const vector<Point>& s, Point& c, double& r){
    int m = s.size() - 1, dim = s[0].size();
    vector<Point> e(m, Point(dim)), a(m, Point(m + 1));
    for (int i = 0; i < m; i++) for (int k = 0; k < dim; k++) e[i][k] = s[i + 1][k] - s[0][k];
    for (int i = 0; i < m; i++){
        for (int j = 0; j <= m; j++){
            int jj = j < m ? j : i;
            double dot = 0;
            for (int k = 0; k < dim; k++) dot += e[i][k] * e[jj][k];
            a[i][j] = j < m ? 2 * dot : dot;
        }
    }
    for (int col = 0; col < m; col++){
        int piv = col;
        for (int r2 = col; r2 < m; r2++) if (abs(a[r2][col]) > abs(a[piv][col])) piv = r2;
        if (abs(a[piv][col]) < 1e-9) return false;
        swap(a[piv], a[col]);
        for (int r2 = 0; r2 < m; r2++){
            if (r2 == col) continue;
            double x = a[r2][col] / a[col][col];
            for (int k = col; k <= m; k++) a[r2][k] -= x * a[col][k];
        }
    }
    c = s[0];
    for (int i = 0; i < m; i++) for (int k = 0; k < dim; k++) c[k] += a[i][m] / a[i][i] * e[i][k];
    r = dist(c, s[0]);
    return true;
}

/// Max of distances against the smallest of all balls through 2, 3 or 4 points that hold every point
void check_enclosing(const vector<Point>& p, const Point& x0, double tol){
    int n = p.size(), dim = x0.size();
    auto f = [&](const Point& x){
        double s = 0;
        for (auto& q : p) s = max(s, dist(x, q));
        return s;
    };

    double fstar = 1e300;
    auto consider = [&](const vector<Point>& s){
        Point c;
        double r;
        if (!ball_through(s, c, r) || r >= fstar) return;
        for (auto& q : p) if (dist(c, q) > r * (1 + 1e-12)) return;
        fstar = r;
    };
    for (int a = 0; a < n; a++){
        for (int b = a + 1; b < n; b++){
            consider({p[a], p[b]});
            for (int c = b + 1; c < n; c++){
                consider({p[a], p[b], p[c]});
                for (int e = c + 1; dim == 3 && e < n; e++) consider({p[a], p[b], p[c], p[e]});
            }
        }
    }
    check_min(x0, f, 20, fstar, tol);
}

/// sum w_i (p_i - t_i)^2 over the probability simplex, against its KKT solution p_i = max(0, t_i - mu / 2w_i)
void check_simplex(const Point& w, const Point& t, bool walls, double tol){
    int d = w.size();
    double lo = -1e6, hi = 1e6;
    for (int it = 0; it < 200; it++){
        double mu = (lo + hi) / 2, s = 0;
        for (int i = 0; i < d; i++) s += max(0.0, t[i] - mu / (2 * w[i]));
        if (s > 1) lo = mu;
        else hi = mu;
    }
    double fstar = 0;
    for (int i = 0; i < d; i++) fstar += w[i] * pow(max(0.0, t[i] - hi / (2 * w[i])) - t[i], 2);

    /// p_i = x_i^2 / sum x^2 has no walls, p_0..p_{d-2} with p_{d-1} = 1 - their sum has the slanted wall sum <= 1
    auto squares = [&](const Point& x){
        double norm = 0, s = 0;
        for (double v : x) norm += v * v;
        for (int i = 0; i < d; i++) s += w[i] * pow(x[i] * x[i] / norm - t[i], 2);
        return s;
    };
    auto walled = [&](const Point& x){
        double s = 0, last = 1;
        for (int i = 0; i + 1 < d; i++){
            if (x[i] < 0) return 1e100;
            s += w[i] * (x[i] - t[i]) * (x[i] - t[i]), last -= x[i];
        }
        return last < 0 ? 1e100 : s + w[d - 1] * (last - t[d - 1]) * (last - t[d - 1]);
    };
    if (walls) check_min(Point(d - 1, 1.0 / d), walled, 1.0 / d, fstar, tol, 1e-9);
    else check_min(Point(d, 1), squares, 0.5, fstar, tol);
}

/// Targets t on the simplex put its minimum inside, targets anywhere often put it on the boundary
void check_random_simplex(int d, bool interior, bool walls){
    Point w = rand_point(d, 0.1, 10), t = rand_point(d, interior ? 0.05 : -0.5, 1);
    if (interior){
        double s = accumulate(t.begin(), t.end(), 0.0);
        for (double& v : t) v /= s;
    }
    check_simplex(w, t, walls, walls ? 1e-2 : 1e-9);
}

/// sum w_i (x_i - t_i)^2 inside the box [0, 1]^d, huge outside: the minimum clamps t to the box, often onto its walls
void check_box(const Point& w, const Point& t, const Point& x0){
    int d = w.size();
    double fstar = 0;
    for (int i = 0; i < d; i++) fstar += w[i] * pow(max(0.0, min(1.0, t[i])) - t[i], 2);

    auto f = [&](const Point& x){
        double s = 0;
        for (int i = 0; i < d; i++){
            if (x[i] < 0 || x[i] > 1) return 1e100;
            s += w[i] * (x[i] - t[i]) * (x[i] - t[i]);
        }
        return s;
    };
    check_min(x0, f, 0.5, fstar);
}

/// beta |x - c|^2 - (beta - 1) (sum (x_i - c_i))^2 / d: a valley along (1, ..., 1), beta times flatter than across
void check_valley(const Point& c, const Point& x0, double beta, double offset, double step){
    int d = c.size();
    auto f = [&](const Point& x){
        double s = 0, t = 0;
        for (int i = 0; i < d; i++) s += (x[i] - c[i]) * (x[i] - c[i]), t += x[i] - c[i];
        return offset + beta * s - (beta - 1) / d * t * t;
    };
    check_min(x0, f, step, offset);
}

/// A bowl stretched by scale, started scale away with step scale / 10
long long bowl_calls(double scale){
    long long calls = 0;
    auto f = [&](const Point& x){
        calls++;
        return pow(x[0] / scale - 1, 2) + pow(x[1] / scale + 0.3, 2) + 7;
    };
    Point x = local_minimum({0, 0}, f, scale / 10);
    assert(f(x) <= 7 + 1e-9 * 7);
    return calls;
}

/// SPOJ NBLINDESC: a walker moves in direction k with probability p_k, blocked moves stay; least expected escape time
struct Maze{
    const int dx[4] = {-1, 0, 1, 0}, dy[4] = {0, 1, 0, -1};
    int cells, src, dst;
    vector<array<int, 4>> next;

    Maze(const vector<string>& g){
        int n = g.size(), m = g[0].size();
        vector<vector<int>> id(n, vector<int>(m, -1));
        cells = 0;
        for (int i = 0; i < n; i++) for (int j = 0; j < m; j++) if (g[i][j] != '#') id[i][j] = cells++;
        next.resize(cells);
        for (int i = 0; i < n; i++){
            for (int j = 0; j < m; j++){
                if (id[i][j] < 0) continue;
                if (g[i][j] == '@') src = id[i][j];
                if (g[i][j] == '*') dst = id[i][j];
                for (int k = 0; k < 4; k++){
                    int a = i + dx[k], b = j + dy[k];
                    bool open = a >= 0 && a < n && b >= 0 && b < m && id[a][b] >= 0;
                    next[id[i][j]][k] = open ? id[a][b] : id[i][j];
                }
            }
        }
    }

    /// Only cells the start reaches get an equation, so a zero probability never leaves a cell without an exit
    double expected_time(const array<double, 4>& p){
        vector<int> idx(cells, -1), order = {src};
        idx[src] = 0;
        for (size_t h = 0; h < order.size(); h++){
            for (int k = 0; k < 4 && order[h] != dst; k++){
                int v = next[order[h]][k];
                if (p[k] > 0 && idx[v] < 0) idx[v] = order.size(), order.push_back(v);
            }
        }
        if (idx[dst] < 0) return 1e100;
        order.erase(find(order.begin(), order.end(), dst));
        int c = order.size();
        for (int r = 0; r < c; r++) idx[order[r]] = r;
        idx[dst] = c;

        /// Elimination by adding probabilities only: a subtraction near a 0 probability cancels to garbage
        vector<Point> a(c, Point(c + 2, 0));
        for (int r = 0; r < c; r++){
            a[r][c + 1] = 1;
            for (int k = 0; k < 4; k++){
                int v = next[order[r]][k];
                if (p[k] > 0 && v != order[r]) a[r][idx[v]] += p[k];
            }
        }
        for (int col = 0; col < c; col++){
            double out = 0;
            for (int k = col + 1; k <= c; k++) out += a[col][k];
            if (out == 0) return 1e100;
            for (int k = col + 1; k <= c + 1; k++) a[col][k] /= out;
            for (int r = col + 1; r < c; r++){
                double x = a[r][col];
                a[r][col] = 0;
                for (int k = col + 1; x > 0 && k <= c + 1; k++) a[r][k] += x * a[col][k];
            }
        }
        Point e(c, 0);
        for (int r = c - 1; r >= 0; r--){
            e[r] = a[r][c + 1];
            for (int k = r + 1; k < c; k++) e[r] += a[r][k] * e[k];
        }
        return e[0];
    }

    /// Directions outside mask get probability 0, the first |mask| - 1 in it get y, the last one 1 - sum y, all > 0
    double on_face(int mask, const Point& y){
        array<double, 4> p = {0, 0, 0, 0};
        double last = 1;
        size_t used = 0;
        for (int k = 0; k < 4; k++){
            if (!(mask >> k & 1)) continue;
            if (used == y.size()){
                p[k] = last;
                break;
            }
            p[k] = y[used++], last -= p[k];
            if (p[k] <= 0) return 1e100;
        }
        return last <= 0 ? 1e100 : expected_time(p);
    }
};

vector<string> random_maze(int n, int m){
    while (true){
        vector<string> g(n, string(m, '.'));
        for (auto& row : g) for (char& ch : row) if (stress::rand_int(0, 3) == 0) ch = '#';
        int s = stress::rand_int(0, n * m - 1), t = stress::rand_int(0, n * m - 1);
        if (s == t) continue;
        g[s / m][s % m] = '@', g[t / m][t % m] = '*';

        vector<int> seen(n * m, 0), stack = {s};
        seen[s] = 1;
        while (!stack.empty()){
            int u = stack.back(), j = u % m;
            stack.pop_back();
            for (int v : {u - m, u + m, j > 0 ? u - 1 : -1, j + 1 < m ? u + 1 : -1}){
                if (v < 0 || v >= n * m || seen[v] || g[v / m][v % m] == '#') continue;
                seen[v] = 1, stack.push_back(v);
            }
        }
        if (seen[t]) return g;
    }
}

/// The expected time jumps where a probability becomes 0, so each face of the simplex is searched on its own
double check_maze(const vector<string>& g){
    Maze maze(g);
    double res = 1e100, reference = 1e100;
    for (int mask = 1; mask < 16; mask++){
        int k = __builtin_popcount(mask);
        auto f = [&](const Point& y){ return maze.on_face(mask, y); };
        Point start(k - 1, 1.0 / k);
        if (k == 1 || f(start) >= 1e100){
            res = min(res, f(start)), reference = min(reference, f(start));
            continue;
        }
        /// A face whose infimum lies on its edge is crawled toward with ever smaller steps, eps = 1e-8 ends that
        long long calls = 0;
        auto counted = [&](const Point& y){ calls++; return f(y); };
        res = min(res, f(local_minimum(start, counted, 0.5 / k, 1e-8)));
        assert(calls <= 2500 * (k - 1) * (k - 1));

        /// Reference: a grid of spacing 1/24 over the open face, then 5^(k-1) grids around the best, each 0.7 as wide
        Point best(k - 1), y(k - 1);
        double fbest = 1e100;
        vector<int> c(k - 1, 1);
        while (true){
            for (int i = 0; i < k - 1; i++) y[i] = c[i] / 24.0;
            if (f(y) < fbest) fbest = f(y), best = y;
            int i = 0;
            while (i < k - 1 && ++c[i] >= 24) c[i] = 1, i++;
            if (i == k - 1) break;
        }
        for (double h = 1 / 24.0; h > 1e-13; h *= 0.7){
            Point center = best;
            vector<int> o(k - 1, -2);
            while (true){
                for (int i = 0; i < k - 1; i++) y[i] = center[i] + o[i] * h;
                if (f(y) < fbest) fbest = f(y), best = y;
                int i = 0;
                while (i < k - 1 && ++o[i] > 2) o[i] = -2, i++;
                if (i == k - 1) break;
            }
        }
        reference = min(reference, fbest);
    }

    /// A face whose infimum lies on its edge has no minimum to agree on, the least over all faces does
    assert(res <= reference * (1 + 1e-9));
    return res;
}

int main(){
    for (long long it = 0; it < stress::scaled(60); it++){
        int dim = 2 + it % 2;
        check_quadratic(it % 6 + 1, pow(10, stress::rand_int(0, 6)));
        check_rosenbrock(dim);
        check_fermat(stress::rand_int(1, 12), dim);
        vector<Point> p(stress::rand_int(2, 10));
        for (auto& q : p) q = rand_point(dim, -50, 50);
        check_enclosing(p, rand_point(dim, -50, 50), dim == 2 ? 1e-9 : 1e-6);
        check_random_simplex(stress::rand_int(2, 7), it % 2, false);
        check_random_simplex(stress::rand_int(2, 7), it % 2, true);
        int d = it % 6 + 1;
        check_box(rand_point(d, 0.1, 10), rand_point(d, -0.5, 1.5), rand_point(d, 0, 1));
    }

    /// The largest size the header claims
    check_quadratic(16, 1e4);

    /// Mazes with a straight path, whose least expected time is its length, reached only with one probability exactly 1
    assert(check_maze({"@..*"}) == 3);
    assert(check_maze({".#.....", "...#.#.", ".@..*..", "....#.."}) == 3);
    /// A maze whose 4 direction face crawls for a million calls when a pass gaining under eps counts as progress
    assert(check_maze({"##...#", ".##@*.", ".....#"}) == 1);
    for (long long it = 0; it < stress::scaled(10); it++){
        check_maze(random_maze(stress::rand_int(1, 3), stress::rand_int(2, 4)));
    }

    /// Ridges where one Nelder-Mead run stalls 2e-6 off, and where restarts without flipping the simplex stall 3e-8 off
    check_enclosing({{-41.8, -1.5, -14.3}, {-17.7, -16.6, 46.5}, {45.3, 44.0, 34.7}, {-41.7, 19.5, -17.2},
                     {38.4, -1.7, -23.2}}, {-38.7, -29.2, -45.0}, 1e-9);
    check_enclosing({{-49.8, 15.5, 6.5}, {39.5, -48.3, 37.5}, {-32.8, 49.4, 31.7}, {-41.5, 16.6, -18.6},
                     {-24.1, 32.3, -1.6}, {7.5, -18.8, 25.2}, {-38.3, 24.4, -36.3}, {35.4, -11.7, -45.1},
                     {27.1, -11.7, -1.3}, {0.1, -47.7, 16.7}}, {-23.7, 2.5, -40.7}, 1e-9);
    /// A corner of the box, where Nelder-Mead left uncapped shrinks for 40000 steps
    check_box({4.9, 3.4}, {-0.2, 1.2}, {0.4, 0.1});
    /// A slanted wall that the compass search, never growing its step, ends 5e-7 short along
    check_simplex({3.4, 6.2, 7.1, 1.2, 7.0, 9.8}, {-0.3, 0.8, 0.8, 0.4, 0.8, -0.3}, true, 1e-7);
    /// A narrow valley the compass search, left uncapped, zigzags down for over 90000 calls
    check_valley({1.6, -2.0, -8.1, 6.7, 3.4, -0.9}, {-5.6, 3.4, 9.6, 9.1, 6.6, 3.3}, 1e6, -97.7, 0.1);

    /// Scale: past 1e6, 1e-12 is below a double's spacing, and an absolute stop would spend 4 times the calls
    assert(bowl_calls(1e6) <= 2 * bowl_calls(1));
    check_min({5.5, 7.5}, [](const Point& x){ return floor(abs(x[0])) + floor(abs(x[1])); }, 1, 0);
    return 0;
}
