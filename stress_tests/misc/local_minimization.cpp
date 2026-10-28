#include "../common.h"

#define main library_main
#include "../../code_library/misc/local_minimization.cpp"
#undef main

template<typename T>
using Point = vector<T>;

/// Every search goes through this type, so local_minimum compiles once per T, not once per lambda: CI's time limit
template<typename T>
using Objective = function<T(const Point<T>&)>;

template<typename T>
T rand_real(T lo, T hi){
    return lo + (hi - lo) * stress::rand_int(0, 1 << 30) / (1 << 30);
}

template<typename T>
Point<T> rand_point(int d, T lo, T hi){
    Point<T> p(d);
    for (T& v : p) v = rand_real(lo, hi);
    return p;
}

template<typename T>
T dist(const Point<T>& a, const Point<T>& b){
    T s = 0;
    for (size_t i = 0; i < a.size(); i++) s += (a[i] - b[i]) * (a[i] - b[i]);
    return sqrt(s);
}

/// Every run: no worse than the start or than tol relative above the reference minimum, in the header's call bound
template<typename T, typename F>
void check_min(Point<T> x0, F f, double step, T fstar, double tol = 1e-9, double eps = 1e-12){
    long long calls = 0;
    Objective<T> counted = [&](const Point<T>& y){ calls++; return f(y); };
    Point<T> x = local_minimum(x0, counted, step, eps);
    assert(calls <= (eps < 1e-12 ? 5000 : 2500) * (long long)(x0.size() * x0.size()));
    assert(x.size() == x0.size());
    assert(f(x) <= f(x0));
    assert(f(x) <= fstar + tol * max(T(1), abs(fstar)));
}

/// (x - c)^T Q (x - c) + offset with Q's eigenvalues spread over 1..cond along random axes
template<typename T>
void check_quadratic(int d, T cond, double tol, double eps){
    vector<Point<T>> axes(d);
    for (int i = 0; i < d; i++){
        axes[i] = rand_point<T>(d, -1, 1);
        for (int j = 0; j < i; j++){
            T dot = 0;
            for (int k = 0; k < d; k++) dot += axes[i][k] * axes[j][k];
            for (int k = 0; k < d; k++) axes[i][k] -= dot * axes[j][k];
        }
        T len = dist(axes[i], Point<T>(d, 0));
        for (T& v : axes[i]) v /= len;
    }

    Point<T> c = rand_point<T>(d, -100, 100), lambda(d);
    for (int i = 0; i < d; i++) lambda[i] = d == 1 ? cond : pow(cond, (T)i / (d - 1));
    T offset = rand_real<T>(-1000, 1000);
    auto f = [&](const Point<T>& x){
        T s = offset;
        for (int i = 0; i < d; i++){
            T y = 0;
            for (int k = 0; k < d; k++) y += axes[i][k] * (x[k] - c[k]);
            s += lambda[i] * y * y;
        }

        return s;
    };
    check_min(rand_point<T>(d, -100, 100), f, pow(10, stress::rand_int(-1, 2)), offset, tol, eps);
}

/// Curved valley, its only minimum is 0 at (1, ..., 1) for 2 and 3 variables
template<typename T>
void check_rosenbrock(int d, double tol, double eps){
    auto f = [&](const Point<T>& x){
        T s = 0;
        for (int i = 0; i + 1 < d; i++) s += 100 * pow(x[i + 1] - x[i] * x[i], 2) + pow(1 - x[i], 2);
        return s;
    };
    check_min(rand_point<T>(d, -2, 2), f, 0.5, T(0), tol, eps);
}

/// Sum of distances against Weiszfeld's iteration, or a data point when that is the minimum
template<typename T>
void check_fermat(int n, int dim, double tol, double eps){
    vector<Point<T>> p(n);
    bool collinear = dim == 2 && stress::rand_int(0, 5) == 0;
    for (int i = 0; i < n; i++){
        p[i] = rand_point<T>(dim, 0, 100);
        if (i > 0 && stress::rand_int(0, 4) == 0) p[i] = p[stress::rand_int(0, i - 1)];
        if (collinear) p[i][1] = p[i][0];  /// the minimum is then a median point, where f has a kink
    }
    auto f = [&](const Point<T>& x){
        T s = 0;
        for (auto& q : p) s += dist(x, q);
        return s;
    };

    /// Stops 10 times finer than the search's eps, so the reference is never the less accurate of the two
    Point<T> y(dim, 0);
    for (auto& q : p) for (int k = 0; k < dim; k++) y[k] += q[k] / n;
    for (int it = 0; it < 100000; it++){
        Point<T> num(dim, 0);
        T den = 0;
        bool on_point = false;
        for (auto& q : p){
            T r = dist(y, q);
            if (r < 1e-12) on_point = true;
            for (int k = 0; k < dim; k++) num[k] += q[k] / max(r, T(1e-12));
            den += 1 / max(r, T(1e-12));
        }
        for (int k = 0; k < dim; k++) num[k] /= den;
        if (on_point || dist(num, y) < eps / 10) break;
        y = num;
    }
    T fstar = f(y);
    for (auto& q : p) fstar = min(fstar, f(q));
    check_min(rand_point<T>(dim, 0, 100), f, 10, fstar, tol, eps);
}

/// Center and radius of the ball with every point of s on its sphere, centered in their affine hull
template<typename T>
bool ball_through(const vector<Point<T>>& s, Point<T>& c, T& r){
    int m = s.size() - 1, dim = s[0].size();
    vector<Point<T>> e(m, Point<T>(dim)), a(m, Point<T>(m + 1));
    for (int i = 0; i < m; i++) for (int k = 0; k < dim; k++) e[i][k] = s[i + 1][k] - s[0][k];
    for (int i = 0; i < m; i++){
        for (int j = 0; j <= m; j++){
            int jj = j < m ? j : i;
            T dot = 0;
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
            T x = a[r2][col] / a[col][col];
            for (int k = col; k <= m; k++) a[r2][k] -= x * a[col][k];
        }
    }
    c = s[0];
    for (int i = 0; i < m; i++) for (int k = 0; k < dim; k++) c[k] += a[i][m] / a[i][i] * e[i][k];
    r = dist(c, s[0]);
    return true;
}

/// Max of distances against the smallest of all balls through 2, 3 or 4 points that hold every point
template<typename T>
void check_enclosing(const vector<Point<T>>& p, const Point<T>& x0, double tol, double eps = 1e-12){
    int n = p.size(), dim = x0.size();
    auto f = [&](const Point<T>& x){
        T s = 0;
        for (auto& q : p) s = max(s, dist(x, q));
        return s;
    };

    T fstar = 1e300;
    auto consider = [&](const vector<Point<T>>& s){
        Point<T> c;
        T r;
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
    check_min(x0, f, 20, fstar, tol, eps);
}

/// sum w_i (p_i - t_i)^2 over the probability simplex, against its KKT solution p_i = max(0, t_i - mu / 2w_i)
template<typename T>
void check_simplex(const Point<T>& w, const Point<T>& t, bool walls, double tol, double eps = 1e-12){
    int d = w.size();
    T lo = -1e6, hi = 1e6;
    for (int it = 0; it < 200; it++){
        T mu = (lo + hi) / 2, s = 0;
        for (int i = 0; i < d; i++) s += max(T(0), t[i] - mu / (2 * w[i]));
        if (s > 1) lo = mu;
        else hi = mu;
    }
    T fstar = 0;
    for (int i = 0; i < d; i++) fstar += w[i] * pow(max(T(0), t[i] - hi / (2 * w[i])) - t[i], 2);

    /// p_i = x_i^2 / sum x^2 has no walls, p_0..p_{d-2} with p_{d-1} = 1 - their sum has the slanted wall sum <= 1
    auto squares = [&](const Point<T>& x){
        T norm = 0, s = 0;
        for (T v : x) norm += v * v;
        for (int i = 0; i < d; i++) s += w[i] * pow(x[i] * x[i] / norm - t[i], 2);
        return s;
    };
    auto walled = [&](const Point<T>& x){
        T s = 0, last = 1;
        for (int i = 0; i + 1 < d; i++){
            if (x[i] < 0) return T(1e100);
            s += w[i] * (x[i] - t[i]) * (x[i] - t[i]), last -= x[i];
        }

        return last < 0 ? T(1e100) : s + w[d - 1] * (last - t[d - 1]) * (last - t[d - 1]);
    };
    if (walls) check_min(Point<T>(d - 1, T(1) / d), walled, 1.0 / d, fstar, tol, 1e-9);
    else check_min(Point<T>(d, 1), squares, 0.5, fstar, tol, eps);
}

/// Targets t on the simplex put its minimum inside, targets anywhere often put it on the boundary
template<typename T>
void check_random_simplex(int d, bool interior, bool walls, double tol, double eps){
    Point<T> w = rand_point<T>(d, 0.1, 10), t = rand_point<T>(d, interior ? 0.05 : -0.5, 1);
    if (interior){
        T s = accumulate(t.begin(), t.end(), T(0));
        for (T& v : t) v /= s;
    }
    check_simplex(w, t, walls, walls ? 1e-2 : tol, eps);
}

/// sum w_i (x_i - t_i)^2 inside the box [0, 1]^d, huge outside: the minimum clamps t to the box, often onto its walls
template<typename T>
void check_box(const Point<T>& w, const Point<T>& t, const Point<T>& x0, double tol = 1e-9, double eps = 1e-12){
    int d = w.size();
    T fstar = 0;
    for (int i = 0; i < d; i++) fstar += w[i] * pow(max(T(0), min(T(1), t[i])) - t[i], 2);

    auto f = [&](const Point<T>& x){
        T s = 0;
        for (int i = 0; i < d; i++){
            if (x[i] < 0 || x[i] > 1) return T(1e100);
            s += w[i] * (x[i] - t[i]) * (x[i] - t[i]);
        }

        return s;
    };
    check_min(x0, f, 0.5, fstar, tol, eps);
}

/// beta |x - c|^2 - (beta - 1) (sum (x_i - c_i))^2 / d: a valley along (1, ..., 1), beta times flatter than across
void check_valley(const Point<double>& c, const Point<double>& x0, double beta, double offset, double step){
    int d = c.size();
    auto f = [&](const Point<double>& x){
        double s = 0, t = 0;
        for (int i = 0; i < d; i++) s += (x[i] - c[i]) * (x[i] - c[i]), t += x[i] - c[i];
        return offset + beta * s - (beta - 1) / d * t * t;
    };
    check_min(x0, f, step, offset);
}

/// A bowl stretched by scale with minimum fmin at scale * center, started scale below on each axis, step scale / 10
long long bowl_calls(const Point<double>& center, double fmin, double scale){
    int d = center.size();
    long long calls = 0;
    Objective<double> f = [&](const Point<double>& x){
        calls++;
        double s = fmin;
        for (int i = 0; i < d; i++) s += (i + 1) * pow(x[i] / scale - center[i], 2);
        return s;
    };
    Point<double> x0(d);
    for (int i = 0; i < d; i++) x0[i] = (center[i] - 1) * scale;
    Point<double> x = local_minimum(x0, f, scale / 10);
    assert(f(x) <= fmin + 1e-9 * max(1.0, fmin));
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
        vector<Point<double>> a(c, Point<double>(c + 2, 0));
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
        Point<double> e(c, 0);
        for (int r = c - 1; r >= 0; r--){
            e[r] = a[r][c + 1];
            for (int k = r + 1; k < c; k++) e[r] += a[r][k] * e[k];
        }

        return e[0];
    }

    /// Directions outside mask get probability 0, the first |mask| - 1 in it get y, the last one 1 - sum y, all > 0
    double on_face(int mask, const Point<double>& y){
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
        auto f = [&](const Point<double>& y){ return maze.on_face(mask, y); };
        Point<double> start(k - 1, 1.0 / k);
        if (k == 1 || f(start) >= 1e100){
            res = min(res, f(start)), reference = min(reference, f(start));
            continue;
        }
        /// A face whose infimum lies on its edge is crawled toward with ever smaller steps, eps = 1e-8 ends that
        long long calls = 0;
        Objective<double> counted = [&](const Point<double>& y){ calls++; return f(y); };
        res = min(res, f(local_minimum(start, counted, 0.5 / k, 1e-8)));
        assert(calls <= 2500 * (k - 1) * (k - 1));

        /// Reference: a grid of spacing 1/24 over the open face, then 5^(k-1) grids around the best, each 0.7 as wide
        Point<double> best(k - 1), y(k - 1);
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
            Point<double> center = best;
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

/// Every random family searched in T with eps: kinks and boxes within tol, smooth f within smooth, 3-dim balls 1e-6
template<typename T>
void check_families(long long iterations, double eps, double tol, double smooth){
    for (long long it = 0; it < iterations; it++){
        int dim = 2 + it % 2;
        check_quadratic<T>(it % 6 + 1, pow(10, stress::rand_int(0, 6)), smooth, eps);
        check_rosenbrock<T>(dim, smooth, eps);
        check_fermat<T>(stress::rand_int(1, 12), dim, tol, eps);
        vector<Point<T>> p(stress::rand_int(2, 10));
        for (auto& q : p) q = rand_point<T>(dim, -50, 50);
        check_enclosing(p, rand_point<T>(dim, -50, 50), dim == 2 ? tol : 1e-6, eps);
        check_random_simplex<T>(stress::rand_int(2, 7), it % 2, false, smooth, eps);
        check_random_simplex<T>(stress::rand_int(2, 7), it % 2, true, tol, eps);
        int d = it % 6 + 1;
        check_box(rand_point<T>(d, 0.1, 10), rand_point<T>(d, -0.5, 1.5), rand_point<T>(d, 0, 1), tol, eps);
    }
}

int main(){
    check_families<double>(stress::scaled(20), 1e-12, 1e-9, 1e-12);
    /// Smooth f in long double within 1e-17, which a value rounded to double inside the search misses
    check_families<long double>(stress::scaled(10), 1e-12, 1e-9, 1e-17);
    /// Below double's spacing kinks and boxes must end 1e3 times closer
    check_families<long double>(stress::scaled(15), 1e-16, 1e-12, 1e-17);

    /// The largest size the header claims, where a too long expansion step passes the call bound
    for (int it = 0; it < 2; it++) check_quadratic<double>(16, 1e4, 1e-9, 1e-12);

    /// Mazes with a straight path, whose least expected time is its length, reached only with one probability exactly 1
    assert(check_maze({"@..*"}) == 3);
    assert(check_maze({".#.....", "...#.#.", ".@..*..", "....#.."}) == 3);
    /// A face whose infimum is on its edge, crawled past the call bound if a pass gaining under eps counts as progress
    assert(check_maze({"##...#", ".##@*.", ".....#"}) == 1);
    for (long long it = 0; it < stress::scaled(5); it++){
        check_maze(random_maze(stress::rand_int(1, 3), stress::rand_int(2, 4)));
    }

    /// Ridges of a max where a single Nelder-Mead run, or restarts that never flip the simplex, stall short
    check_enclosing<double>({{-41.8, -1.5, -14.3}, {-17.7, -16.6, 46.5}, {45.3, 44.0, 34.7}, {-41.7, 19.5, -17.2},
                             {38.4, -1.7, -23.2}}, {-38.7, -29.2, -45.0}, 1e-9);
    check_enclosing<double>({{-49.8, 15.5, 6.5}, {39.5, -48.3, 37.5}, {-32.8, 49.4, 31.7}, {-41.5, 16.6, -18.6},
                             {-24.1, 32.3, -1.6}, {7.5, -18.8, 25.2}, {-38.3, 24.4, -36.3}, {35.4, -11.7, -45.1},
                             {27.1, -11.7, -1.3}, {0.1, -47.7, 16.7}}, {-23.7, 2.5, -40.7}, 1e-9);
    /// Ridges where restarts that stop shrinking, or that end after 2 idle passes instead of 3, stall short
    check_enclosing<double>({{17.0, 15.4, -15.7}, {39.1, -44.4, 12.8}, {13.9, -46.4, 22.8}, {-8.2, 0.7, 49.6},
                             {6.3, -37.3, 5.6}, {-38.2, -31.1, 11.7}, {-26.0, 23.1, -3.7}, {-4.4, -8.6, 38.5}},
                            {32.1, 4.8, 35.4}, 1e-9);
    /// A ridge where reflecting whenever the point beats the worst vertex, not the second worst, stalls short
    check_enclosing<double>({{34.9, -33.9, 33.9}, {-46.7, -41.5, 47.6}, {-37.6, -18.3, -12.7}, {-17.6, -19.1, 24.1},
                             {40.7, -46.7, -24.0}, {9.9, -44.4, -10.0}, {-28.7, 41.8, -42.1}, {36.9, 23.8, -44.3}},
                            {-11.3, -24.6, -47.1}, 1e-9);
    /// A corner of the box, where Nelder-Mead without its iteration cap shrinks past the call bound
    check_box<double>({4.9, 3.4}, {-0.2, 1.2}, {0.4, 0.1});
    /// A slanted wall the compass search ends short along when it never grows its step
    check_simplex<double>({3.4, 6.2, 7.1, 1.2, 7.0, 9.8}, {-0.3, 0.8, 0.8, 0.4, 0.8, -0.3}, true, 1e-7);
    /// A narrow valley the compass search zigzags down past the call bound when its rounds are not capped
    check_valley({1.6, -2.0, -8.1, 6.7, 3.4, -0.9}, {-5.6, 3.4, 9.6, 9.1, 6.6, 3.3}, 1e6, -97.7, 0.1);

    /// Scale: past 1e6, 1e-12 is below a double's spacing, so an absolute stop would spend far more calls
    assert(bowl_calls({1, -0.3}, 7, 1e6) <= 2 * bowl_calls({1, -0.3}, 7, 1));
    check_min<double>({5.5, 7.5}, [](const Point<double>& x){ return floor(abs(x[0])) + floor(abs(x[1])); }, 1, 0);
    /// Stop tests floored at 1: a minimum of 0 at x = 0 costs about as many calls as one away from both zeros
    assert(bowl_calls(Point<double>(4, 0), 0, 1) <= 1.5 * bowl_calls(Point<double>(4, 1), 7, 1));

    /// A step up to eps * max(1, max |x_i|) is under the stop test from the start, and could return x unmoved
    Objective<double> bowl = [](const Point<double>& x){ return pow(x[0] - 1005, 2) + pow(x[1] - 1000, 2); };
    stress::expect_assert("misc/local_minimization.cpp", [&]{ local_minimum({1000, 1000}, bowl, 1e-9); });
    stress::expect_assert("misc/local_minimization.cpp", [&]{ local_minimum({0, 0}, bowl, 1e-12); });
    check_min<double>({1000, 1000}, bowl, 2e-9, 0);

    return 0;
}
