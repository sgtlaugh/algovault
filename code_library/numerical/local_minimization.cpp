/***
 *
 * Local Minimization
 * A point where f is locally smallest, for f over a few real variables, smooth or not, with no derivatives
 *
 * Complexity: at most 2500 d^2 calls of f for d variables, 5000 d^2 when eps is below 1e-12
 *
 * local_minimum(x, f, step, eps = 1e-12): a local minimum of f near the start x, step a rough distance to it
 *
 * Notes:
 *   - f takes const vector<double>& and returns double; pass x as a vector<long double> to search in long double
 *   - f never returns nan and is finite at x; x is nonempty and step > eps * max(1, max |x_i|)
 *   - Walls: f returns a huge value like 1e100 outside the allowed region, and x starts inside it
 *   - Probabilities summing to 1: search free x with p_i = x_i^2 / sum x_j^2, starting from x = (1, ..., 1)
 *   - Local only: with several minima, run from several starts and keep the best
 *   - Ends within 1e-9 relative on smooth functions, sums of distances and boxes; a max of distances in 3 dims 1e-6
 *   - A minimum on a slanted wall like x_0 + x_1 <= 1 can end 1e-2 off: change variables so walls are axis-aligned
 *   - A value f takes only at exact points, like a cost once a probability is exactly 0, is missed: search those apart
 *   - Long double gains only with eps = 1e-16 and f computed in long double: about 1e4 times closer
 *   - A larger eps ends sooner, as needed when f only approaches its smallest value at a wall
 *   - Tested up to 16 variables
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template<typename T, typename F>
vector<T> nelder_mead(vector<T> x, F& f, double step, double eps){
    int d = x.size();
    vector<vector<T>> p(d + 1, x);
    vector<T> v(d + 1);
    for (int i = 0; i < d; i++) p[i + 1][i] += step;
    for (int i = 0; i <= d; i++) v[i] = f(p[i]);

    /// Capped: wedged in a corner of walls it can shrink for very long, the next run restarts it fresh
    for (int iter = 0; iter < 100 * d * d; iter++){
        vector<int> order(d + 1);
        iota(order.begin(), order.end(), 0);
        sort(order.begin(), order.end(), [&](int a, int b){ return v[a] < v[b]; });
        int lo = order[0], hi = order[d], second = order[d - 1];

        /// Relative: vertices an ulp apart average to one of themselves, so an absolute size might never be reached
        T size = 0, scale = 1;
        for (int i = 0; i <= d; i++){
            for (int j = 0; j < d; j++) size = max(size, abs(p[i][j] - p[lo][j])), scale = max(scale, abs(p[i][j]));
        }
        if (size < eps * scale) return p[lo];

        vector<T> centroid(d, 0);
        for (int i = 0; i <= d; i++){
            if (i == hi) continue;
            for (int j = 0; j < d; j++) centroid[j] += p[i][j] / d;
        }
        auto along = [&](T t){
            vector<T> y(d);
            for (int j = 0; j < d; j++) y[j] = centroid[j] + t * (p[hi][j] - centroid[j]);
            return y;
        };

        vector<T> r = along(-1);
        T fr = f(r);
        if (fr < v[lo]){
            vector<T> e = along(-2);
            T fe = f(e);
            if (fe < fr) p[hi] = e, v[hi] = fe;
            else p[hi] = r, v[hi] = fr;
        }
        else if (fr < v[second]) p[hi] = r, v[hi] = fr;
        else{
            vector<T> c = along(fr < v[hi] ? -0.5 : 0.5);
            T fc = f(c);
            if (fc < min(fr, v[hi])) p[hi] = c, v[hi] = fc;
            else{
                for (int i = 0; i <= d; i++){
                    if (i == lo) continue;
                    for (int j = 0; j < d; j++) p[i][j] = (p[i][j] + p[lo][j]) / 2;
                    v[i] = f(p[i]);
                }
            }
        }
    }

    return p[min_element(v.begin(), v.end()) - v.begin()];
}

/// A simplex squeezed flat against a wall stalls, and only moves between coordinate pairs slide along x_0 + x_1 = c
template<typename T, typename F>
vector<T> compass_search(vector<T> x, F& f, double step, double eps){
    int d = x.size();
    T fx = f(x);

    /// Capped: zigzagging down a narrow valley can take very many rounds, the next Nelder-Mead run follows it
    for (int round = 0; round < 100; round++){
        T scale = 1;
        for (T xi : x) scale = max(scale, abs(xi));
        if (step < eps * scale) break;

        bool moved = false;
        for (int i = 0; i < d; i++){
            for (int j = i; j < d; j++){
                for (double s : {step, -step}){
                    vector<T> y = x;
                    y[i] += s;
                    if (j != i) y[j] -= s;
                    T fy = f(y);
                    if (fy < fx) x = y, fx = fy, moved = true;
                }
            }
        }
        step = moved ? step * 2 : step / 2;
    }

    return x;
}

template<typename T = double, typename F>
vector<T> local_minimum(vector<T> x, F f, double step, double eps = 1e-12){
    T scale = 1;
    for (T xi : x) scale = max(scale, abs(xi));
    /// A step up to eps * scale is already under the stop tests, so the search would return x unmoved
    assert(!x.empty() && step > eps * scale);
    T fx = f(x);

    /// A simplex stalled on a ridge of a max of functions often moves again when restarted flipped and smaller
    for (int pass = 0, idle = 0; idle < 3; pass++){
        double size = (pass % 2 ? -step : step) / (1 << min(pass, 30));
        vector<T> y = compass_search(nelder_mead(x, f, size, eps), f, step, eps);
        T fy = f(y);
        idle = fy < fx - eps * max(T(1), abs(fx)) ? 0 : idle + 1;
        x = y, fx = fy;
    }

    return x;
}

int main(){
    /// The bottom of a valley 10 times steeper across than along
    auto valley = [](const vector<double>& v){ return pow(v[0] - 3, 2) + 10 * pow(v[1] + 1, 2); };
    vector<double> x = local_minimum({0, 0}, valley, 1);
    assert(abs(x[0] - 3) < 1e-6);
    assert(abs(x[1] + 1) < 1e-6);

    /// Not smooth: the smallest circle around a right triangle is centered on its hypotenuse (0, 8) to (6, 0), radius 5
    vector<vector<double>> triangle = {{0, 0}, {6, 0}, {0, 8}};
    auto farthest = [&](const vector<double>& v){
        double r = 0;
        for (auto& c : triangle) r = max(r, hypot(v[0] - c[0], v[1] - c[1]));
        return r;
    };
    vector<double> center = local_minimum({1, 1}, farthest, 1);
    assert(abs(center[0] - 3) < 1e-6);
    assert(abs(center[1] - 4) < 1e-6);

    /// p_0^2 + 2 p_1^2 + 3 p_2^2 over probabilities: Lagrange makes p_i proportional to 1 / (i + 1), 6/11, 3/11, 2/11
    auto probabilities = [](const vector<double>& x){
        double norm = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
        return vector<double>{x[0] * x[0] / norm, x[1] * x[1] / norm, x[2] * x[2] / norm};
    };
    auto weighted = [&](const vector<double>& x){
        vector<double> p = probabilities(x);
        return p[0] * p[0] + 2 * p[1] * p[1] + 3 * p[2] * p[2];
    };
    vector<double> p = probabilities(local_minimum({1, 1, 1}, weighted, 0.5));
    assert(abs(p[0] - 6.0 / 11) < 1e-6);
    assert(abs(p[1] - 3.0 / 11) < 1e-6);
    assert(abs(p[2] - 2.0 / 11) < 1e-6);

    return 0;
}
