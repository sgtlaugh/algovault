// LINK: -lquadmath
#include "../common.h"
#include <quadmath.h>

#define main library_main
#include "../../code_library/linear_algebra/gauss_band_matrix.cpp"
#undef main

typedef GaussBand<long double> Band;

/// Dense elimination with partial pivoting in quad precision
vector<__float128> dense_solve(vector<vector<__float128>> a, vector<__float128> b){
    int n = b.size();
    for (int c = 0; c < n; c++){
        int p = c;
        for (int r = c + 1; r < n; r++) if (fabsq(a[r][c]) > fabsq(a[p][c])) p = r;
        swap(a[p], a[c]), swap(b[p], b[c]);
        for (int r = c + 1; r < n; r++){
            __float128 f = a[r][c] / a[c][c];
            if (f == 0) continue;
            for (int k = c; k < n; k++) a[r][k] -= f * a[c][k];
            b[r] -= f * b[c];
        }
    }

    vector<__float128> x(n);
    for (int r = n - 1; r >= 0; r--){
        __float128 s = b[r];
        for (int k = r + 1; k < n; k++) s -= a[r][k] * x[k];
        x[r] = s / a[r][r];
    }

    return x;
}

/// Walk x(u) = rhs(u) + sum of flow(u, v) x(v) stopping at u with chance stop(u): dense, quad, random elimination order
vector<__float128> dense_walk(vector<vector<__float128>> flow, vector<__float128> stop, vector<__float128> rhs){
    int n = rhs.size();
    vector<int> order(n);
    iota(order.begin(), order.end(), 0);
    shuffle(order.begin(), order.end(), stress::rng());

    vector<bool> done(n, false);
    for (int k : order){
        done[k] = true;
        __float128 pivot = stop[k];
        for (int v = 0; v < n; v++) if (!done[v]) pivot += flow[k][v];
        assert(pivot > 0);
        for (int v = 0; v < n; v++) if (!done[v]) flow[k][v] /= pivot;
        stop[k] /= pivot, rhs[k] /= pivot;

        for (int i = 0; i < n; i++){
            __float128 f = flow[i][k];
            if (done[i] || f == 0) continue;
            flow[i][k] = 0;
            for (int v = 0; v < n; v++) if (!done[v] && v != i) flow[i][v] += f * flow[k][v];
            stop[i] += f * stop[k], rhs[i] += f * rhs[k];
        }
    }

    vector<__float128> x(n, 0);
    for (int t = n - 1; t >= 0; t--){
        int k = order[t];
        x[k] = rhs[k];
        for (int v = 0; v < n; v++) if (v != k) x[k] += flow[k][v] * x[v];
    }

    return x;
}

/// Strictly diagonally dominant systems: signed weights take the plain regime, weights in [0, 1] the probability one
void check_dominant(){
    for (long long it = 0; it < stress::scaled(400); it++){
        int n = stress::rand_int(1, it % 10 ? 12 : 400), m = stress::rand_int(1, it % 10 ? 8 : 6);
        if (it % 50 == 0) m = stress::rand_int(31, 60), n = min(n, 30);  /// widest band, kept small: cells times band
        if (stress::rand_int(0, 1)) swap(n, m);  /// so that numbering by column is the narrower band as often as by row

        /// Generated as a band of half width h in its own order, the library must find an order at least as narrow
        bool by_column = stress::rand_int(0, 1);
        int cells = n * m, h = stress::rand_int(0, min(34, cells - 1));
        auto cell = [&](int u){ return by_column ? make_pair(u % n, u / n) : make_pair(u / m, u % m); };

        vector<map<int, long double>> coef(cells);
        vector<long double> rhs(cells, 0);
        vector<bool> untouched(cells, false);
        bool signed_weights = it % 3, signed_parts = signed_weights || stress::rand_int(0, 3) == 0;
        bool heavy_self = !signed_weights && stress::rand_int(0, 3) == 0;
        for (int u = 0; u < cells; u++){
            /// A cell without calls stands for x(u) = 0, other equations may still use it
            if (it % 4 == 0 && stress::rand_int(0, 3) == 0){
                coef[u][u] = 1, untouched[u] = true;
                continue;
            }

            long double off = 0;
            for (int v = max(0, u - h); v <= min(cells - 1, u + h); v++){
                if (v == u || stress::rand_int(0, 2)) continue;
                long double c = stress::rand_int(-1000, 1000) / 250.0L;
                if (!signed_weights) c = -stress::rand_int(0, 1000) / 1e3L / (2 * h);
                coef[u][v] = c, off += fabsl(c);
            }
            /// Strict diagonal dominance, the library pivots on the diagonal only; unsigned weights stay in [0, 1]
            coef[u][u] = (off + stress::rand_int(1, 100) / 10.0L) * (stress::rand_int(0, 1) ? 1 : -1);
            if (!signed_weights) coef[u][u] = off + stress::rand_int(1, 100) / 100.0L * (1 - off);
            /// A weight above 1 on itself: every weight >= 0 yet the sum passes 1, so the regime must be plain
            if (heavy_self && stress::rand_int(0, 15) == 0) coef[u][u] = -coef[u][u] - 1;
            if (stress::rand_int(0, 1)) rhs[u] = stress::rand_int(-1000, 1000) / 8.0L;
        }

        /// sum of a(u, v) x(v) = b(u) is x(u) = b(u) + (1 - a(u, u)) x(u) - sum over v != u of a(u, v) x(v)
        Band gauss(n, m);
        int row_reach = 0, column_reach = 0;
        bool probabilities = true;
        vector<long double> stop(cells, 1);
        auto add_term = [&](int i, int j, int k, int l, long double weight){
            gauss.add_term(i, j, k, l, weight);
            row_reach = max(row_reach, abs((i * m + j) - (k * m + l)));
            column_reach = max(column_reach, abs((j * n + i) - (l * n + k)));
            if (weight < 0) probabilities = false;
            stop[by_column ? j * n + i : i * m + j] -= weight;
        };

        for (int u = 0; u < cells; u++){
            if (untouched[u]) continue;

            auto [i, j] = cell(u);
            for (auto [v, c] : coef[u]){
                long double weight = v == u ? 1 - c : -c;
                if (weight == 0 && stress::rand_int(0, 1)) continue;

                auto [k, l] = cell(v);
                /// Two calls must accumulate, a negative part makes the regime plain even if the sum is not
                long double part = stress::rand_int(0, 1) ? stress::rand_int(-1000, 1000) / 64.0L : 0;
                if (!signed_parts) part = weight * stress::rand_int(0, 4) / 4;
                add_term(i, j, k, l, part);
                add_term(i, j, k, l, weight - part);
            }

            long double part = stress::rand_int(0, 1) ? stress::rand_int(-1000, 1000) / 64.0L : 0;
            if (rhs[u] != 0 || stress::rand_int(0, 1)){
                gauss.add_constant(i, j, part);
                gauss.add_constant(i, j, rhs[u] - part);
            }
        }
        for (int u = 0; u < cells; u++) if (stop[u] < -1e-14L) probabilities = false;

        auto res = gauss.solve();
        assert(gauss.reach == min(row_reach, column_reach));
        if (row_reach != column_reach) assert(gauss.transposed == (column_reach < row_reach));
        assert(gauss.probabilities == probabilities);

        assert((int)res.size() == n);
        vector<long double> x(cells);
        for (int u = 0; u < cells; u++){
            auto [i, j] = cell(u);
            assert((int)res[i].size() == m);
            x[u] = res[i][j];
        }

        long double norm = 1;
        for (auto v : x) norm = max(norm, fabsl(v));
        for (int u = 0; u < cells; u++){
            long double s = 0;
            for (auto [v, c] : coef[u]) s += c * x[v];
            assert(fabsl(s - rhs[u]) <= 1e-12L * norm * (1 + fabsl(coef[u][u])));
            if (untouched[u]) assert(x[u] == 0);
        }

        if (cells <= 80){
            vector<vector<__float128>> a(cells, vector<__float128>(cells, 0));
            vector<__float128> b(cells);
            for (int u = 0; u < cells; u++){
                for (auto [v, c] : coef[u]) a[u][v] = c;
                b[u] = rhs[u];
            }
            auto expected = dense_solve(a, b);
            for (int u = 0; u < cells; u++) assert(fabsq(expected[u] - x[u]) <= (__float128)1e-13L * norm);
        }
    }
}

/// Random walks in double, as contestants call it, with move chances down to 1e-15 where 1 - (stay) cancels
void check_walks(){
    const int dx[] = {-1, 0, 1, 0, -1, -1, 1, 1}, dy[] = {0, 1, 0, -1, -1, 1, -1, 1};

    for (long long it = 0; it < stress::scaled(1500); it++){
        int n = stress::rand_int(1, it % 10 ? 8 : 12), m = stress::rand_int(1, it % 10 ? 8 : 10);
        int moves = stress::rand_int(0, 3) ? 4 : 8, tiny = stress::rand_int(0, 15);
        int walls = stress::rand_int(0, 2) * 15, encoding = stress::rand_int(0, 2);
        bool same_everywhere = stress::rand_int(0, 1);
        int cells = n * m;

        auto chances = [&](){
            vector<double> p(moves);
            double total = 0;
            while (total == 0){
                for (auto& q : p){
                    int kind = stress::rand_int(0, 7);
                    q = kind == 0 ? 0 : stress::rand_int(1, 100);
                    if (kind && kind <= 3 && tiny) q = stress::rand_int(1, 9) * pow(10.0, -tiny);
                    total += q;
                }
            }
            for (auto& q : p) q /= total;
            return p;
        };

        /// Encodings, 0: a step into a wall is a self term, 1: divided by the chance to move, 2: a chance to stop
        vector<bool> wall(cells), target(cells, false);
        for (int u = 0; u < cells; u++) wall[u] = stress::rand_int(0, 99) < walls;
        for (int t = stress::rand_int(1, 3); t > 0; t--) target[stress::rand_int(0, cells - 1)] = true;
        for (int u = 0; u < cells; u++) if (target[u]) wall[u] = false;

        vector<vector<double>> flow(cells, vector<double>(cells, 0));
        vector<double> self(cells, 0), rhs(cells, 0);
        vector<bool> moving(cells, false), dropped(cells, false);
        vector<double> shared = chances();
        for (int u = 0; u < cells; u++){
            if (wall[u] || target[u]) continue;

            int i = u / m, j = u % m;
            vector<double> p = same_everywhere ? shared : chances();
            /// Chances to stop down to 1e-13, above the rounding cutoff, must stay chances to stop
            double keep = 1 - stress::rand_int(1, 5000) / 1e4;
            if (stress::rand_int(0, 1)) keep = 1 - stress::rand_int(1, 9) * pow(10.0, -stress::rand_int(4, 13));
            if (encoding != 2 || stress::rand_int(0, 1)) keep = 1;
            double out = 0, stay = 0;
            for (int d = 0; d < moves; d++){
                int k = i + dx[d], l = j + dy[d];
                if (k < 0 || l < 0 || k >= n || l >= m || wall[k * m + l]) stay += p[d];
                else out += p[d];
            }
            if (encoding == 1 && out == 0) continue;

            double c = stress::rand_int(0, 2) ? 1 : stress::rand_int(0, 10);
            rhs[u] = encoding == 1 ? c / out : c;
            for (int d = 0; d < moves; d++){
                int k = i + dx[d], l = j + dy[d];
                if (k < 0 || l < 0 || k >= n || l >= m || wall[k * m + l] || p[d] == 0) continue;
                flow[u][k * m + l] += encoding == 1 ? p[d] / out : p[d] * keep;
            }
            if (encoding != 1) self[u] = stay * keep;
            dropped[u] = keep < 1;
            moving[u] = true;
        }

        /// A cell that cannot reach a stop is left without calls, so it stops too and the system stays regular
        vector<bool> escapes(cells, false);
        for (int u = 0; u < cells; u++) escapes[u] = !moving[u] || dropped[u];
        for (bool grew = true; grew;){
            grew = false;
            for (int u = 0; u < cells; u++){
                for (int v = 0; v < cells && !escapes[u]; v++){
                    if (flow[u][v] > 0 && escapes[v]) escapes[u] = grew = true;
                }
            }
        }

        GaussBand<double> gauss(n, m);
        vector<vector<__float128>> exact_flow(cells, vector<__float128>(cells, 0));
        vector<__float128> exact_stop(cells, 1), exact_rhs(cells, 0);
        for (int u = 0; u < cells; u++){
            if (!moving[u] || !escapes[u]) continue;

            /// The chance of stopping is exactly 1 - (weights) for 2, and exactly 0 otherwise however the weights round
            exact_stop[u] = 0, exact_rhs[u] = rhs[u];
            if (dropped[u]) exact_stop[u] = 1 - (__float128)self[u];
            gauss.add_constant(u / m, u % m, rhs[u]);
            if (self[u] > 0) gauss.add_term(u / m, u % m, u / m, u % m, self[u]);
            for (int v = 0; v < cells; v++){
                if (flow[u][v] == 0) continue;
                exact_flow[u][v] = flow[u][v];
                if (dropped[u]) exact_stop[u] -= flow[u][v];
                /// Splitting a weight in two nonnegative calls keeps the probability regime
                double part = stress::rand_int(0, 3) ? 0 : flow[u][v] / 3;
                if (part > 0) gauss.add_term(u / m, u % m, v / m, v % m, part);
                gauss.add_term(u / m, u % m, v / m, v % m, flow[u][v] - part);
            }
        }

        auto res = gauss.solve();
        assert(gauss.probabilities);
        auto expected = dense_walk(exact_flow, exact_stop, exact_rhs);
        /// A chance s to stop is read as 1 - (sum of weights) in double, so it carries about 1e-16 / s relative error
        __float128 tolerance = 1e-13;
        for (int u = 0; u < cells; u++) if (exact_stop[u] > 0) tolerance = max(tolerance, 1e-14 / exact_stop[u]);
        for (int u = 0; u < cells; u++){
            assert(res[u / m][u % m] >= 0);
            assert(fabsq(expected[u] - res[u / m][u % m]) <= tolerance * expected[u]);
        }
    }
}

/// SPOJ NBLINDESC's maze, where the old pivots, 1 - (chance of coming back), cancelled to -1.5e20 as small -> 0
void check_vanishing_direction(){
    const vector<string> maze = {".#.....", "...#.#.", ".@..*..", "....#.."};
    const int dx[] = {-1, 0, 1, 0}, dy[] = {0, 1, 0, -1};
    int n = maze.size(), m = maze[0].size(), cells = n * m;
    auto open = [&](int i, int j){ return i >= 0 && j >= 0 && i < n && j < m && maze[i][j] != '#'; };

    for (int e = 1; e <= 16; e++){
        double small = e == 16 ? 0 : pow(10.0, -e), p[4] = {small, 1 - 3 * small, small, small};

        /// Only cells the start reaches get an equation, as a cell that cannot move has none
        vector<bool> reached(cells, false);
        vector<int> queue = {2 * m + 1};
        reached[2 * m + 1] = true;
        for (int t = 0; t < (int)queue.size(); t++){
            int i = queue[t] / m, j = queue[t] % m;
            for (int d = 0; d < 4 && maze[i][j] != '*'; d++){
                int k = i + dx[d], l = j + dy[d];
                if (p[d] == 0 || !open(k, l) || reached[k * m + l]) continue;
                reached[k * m + l] = true, queue.push_back(k * m + l);
            }
        }

        GaussBand<double> gauss(n, m);
        vector<vector<__float128>> flow(cells, vector<__float128>(cells, 0));
        vector<__float128> stop(cells, 1), rhs(cells, 0);
        for (int u : queue){
            int i = u / m, j = u % m;
            if (maze[i][j] == '*') continue;

            double out = 0;
            for (int d = 0; d < 4; d++) if (open(i + dx[d], j + dy[d])) out += p[d];
            gauss.add_constant(i, j, 1 / out);
            stop[u] = 0, rhs[u] = 1 / out;
            for (int d = 0; d < 4; d++){
                int k = i + dx[d], l = j + dy[d];
                if (!open(k, l) || p[d] == 0) continue;
                gauss.add_term(i, j, k, l, p[d] / out);
                flow[u][k * m + l] = p[d] / out;
            }
        }

        double res = gauss.solve()[2][1];
        assert(gauss.probabilities);
        assert(res >= 3);
        if (small == 0) assert(res == 3);
        __float128 expected = dense_walk(flow, stop, rhs)[2 * m + 1];
        assert(fabsq(expected - res) <= (__float128)1e-13 * expected);
    }
}

int main(){
    check_dominant();
    check_walks();
    check_vanishing_direction();

    return 0;
}
