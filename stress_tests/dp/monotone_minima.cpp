#include "../common.h"

#define main library_main
#include "../../code_library/dp/monotone_minima.cpp"
#undef main

using Matrix = vector<vector<long long>>;

vector<int> brute(const Matrix& a){
    vector<int> res;
    for (auto& row : a) res.push_back(min_element(row.begin(), row.end()) - row.begin());

    return res;
}

bool monotone_argmin(const Matrix& a){
    vector<int> best = brute(a);

    return is_sorted(best.begin(), best.end());
}

bool totally_monotone(const Matrix& a){
    int n = a.size(), m = a[0].size();
    for (int i = 0; i < n; i++){
        for (int i2 = i + 1; i2 < n; i2++){
            for (int j = 0; j < m; j++){
                for (int k = j + 1; k < m; k++){
                    if (a[i][k] < a[i][j] && !(a[i2][k] < a[i2][j])) return false;
                }
            }
        }
    }

    return true;
}

/// Mixed differences -d[i][j] <= 0 over random row and column terms give every Monge matrix; small d makes ties
Matrix random_monge(int n, int m, int spread){
    Matrix a(n, vector<long long>(m));
    vector<long long> row(n), col(m);
    for (auto& x : row) x = stress::rand_int(-spread, spread);
    for (auto& x : col) x = stress::rand_int(-spread, spread);

    vector<vector<long long>> drop(n + 1, vector<long long>(m + 1, 0));
    for (int i = 1; i <= n; i++){
        for (int j = 1; j <= m; j++){
            long long d = stress::rand_int(0, 2) ? 0 : stress::rand_int(0, spread);
            drop[i][j] = drop[i - 1][j] + drop[i][j - 1] - drop[i - 1][j - 1] + d;
        }
    }
    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++) a[i][j] = row[i] + col[j] - drop[i][j];
    }

    return a;
}

/// g(b[j] - a[i]) for a convex g given by non-decreasing slopes, with sorted a and b
Matrix random_convex(int n, int m, int range){
    vector<long long> ra(n), rb(m), slope(4 * range + 1);
    for (auto& x : ra) x = stress::rand_int(-range, range);
    for (auto& x : rb) x = stress::rand_int(-range, range);
    for (auto& x : slope) x = stress::rand_int(-range, range);
    sort(ra.begin(), ra.end());
    sort(rb.begin(), rb.end());
    sort(slope.begin(), slope.end());

    vector<long long> g(4 * range + 2, 0);
    for (int x = 1; x < (int)g.size(); x++) g[x] = g[x - 1] + slope[x - 1];

    Matrix a(n, vector<long long>(m));
    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++) a[i][j] = g[rb[j] - ra[i] + 2 * range];
    }

    return a;
}

/// Picks a non-decreasing argmin per row, then fills each row so that column is its leftmost minimum
Matrix random_monotone(int n, int m, int spread){
    vector<int> best(n);
    for (auto& x : best) x = stress::rand_int(0, m - 1);
    sort(best.begin(), best.end());

    Matrix a(n, vector<long long>(m));
    for (int i = 0; i < n; i++){
        long long low = stress::rand_int(-spread, spread);
        for (int j = 0; j < m; j++) a[i][j] = low + stress::rand_int(j < best[i], spread);
        a[i][best[i]] = low;
    }

    return a;
}

void check(const Matrix& a, bool monge){
    int n = a.size(), m = a[0].size();
    auto f = [&](int i, int j){
        assert(0 <= i && i < n && 0 <= j && j < m);
        return a[i][j];
    };

    vector<int> expected = brute(a);
    assert(monotone_minima(n, m, f) == expected);
    if (monge) assert(smawk(n, m, f) == expected);
}

/// Counts calls of f on an n x m Monge matrix so the complexity bounds are checked
pair<long long, long long> calls(int n, int m){
    long long by_minima = 0, by_smawk = 0;
    auto counted = [n, m](long long* counter){
        return [counter, n, m](int i, int j){
            (*counter)++;
            long long centre = 1LL * i * m / n;
            return (j - centre) * (j - centre);
        };
    };
    vector<int> x = monotone_minima(n, m, counted(&by_minima)), y = smawk(n, m, counted(&by_smawk));
    assert(x == y);

    return {by_minima, by_smawk};
}

int main(){
    /// Every matrix up to 3 x 3 with entries 0..2: smawk on the totally monotone ones, monotone_minima on the rest too
    for (int n = 1; n <= 3; n++){
        for (int m = 1; m <= 3; m++){
            int total = 1;
            for (int t = 0; t < n * m; t++) total *= 3;
            for (int code = 0; code < total; code++){
                Matrix a(n, vector<long long>(m));
                for (int t = 0, c = code; t < n * m; t++, c /= 3) a[t / m][t % m] = c % 3;
                bool tm = totally_monotone(a);
                if (tm || monotone_argmin(a)) check(a, tm);
            }
        }
    }

    assert(monotone_minima(0, 3, [](int, int){ return 0; }).empty());
    assert(smawk(0, 3, [](int, int){ return 0; }).empty());

    for (int it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(1, 40), m = stress::rand_int(1, 40), spread = stress::rand_int(1, 1000);
        if (it % 3 == 0) check(random_monge(n, m, spread), true);
        else if (it % 3 == 1) check(random_convex(n, m, stress::rand_int(1, 30)), true);
        else check(random_monotone(n, m, spread), false);
    }

    /// Wide, tall and square shapes up to 5e5 cells, every row against the brute force
    for (auto [n, m] : vector<pair<int, int>>{{1, 2000}, {2000, 1}, {700, 700}, {3, 1500}, {1500, 3}}){
        check(random_monge(n, m, 1000000), true);
        check(random_convex(n, m, 1000), true);
        check(random_monotone(n, m, 5), false);
    }

    /// Complexity: smawk within 8 (n + m) calls, monotone_minima within (n + m)(log2 n + 2)
    for (auto [n, m] : vector<pair<int, int>>{{200000, 200000}, {1, 300000}, {300000, 1}, {100000, 7}, {7, 100000}}){
        auto [by_minima, by_smawk] = calls(n, m);
        assert(by_smawk <= 8LL * (n + m));
        assert(by_minima <= 1LL * (n + m) * (__lg(n) + 2));
    }

    stress::expect_assert("dp/monotone_minima.cpp", []{ monotone_minima(2, 0, [](int, int){ return 0; }); });
    stress::expect_assert("dp/monotone_minima.cpp", []{ smawk(2, 0, [](int, int){ return 0; }); });

    return 0;
}
