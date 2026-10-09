#include "../common.h"

#define main library_main
#include "../../code_library/misc/dancing_links.cpp"
#undef main

/// Exact cover by always covering the lowest uncovered column, memoized on the covered set
bool brute(int covered, int full, const vector<int>& masks, vector<signed char>& memo){
    if (covered == full) return true;
    if (memo[covered] != -1) return memo[covered];

    int c = __builtin_ctz(~covered);
    bool res = false;
    for (int m : masks){
        if ((m >> c & 1) && !(m & covered) && brute(covered | m, full, masks, memo)){
            res = true;
            break;
        }
    }

    memo[covered] = res;
    return res;
}

bool valid_sudoku(const vector<vector<int>>& g, const vector<vector<int>>& givens){
    int n = g.size(), m = sqrt(n + 0.5);

    for (int i = 0; i < n; i++){
        vector<int> row(n + 1), col(n + 1), box(n + 1);
        for (int j = 0; j < n; j++){
            if (givens[i][j] && givens[i][j] != g[i][j]) return false;
            if (g[i][j] < 1 || g[i][j] > n) return false;
            if (row[g[i][j]]++ || col[g[j][i]]++ || box[g[(i / m) * m + j / m][(i % m) * m + j % m]]++) return false;
        }
    }
    return true;
}

/// A solved grid shuffled by digit, band, stack, row-in-band and column-in-stack permutations
vector<vector<int>> random_solved(int m){
    int n = m * m;
    vector<int> digits(n), rows, cols;
    iota(digits.begin(), digits.end(), 1);
    shuffle(digits.begin(), digits.end(), stress::rng());

    for (auto* order : {&rows, &cols}){
        vector<int> bands(m);
        iota(bands.begin(), bands.end(), 0);
        shuffle(bands.begin(), bands.end(), stress::rng());
        for (int b : bands){
            vector<int> inner(m);
            iota(inner.begin(), inner.end(), 0);
            shuffle(inner.begin(), inner.end(), stress::rng());
            for (int x : inner) order->push_back(b * m + x);
        }
    }

    vector<vector<int>> g(n, vector<int>(n));
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++) g[i][j] = digits[(rows[i] * m + rows[i] / m + cols[j]) % n];
    }
    return g;
}

struct Instance{
    int nc;
    vector<int> masks, ids;
    vector<vector<int>> columns;
};

Instance random_instance(){
    Instance in;
    in.nc = stress::rand_int(0, 12);
    int nr = stress::rand_int(0, 14), full = (1 << in.nc) - 1;
    if (in.nc && stress::rand_int(0, 1)){
        /// Plant a solution: a random partition of the columns
        vector<int> part(in.nc);
        for (int c = 0; c < in.nc; c++) part[c] = stress::rand_int(0, c);
        for (int p = 0; p < in.nc; p++){
            int m = 0;
            for (int c = 0; c < in.nc; c++) if (part[c] == p) m |= 1 << c;
            if (m) in.masks.push_back(m);
        }
    }
    while ((int)in.masks.size() < nr) in.masks.push_back(stress::rand_int(0, 3) ? stress::rng()() & full : 0);
    shuffle(in.masks.begin(), in.masks.end(), stress::rng());

    in.ids.resize(in.masks.size());
    iota(in.ids.begin(), in.ids.end(), 1);
    shuffle(in.ids.begin(), in.ids.end(), stress::rng());

    for (int m : in.masks){
        in.columns.emplace_back();
        for (int c = 0; c < in.nc; c++) if (m >> c & 1) in.columns.back().push_back(c + 1);
        shuffle(in.columns.back().begin(), in.columns.back().end(), stress::rng());
    }
    return in;
}

void check(const Instance& in, bool found, vector<int> rows){
    int full = (1 << in.nc) - 1;
    vector<signed char> memo(1 << in.nc, -1);
    assert(found == brute(0, full, in.masks, memo));
    if (!found) return;

    int covered = 0;
    sort(rows.begin(), rows.end());
    assert(unique(rows.begin(), rows.end()) == rows.end());
    for (int id : rows){
        int m = in.masks[find(in.ids.begin(), in.ids.end(), id) - in.ids.begin()];
        assert(!(covered & m));
        covered |= m;
    }
    assert(covered == full);
}

/// Two live instances with interleaved add_row calls, each checked against brute force, then re-solved to check the links were restored,
/// then one more row added after the solve and checked against brute force again
int main(){
    for (long long it = 0; it < stress::scaled(3000); it++){
        Instance a = random_instance(), b = random_instance();
        DancingLinks da(a.nc), db(b.nc, stress::rand_int(0, 20));
        for (size_t i = 0; i < max(a.masks.size(), b.masks.size()); i++){
            if (i < a.masks.size()) da.add_row(a.ids[i], a.columns[i]);
            if (i < b.masks.size()) db.add_row(b.ids[i], b.columns[i]);
        }

        vector<int> rows_a, rows_b, again;
        bool found_a = da.exact_cover(rows_a), found_b = db.exact_cover(rows_b);
        check(a, found_a, rows_a);
        check(b, found_b, rows_b);
        assert(da.exact_cover(again) == found_a);
        if (found_a) assert(again == rows_a);

        int extra = stress::rng()() & ((1 << a.nc) - 1);
        a.masks.push_back(extra);
        a.ids.push_back(a.ids.size() + 1);
        a.columns.emplace_back();
        for (int c = 0; c < a.nc; c++) if (extra >> c & 1) a.columns.back().push_back(c + 1);
        shuffle(a.columns.back().begin(), a.columns.back().end(), stress::rng());
        da.add_row(a.ids.back(), a.columns.back());
        found_a = da.exact_cover(rows_a);
        check(a, found_a, rows_a);
    }

    for (long long it = 0; it < stress::scaled(150); it++){
        int m = it % 10 ? stress::rand_int(1, 3) : 4, n = m * m;
        auto givens = random_solved(m);
        int blank = m < 4 ? stress::rand_int(0, 100) : stress::rand_int(0, 1) ? stress::rand_int(0, 40) : 100;  /// sparse 16 x 16 givens can make search exponential
        for (auto& row : givens) for (auto& x : row) if (stress::rand_int(0, 99) < blank) x = 0;
        bool unsolvable = n > 1 && stress::rand_int(0, 3) == 0;
        if (unsolvable){
            int i = stress::rand_int(0, n - 1), j = stress::rand_int(0, n - 1), k = (j + stress::rand_int(1, n - 1)) % n;
            givens[i][j] = givens[i][k] = stress::rand_int(1, n);  /// a row repeating a digit
        }

        auto g = givens;
        assert(sudoku::solve(g) == !unsolvable);
        if (!unsolvable) assert(valid_sudoku(g, givens));
    }

    return 0;
}
