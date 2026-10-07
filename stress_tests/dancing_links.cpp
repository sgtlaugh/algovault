#include "common.h"

#define main library_main
#include "../code_library/dancing_links.cpp"
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

int main(){
    for (long long it = 0; it < stress::scaled(3000); it++){
        int nc = stress::rand_int(0, 12), nr = stress::rand_int(0, 14), full = (1 << nc) - 1;
        vector<int> masks;
        if (nc && stress::rand_int(0, 1)){
            /// Plant a solution: a random partition of the columns
            vector<int> part(nc);
            for (int c = 0; c < nc; c++) part[c] = stress::rand_int(0, c);
            for (int p = 0; p < nc; p++){
                int m = 0;
                for (int c = 0; c < nc; c++) if (part[c] == p) m |= 1 << c;
                if (m) masks.push_back(m);
            }
        }
        while ((int)masks.size() < nr) masks.push_back(stress::rand_int(0, 3) ? stress::rng()() & full : 0);
        shuffle(masks.begin(), masks.end(), stress::rng());

        vector<int> ids(masks.size());
        iota(ids.begin(), ids.end(), 1);
        shuffle(ids.begin(), ids.end(), stress::rng());

        dlx::init(nc);
        for (int i = 0; i < (int)masks.size(); i++){
            vector<int> columns;
            for (int c = 0; c < nc; c++) if (masks[i] >> c & 1) columns.push_back(c + 1);
            shuffle(columns.begin(), columns.end(), stress::rng());
            dlx::addrow(ids[i], columns);
        }

        vector<signed char> memo(1 << nc, -1);
        vector<int> rows;
        bool found = dlx::exact_cover(rows);
        assert(found == brute(0, full, masks, memo));

        if (found){
            int covered = 0;
            sort(rows.begin(), rows.end());
            assert(unique(rows.begin(), rows.end()) == rows.end());
            for (int id : rows){
                int m = masks[find(ids.begin(), ids.end(), id) - ids.begin()];
                assert(!(covered & m));
                covered |= m;
            }
            assert(covered == full);
        }
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
