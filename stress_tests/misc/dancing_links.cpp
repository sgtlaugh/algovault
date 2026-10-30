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

/// Number of exact covers over every subset of rows from i on, skipping empty rows since add_row ignores them
long long brute_count(int i, int covered, int full, const vector<int>& masks){
    if (i == (int)masks.size()) return covered == full;

    long long res = brute_count(i + 1, covered, full, masks);
    if (masks[i] && !(masks[i] & covered)) res += brute_count(i + 1, covered | masks[i], full, masks);

    return res;
}

/// Placements of queens in rows r..n-1 attacking no queen in an earlier row
long long brute_queens(int n, int r, vector<int>& queen_col){
    if (r == n) return 1;

    long long res = 0;
    for (int c = 0; c < n; c++){
        bool free = true;
        for (int pr = 0; pr < r; pr++) free = free && queen_col[pr] != c && abs(queen_col[pr] - c) != r - pr;
        if (!free) continue;
        queen_col[r] = c;
        res += brute_queens(n, r + 1, queen_col);
    }

    return res;
}

/// N queens as exact cover: rows and columns exactly once, each diagonal at most once through a slack row
DancingLinks queens_cover(int n){
    int diagonals = 2 * n - 1;
    DancingLinks dlx(2 * n + 2 * diagonals);
    for (int r = 0; r < n; r++){
        for (int c = 0; c < n; c++){
            dlx.add_row(r * n + c, {r, n + c, 2 * n + r + c, 2 * n + diagonals + r - c + n - 1});
        }
    }
    for (int d = 2 * n; d < 2 * n + 2 * diagonals; d++) dlx.add_row(n * n + d, {d});

    return dlx;
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
        for (int c = 0; c < in.nc; c++) if (m >> c & 1) in.columns.back().push_back(c);
        shuffle(in.columns.back().begin(), in.columns.back().end(), stress::rng());
    }

    return in;
}

void check(const Instance& in, optional<vector<int>> rows){
    int full = (1 << in.nc) - 1;
    vector<signed char> memo(1 << in.nc, -1);
    assert(rows.has_value() == brute(0, full, in.masks, memo));
    if (!rows) return;

    int covered = 0;
    sort(rows->begin(), rows->end());
    assert(unique(rows->begin(), rows->end()) == rows->end());
    for (int id : *rows){
        int m = in.masks[find(in.ids.begin(), in.ids.end(), id) - in.ids.begin()];
        assert(!(covered & m));
        covered |= m;
    }
    assert(covered == full);
}

/// Two live instances with interleaved add_row calls, each checked against brute force
int main(){
    for (long long it = 0; it < stress::scaled(3000); it++){
        Instance a = random_instance(), b = random_instance();
        DancingLinks da(a.nc), db(b.nc, stress::rand_int(0, 20));
        for (size_t i = 0; i < max(a.masks.size(), b.masks.size()); i++){
            if (i < a.masks.size()) da.add_row(a.ids[i], a.columns[i]);
            if (i < b.masks.size()) db.add_row(b.ids[i], b.columns[i]);
        }

        auto rows_a = da.exact_cover(), rows_b = db.exact_cover();
        check(a, rows_a);
        check(b, rows_b);

        long long covers_a = brute_count(0, 0, (1 << a.nc) - 1, a.masks), big = stress::rand_int(3, 1 << 14);
        assert(da.count(1) == min(covers_a, 1LL));
        assert(da.count(2) == min(covers_a, 2LL));
        assert(da.count(big) == min(covers_a, big));
        assert(da.count() == covers_a);
        assert(da.count(0) == 0);
        assert(db.count() == brute_count(0, 0, (1 << b.nc) - 1, b.masks));

        auto again = da.exact_cover();
        assert(again == rows_a);  /// a solve after counts finds the same rows only if every call restored the links

        /// A row added after a solve must link into the restored structure
        int extra = stress::rng()() & ((1 << a.nc) - 1);
        a.masks.push_back(extra);
        a.ids.push_back(a.ids.size() + 1);
        a.columns.emplace_back();
        for (int c = 0; c < a.nc; c++) if (extra >> c & 1) a.columns.back().push_back(c);
        shuffle(a.columns.back().begin(), a.columns.back().end(), stress::rng());
        da.add_row(a.ids.back(), a.columns.back());
        rows_a = da.exact_cover();
        check(a, rows_a);
        assert(da.count() == brute_count(0, 0, (1 << a.nc) - 1, a.masks));
    }

    /// Limits below, at and above the true count, where counting must stop exactly at the limit
    for (int n = 1; n <= 9; n++){
        vector<int> queen_col(n);
        long long covers = brute_queens(n, 0, queen_col);
        DancingLinks dlx = queens_cover(n);
        assert(dlx.count() == covers);
        for (int it = 0; it < 20; it++){
            long long limit = stress::rand_int(0, covers + 2);
            assert(dlx.count(limit) == min(covers, limit));
        }
    }

    /// The empty matrix has exactly one cover, the empty one
    DancingLinks nothing(0);
    assert(nothing.count() == 1);
    assert(nothing.count(2) == 1);

    stress::expect_assert("misc/dancing_links.cpp", []{
        DancingLinks dlx(1);
        dlx.count(-1);
    });

    for (long long it = 0; it < stress::scaled(150); it++){
        int m = it % 10 ? stress::rand_int(1, 3) : 4, n = m * m;
        auto givens = random_solved(m);
        /// Sparse 16 x 16 givens can make the search exponential
        int blank = m < 4 ? stress::rand_int(0, 100) : stress::rand_int(0, 1) ? stress::rand_int(0, 40) : 100;
        for (auto& row : givens) for (auto& x : row) if (stress::rand_int(0, 99) < blank) x = 0;
        bool unsolvable = n > 1 && stress::rand_int(0, 3) == 0;
        if (unsolvable){
            int i = stress::rand_int(0, n - 1), j = stress::rand_int(0, n - 1), k = (j + stress::rand_int(1, n - 1)) % n;
            givens[i][j] = givens[i][k] = stress::rand_int(1, n);  /// a row repeating a digit
        }

        auto before = givens;
        auto solved = sudoku::solution(givens);
        assert(solved.has_value() == !unsolvable);
        assert(givens == before);
        if (!unsolvable) assert(valid_sudoku(*solved, givens));
    }

    /// Guards the header's stack claim: one frame per selected row, CI's flags overflow 8 MB a little past depth 4e4
    int depth = 30000;
    DancingLinks chain(depth, depth);
    for (int i = 1; i <= depth; i++) chain.add_row(i, {i - 1});
    auto rows = chain.exact_cover();
    assert(rows && (int)rows->size() == depth);

    return 0;
}
