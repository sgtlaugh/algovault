#include "../common.h"

#define main library_main
#include "../../code_library/graphs/hungarian_algorithm.cpp"
#undef main

/// Best sum over every way to give each row of the smaller side a distinct column
template <class T>
T brute(const vector<vector<T>>& mat, bool minimize){
    int row = mat.size(), col = mat[0].size();
    if (row > col){
        vector<vector<T>> t(col, vector<T>(row));
        for (int i = 0; i < row; i++) for (int j = 0; j < col; j++) t[j][i] = mat[i][j];
        return brute(t, minimize);
    }

    vector<int> perm(col);
    iota(perm.begin(), perm.end(), 0);
    T best = minimize ? numeric_limits<T>::max() : numeric_limits<T>::lowest();
    do{
        T total = 0;
        for (int i = 0; i < row; i++) total += mat[i][perm[i]];
        best = minimize ? min(best, total) : max(best, total);
    } while (next_permutation(perm.begin(), perm.end()));
    return best;
}

template <class T>
void check_matching(const vector<vector<T>>& mat, const pair<T, vector<pair<int, int>>>& res){
    int row = mat.size(), col = mat[0].size();
    auto& matches = res.second;
    assert((int)matches.size() == min(row, col));

    set<int> rows, cols;
    T total = 0;
    for (auto [i, j] : matches){
        assert(0 <= i && i < row && 0 <= j && j < col);
        rows.insert(i), cols.insert(j);
        total += mat[i][j];
    }
    assert(rows.size() == matches.size() && cols.size() == matches.size());
    if constexpr (is_floating_point_v<T>) assert(fabs(total - res.first) < 1e-6);
    else assert(total == res.first);
}

template <class T>
vector<vector<T>> random_matrix(int row, int col, long long lo, long long hi){
    vector<vector<T>> mat(row, vector<T>(col));
    for (auto& r : mat) for (auto& x : r) x = stress::rand_int(lo, hi);
    return mat;
}

int main(){
    for (long long it = 0; it < stress::scaled(3000); it++){
        int row = stress::rand_int(1, 6), col = stress::rand_int(1, 6);
        bool minimize = stress::rand_int(0, 1);

        auto mat = random_matrix<long long>(row, col, it % 3 ? 0 : -1000000000000LL, it % 3 ? 9 : 1000000000000LL);
        auto res = hungarian(mat, minimize);
        assert(res.first == brute(mat, minimize));
        check_matching(mat, res);

        auto small = random_matrix<int>(row, col, -50, 50);
        auto res_int = hungarian(small, minimize);
        assert(res_int.first == brute(small, minimize));
        check_matching(small, res_int);

        vector<vector<double>> real(row, vector<double>(col));
        for (auto& r : real) for (auto& x : r) x = stress::rand_int(-1000000, 1000000) / 1024.0;
        auto res_real = hungarian(real, minimize);
        assert(fabs(res_real.first - brute(real, minimize)) < 1e-6);
        check_matching(real, res_real);
    }

    /// Too large to enumerate: transposing must not change the optimum, and the two take different paths
    for (long long it = 0; it < stress::scaled(40); it++){
        int row = stress::rand_int(1, 60), col = stress::rand_int(1, 60);
        /// Costs in [0, 2] tie everywhere, so free columns often tie with matched ones at the closest distance
        auto mat = random_matrix<long long>(row, col, it % 2 ? -1000 : 0, it % 2 ? 1000 : 2);
        vector<vector<long long>> t(col, vector<long long>(row));
        for (int i = 0; i < row; i++) for (int j = 0; j < col; j++) t[j][i] = mat[i][j];

        bool minimize = stress::rand_int(0, 1);
        auto res = hungarian(mat, minimize);
        check_matching(mat, res);
        assert(res.first == hungarian(t, minimize).first);
    }

    /// A tall input costs O(row^2 * col) unless the solver transposes it internally
    for (int minimize = 0; minimize < 2; minimize++){
        auto tall = random_matrix<long long>(3000, 50, -1000000, 1000000);
        vector<vector<long long>> wide(50, vector<long long>(3000));
        for (int i = 0; i < 3000; i++) for (int j = 0; j < 50; j++) wide[j][i] = tall[i][j];

        auto res = hungarian(tall, minimize);
        check_matching(tall, res);
        assert(is_sorted(res.second.begin(), res.second.end()));
        assert(res.first == hungarian(wide, minimize).first);
    }

    /// At the documented limit, 2 max(n, m) max |mat| fits in T, int cells must agree with the same cells in long long
    for (long long it = 0; it < stress::scaled(20000); it++){
        int row = stress::rand_int(1, 8), col = stress::rand_int(1, 8), limit = INT_MAX / (2 * max(row, col));
        vector<vector<int>> mat(row, vector<int>(col));
        vector<vector<long long>> wide(row, vector<long long>(col));
        for (int i = 0; i < row; i++){
            for (int j = 0; j < col; j++){
                int pick = stress::rand_int(0, 2);
                mat[i][j] = wide[i][j] = pick == 0 ? limit : pick == 1 ? -limit : stress::rand_int(-limit, limit);
            }
        }

        bool minimize = stress::rand_int(0, 1);
        auto res = hungarian(mat, minimize);
        check_matching(mat, res);
        assert(res.first == hungarian(wide, minimize).first);
    }

    /// 1000 x 1000 with a planted zero-cost permutation and every other cell positive: it is the unique optimum
    for (int top : {1000000, 2}){
        int n = 1000;
        vector<int> planted(n);
        iota(planted.begin(), planted.end(), 0);
        shuffle(planted.begin(), planted.end(), stress::rng());
        auto mat = random_matrix<int>(n, n, 1, top);
        for (int i = 0; i < n; i++) mat[i][planted[i]] = 0;

        auto [total, pairs] = hungarian(mat);
        assert(total == 0);
        for (int i = 0; i < n; i++) assert((pairs[i] == make_pair(i, planted[i])));
    }

    assert(hungarian(vector<vector<int>>{}).first == 0);
    return 0;
}
