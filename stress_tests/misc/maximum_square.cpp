#include "../common.h"

#define main library_main
#include "../../code_library/misc/maximum_square.cpp"
#undef main

/// Grow each shape from size 1 while every one of its cells is filled
int main(){
    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(1, 12), m = stress::rand_int(1, 12), fill = stress::rand_int(30, 100);
        vector<vector<int>> g(n, vector<int>(m));
        for (auto& row : g){
            for (auto& x : row) x = stress::rand_int(1, 100) <= fill;
        }
        auto sq = square_sizes(g), dm = diamond_sizes(g);
        auto filled = [&](int i, int j){ return i >= 0 && j >= 0 && i < n && j < m && g[i][j]; };

        for (int i = 0; i < n; i++){
            for (int j = 0; j < m; j++){
                int k = 0;
                while (true){
                    bool ok = true;
                    for (int a = 0; a <= k && ok; a++){
                        for (int b = 0; b <= k && ok; b++) ok = filled(i - a, j - b);
                    }
                    if (!ok) break;
                    k++;
                }
                assert(sq[i][j] == k);

                k = 0;
                while (true){
                    int r = k, ci = i - r;
                    bool ok = true;
                    for (int a = ci - r; a <= ci + r && ok; a++){
                        for (int b = j - r; b <= j + r && ok; b++){
                            if (abs(a - ci) + abs(b - j) <= r) ok = filled(a, b);
                        }
                    }
                    if (!ok) break;
                    k++;
                }
                assert(dm[i][j] == k);
            }
        }
    }
    return 0;
}
