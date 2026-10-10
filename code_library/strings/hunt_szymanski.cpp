/***
 * Hunt-Szymanski Algorithm for LCS
 *
 * Complexity: O(R + N) log N
 * R = numbered of ordered pairs of positions where the two strings match (worst case, R = N^2)
 * For dense inputs (small alphabet, large R) bit_string_lcs.cpp is faster at O(N * M / 64)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

int lcs(const char* A, const char* B){
    vector<int> adj[256], ar(1, -1);
    int i, j, k, n = strlen(A), m = strlen(B);
    for (i = 0; i < m; i++) adj[(unsigned char)B[i]].push_back(i);

    for (i = 0; i < n; i++){
        const auto& pos = adj[(unsigned char)A[i]];
        /// Positions come in decreasing order, so each one lands at or before the slot of the previous
        for (j = (int)pos.size() - 1, k = ar.size(); j >= 0; j--){
            int x = pos[j];
            if (x > ar.back()) ar.push_back(x), k = ar.size();
            else{
                k = lower_bound(ar.begin(), ar.begin() + k, x) - ar.begin();
                ar[k++] = x;
            }
        }
    }

    return ar.size() - 1;
}

int main(){
    mt19937 rng(666);

    int n = 50000, m = 50000;
    string A(n, 0), B(m, 0);
    for (int i = 0; i < n; i++) A[i] = (rng() % 26) + 'A';
    for (int i = 0; i < m; i++) B[i] = (rng() % 26) + 'A';

    clock_t start = clock();
    assert(lcs(A.c_str(), B.c_str()) == 16259);

    fprintf(stderr, "\nTime taken = %0.6f\n", (clock() - start) / (1.0 * CLOCKS_PER_SEC)); /// took 3.96 s locally
    return 0;
}
