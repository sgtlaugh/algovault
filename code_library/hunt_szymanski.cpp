/***
 * Hunt-Szymanski Algorithm for LCS
 *
 * Complexity: O(R + N) log N
 * R = numbered of ordered pairs of positions where the two strings match (worst case, R = N^2)
 *
***/

#include <bits/stdc++.h>

#define MAX 50010

using namespace std;

char A[MAX], B[MAX];

int lcs(char* A, char* B){
    vector <int> adj[256], ar(1, -1);
    int i, j, n = strlen(A), m = strlen(B);
    for (i = 0; i < m; i++) adj[(unsigned char)B[i]].push_back(i);

    for (i = 0; i < n; i++){
        const auto& pos = adj[(unsigned char)A[i]];
        for (j = (int)pos.size() - 1; j >= 0; j--){
            int x = pos[j];
            if (x > ar.back()) ar.push_back(x);
            else *lower_bound(ar.begin(), ar.end(), x) = x;
        }
    }
    return ar.size() - 1;
}

int main(){
    mt19937 rng(666);

    int n = MAX - 10, m = MAX - 10;
    for (int i = 0; i < n; i++) A[i] = (rng() % 26) + 'A';
    for (int i = 0; i < m; i++) B[i] = (rng() % 26) + 'A';
    A[n] = B[m] = 0;

    clock_t start = clock();
    assert(lcs(A, B) == 16259);

    fprintf(stderr, "\nTime taken = %0.6f\n", (clock() - start) / (1.0 * CLOCKS_PER_SEC)); /// took 3.96 s locally
    return 0;
}
