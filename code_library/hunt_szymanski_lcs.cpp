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

int ar[MAX];
char A[MAX], B[MAX];

int lcs(char* A, char* B){
    vector <int> adj[256];
    int i, j, l = 0, n = strlen(A), m = strlen(B);
    for (i = 0; i < m; i++) adj[(int)B[i]].push_back(i);

    ar[l++] = -1;
    for (i = 0; i < n; i++){
        for (j = (int)adj[(int)A[i]].size() - 1; j >= 0; j--){
            int x = adj[(int)A[i]][j];
            if (x > ar[l - 1]) ar[l++] = x;
            else ar[lower_bound(ar, ar + l, x) - ar] = x;
        }
    }
    return l - 1;
}

int main(){
    srand(666);

    int n = MAX - 10, m = MAX - 10;
    for (int i = 0; i < n; i++) A[i] = (rand() % 26) + 'A';
    for (int i = 0; i < m; i++) B[i] = (rand() % 26) + 'A';
    A[n] = B[m] = 0;

    clock_t start = clock();
    assert(lcs(A, B) == 16259);

    fprintf(stderr, "\nTime taken = %0.6f\n", (clock() - start) / (1.0 * CLOCKS_PER_SEC)); /// took 3.96 s locally
    return 0;
}
