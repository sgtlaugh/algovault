#include <stdio.h>
#include <bits/stdtr1c++.h>


using namespace std;

/***
 *
 * returns the rank of the permutation
 * rank and permutation elements follows 1 based indexing
 *
***/

long long find_rank(vector <int> permutation){
    long long res = 1;
    int i, j, n = permutation.size();
    vector <bool> visited(n + 1, false);
    vector <long long> factorial(n + 1, 1);
    for (i = 1; i <= n; i++) factorial[i] = factorial[i - 1] * i;

    for (i = 1; i <= n; i++){
        for (j = 1; j <= n; j++){
            if (!visited[j]){
                if (permutation[i - 1] == j){
                    visited[j] = true;
                    break;
                }
                res += factorial[n - i];
            }
        }
    }
    return res;
}


/***
 *
 * returns the k'th permutation from 1 to n
 * k must be between 1 and n!
 *
***/

vector <int> find_permutation(int n, long long k){
    int i, j;
    vector <int> res(n, 0);
    vector <bool> visited(n + 1, false);
    vector <long long> factorial(n + 1, 1);
    for (i = 1; i <= n; i++) factorial[i] = factorial[i - 1] * i;

    for (i = 1; i <= n; i++){
        for (j = 1; j <= n; j++){
            if (!visited[j]){
                if (factorial[n - i] >= k){
                    visited[j] = true;
                    res[i - 1] = j;
                    break;
                }
                k -= factorial[n - i];
            }
        }
    }
    return res;
}

int main(){
    printf("%lld\n", find_rank({1, 2, 3}));      /// 1
    printf("%lld\n", find_rank({1, 3, 2}));      /// 2
    printf("%lld\n", find_rank({2, 1, 3}));      /// 3
    printf("%lld\n", find_rank({2, 3, 1}));      /// 4
    printf("%lld\n", find_rank({3, 1, 2}));      /// 5
    printf("%lld\n", find_rank({3, 2, 1}));      /// 6

    vector <int> permutation = find_permutation(15, 1000000000007LL);
    for (auto v: permutation) printf("%d ", v);  /// 12 7 9 11 2 3 15 5 6 14 4 13 10 1 8
    puts("");

    return 0;
}
