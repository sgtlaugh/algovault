/***
 *
 * Josephus problem, n people numbered from 1 to n stand in a circle.
 * Counting starts from 1 and every k'th person dies
 * Returns the id of the m'th killed person
 * For example if n = 10 and k = 3, then the people killed are 3, 6, 9, 2, 7, 1, 8, 5, 10, 4 respectively
 *
***/

#include <bits/stdc++.h>

using namespace std;

/// O(n)
int josephus1(int n, int k, int m){
    int i;
    for (m = n - m, i = m + 1; i <= n; i++){
        m += k;
        if (m >= i) m %= i;
    }

    return m + 1;
}

/// O(k log(n))
long long josephus2(long long n, long long k, long long m){
    m = n - m;
    if (k <= 1) return n - m;

    long long i = m;
    while (i < n){
        long long r = (i - m + k - 2) / (k - 1);
        if ((i + r) > n) r = n - i;
        else if (!r) r = 1;
        i += r;
        m = (m + (r * k)) % i;
    }

    return m + 1;
}

int main(){
    assert(josephus1(100, 13, 50) == 61);
    assert(josephus2(100, 13, 50) == 61);

    for (int n = 1; n <= 50; n++){
        for (int k = 1; k <= 50; k++){
            for (int m = 1; m <= n; m++){
                assert(josephus1(n, k, m) == josephus2(n, k, m));
            }
        }
    }

    return 0;
}
