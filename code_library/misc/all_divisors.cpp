/***
 *
 * Generates all divisors from 1 to MAX efficiently
 *
 * The idea is to generate lp[x], the largest prime factor of x using sieve
 * Then generate the divisors iteratively
 * To get divisors of x, use lp[x] and the divisors of [x / lp[x]]
 * This is order of magnitudes faster than the naive approach because of reduced cache miss
 * Takes 0.8 seconds when MAX=10^7 locally, where as the naive approach takes 11.5 seconds
 *
 * Naive approach:
 *
 * vector <int> divisors[MAX];
 *
 * void generate(){
 *     for (int i = 1; i < MAX; i++){
 *         for (int j = i; j < MAX; j += i){
 *              divisors[j].push_back(i);
 *         }
 *     }
 * }
 *
 * Note, with MAX=10^7 it might crash if not sufficient memory available to store the divisors
 *
***/

#include <bits/stdc++.h>

#define MAX 10000010

#pragma GCC optimize("Ofast,no-stack-protector")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,tune=native")
#pragma GCC optimize("unroll-loops")

using namespace std;

vector <int> divisors[MAX];

void gen(bool make_sorted=false){ /// sorting makes it 2x slower
    short sp[MAX];
    int i, j, k, c, l, v, lp[MAX];

    sp[1] = lp[1] = 1;
    for (i = 4; i < MAX; i++, i++) sp[i] = 2;

    for (i = 3; i * i < MAX; i += 2){
        for (j = i * i; j < MAX && !sp[i]; j += 2 * i){
            if (!sp[j]) sp[j] = i;
        }
    }

    for (i = 2; i < MAX; i++){
        if (!sp[i]) lp[i] = i;
        else lp[i] = max((int)sp[i], lp[i /sp[i]]);
    }

    divisors[1].push_back(1);
    for (i = 2; i < MAX; i++){
        c = 1, l = 0;
        for (k = i; k > 1 && lp[k] == lp[i]; c++) k /= lp[k];
        divisors[i].resize(c * divisors[k].size());

        for (v = 1, j = 0; j < c; j++, v *= lp[i]){
            for (const auto d: divisors[k]){
                divisors[i][l++] = d * v;
            }
        }
        if (make_sorted) sort(divisors[i].begin(), divisors[i].end());
    }
}

int main(){
    clock_t start = clock();
    gen();

    assert(divisors[1] == vector<int>({1}));
    assert(divisors[2] == vector<int>({1, 2}));
    assert(divisors[10007] == vector<int>({1, 10007}));

    assert(divisors[8].size() == 4);
    assert(divisors[24].size() == 8);
    assert(divisors[840000].size() == 140);

    fprintf(stderr, "\nTime taken = %0.6f\n", (clock() - start) / (1.0 * CLOCKS_PER_SEC));  /// 0.831455
    return 0;
}
