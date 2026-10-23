/***
 *
 * All Divisors
 * Generates the divisors of every number from 1 to n
 *
 * Complexity:
 *   - O(n log n) time and memory, the lists hold the sum of d(x) over x <= n, about n ln n divisors
 *   - make_sorted adds O(d log d) per list of d divisors
 *   - about 1.1 GB peak memory at n = 10^7
 *
 * all_divisors(n, make_sorted = false)[x] = the divisors of x for 1 <= x <= n, index 0 is an empty list
 *
 * Notes:
 *   - 0.8 s at n = 10^7, the naive loop over multiples misses the cache far more and takes 11.5 s
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

vector<vector<int>> all_divisors(int n, bool make_sorted=false){ /// sorting makes it 2x slower
    vector<vector<int>> divisors(max(n, 0) + 1);
    if (n < 1) return divisors;

    vector<short> sp(n + 1, 0);
    vector<int> lp(n + 1, 0);
    int i, j, k, c, l;
    long long v;

    sp[1] = lp[1] = 1;
    for (i = 4; i <= n; i++, i++) sp[i] = 2;

    for (i = 3; (long long)i * i <= n; i += 2){
        for (j = i * i; j <= n && !sp[i]; j += 2 * i){
            if (!sp[j]) sp[j] = i;
        }
    }

    for (i = 2; i <= n; i++){
        if (!sp[i]) lp[i] = i;
        else lp[i] = max((int)sp[i], lp[i /sp[i]]);
    }

    divisors[1].push_back(1);
    for (i = 2; i <= n; i++){
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

    return divisors;
}

int main(){
    clock_t start = clock();
    auto divisors = all_divisors(10000000);

    assert(divisors[1] == vector<int>({1}));
    assert(divisors[2] == vector<int>({1, 2}));
    assert(divisors[10007] == vector<int>({1, 10007}));

    assert(divisors[8].size() == 4);
    assert(divisors[24].size() == 8);
    assert(divisors[840000].size() == 140);

    assert(all_divisors(12, true)[12] == vector<int>({1, 2, 3, 4, 6, 12}));
    assert(all_divisors(1).size() == 2 && all_divisors(0).size() == 1);

    fprintf(stderr, "\nTime taken = %0.6f\n", (clock() - start) / (1.0 * CLOCKS_PER_SEC));  /// 0.831455
    return 0;
}
