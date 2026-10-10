/***
 * Fast radix sort in O(w * n) with loop unrolling
 * w is the key length, and is fixed to 8 bits
 * This enables sorting any list of integers in O(32/8 * n) or O(4 * n) passes
 *
 * Complexity: O(n + 4 * 256), one counting pass and 4 scatter passes over the array
 * Memory: a scratch vector of n unsigned ints, about 40 MB for n = 1e7
 *
 * radix_sort(ar, n) sorts ar[0, n) with a scratch vector allocated per call
 * radix_sort(ar, n, tmp) uses the caller's tmp (grown to n if smaller), pass the same tmp to repeated sorts:
 * at -O2 on n = 1e7 a fresh scratch per call measured 98 ms per sort vs 78 ms reused, its pages fault in on every call
 *
 * Most useful in scenarios where you want to sort a large list of items fast, usually in sub-optimal solutions
 * Can be generalized with templates but that usually makes it 2-3 x slower
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

void radix_sort(unsigned int ar[], int n, vector<unsigned int>& tmp){
    assert(sizeof(unsigned int) == 4);

    unsigned int cnt[4][256] = {};
    if ((int)tmp.size() < n) tmp.resize(n);
    for (int i = 0; i < n; i++){
        cnt[0][ar[i] & 255]++;
        cnt[1][(ar[i] >> 8) & 255]++;
        cnt[2][(ar[i] >> 16) & 255]++;
        cnt[3][(ar[i] >> 24) & 255]++;
    }

    for (int j = 0; j < 4; j++){
        for (int i = 1; i < 256; i++){
            cnt[j][i] += cnt[j][i - 1];
        }
    }

    for (int i = n - 1; i >= 0; i--) tmp[--cnt[0][ar[i] & 255]] = ar[i];
    for (int i = n - 1; i >= 0; i--) ar[--cnt[1][(tmp[i] >> 8) & 255]] = tmp[i];
    for (int i = n - 1; i >= 0; i--) tmp[--cnt[2][(ar[i] >> 16) & 255]] = ar[i];
    for (int i = n - 1; i >= 0; i--) ar[--cnt[3][(tmp[i] >> 24) & 255]] = tmp[i];
}

void radix_sort(unsigned int ar[], int n){
    vector<unsigned int> tmp;
    radix_sort(ar, n, tmp);
}

int main(){
    mt19937 rng(42);
    int i, n = 100000000;
    vector<unsigned int> ar(n);

    puts("Generating array");
    for (i = 0; i < n; i++) ar[i] = rng();

    clock_t start = clock();
    radix_sort(ar.data(), n);
    printf("Time taken to sort = %0.6f s\n", (clock() - start) / (double)CLOCKS_PER_SEC);  /// Time taken = 0.542000

    for (i = 0; (i + 1) < n; i++){
        assert(ar[i] <= ar[i + 1]);
    }

    vector<unsigned int> tmp, big = {7, 4294967295u, 0, 256, 7, 65536}, small = {3, 1, 2};
    radix_sort(big.data(), big.size(), tmp);
    radix_sort(small.data(), small.size(), tmp);
    assert((big == vector<unsigned int>{0, 7, 7, 256, 65536, 4294967295u}));
    assert((small == vector<unsigned int>{1, 2, 3}));

    return 0;
}
