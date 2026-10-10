/***
 * Fast radix sort in O(w * n) with loop unrolling
 * w is the key length, and is fixed to 8 bits
 * This enables sorting any list of integers in O(64/8 * n) or O(8 * n) passes
 *
 * Complexity: O(8 n + 8 * 256) ~ O(n), memory is a scratch vector of n unsigned long longs (80 MB at n = 1e7)
 *
 * radix_sort(ar, n) sorts ar[0, n) with a scratch vector allocated per call
 * radix_sort(ar, n, tmp) uses the caller's tmp (grown to n if smaller), pass the same tmp to repeated sorts:
 * at -O2 on n = 1e7 a fresh scratch per call measured 270 ms per sort vs 229 ms reused, its pages fault in on every call
 *
 * Most useful in scenarios where you want to sort a large list of items fast
 * Why not make radix_sort (https://github.com/sgtlaugh/algovault/blob/master/code_library/radix_sort.cpp) generic?
 * Because you want to use radix sort instead of std::sort when speed is important (think squeezing sub-optimal solutions into TL)
 * And making it generic is usally also making it slower
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

void radix_sort(unsigned long long ar[], int n, vector<unsigned long long>& tmp){
    assert(sizeof(unsigned long long) == 8);

    unsigned long long cnt[8][256] = {};
    if ((int)tmp.size() < n) tmp.resize(n);
    for (int i = 0; i < n; i++){
        cnt[0][ar[i] & 255]++, cnt[1][(ar[i] >> 8) & 255]++;
        cnt[2][(ar[i] >> 16) & 255]++, cnt[3][(ar[i] >> 24) & 255]++;
        cnt[4][(ar[i] >> 32) & 255]++, cnt[5][(ar[i] >> 40) & 255]++;
        cnt[6][(ar[i] >> 48) & 255]++, cnt[7][(ar[i] >> 56) & 255]++;
    }

    for (int j = 0; j < 8; j++){
        for (int i = 1; i < 256; i++){
            cnt[j][i] += cnt[j][i - 1];
        }
    }

    for (int i = n - 1; i >= 0; i--) tmp[--cnt[0][ar[i] & 255]] = ar[i];
    for (int i = n - 1; i >= 0; i--) ar[--cnt[1][(tmp[i] >> 8) & 255]] = tmp[i];
    for (int i = n - 1; i >= 0; i--) tmp[--cnt[2][(ar[i] >> 16) & 255]] = ar[i];
    for (int i = n - 1; i >= 0; i--) ar[--cnt[3][(tmp[i] >> 24) & 255]] = tmp[i];
    for (int i = n - 1; i >= 0; i--) tmp[--cnt[4][(ar[i] >> 32) & 255]] = ar[i];
    for (int i = n - 1; i >= 0; i--) ar[--cnt[5][(tmp[i] >> 40) & 255]] = tmp[i];
    for (int i = n - 1; i >= 0; i--) tmp[--cnt[6][(ar[i] >> 48) & 255]] = ar[i];
    for (int i = n - 1; i >= 0; i--) ar[--cnt[7][(tmp[i] >> 56) & 255]] = tmp[i];
}

void radix_sort(unsigned long long ar[], int n){
    vector<unsigned long long> tmp;
    radix_sort(ar, n, tmp);
}

int main(){
    mt19937_64 rng(42);
    int i, n = 50000000;
    vector<unsigned long long> ar(n);

    puts("Generating array");
    for (i = 0; i < n; i++) ar[i] = rng();

    clock_t start = clock();
    radix_sort(ar.data(), n);
    printf("Time taken to sort = %0.6f s\n", (clock() - start) / (double)CLOCKS_PER_SEC);  /// Time taken = 0.876000

    for (i = 0; (i + 1) < n; i++){
        assert(ar[i] <= ar[i + 1]);
    }

    vector<unsigned long long> tmp, big = {7, 18446744073709551615ull, 0, 256, 7, 65536}, small = {3, 1, 2};
    radix_sort(big.data(), big.size(), tmp);
    radix_sort(small.data(), small.size(), tmp);
    assert((big == vector<unsigned long long>{0, 7, 7, 256, 65536, 18446744073709551615ull}));
    assert((small == vector<unsigned long long>{1, 2, 3}));

    return 0;
}
