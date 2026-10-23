/***
 *
 * Radix Sort
 * LSD radix sort with 8-bit digits, 4 counting-sort passes over 32-bit integers, loop unrolled
 *
 * Complexity: O(n + 1024), a scratch vector of n unsigned ints (40 MB at n = 1e7)
 *
 * radix_sort(ar, n): sorts unsigned ints ar[0, n) with a scratch vector allocated per call
 * radix_sort(v): sorts a vector<int>, negatives first
 * radix_sort(ar, n, tmp): sorts with the caller's tmp as scratch, grown to n if smaller
 *
 * Notes:
 *   - Pass the same tmp to repeated sorts: at n = 1e7 a fresh scratch measured 98 ms per sort, a reused one 78 ms
 *   - Use radix_sort_64.cpp for 64-bit keys, a template over the key type made it 2-3x slower
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

/// Flipping the sign bit maps int order onto unsigned order
void radix_sort(vector<int>& v){
    for (int& x : v) x ^= INT_MIN;
    radix_sort((unsigned int*)v.data(), v.size());
    for (int& x : v) x ^= INT_MIN;
}

int main(){
    mt19937 rng(42);
    int i, n = 10000000;
    vector<unsigned int> ar(n);

    puts("Generating array");
    for (i = 0; i < n; i++) ar[i] = rng();

    clock_t start = clock();
    radix_sort(ar.data(), n);
    /// Time taken = 0.107421 at -O2
    printf("Time taken to sort = %0.6f s\n", (clock() - start) / (double)CLOCKS_PER_SEC);

    for (i = 0; (i + 1) < n; i++){
        assert(ar[i] <= ar[i + 1]);
    }

    vector<unsigned int> tmp, big = {7, 4294967295u, 0, 256, 7, 65536}, small = {3, 1, 2};
    radix_sort(big.data(), big.size(), tmp);
    radix_sort(small.data(), small.size(), tmp);
    assert((big == vector<unsigned int>{0, 7, 7, 256, 65536, 4294967295u}));
    assert((small == vector<unsigned int>{1, 2, 3}));

    vector<int> signed_values = {3, -5, INT_MAX, 0, -1, INT_MIN, 7, -1};
    radix_sort(signed_values);
    assert((signed_values == vector<int>{INT_MIN, -5, -1, -1, 0, 3, 7, INT_MAX}));

    vector<int> none;
    radix_sort(none);
    assert(none.empty());

    return 0;
}
