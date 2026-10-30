/***
 *
 * Radix Sort (64-bit)
 * LSD radix sort with 8-bit digits, 8 counting-sort passes over 64-bit integers, loop unrolled
 *
 * Complexity: O(8 n + 2048), a scratch vector of n unsigned long longs (80 MB at n = 1e7)
 *
 * radix_sort(ar, n): sorts unsigned long longs ar[0, n) with a scratch vector allocated per call
 * radix_sort(v): sorts a vector<long long>, negatives first
 * radix_sort(ar, n, tmp): sorts with the caller's tmp as scratch, grown to n if smaller
 *
 * Notes:
 *   - Pass the same tmp to repeated sorts: at n = 1e7 a fresh scratch measured 270 ms per sort, a reused one 229 ms
 *   - Use radix_sort.cpp for 32-bit keys, a template over the key type made it 2-3x slower
 *   - Packing a record's fields into one key, x << 18 | rank << 1 | sign, sorts 2x faster than std::sort of tuples
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

/// Flipping the sign bit maps long long order onto unsigned order
void radix_sort(vector<long long>& v){
    for (auto& x : v) x ^= LLONG_MIN;
    radix_sort((unsigned long long*)v.data(), v.size());
    for (auto& x : v) x ^= LLONG_MIN;
}

int main(){
    mt19937_64 rng(42);
    int i, n = 10000000;
    vector<unsigned long long> ar(n);

    puts("Generating array");
    for (i = 0; i < n; i++) ar[i] = rng();

    clock_t start = clock();
    radix_sort(ar.data(), n);
    /// Time taken = 0.302028 at -O2
    printf("Time taken to sort = %0.6f s\n", (clock() - start) / (double)CLOCKS_PER_SEC);

    for (i = 0; (i + 1) < n; i++){
        assert(ar[i] <= ar[i + 1]);
    }

    vector<unsigned long long> tmp, big = {7, 18446744073709551615ull, 0, 256, 7, 65536}, small = {3, 1, 2};
    radix_sort(big.data(), big.size(), tmp);
    radix_sort(small.data(), small.size(), tmp);
    assert((big == vector<unsigned long long>{0, 7, 7, 256, 65536, 18446744073709551615ull}));
    assert((small == vector<unsigned long long>{1, 2, 3}));

    vector<long long> signed_values = {3, -5, LLONG_MAX, 0, -1, LLONG_MIN, 1LL << 40, -(1LL << 40), -1};
    radix_sort(signed_values);
    assert((signed_values == vector<long long>{LLONG_MIN, -(1LL << 40), -5, -1, -1, 0, 3, 1LL << 40, LLONG_MAX}));

    vector<long long> none;
    radix_sort(none);
    assert(none.empty());

    return 0;
}
