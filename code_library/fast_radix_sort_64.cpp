/***
 * Fast radix sort in O(w * n) with loop unrolling
 * w is the key length, and is fixed to 8 bits
 * This enables sorting any list of integers in O(64/8 * n) or O(8 * n) passes
 *
 * Overall complexity: O(8 * n) ~ O(n)
 *
 * Most useful in scenarios where you want to sort a large list of items fast, usually in sub-optimal solutions
 * Can be generalized with templates and vectors but that usually makes it 2-3 x slower
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

const int MAXN = 50000010;

unsigned long long cnt[8][256], tmp[MAXN + 5];

void radix_sort(unsigned long long ar[], int n){
    assert(sizeof(unsigned long long) == 8);

    memset(cnt, 0, sizeof(cnt));
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


int main(){
    mt19937_64 rng(42);
    array <unsigned long long, MAXN> ar;

    puts("Generating array");
    int i, n = ar.size();
    for (i = 0; i < n; i++) ar[i] = rng();

    clock_t start = clock();
    radix_sort(ar.data(), n);
    printf("Time taken to sort = %0.6f s\n", (clock() - start) / (double)CLOCKS_PER_SEC);  /// Time taken = 0.876000

    for (i = 0; (i + 1) < n; i++){
        assert(ar[i] <= ar[i + 1]);
    }
    return 0;
}
