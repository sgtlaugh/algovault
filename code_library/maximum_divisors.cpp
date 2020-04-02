#include <stdio.h>


const int primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71};
const int p = sizeof(primes) / sizeof(int);


unsigned long long n, res, idx;


unsigned long long mul(unsigned long long a, unsigned long long b){
    unsigned long long res = 0;

    while (b){
        if (b & 1LL) res = (res + a);
        if (res > n) return 0;
        a = (a << 1LL);
        b >>= 1LL;
    }

    return res;
}

void backtrack(int i, int lim, unsigned long long val, unsigned long long r){
    if ((r > res) || (r == res && val < idx)) res = r, idx = val;
    if (i == p) return;

    int d;
    unsigned long long x = val;

    for (d = 1; d <= lim; d++){
        x = mul(x, primes[i]);
        if (x == 0) return;
        backtrack(i + 1, d, x, r * (d + 1));
    }
}

void print_smallest_num_with_max_div(long long limit) {
    /***
     * Prints the smallest number not exceeding limit with maximum number of divisors along with its count
     * Tested for n <= 10^18
    ***/

    res = 0;
    n = limit;
    backtrack(0, 100, 1, 1);
    printf("Smallest number is %llu with %llu divisors\n", idx, res);
}

int main(){
    print_smallest_num_with_max_div(1000000);               /// Smallest number is 720720 with 240 divisors
    print_smallest_num_with_max_div(1000000000000LL);       /// Smallest number is 963761198400 with 6720 divisors
    print_smallest_num_with_max_div(1000000000000000000LL); /// Smallest number is 897612484786617600 with 103680 divisors

    return 0;
}
