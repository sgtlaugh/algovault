#include <bits/stdc++.h>

using namespace std;

unsigned long long n, max_cnt, number;
const vector <int> primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71};

unsigned long long multiply(unsigned long long a, unsigned long long b){
    if ((long double)a * b > (n + 0.5)) return 0;
    return a * b;
}

void backtrack(int i, int lim, unsigned long long val, unsigned long long r){
    if (r > max_cnt || (r == max_cnt && val < number)) max_cnt = r, number = val;
    if (i == (int)primes.size()) return;

    unsigned long long x = val;
    for (int d = 1; d <= lim; d++){
        x = multiply(x, primes[i]);
        if (x == 0) return;
        backtrack(i + 1, d, x, r * (d + 1));
    }
}

void print_smallest_num_with_max_div(long long limit) {
    /***
     * Prints the smallest number not exceeding limit with maximum number of divisors along with its count
     * Tested for n <= 10^18
    ***/

    max_cnt = 0, n = limit;
    backtrack(0, 100, 1, 1);
    printf("Smallest number is %llu with %llu divisors\n", number, max_cnt);
}

int main(){
    print_smallest_num_with_max_div(1000000);               /// Smallest number is 720720 with 240 divisors
    print_smallest_num_with_max_div(1000000000000LL);       /// Smallest number is 963761198400 with 6720 divisors
    print_smallest_num_with_max_div(1000000000000000000LL); /// Smallest number is 897612484786617600 with 103680 divisors

    return 0;
}
