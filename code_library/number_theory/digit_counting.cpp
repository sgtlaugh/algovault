/***
 *
 * Digit Counting
 * How often each digit 0-9 appears when every integer of a range is written in decimal
 *
 * Complexity: O(log r) per call
 *
 * digit_counts(n): count of each digit 0-9 over the integers 0..n, all zero for n < 0
 * digit_counts(l, r): count of each digit 0-9 over the integers l..r, 0 <= l <= r
 *
 * Notes:
 *   - Leading zeros are not written, and 0 itself is the single digit 0
 *   - r <= 4e18, so every count fits long long
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

array<long long, 10> digit_counts(long long n){
    assert(n <= 4000000000000000000LL);
    array<long long, 10> res{};
    if (n < 0) return res;

    /// n / p / 10, not n / (10 * p): 10 * p overflows long long when n has 19 digits
    for (long long p = 1; ; p *= 10){
        long long high = n / p / 10, cur = n / p % 10, low = n % p;
        for (int d = 0; d < 10; d++){
            res[d] += high * p;
            if (d < cur) res[d] += p;
            if (d == cur) res[d] += low + 1;
        }
        if (p > 1) res[0] -= p;  /// the numbers below p have a leading zero here, except the units digit of 0

        if (n / p < 10) break;
    }

    return res;
}

array<long long, 10> digit_counts(long long l, long long r){
    assert(0 <= l && l <= r);
    array<long long, 10> res = digit_counts(r), below = digit_counts(l - 1);
    for (int d = 0; d < 10; d++) res[d] -= below[d];

    return res;
}

int main(){
    /// 1..13 writes 1 2 3 4 5 6 7 8 9 10 11 12 13
    array<long long, 10> c = digit_counts(1, 13);
    assert(c[0] == 1);
    assert(c[1] == 6);
    assert(c[2] == 2);
    assert(c[3] == 2);
    assert(c[9] == 1);

    /// 0..100: zeros in 0, 10, 20, ..., 90 and two in 100; ones in 10 units, 10 tens and 100
    array<long long, 10> upto = digit_counts(100);
    assert(upto[0] == 12);
    assert(upto[1] == 21);
    assert(upto[5] == 20);

    assert(digit_counts(0)[0] == 1);
    return 0;
}
