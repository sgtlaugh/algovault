/***
 *
 * Finds a number with maximum divisor count not exceeding the given limit
 * If there are multiple candidates, picks the smallest one
 * Returns the number as well as the maximum divisor count
 *
 * Tested and verified up to limit <= 10^18
 *
 * Table for n=10^d
 *
 *    1 =                      1              1
 *   10 =                      6              4
 *  100 =                     60             12
 *  1e3 =                    840             32
 *  1e4 =                   7560             64
 *  1e5 =                  83160            128
 *  1e6 =                 720720            240
 *  1e7 =                8648640            448
 *  1e8 =               73513440            768
 *  1e9 =              735134400           1344
 * 1e10 =             6983776800           2304
 * 1e11 =            97772875200           4032
 * 1e12 =           963761198400           6720
 * 1e13 =          9316358251200          10752
 * 1e14 =         97821761637600          17280
 * 1e15 =        866421317361600          26880
 * 1e16 =       8086598962041600          41472
 * 1e17 =      74801040398884800          64512
 * 1e18 =     897612484786617600         103680
 *
 * Check https://oeis.org/A066151 and https://oeis.org/A002182 for more info
 *
***/

#include <bits/stdc++.h>

using namespace std;

uint64_t limit, max_cnt, number;
const vector<int> primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71};

uint64_t multiply(uint64_t a, uint64_t b){
    if ((long double)a * b > limit + 0.5) return 0;
    return a * b;
}

void backtrack(int i, int lim, uint64_t cur, uint64_t cnt){
    if (cnt > max_cnt || (cnt == max_cnt && cur < number)) max_cnt = cnt, number = cur;
    if (i == (int)primes.size()) return;

    uint64_t x = cur;
    for (int c = 1; c <= lim && x; c++){
        x = multiply(x, primes[i]);
        if (x) backtrack(i + 1, c, x, cnt * (c + 1));
    }
}

pair<uint64_t, uint64_t> solve(long long n) {
    limit = n;
    max_cnt = 0;
    backtrack(0, 64, 1, 1);
    return {number, max_cnt};
}

int main(){
    typedef pair<uint64_t, uint64_t> Pair;

    assert(solve(1000000) == Pair(720720, 240));
    assert(solve(1000000000000LL) == Pair(963761198400LL, 6720));
    assert(solve(1000000000000000000LL) == Pair(897612484786617600LL, 103680));

    return 0;
}
