/***
 *
 * Maximum Divisors
 * Smallest number not exceeding a limit with the most divisors, and its divisor count
 *
 * Complexity: O(H) for H the count of x <= limit with non-increasing exponents over consecutive primes (OEIS A025487)
 *   32749 calls and a recursion depth of 16 at limit = 10^18
 *
 * max_divisors(limit): {smallest x <= limit with the most divisors, its divisor count}, 1 <= limit <= 1e18
 *
 * Notes:
 *   - The answers for limit = 10^d, from https://oeis.org/A066151 and https://oeis.org/A002182:
 *
 *        1 =                      1              1
 *       10 =                      6              4
 *      100 =                     60             12
 *      1e3 =                    840             32
 *      1e4 =                   7560             64
 *      1e5 =                  83160            128
 *      1e6 =                 720720            240
 *      1e7 =                8648640            448
 *      1e8 =               73513440            768
 *      1e9 =              735134400           1344
 *     1e10 =             6983776800           2304
 *     1e11 =            97772875200           4032
 *     1e12 =           963761198400           6720
 *     1e13 =          9316358251200          10752
 *     1e14 =         97821761637600          17280
 *     1e15 =        866421317361600          26880
 *     1e16 =       8086598962041600          41472
 *     1e17 =      74801040398884800          64512
 *     1e18 =     897612484786617600         103680
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

pair<long long, long long> max_divisors(long long limit){
    assert(limit >= 1);
    static const int primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71};
    long long best = 1, best_cnt = 1;

    /// exponents never increase from one prime to the next, any other shape has a smaller number with as many divisors
    auto search = [&](auto& self, int i, int max_e, long long x, long long cnt) -> void{
        if (cnt > best_cnt || (cnt == best_cnt && x < best)) best = x, best_cnt = cnt;
        for (int e = 1; i < 20 && e <= max_e && x <= limit / primes[i]; e++) self(self, i + 1, e, x *= primes[i], cnt * (e + 1));
    };
    search(search, 0, 64, 1, 1);
    return {best, best_cnt};
}

int main(){
    /// The smallest number with the most divisors up to each limit, with its divisor count
    typedef pair<long long, long long> Pair;
    assert(max_divisors(1) == Pair(1, 1));
    assert(max_divisors(2) == Pair(2, 2));
    assert(max_divisors(12) == Pair(12, 6));
    assert(max_divisors(1000000) == Pair(720720, 240));
    assert(max_divisors(1000000000000LL) == Pair(963761198400LL, 6720));
    assert(max_divisors(1000000000000000000LL) == Pair(897612484786617600LL, 103680));

    return 0;
}
