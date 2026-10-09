/***
 *
 * Partition Numbers
 * p(0..n) modulo m, p(i) is the number of ways to write i as a sum of positive integers ignoring order
 *
 * Complexity: O(n sqrt n), O(n) memory
 *
 * partition_numbers(n, m) returns p(0), ..., p(n) modulo m, any modulus 1 <= m <= 2^62
 * Euler's pentagonal number theorem gives p(i) = sum over k >= 1 of (-1)^(k + 1) (p(i - k (3k - 1) / 2) + p(i - k (3k + 1) / 2)),
 * and only O(sqrt i) generalized pentagonal numbers are at most i
 *
 * Only additions and subtractions modulo m, so the modulus need not be prime
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

vector<long long> partition_numbers(int n, long long m){
    assert(n >= 0 && 1 <= m && m <= (1LL << 62));

    /// generalized pentagonal numbers k (3k - 1) / 2 and k (3k + 1) / 2 in increasing order, with the sign of their k
    vector<pair<long long, bool>> pentagonal;
    for (long long k = 1; k * (3 * k - 1) / 2 <= n; k++){
        pentagonal.push_back({k * (3 * k - 1) / 2, k & 1});
        pentagonal.push_back({k * (3 * k + 1) / 2, k & 1});
    }

    vector<long long> p(n + 1, 0);
    p[0] = 1 % m;
    for (int i = 1; i <= n; i++){
        long long sum = 0;
        for (auto [g, positive] : pentagonal){
            if (g > i) break;
            long long x = p[i - g];
            if (positive) sum = sum + x >= m ? sum + x - m : sum + x;
            else sum = sum >= x ? sum - x : sum - x + m;
        }
        p[i] = sum;
    }

    return p;
}

int main(){
    assert((partition_numbers(0, 1000000007) == vector<long long>{1}));
    assert((partition_numbers(1, 1000000007) == vector<long long>{1, 1}));
    assert((partition_numbers(20, 1000000007) == vector<long long>{1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42, 56, 77, 101, 135, 176, 231, 297, 385, 490, 627}));
    assert((partition_numbers(11, 2) == vector<long long>{1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0}));
    assert((partition_numbers(3, 1) == vector<long long>{0, 0, 0, 0}));

    assert(partition_numbers(100, 1000000007)[100] == 190569292);
    assert(partition_numbers(100, 7)[100] == 4);
    assert(partition_numbers(200, 1LL << 62)[200] == 3972999029388LL);
    assert(partition_numbers(1000, 1000000007)[1000] == 709496666);
    assert(partition_numbers(1000, 998244353)[1000] == 627356119);
    assert(partition_numbers(1000, 1LL << 62)[1000] == 3983852512315458295LL);
    assert(partition_numbers(1000, (1LL << 62) - 1)[1000] == 3983857729814604569LL);

    return 0;
}
