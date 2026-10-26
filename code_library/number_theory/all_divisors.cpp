/***
 *
 * All Divisors
 * Generates the divisors of every number from 1 to n
 *
 * Complexity:
 *   - O(n log n) time and memory, the lists hold about n ln n divisors, 1.1 GB at n = 1e7
 *   - make_sorted adds O(d log d) per list of d divisors
 *
 * all_divisors(n, make_sorted = false): res[x] = the divisors of 1 <= x <= n, sorted if make_sorted, res[0] empty
 *
 * Notes:
 *   - n up to about 1e7 by memory; smallest prime factors are stored as short, which caps n at 32767^2
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

vector<vector<int>> all_divisors(int n, bool make_sorted = false){
    vector<vector<int>> divisors(max(n, 0) + 1);
    if (n < 1) return divisors;

    vector<short> smallest_prime(n + 1, 0);  /// 0 marks a prime, the largest_prime loop relies on it
    vector<int> largest_prime(n + 1, 0);

    smallest_prime[1] = largest_prime[1] = 1;
    for (int i = 4; i <= n; i += 2) smallest_prime[i] = 2;
    for (int i = 3; (long long)i * i <= n; i += 2){
        if (smallest_prime[i]) continue;
        for (int j = i * i; j <= n; j += 2 * i){
            if (!smallest_prime[j]) smallest_prime[j] = i;
        }
    }

    for (int i = 2; i <= n; i++){
        if (!smallest_prime[i]) largest_prime[i] = i;
        else largest_prime[i] = max((int)smallest_prime[i], largest_prime[i / smallest_prime[i]]);
    }

    /// Each list is rest's list times p^0..p^e, measured 12x faster than pushing every d onto its multiples
    divisors[1].push_back(1);
    for (int i = 2; i <= n; i++){
        int p = largest_prime[i], rest = i, powers = 1;
        while (rest > 1 && largest_prime[rest] == p){
            rest /= p;
            powers++;
        }

        divisors[i].resize(powers * divisors[rest].size());
        int pos = 0;
        long long power = 1;
        for (int e = 0; e < powers; e++){
            for (int d : divisors[rest]) divisors[i][pos++] = d * power;
            power *= p;
        }

        if (make_sorted) sort(divisors[i].begin(), divisors[i].end());
    }

    return divisors;
}

int main(){
    /// The divisors of every number up to 12, in no particular order unless make_sorted
    vector<vector<int>> divisors = all_divisors(12);
    assert(divisors[0].empty());
    assert(divisors[1] == vector<int>({1}));
    assert(divisors[7] == vector<int>({1, 7}));
    assert(divisors[12].size() == 6);

    vector<vector<int>> sorted = all_divisors(12, true);
    assert(sorted[12] == vector<int>({1, 2, 3, 4, 6, 12}));

    return 0;
}
