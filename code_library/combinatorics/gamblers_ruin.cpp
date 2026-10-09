/***
 *
 * Gambler's ruin problem (https://en.wikipedia.org/wiki/Gambler%27s_ruin)
 * Probability that the first player goes broke
 *
 * Complexity: O(1)
 *
 * First player has n1 coins, second player has n2 coins
 * After each move, first player wins with probability p and second player wins with probability q and (p + q = 1)
 * The loser gives 1 coin to the winner
 * When number of coins reaches 0, a player loses
 *
 * gamblers_ruin(n1, n2, p) returns the probability of the first player losing, for 0 <= n1, n2 <= INT_MAX, n1 + n2 > 0
 * Relative error stays below ~1e-15 for any n1, n2 and p in [0, 1], including p very close to 0.5 (results under ~1e-4900 underflow)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

long double gamblers_ruin(int n1, int n2, long double p){
    long double q = 1 - p, d = 2 * (p - 0.5L);  /// d = p - q, computed exactly
    if (n1 == 0) return 1;
    if (n2 == 0) return 0;
    if (d == 0) return (long double)n2 / ((long double)n1 + n2);

    /// With L = log of a ratio below 1, 1 - ratio^k = -expm1l(k * L) avoids cancellation when the ratio is near 1
    /// log1pl is accurate for a ratio near 1, logl for a ratio far from 1
    if (p < q){
        long double r = p / q, L = r < 0.5 ? logl(r) : log1pl(d / q);
        return expm1l(n2 * L) / expm1l(((long double)n1 + n2) * L);
    }

    long double s = q / p, L = s < 0.5 ? logl(s) : log1pl(-d / p);
    return expl(n1 * L) * expm1l(n2 * L) / expm1l(((long double)n1 + n2) * L);
}

int main(){
    auto close = [](long double got, long double expected){ return fabsl(got - expected) <= 1e-15L * expected; };

    assert(close(gamblers_ruin(3, 50, 0.49L), 0.9826119844357831L));
    assert(close(gamblers_ruin(1, 1, 0.9L), 0.1L));
    assert(close(gamblers_ruin(2, 2, 2.0L / 3), 0.2L));
    assert(close(gamblers_ruin(1, 2, 0.6L), 10.0L / 19));
    assert(close(gamblers_ruin(10, 10, 0.3L), 0.9997910023653128L));
    assert(close(gamblers_ruin(100, 1, 0.51L), 0.0007307293496802705L));
    assert(close(gamblers_ruin(3, 9, 0.5L), 0.75L));
    assert(close(gamblers_ruin(INT_MAX, INT_MAX, 0.5L), 0.5L));

    assert(gamblers_ruin(0, 7, 0.3L) == 1 && gamblers_ruin(5, 0, 0.7L) == 0);
    assert(gamblers_ruin(5, 7, 0.0L) == 1 && gamblers_ruin(5, 7, 1.0L) == 0);

    return 0;
}
