/***
 * 
 * Gambler's ruin problem (https://en.wikipedia.org/wiki/Gambler%27s_ruin)
 * 
 * First player has n1 coins, second player has n2 coins
 * After each move, first player wins with probability p and second player wins with probability q and (p + q = 1)
 * The loser gives 1 coin to the winner
 * When number of coins reaches 0, a player loses
 * 
***/

#include <stdio.h>
#include <math.h>
#include <assert.h>

/***
 *
 * Returns the probability of first player losing
 * Relative error stays around 1e-16 for any n1, n2 and p, including p very close to 0.5
 *
***/

long double gamblers_ruin(int n1, int n2, long double p){
    long double q = 1 - p, d = 2 * (p - 0.5L);  /// d = p - q, computed exactly
    if (n1 == 0) return 1;
    if (n2 == 0) return 0;
    if (d == 0) return (long double)n2 / (n1 + n2);

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
    const int fp_coins = 3;
    const int sp_coins = 50;
    const double fp_win_prob = 0.49;

    assert(fabsl(gamblers_ruin(fp_coins, sp_coins, fp_win_prob) - 0.98261198444) < 1e-9);
    assert(fabsl(gamblers_ruin(1, 1, 0.9) - 0.1) < 1e-9);

    return 0;
}
