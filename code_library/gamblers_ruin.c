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

#define EPS 1e-12

long double expo(long double x, int n){
    if (n == 0) return 1;
    if (n & 1) return expo(x, n - 1) * x;
    long double v = expo(x, n >> 1);
    return v * v;

}

/***
 * 
 * Returns the probability of first player losing
 * Watch out for overflow/precission issues for bigger numbers
 * 
***/

long double gamblers_ruin(int n1, int n2, long double p){
    long double q = 1 - p;
    if (fabsl(p - q) < EPS){
        return (long double)n2 / (n1 + n2);
    }

    /// Always exponentiate a ratio <= 1 so large n1 + n2 underflows to 0 instead of overflowing to inf / inf
    if (p < q){
        long double r = p / q;
        return (1 - expo(r, n2)) / (1 - expo(r, n1 + n2));
    }

    long double s = q / p;
    return (expo(s, n1) - expo(s, n1 + n2)) / (1 - expo(s, n1 + n2));
}

int main(){
    const int fp_coins = 3;
    const int sp_coins = 50;
    const double fp_win_prob = 0.49;

    assert(fabsl(gamblers_ruin(fp_coins, sp_coins, fp_win_prob) - 0.98261198444) < 1e-9);
    assert(fabsl(gamblers_ruin(1, 1, 0.9) - 0.1) < 1e-9);

    return 0;
}
