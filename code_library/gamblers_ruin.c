#include <stdio.h>
#include <math.h>

#define EPS 1e-9

long double expo(long double x, int n){
    if (n == 0) return 1;
    if (n & 1) return expo(x, n - 1) * x;
    else{
        long double res = expo(x, n >> 1);
        return res * res;
    }
}

/***
 * Gambler's ruin problem
 * First player has n1 coins, second player has n2 coins
 * After each move, first player wins with probability p and second player wins with probability q and (p + q = 1)
 * The loser gives 1 coin to the winner
 * When number of coins reaches 0, a player loses
 * Returns the probability of first player losing
***/

long double gamblers_ruin(int n1, int n2, long double p, long double q){
    if (fabs(p - q) < EPS){
        return (long double)n2 / (n1 + n2);
    }

    long double x = 1.0 - expo(q / p, n2);
    long double y = 1.0 - expo(q / p, n1 + n2);
    return (x / y);
}

int main(){
    printf("%0.5f\n", (double)gamblers_ruin(3, 50, 0.49, 0.51));   /// 0.87148
    return 0;
}
