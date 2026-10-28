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

long double p, keep;

long double profit(int w, int l){
    long double win = gamblers_ruin(w, l, 1 - p);
    return win * w - (1 - win) * l * keep;
}

int main(){
    double x, pct;

    while (scanf("%lf %lf", &x, &pct) == 2){
        p = pct / 100, keep = 1 - x / 100;

        long double res = 0;
        if (p > 0){
            auto best = [](int w){
                int lo = 0, hi = 20540;
                while (hi - lo > 3){
                    int third = (hi - lo) / 3, a = lo + third, b = hi - third;
                    if (profit(w, a) < profit(w, b)) lo = a;
                    else hi = b;
                }
                long double res = 0;
                for (int l = lo; l <= hi; l++) res = max(res, profit(w, l));
                return res;
            };

            int lo = 0, hi = 2500;
            while (hi - lo > 3){
                int third = (hi - lo) / 3, a = lo + third, b = hi - third;
                if (best(a) < best(b)) lo = a;
                else hi = b;
            }
            for (int w = lo; w <= hi; w++) res = max(res, best(w));
        }
        printf("%0.9f\n", (double)res + 1e-12);
    }

    return 0;
}
