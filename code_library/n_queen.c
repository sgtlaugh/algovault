#include <stdio.h>

int n, lim, counter;

void backtrack(int i, int c, int l, int r){
    if (!i){
        counter++;
        return;
    }

    int bitmask, x;
    --i, bitmask = lim & ~(l | r | c);

    while (bitmask){
        x = (-bitmask & bitmask);
        if (!x) return;
        bitmask ^= x;
        backtrack(i, c | x, (l | x) << 1, (r | x) >> 1);
    }
}

int count_ways(int dimension){
	n = dimension;
	counter = 0, lim = (1 << n) - 1;
	backtrack(n, 0, 0, 0);

	return counter;
}

int main(){
    printf("%d\n", count_ways(8));   /// 92
    printf("%d\n", count_ways(13));  /// 73712
    printf("%d\n", count_ways(15));  /// 2279184

    return 0;
}
