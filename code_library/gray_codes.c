/***
 * 
 * https://en.wikipedia.org/wiki/Gray_code
 * 
 * Gray code is an ordering of the binary numeral system such that two successive values differ in only one bit
 * 
***/

#include <stdio.h>
#include <assert.h>

unsigned long long gray_code(unsigned long long x){
    return x ^ (x >> 1);
}

unsigned long long inverse_gray_code(unsigned long long x){
    unsigned long long inv = 0;
    for (inv = 0; x; x >>= 1) inv ^= x;
    return inv;
}

int main(){
	int n;
	for (n = 0; n < 1048576; n++){
		assert(inverse_gray_code(gray_code(n)) == n);
	}

    assert(gray_code(1000000007) == 643280644);
    assert(gray_code(1000000000000000003LL) == 797398725282889730LL);
    
    assert(inverse_gray_code(643280644) == 1000000007);
    assert(inverse_gray_code(797398725282889730LL) == 1000000000000000003LL);

    return 0;
}
