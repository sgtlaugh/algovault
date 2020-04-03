#include <stdio.h>


//// https://en.wikipedia.org/wiki/Gray_code
///  Gray code is an ordering of the binary numeral system such that two successive values differ in only one bit (binary digit)

long long gray_code(long long x){
    return x ^ (x >> 1);
}

long long inverse_gray_code(long long x){
    long long h = 1, res = 0;
    do{
        if (x & 1) res ^= h;
        x >>= 1, h = (h << 1) + 1;
    } while (x);
    return res;
}

int main(){
    printf("%lld\n", gray_code(100));                /// 86
    printf("%lld\n", gray_code(1000000007));         /// 643280644

    printf("%lld\n", inverse_gray_code(86));         /// 100
    printf("%lld\n", inverse_gray_code(643280644));  /// 1000000007

    return 0;
}
