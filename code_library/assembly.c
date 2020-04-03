#include <stdio.h>


int popcount(int x){
    int counter = 0;
    __asm__ volatile("POPCNT %1, %0;"
        :"=r"(counter)
        :"r"(x)
        :
    );
    return counter;
}

int lzcount(int x){
    int counter = 0;
    __asm__ volatile("LZCNT %1, %0;"
        :"=r"(counter)
        :"r"(x)
        :
    );
    return counter;
}

int bsf(unsigned int x){
    __asm__ volatile("bsf %0, %0" : "=r" (x) : "0" (x));
    return x;
}

int gcd(int a, int b){
    int res;
    __asm__ volatile(
                      "movl %1, %%eax;"
                      "movl %2, %%ebx;"
                      "repeat_%=:\n"
                      "cmpl $0, %%ebx;"
                      "je terminate_%=\n;"
                      "xorl %%edx, %%edx;"
                      "idivl %%ebx;"
                      "movl %%ebx, %%eax;"
                      "movl %%edx, %%ebx;"
                      "jmp repeat_%=\n;"
                      "terminate_%=:\n"
                      "movl %%eax, %0;"

                      : "=g"(res)
                      : "g"(a), "g"(b)
                      : "eax", "ebx", "edx"

    );

    return res;
}

double fsqrt(double x){ /// Also works for long double
    __asm__ volatile("fsqrt" : "+t" (x));
    return x;
}

int main(){
    printf("%d\n", popcount(0));   /// 0
    printf("%d\n", popcount(13));  /// 3

    printf("%d\n", lzcount(0));    /// 32
    printf("%d\n", lzcount(100));  /// 25

    printf("%d\n", bsf(100));      /// 2
    printf("%d\n", bsf(15));       /// 0

    printf("%d\n", gcd(0, 10));    /// 10
    printf("%d\n", gcd(120, 180)); /// 60

    printf("%0.3f\n", fsqrt(1));   /// 1.000
    printf("%0.3f\n", fsqrt(50));  /// 7.071

    return 0;
}
