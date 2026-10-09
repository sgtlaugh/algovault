#include "../common.h"

#define main library_main
#include "../../code_library/algebra/subset_convolution.cpp"
#undef main

long long brute_reduce(long long x, int mod){
    return ((x % mod) + mod) % mod;
}

/// Every pair (S, T subset of S) through submask enumeration, O(3^n)
void brute_all(const vector<int>& a, const vector<int>& b, int mod, vector<int>& conv, vector<int>& sub, vector<int>& sub_inv, vector<int>& sup, vector<int>& sup_inv){
    int size = a.size();
    conv.assign(size, 0), sub.assign(size, 0), sub_inv.assign(size, 0), sup.assign(size, 0), sup_inv.assign(size, 0);

    for (int s = 0; s < size; s++){
        for (int t = s; ; t = (t - 1) & s){
            long long sign = __builtin_popcount(s ^ t) & 1 ? -1 : 1;
            conv[s] = brute_reduce(conv[s] + brute_reduce(a[t], mod) * brute_reduce(b[s ^ t], mod), mod);
            sub[s] = brute_reduce(sub[s] + brute_reduce(a[t], mod), mod);
            sub_inv[s] = brute_reduce(sub_inv[s] + sign * brute_reduce(a[t], mod), mod);
            sup[t] = brute_reduce(sup[t] + brute_reduce(a[s], mod), mod);
            sup_inv[t] = brute_reduce(sup_inv[t] + sign * brute_reduce(a[s], mod), mod);
            if (t == 0) break;
        }
    }
}

int random_value(int mod){
    int kind = stress::rand_int(0, 5);
    if (kind == 0) return stress::rand_int(0, 1) ? INT_MIN : INT_MAX;
    if (kind == 1) return stress::rand_int(-3, 3);
    if (kind == 2) return mod - stress::rand_int(1, min(mod, 3));
    return stress::rand_int(INT_MIN, INT_MAX);
}

void check(int n, int mod){
    int size = 1 << n;
    vector<int> a(size), b(size);
    for (auto& x : a) x = random_value(mod);
    for (auto& x : b) x = random_value(mod);

    vector<int> conv, sub, sub_inv, sup, sup_inv;
    brute_all(a, b, mod, conv, sub, sub_inv, sup, sup_inv);

    assert(subset_convolution(a, b, mod) == conv);

    vector<int> f = a;
    subset_zeta(f, mod);
    assert(f == sub);
    subset_mobius(f, mod);
    for (int i = 0; i < size; i++) assert(f[i] == brute_reduce(a[i], mod));

    f = a;
    subset_mobius(f, mod);
    assert(f == sub_inv);

    f = a;
    superset_zeta(f, mod);
    assert(f == sup);
    superset_mobius(f, mod);
    for (int i = 0; i < size; i++) assert(f[i] == brute_reduce(a[i], mod));

    f = a;
    superset_mobius(f, mod);
    assert(f == sup_inv);
}

/// Every transform and the subset convolution against direct sums over all submasks
int main(){
    const int mods[] = {1, 2, 3, 6, 1024, 998244353, 1000000007, 2147483646, 2147483647};

    for (long long it = 0; it < stress::scaled(1500); it++){
        int n = it < 300 ? it % 11 : stress::rand_int(0, 10);
        int mod = stress::rand_int(0, 4) ? mods[stress::rand_int(0, 8)] : stress::rand_int(1, INT_MAX);
        check(n, mod);
    }

    for (int mod : {998244353, 2147483647}) check(14, mod);

    return 0;
}
