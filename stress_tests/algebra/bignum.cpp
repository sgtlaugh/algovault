// LINK: -lgmpxx -lgmp
#include "../common.h"
#include <gmpxx.h>

#define main library_main
#include "../../code_library/algebra/bignum.cpp"
#undef main

/// Random decimal strings with runs of 0s and 9s, which stress carries, borrows and Knuth D's correction steps
string random_number(int digits){
    if (digits <= 0) return "0";
    int mode = stress::rand_int(0, 3);
    string s(digits, '0');

    for (auto& c : s){
        if (mode == 0) c = '0' + stress::rand_int(0, 9);
        else if (mode == 1) c = stress::rand_int(0, 7) ? '9' : '0' + stress::rand_int(0, 9);
        else if (mode == 2) c = stress::rand_int(0, 7) ? '0' : '0' + stress::rand_int(0, 9);
        else c = stress::rand_int(0, 1) ? '9' : '0';
    }

    if (s[0] == '0') s[0] = '1' + stress::rand_int(0, 8);
    if (stress::rand_int(0, 3) == 0) s = string(stress::rand_int(1, 3), '0') + s;

    int sign = stress::rand_int(0, 7);
    return (sign < 4 ? "-" : (sign == 4 ? "+" : "")) + s;
}

mpz_class to_gmp(const string& s){
    return mpz_class(s[0] == '+' ? s.substr(1) : s, 10);
}

void check(const Bignum& got, const mpz_class& want){
    assert(got.to_string() == want.get_str());
}

int main(){
    /// Digit counts around the 18 row carry fold (162 digits) and the Karatsuba cutoff (96 limbs = 864 digits)
    const vector<int> sizes = {0, 1, 8, 9, 10, 18, 19, 153, 162, 171, 855, 864, 873, 1728, 1729, 3000};

    for (long long it = 0; it < stress::scaled(2000); it++){
        int da = sizes[stress::rand_int(0, sizes.size() - 1)] + stress::rand_int(0, 3);
        int db = stress::rand_int(0, 2) ? sizes[stress::rand_int(0, sizes.size() - 1)] : stress::rand_int(0, da + 1);
        if (it % 100 == 99) da = stress::rand_int(5000, 20000), db = stress::rand_int(1, da);

        string sa = random_number(da), sb = random_number(db);
        Bignum a(sa), b(sb);
        mpz_class x = to_gmp(sa), y = to_gmp(sb), q, r;

        check(a, x);
        check(-a, -x);
        assert(-(a - a) == 0 && !(-(a - a) < 0));  /// to_string prints "0" even for a negative zero, comparisons expose it
        check(a + b, x + y);
        check(a - b, x - y);
        check(a * b, x * y);
        assert((a < b) == (x < y) && (a == b) == (x == y) && (a > b) == (x > y));
        if (y != 0){
            mpz_tdiv_qr(q.get_mpz_t(), r.get_mpz_t(), x.get_mpz_t(), y.get_mpz_t());
            check(a / b, q);
            check(a % b, r);
        }

        long long s = (long long)(stress::rng()() >> stress::rand_int(0, 63)) * (stress::rand_int(0, 1) ? 1 : -1);
        if (it % 50 == 0) s = stress::rand_int(0, 1) ? LLONG_MIN : LLONG_MAX;
        mpz_class z;
        mpz_set_si(z.get_mpz_t(), s);
        check(a * s, x * z);
        check(s - a, z - x);
        if (s){
            mpz_tdiv_qr(q.get_mpz_t(), r.get_mpz_t(), x.get_mpz_t(), z.get_mpz_t());
            check(a / s, q);
            check(a % s, r);
        }
    }

    return 0;
}
